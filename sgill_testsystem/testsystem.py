
# coding: utf-8

# In[5]:

import argparse
from simtk import unit, openmm
from openmmtools import cache, alchemy
from openmmtools.states import SamplerState, ThermodynamicState, CompoundThermodynamicState
from openmmtools import storage

from blues import utils
import parmed
import logging
import os, sys, copy
import numpy as np
from blues.moves import RandomLigandRotationMove, ModLangevinDynamicsMove
from blues.simulation import NCMCSampler, SystemFactory
from blues.reporters import NetCDF4Reporter, NetCDF4Storage, BLUESStateDataStorage

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
#from ncmc_switching import *
import mdtraj as md
from openmmtools import testsystems
import math
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalState
import numpy as np
from mdtraj.reporters import HDF5Reporter
from blues.reporters import *
import parmed, sys
from simtk.openmm.app.simulation import Simulation
import mdtraj as md
import matplotlib.pyplot as plt


finfo = np.finfo(np.float32)
rtol = finfo.precision
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logging.getLogger("parmed").setLevel(logging.ERROR)
logging.getLogger("openmmtools.alchemy").setLevel(logging.ERROR)
np.random.RandomState(seed=3134)

from blues.reporters import init_logger
#logging.basicConfig(format='%(asctime)s | %(levelname)s : %(message)s', level=logging.INFO, stream=sys.stdout)


# In[6]:


def runTestSystem(outfname, nsteps, nIter, reportInterval):
    # Define parameters
    temperature = 200 * unit.kelvin
    collision_rate = 1 / unit.picoseconds
    timestep = 1.0 * unit.femtoseconds
    n_steps = nsteps
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    logging.getLogger("parmed").setLevel(logging.ERROR)
    logging.getLogger("openmmtools.alchemy").setLevel(logging.ERROR)
    logger = init_logger(logger, outfname=outfname)
    #reportInterval = 25
    #nIter=100

    context_cache = cache.ContextCache(capacity=4)
    struct1 = parmed.load_file('/export/home/limn1/beegfs-data/molssi-project/blues/sgill_testsystem/sqB.pdb')
    struct2 = parmed.load_file('/export/home/limn1/beegfs-data/molssi-project/blues/sgill_testsystem/eth.prmtop', xyz='/export/home/limn1/beegfs-data/molssi-project/blues/sgill_testsystem/eth.inpcrd')
    struct = struct1 + struct2
    system = struct.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds, removeCMMotion=True)
    system.removeForce(4)

    nonbonded = system.getForce(3)
    particles = system.getNumParticles()
    parameter_list = [nonbonded.getParticleParameters(i) for i in range(particles)]
    #for i in range(particles):
    #    print(nonbonded.getParticleParameters(i))
    #print(parameter_list)
    system.removeForce(3)

    #if you want to make the system use lambda parameters
    pairwiseForce = CustomNonbondedForce("q/(r^2) + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2)*lambda_sterics; epsilon=sqrt(epsilon1*epsilon2)*lambda_electrostatics; q = lambda_charge*(q1*q2)")
    #else just use normal parameters
    #pairwiseForce = CustomNonbondedForce("q/(r^2) + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2); q = q1*q2"
    pairwiseForce.addPerParticleParameter("sigma")
    pairwiseForce.addPerParticleParameter("epsilon")
    pairwiseForce.addPerParticleParameter("q")
    pairwiseForce.addPerParticleParameter("lambda_on")
    #pairwiseForce.addPerParticleParameter("lambda_on")
    pairwiseForce.addGlobalParameter("lambda_sterics", 1)
    pairwiseForce.addGlobalParameter("lambda_electrostatics", 1)
    pairwiseForce.addGlobalParameter("lambda_charge", 1)

    for i, parameter in enumerate(parameter_list):
        new_param = [parameter[1]*1.2] + [parameter[2]]+[parameter[0]]
        pairwiseForce.addParticle()
        pairwiseForce.setParticleParameters(i,new_param+[0])
    pairwiseForce.setParticleParameters(0, [0.324999852378,0.71128, -0.2, 10])
    pairwiseForce.setParticleParameters(1, [0.324999852378,0.71128, -0.5, 10])
    numParticles = system.getNumParticles()
    rangeparticles = range(numParticles)
    pairwiseForce.addInteractionGroup([0,1], rangeparticles[2:])
    num_params = pairwiseForce.getNumPerParticleParameters()
    system.addForce(pairwiseForce)
#     for i, parameter in enumerate(parameter_list):
#         print('param', pairwiseForce.getParticleParameters(i))
#     for i in range(pairwiseForce.getNumPerParticleParameters()):
#         print(i, pairwiseForce.getPerParticleParameterName(i))


    centoridForce = CustomCentroidBondForce(2, ('0.5*k*distance(g1,g2)^2'))
    centoridForce.addPerBondParameter('k')
    centoridForce.addGroup([0,1], [1,1])
    centoridForce.addGroup([2,3])
    bondGroups = [0,1]
    bondParameters=[100000]
    centoridForce.addBond(bondGroups, bondParameters)
    system.addForce(centoridForce)
    system.setParticleMass(0, 0)
    system.setParticleMass(1, 0)

    # Create our State objects
    sampler_state = SamplerState(positions=struct.positions)
    thermodynamic_state = ThermodynamicState(system=system, temperature=temperature)

    # Create our AlchemicalState
    alchemical_atoms = utils.atomIndexfromTop('TMP',struct.topology)
    toluene_alchemical_system = SystemFactory.generateAlchSystem(system, alchemical_atoms)
    alchemical_state = alchemy.AlchemicalState.from_system(toluene_alchemical_system)
    alch_thermodynamic_state = ThermodynamicState(system=toluene_alchemical_system, temperature=temperature)
    alch_thermodynamic_state = CompoundThermodynamicState(alch_thermodynamic_state, composable_states=[alchemical_state])
    alch_thermodynamic_state.topology = struct.topology

    # Iniitialize our Move set
    ncmc_move = RandomLigandRotationMove(atom_subset=alchemical_atoms, timestep=timestep, context_cache=context_cache, n_steps=n_steps)
    langevin_move = ModLangevinDynamicsMove(timestep, collision_rate, n_steps,
                                         reassign_velocities=True,
                                        context_cache=context_cache)

    # Create a reporter
    with open(outfname+'.pdb', 'w') as pdb:
        pdbfile.PDBFile.writeFile(struct.topology, struct.positions, pdb)
    filename = outfname+'.nc'
    if os.path.exists(filename):
        os.remove(filename)
    else:
        print("Sorry, I can not remove %s file." % filename)
    nc_reporter = NetCDF4Storage(filename, reportInterval)

    state_reporter = BLUESStateDataStorage(file=logger, step=True, reportInterval=reportInterval, progress=True, remainingTime=True, speed=True, totalSteps=n_steps, currentIter=False)

    sampler = NCMCSampler(alchemical_atoms, thermodynamic_state, alch_thermodynamic_state,
                     sampler_state, move=[langevin_move, ncmc_move], platform=None,
                     reporter=[nc_reporter, state_reporter], topology=struct.topology)




    sampler.minimize()
    sampler.run(nIter)

    acceptance = sampler.n_accepted / sampler.n_proposed


    traj = md.load(outfname+'.nc', top=outfname+'.pdb')
    print(traj)
    dist = md.compute_distances(traj, [[0,2]])
    count_list = [0 if i <=0.49 else 1 for i in dist]
    counts = (count_list.count(0), count_list.count(1))
    #plt.figure(1)
    #plt.subplot(211)
    #plt.plot(dist)
    #plt.subplot(212)
    #plt.bar([0,1], counts)
    total = counts[0] + counts[1]
    counts0 = counts[0]/total
    counts1 = counts[1]/total
    #plt.legend(['{:.2f} {:.2f}'.format(counts0,counts1)])
    #plt.show()

    return (acceptance, counts0, counts1)


# In[7]:

parser = argparse.ArgumentParser(description='Restart file name')
parser.add_argument('-j', '--jobname', type=str, help="store jobname")
parser.add_argument('-n', '--nIter', default=1, type=int, help="number of Iterations")
parser.add_argument('-s', '--nsteps', default=100, type=int, help="number of steps")
parser.add_argument('-r', '--reportInterval', default=10, type=int, help="reportInterval")
args = parser.parse_args()

outfname = '%s' % args.jobname
counts = runTestSystem(outfname, args.nsteps, args.nIter, args.reportInterval)
print("{},{},{}".format(counts[0],counts[1],counts[2]))
