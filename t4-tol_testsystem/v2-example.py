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

finfo = np.finfo(np.float32)
rtol = finfo.precision
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logging.getLogger("parmed").setLevel(logging.ERROR)
logging.getLogger("openmmtools.alchemy").setLevel(logging.ERROR)
np.random.RandomState(seed=3134)
#logging.basicConfig(format='%(asctime)s | %(levelname)s : %(message)s', level=logging.INFO, stream=sys.stdout)


from blues import utils
parser = argparse.ArgumentParser(description='Restart file name')
parser.add_argument('-j', '--jobname', type=str, help="store jobname")
parser.add_argument('-n', '--nIter', default=1, type=int, help="number of Iterations")
parser.add_argument('-s', '--nsteps', default=100, type=int, help="number of steps")
parser.add_argument('-r', '--reportInterval', default=10, type=int, help="reportInterval")
args = parser.parse_args()

# Define parameters
outfname = args.jobname
temperature = 300 * unit.kelvin
collision_rate = 1 / unit.picoseconds
timestep = 4.0 * unit.femtoseconds
n_steps = args.nsteps
reportInterval = args.reportInterval
nIter = args.nIter

context_cache = cache.ContextCache(capacity=4)
prmtop = utils.get_data_filename('blues', 'tests/data/eqToluene.prmtop') #TOL-parm
inpcrd = utils.get_data_filename('blues', 'tests/data/eqToluene.inpcrd')
tol = parmed.load_file(prmtop, xyz=inpcrd)
tol.system = tol.createSystem(nonbondedMethod=openmm.app.PME,
                             nonbondedCutoff=10*unit.angstrom,
                             constraints=openmm.app.HBonds,
                             hydrogenMass=3.024*unit.dalton,
                             rigidWater=True,removeCMMotion=True,
                             flexibleConstraints=True, splitDihedrals=False
                             )

# Create our State objects
sampler_state = SamplerState(positions=tol.positions)
thermodynamic_state = ThermodynamicState(system=tol.system, temperature=temperature)

# Create our AlchemicalState
alchemical_atoms = utils.atomIndexfromTop('LIG',tol.topology)
toluene_alchemical_system = SystemFactory.generateAlchSystem(tol.system, alchemical_atoms)
alchemical_state = alchemy.AlchemicalState.from_system(toluene_alchemical_system)
alch_thermodynamic_state = ThermodynamicState(system=toluene_alchemical_system, temperature=temperature)
alch_thermodynamic_state = CompoundThermodynamicState(alch_thermodynamic_state, composable_states=[alchemical_state])
alch_thermodynamic_state.topology = tol.topology

# Iniitialize our Move set
ncmc_move = RandomLigandRotationMove(atom_subset=alchemical_atoms, timestep=timestep, context_cache=context_cache, n_steps=n_steps)
langevin_move = ModLangevinDynamicsMove(timestep, collision_rate, n_steps,
                                     reassign_velocities=True,
                                    context_cache=context_cache)

# Create a reporter
with open(outfname+'.pdb', 'w') as pdb:
    openmm.app.pdbfile.PDBFile.writeFile(tol.topology, tol.positions, pdb)
filename = outfname+'.nc'
if os.path.exists(filename):
    os.remove(filename)
else:
    print("Sorry, I can not remove %s file." % filename)
nc_reporter = NetCDF4Storage(filename, reportInterval)
from blues.reporters import init_logger
logger = init_logger(logger, outfname=outfname)
state_reporter = BLUESStateDataStorage(file=logger, step=True, reportInterval=reportInterval, progress=True, remainingTime=True, speed=True, totalSteps=n_steps, currentIter=False)

sampler = NCMCSampler(alchemical_atoms, thermodynamic_state, alch_thermodynamic_state,
                 sampler_state, move=[langevin_move, ncmc_move], platform=None,
                 reporter=[nc_reporter, state_reporter], topology=tol.topology)


sampler.minimize()
sampler.run(nIter)
