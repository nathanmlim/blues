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
#new

platform = Platform.getPlatformByName('CUDA')
struct1 = parmed.load_file('sqB.pdb')
struct2 = parmed.load_file('eth.prmtop', xyz='eth.inpcrd')
struct = struct1 + struct2
system = struct.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds, removeCMMotion=True)
system.removeForce(4)
nonbonded = system.getForce(3)
particles = system.getNumParticles()
parameter_list = [nonbonded.getParticleParameters(i) for i in range(particles)]
#for i in range(particles):
#    print(nonbonded.getParticleParameters(i))
print(parameter_list)
system.removeForce(3)

#if you want to make the system use lambda parameters
pairwiseForce = CustomNonbondedForce("q/(r^2) + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2)*lambda_sterics; epsilon=sqrt(epsilon1*epsilon2)*lambda_electrostatics; q = lambda_charge*(q1*q2)")
#else just use normal parameters
#pairwiseForce = CustomNonbondedForce("q/(r^2) + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2); q = q1*q2")


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
for i, parameter in enumerate(parameter_list):
    print('param', pairwiseForce.getParticleParameters(i))
for i in range(pairwiseForce.getNumPerParticleParameters()):
    print(i, pairwiseForce.getPerParticleParameterName(i))


#centroid force

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
xml = XmlSerializer.serialize(system)
with open('test_system.xml', 'w') as f:
    f.write(xml)
#temp at 200 gives good sampling
temp = 200*kelvin
print(list(enumerate(system.getForces())))
integrator = LangevinIntegrator(temp, 1.0*1/picoseconds, 0.001*picoseconds)
#context = Context(system, integrator)
sim = Simulation(struct.topology, system, integrator, platform)
sim.context.setPositions(struct.positions)
print(sim.context.getState(getPositions=True).getPositions())

sim.minimizeEnergy()
totalSteps = 100000000
min_pos = sim.context.getState(getPositions=True).getPositions()
sim.context.setVelocitiesToTemperature(temp)
#reporter = DCDReporter('output.dcd', 1000)
#nc_reporter = NetCDF4Reporter('output.nc', 25000)
state_reporter = statedatareporter.StateDataReporter(sys.stdout, 250000, step=True, speed=True, progress=True, totalSteps=totalSteps, remainingTime=True, separator='\t')
#sim.reporters.append(nc_reporter)
sim.reporters.append(state_reporter)
#if running test
sim.step(totalSteps)
struct.save('min.pdb')

print('done')
print(struct)
