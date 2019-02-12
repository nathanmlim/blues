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
from blues.reporters import NetCDF4Reporter, NetCDF4Storage

finfo = np.finfo(np.float32)
rtol = finfo.precision
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logging.getLogger("parmed").setLevel(logging.ERROR)
logging.getLogger("openmmtools.alchemy").setLevel(logging.ERROR)
np.random.RandomState(seed=3134)
logging.basicConfig(format='%(asctime)s | %(levelname)s : %(message)s', level=logging.INFO, stream=sys.stdout)

# Define parameters
temperature = 300 * unit.kelvin
collision_rate = 1 / unit.picoseconds
timestep = 4.0 * unit.femtoseconds
n_steps = 1000
reportInterval = 250
nIter=2
context_cache = cache.ContextCache(capacity=1)
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
ncmc_move = RandomLigandRotationMove(atom_subset=alchemical_atoms, context_cache=context_cache, n_steps=n_steps)
langevin_move = ModLangevinDynamicsMove(timestep, collision_rate, n_steps,
                                     reassign_velocities=True,
                                    context_cache=context_cache)

# Create a reporter
with open('test.pdb', 'w') as pdb:
    openmm.app.pdbfile.PDBFile.writeFile(tol.topology, tol.positions, pdb)
filename = 'test.nc'
if os.path.exists(filename):
    os.remove(filename)
else:
    print("Sorry, I can not remove %s file." % filename)
#nc_reporter = NetCDF4Reporter(filename, reportInterval)
nc_reporter = NetCDF4Storage(filename, reportInterval)

sampler = NCMCSampler(alchemical_atoms, thermodynamic_state, alch_thermodynamic_state,
                 sampler_state, move=[langevin_move, ncmc_move], platform=None,
                 reporter=nc_reporter, topology=tol.topology)


sampler.equil()
sampler.run(nIter)
