from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Function to add backbone position restraints
def add_backbone_posres(system, positions, atoms, restraint_force):
  force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
  force_amount = restraint_force * kilocalories_per_mole/angstroms**2
  force.addGlobalParameter("k", force_amount)
  force.addPerParticleParameter("x0")
  force.addPerParticleParameter("y0")
  force.addPerParticleParameter("z0")
  for i, (atom_crd, atom) in enumerate(zip(positions, atoms)):
    if atom.name in  ('CA', 'C', 'N'):
      force.addParticle(i, atom_crd.value_in_unit(nanometers))
  posres_sys = copy.deepcopy(system)
  posres_sys.addForce(force)
  return posres_sys

# Set up 
job_name = 'WAT'
prmtop = AmberPrmtopFile('water_box.prmtop')
inpcrd = AmberInpcrdFile('water_box.inpcrd')
platform = Platform.getPlatformByName('CUDA')

# Simulation Parameters
temperature = 300 * kelvin # Final Room Temcperature
friction_coefficient = 2.0 / picosecond # Collision Frequency
step_size = 1 * femtoseconds # Slower the better but too small too slow sims
pressure = 1 * bar # STP Conditions

#simulation options minimization, qht: quick heating, peq: pressure eq, mdr: production temp eq
# Minimization Parameters
min_length = 10000 # = Amber's maxcyc = 10000, ! number of minimization steps

# NVT: Quick Heating Parameters
initial_temp = 5.0 * kelvin # Initial Box Temperature
qht_run = 50000 # Quick Heating for 50000fs = 50ps
qht_log_int = 1000 # Log Simulation every 1000ps = 1ps

# NPT: Pressure Equilibrations Parameters
peq_run = 500000 # = Amber's nstlim   = 500000, ! number of steps, : 500 ps NPT Run 
peq_trj_int = 10000 # = Amber's ntwx = 10000, ! write coordinates every 10th step
peq_log_int = 10000 # = Amber's ntpr = 10000, ! logfile print frequency
peq_res_int = 10000 # = Amber's ntwr   = 10000,  ! write restart file at last step

# NVT: Production Run Parameters
prod_run=10000000 # 1E+fs = 10000fs = 10ns : Total Production Run
prod_trj_int=1000 # Save Trajectory every 1000fs = 1ps: 
prod_log_int=1000 # Log Simulations every 1000fs = 1ps: 
prod_res_int=10000 # Save Restart file every 10000 = 10ps: 
print('This simulation will minimise the system, before pressure equilibration (npt) for', peq_run*step_size,\
      'and then a production run for',prod_run*step_size)

# Output files
# Minimization
min_data_filename = job_name + '_min.csv'
min_trajectory_filename = job_name + '_min_traj.pdb'


# NVT: Quick Heating
qht_data_filename = job_name + '_qht.csv'
qht_trajectory_filename = job_name + '_qht_traj.pdb'
qht_state_filename = job_name + '_qht_traj.state'
qht_restart_filename = job_name + '_qht_restart.chk'

# NPT: pressure equilibrations
peq_data_filename = job_name + '_peq.csv'
peq_trajectory_filename = job_name + '_peq_traj.pdb'
peq_state_filename = job_name + '_peq_traj.state'
peq_restart_filename = job_name + '_peq_restart.chk'

# NVT: production run output
prod_data_filename = job_name + '_prod_traj.csv'
prod_trajectory_filename = job_name + '_prod_traj'
prod_state_filename = job_name + '_prod_traj.state'
prod_restart_filename = job_name + '_prod_traj.chk'

#System preperation 
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)

posres_sys = add_backbone_posres(system, inpcrd.positions, prmtop.topology.atoms(), 100)
integrator = LangevinIntegrator(temperature, friction_coefficient, step_size)
integrator.setConstraintTolerance(0.00001)
simulation = Simulation(prmtop.topology, posres_sys, integrator, platform)
simulation.reporters.append(
    StateDataReporter(min_data_filename,
                      10000,
                      step=True,
                      potentialEnergy=True,
                      temperature=True,
                      volume=True,
                      density=True,
                      separator=','))

# Minimize
print('Minimizing...')
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
position = simulation.context.getState(getPositions=True).getPositions()
energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
PDBFile.writeFile(simulation.topology, position,
                  open(min_trajectory_filename, 'w'))
print ('Energy at Minimum is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
simulation.step(min_length)

# WarmUp with a NVT run.  Slowly warm up temperature - every 1000 steps raise the temperature by 5 K, ending at 300 K
simulation.context.setVelocitiesToTemperature(initial_temp*kelvin)
print('Warming up the system...')
simulation.reporters.append(PDBReporter((qht_trajectory_filename), reportInterval = qht_log_int))

for i in range(60):
  simulation.step(int(qht_run/60) )
  temperature = (initial_temp+(i*initial_temp))*kelvin 
  integrator.setTemperature(temperature)

simulation.saveState((qht_state_filename))
simulation.saveCheckpoint((qht_restart_filename))

# NPT equilibration, reducing backbone constraints
barostat = system.addForce(MonteCarloBarostat(pressure, 300*kelvin, 1))
simulation.context.reinitialize(True)
print('Running NPT equilibration...')
simulation.context.setParameter('k', 10*kilocalories_per_mole/angstroms**2)
simulation.step(peq_run)

# save the equilibration results to file : state is platform independent but less precise, checkpoint file
simulation.saveState((peq_state_filename))
simulation.saveCheckpoint((peq_restart_filename))

# Load checkpoint
# reset step and time counters
simulation.context.setParameter('k', 0)
simulation.loadCheckpoint((peq_restart_filename))
eq_state = simulation.context.getState(getVelocities=True, getPositions=True)
positions = eq_state.getPositions()
velocities = eq_state.getVelocities()
integrator = LangevinIntegrator(300*kelvin, friction_coefficient, step_size)
integrator.setConstraintTolerance(0.00001)
simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

# append reporters
simulation.reporters.append(
    DCDReporter((prod_trajectory_filename + '.dcd'), peq_trj_int))
simulation.reporters.append(
    StateDataReporter(prod_data_filename,
                      prod_log_int,
                      step=True,
                      potentialEnergy=True,
                      temperature=True,
                      progress=True,
                      remainingTime=True,
                      speed=True,
                      totalSteps=prod_run,
                      separator=','))
simulation.reporters.append(
    PDBReporter((prod_trajectory_filename + '.pdb'), reportInterval = peq_trj_int))

# run production simulation fpr 10ns
print('Running Production...')
simulation.step(prod_run)
simulation.saveState(prod_state_filename)
simulation.saveCheckpoint(prod_restart_filename)
print('Simulation Finished!, Buy Developer A Beer!')
