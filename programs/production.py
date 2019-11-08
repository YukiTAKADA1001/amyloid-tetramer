# -*- coding: utf-8 -*-
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


# Input
PRMTOP_FILENAME = 'initial_wat_ion.prmtop'
INPCRD_FILENAME = 'initial_wat_ion.inpcrd'
# Simulation parameters
PRODUCTION_STEP = int(250e6)     # 500ns
RECORD_STEP = int(50e3)          # every 100ps
# Output
TRAJ_FILENAME = 'traj.dcd'
STATE_FILENAME = 'state.csv'
RESTART_FILENAME = 'restart.xml'


# System settings
prmtop = AmberPrmtopFile(PRMTOP_FILENAME)
inpcrd = AmberInpcrdFile(INPCRD_FILENAME)
system_params = {
        'nonbondedMethod': PME,
        'nonbondedCutoff': 12*angstroms,
        'constraints': HBonds
}
system = prmtop.createSystem(**system_params)

# Set random seeds const to reproduce simulation later.
# Pressure coupling
barostat = MonteCarloBarostat(1*bar, 310*kelvin)
barostat.setRandomNumberSeed(seed=42)
system.addForce(barostat)
# Temperature coupling
integrator = LangevinIntegrator(310*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setRandomNumberSeed(seed=42)

simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)


# Energy minimization
simulation.minimizeEnergy()


# Trajectory file setting
simulation.reporters.append(DCDReporter(TRAJ_FILENAME, RECORD_STEP))
# State data file setting
reporter_params = {
        'step': True,
        'time': True,
        'progress': True,
        'remainingTime': True,
        'potentialEnergy': True,
        'kineticEnergy': True,
        'totalEnergy': True,
        'temperature': True,
        'volume': True,
        'density': True,
        'speed': True,
        'totalSteps': PRODUCTION_STEP,
        'separator': ','
}
simulation.reporters.append( StateDataReporter(STATE_FILENAME, RECORD_STEP, **reporter_params) )


# Heating and equilibration of solvent 
simulation.context.setVelocitiesToTemperature(310*kelvin, randomSeed=42)


# Production run
# Solvent equilibration step is contained in production run steps.
# These steps are not used in analysis (after aggregation will be used), so that's not the issue,'
simulation.step(PRODUCTION_STEP)


# Create restart file
simulation.saveState(RESTART_FILENAME)
