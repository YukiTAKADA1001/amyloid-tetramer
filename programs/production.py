# -*- coding: utf-8 -*-
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


# Input
PRMTOP_FILENAME = 'initial_wat_ion.prmtop'
INPCRD_FILENAME = 'initial_wat_ion.inpcrd'
# Simulation parameters
EQUIL_STEP = 500e3          # 1ns
PRODUCTION_STEP = 500e4     # 10ns
RECORD_STEP = 50e3          # every 100ps
# Output
TRAJ_FILENAME = 'traj.dcd'
STATE_FILENAME = 'state.csv'
RESTART_FILENAME = 'restart.xml'


# System settings
prmtop = AmberPrmtopFile(PRMTOP_FILENAME)
inpcrd = AmberInpcrdFile(INPCRD_FILENAME)
system = prmtop.createSystem(nonbondedMethod=PME,
                             nonbondedCutoff=12*angstroms,
                             constraints=HBonds)
system.addForce(MonteCarloBarostat(1*bar, 310*kelvin))
integrator = LangevinIntegrator(310*kelvin, 1/picosecond, 2*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)


# Heating and equilibration of solvent 
simulation.context.setVelocitiesToTemperature(310*kelvin)
simulation.step(EQUIL_STEP)


# Production run
# Trajectory file setting
simulation.reporters.append(DCDReporter(TRAJ_FILENAME, RECORD_STEP))
# State data file setting
simulation.reporters.append(StateDataReporter(STATE_FILENAME,
                                              RECORD_STEP,
                                              step=True,
                                              time=True,
                                              progress=True,
                                              remainingTime=True,
                                              potentialEnergy=True,
                                              kineticEnergy=True,
                                              totalEnergy=True,
                                              temperature=True,
                                              speed=True,
                                              totalSteps=PRODUCTION_STEP,
                                              separator='▸----'))
simulation.step(PRODUCTION_STEP)


# Create restart file
simulation.saveState(RESTART_FILENAME)
