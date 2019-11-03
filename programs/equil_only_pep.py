# -*- encoding: utf-8 -*-
import os
from sys import stdout
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


## parameters
allstep = 2500000


def calcSinglePep(label):
    ## 系の設定
    prmtop = AmberPrmtopFile('straight.prmtop')
    inpcrd = AmberInpcrdFile('straight.inpcrd')
    system = prmtop.createSystem(implicitSolvent=GBn2,
                 nonbondedMethod=NoCutoff, constraints=HBonds)
    integrator = LangevinIntegrator(310*kelvin, 1/picosecond, 2*femtoseconds)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    
    ## エネルギー極小化
    simulation.minimizeEnergy()
    
    ## 平衡化MD計算(1ns)
    simulation.reporters.append(StateDataReporter('data.tsv', 5000, step=True,
            time=True, progress=True, remainingTime=True, potentialEnergy=True,
            kineticEnergy=True, totalEnergy=True, temperature=True, speed=True,
            totalSteps=allstep, separator='	'))
    simulation.step(allstep)

    ## Saving result as PDB file format
    positions = simulation.context.getState(getPositions=True).getPositions()
    filename = f'amyloid{label}.pdb'
    PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))
    

print(os.getcwd())
for label in range(4):
    calcSinglePep(label)
