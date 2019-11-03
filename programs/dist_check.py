# -*- coding: utf-8 -*-
# This script check the minimum distance of all peptide pairs in the system.
# If any pair is too close (5 angstrom) exit with error code.

import itertools
import sys
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import pandas as pd


TOP_FILENAME = 'initial_wat_ion.pdb'
TRAJ_FILENAME = 'initial_wat_ion.pdb'


universe = mda.Universe(TOP_FILENAME, TRAJ_FILENAME)


for i, j in itertools.combinations( range(4), 2 ):
    # ref: reference, conf: configure
    segment_ref = universe.select_atoms(f'segid {i}')
    segment_conf = universe.select_atoms(f'segid {j}')
    minimum_dist = distance_array(segment_ref.atoms.positions, segment_conf.atoms.positions).min()

    if (minimum_dist < 5):
        sys.exit(1)
