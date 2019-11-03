# coding: utf-8
import string


PDB_FILE = 'initial_wat_ion.pdb'


with open(PDB_FILE) as f:
    lines = [l.strip() for l in f.readlines()]


newlines=[]
pep_label = 0
for line in lines:
    newline = line

    if line[:4] == 'ATOM' and pep_label < 4:
        newline = newline[:21] + str(pep_label) + newline[22:]
    elif line[:3] == 'TER':
        pep_label += 1

    newlines.append(newline)


with open(PDB_FILE, mode='w') as f:
    f.write('\n'.join(newlines))
