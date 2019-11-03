# coding: utf-8
import string


PATH_ORIG = 'initial_wat_ion.pdb'


with open(PATH_ORIG) as f:
    lines = [l.strip() for l in f.readlines()]


newlines = []
pep_label = 0
for line in lines:
    newline = line

    if line[:4] == 'ATOM' and pep_label < 4:
        newline = newline[:21] + str(pep_label) + newline[22:]
    elif line[:3] == 'TER':
        pep_label += 1

    newlines.append(newline)
    print(newline)
        

