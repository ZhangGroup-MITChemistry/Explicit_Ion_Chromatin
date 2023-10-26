#!/usr/bin/env python

import sys
sys.path.append('/home/binz/Program/lmp_tools/lib/')
import lammps_tools as lmp

from numpy import *

dna = lmp.Data()
dna.read_from_file('dna_from_xtal.in', peptideFlag=1)

fo = open('dna_bond.psf', 'w')
fo.write('%8d !NBONDS\n'%len(dna.bonds))

ib = 0
for bond in dna.bonds:
    ib += 1 
    fo.write('%8d%8d'%(bond['atom1']['i'], bond['atom2']['i']))
    if mod(ib,4) == 0:
        fo.write('\n')

fo.close()
