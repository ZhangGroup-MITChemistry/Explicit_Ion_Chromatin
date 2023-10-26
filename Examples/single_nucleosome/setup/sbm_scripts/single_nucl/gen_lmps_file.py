#!/usr/bin/env python

import sys
sys.path.append('../')
from convert_gromacs_to_lammps import *
from numpy import *

lmps_topo = build_sbm_4_lmps('smog.top', 'smog.gro', energy_scale = 2.5)

chainLen = array([135, 102, 128, 122, 135, 102, 128, 122])
nc = shape(chainLen)[0]
chainSta = ones(nc, dtype=int)
chainSta[1::] = cumsum(chainLen)[0:-1] + 1

# fix molecule ID
for atom in lmps_topo.atoms:
    for ic in range(nc-1):
        if (atom['res_i'] >= chainSta[ic]) and (atom['res_i'] < chainSta[ic+1]): 
            atom['mol_i'] = ic+1
            continue
    if (atom['res_i'] >= chainSta[nc-1]):
        atom['mol_i'] = nc

lmps_topo.write_to_file('data.histone', peptideFlag=1)
lmps_topo.write_to_file('data.histone.4vmd', peptideFlag=0)
