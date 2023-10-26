#!/usr/bin/env python

from numpy import *
from numpy import linalg as LA
import sys
sys.path.append('/Users/xl23/bin/lmp_tools/lib/')
import lammps_tools as lmp

orderedSegment = arange(1, 184)

print orderedSegment

topo = lmp.Data()
topo.read_from_file('../data.prot_dna', peptideFlag=1)
n_atom = len(topo.atoms)

atoms = zeros((n_atom,4), dtype=int)
for ia, atom in enumerate(topo.atoms):
    atoms[ia,0] = atom['i']
    atoms[ia,1] = atom['mol_i']
    atoms[ia,2] = atom['res_i']
    atoms[ia,3] = atom['atom_type_i']

#   ----------------------------------------    #
#   selections
#   ----------------------------------------    #

#
# get alpha carbon coordinates

mask        = atoms[:, 1]<=3
prc2Ind    = squeeze( atoms[where(mask), 0] )

#   ----------------------------------------    #
#   output
#   ----------------------------------------    #

fh = open('group_rigid.txt', 'w')

for inucl, atomInd in enumerate([prc2Ind]):
    fh.write('group 1apl id ')
    natom = 0

    for id in range(len(atomInd)):
        # Only output those residues in the ordered segment as defined above;
        if topo.atoms[atomInd[id]-1]['res_i'] in orderedSegment:
            fh.write('%d '%atomInd[id])
            natom+=1
    fh.write('\n')
    print natom

# Output the core region;
fh.write('group core id ')
natom = 0
for inucl, atomInd in enumerate([prc2Ind]):

    for id in range(len(atomInd)):
        # Only output those residues in the ordered segment as defined above;
        if topo.atoms[atomInd[id]-1]['res_i'] in orderedSegment:
            fh.write('%d '%atomInd[id])
            natom+=1

fh.write('\n')
print natom

fh.close()
