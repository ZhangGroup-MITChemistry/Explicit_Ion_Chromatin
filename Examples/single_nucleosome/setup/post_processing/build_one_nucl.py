#!/usr/bin/env python

from merge_mols import *

prot = lmp.Data()
prot.read_from_file('../sbm_scripts/single_nucl/data.histone', peptideFlag=1)

dna = lmp.Data()
dna.read_from_file('./dnaCoord_from_xtal/dna_from_xtal.in', peptideFlag=1)

dna.atom_types = dna.atom_types[0:14]

res_mol1 = 0
for atoms in dna.atoms:
    if atoms['mol_i'] == 1:
        res_mol1 = atoms['res_i']

for atoms in dna.atoms:
    if atoms['mol_i'] == 2:
        atoms['res_i'] += res_mol1 

#update_dna_topology(len(prot.atoms), len(prot.atoms)+len(dna.atoms))

print 'start merging data files ...'
topo = merge(prot, dna)
topo.write_to_file('data.prot_dna', peptideFlag=1)
topo.write_to_file('data.prot_dna.4vmd', peptideFlag=0)
