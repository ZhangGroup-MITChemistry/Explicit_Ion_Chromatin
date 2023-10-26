#!/usr/bin/env python

import sys
sys.path.append('/Users/xl23/bin/lmp_tools/lib/')
import lammps_tools as lmp

def update_dna_topology(n_prot_atoms, n_atom_per_dna):
    '''
    modify atom index in the list files
    '''

    # update bond
    bfile = '../../buildDna/in00_bond.list'
    fin = open(bfile, 'r')
    rawdata = fin.readlines(); fin.close()

    fo = open('in00_bond.list', 'w')
    for line in rawdata:
        items = line.split()
        atom1 = int(items[0]) + n_prot_atoms
        atom2 = int(items[1]) + n_prot_atoms
        fo.write("%ld\t%ld\t%lf\t%lf\t%lf\t%lf\n" % (atom1,atom2,float(items[2]),float(items[3]),0.0,float(items[5])))

    fo.close()

    # update angle
    bfile = '../../buildDna/in00_angl.list'
    fin = open(bfile, 'r')
    rawdata = fin.readlines(); fin.close()

    fo = open('in00_angl.list', 'w')
    for line in rawdata:
        items = line.split()
        atom1 = int(items[0]) + n_prot_atoms
        atom2 = int(items[1]) + n_prot_atoms
        atom3 = int(items[2]) + n_prot_atoms
        fo.write("%ld\t%ld\t%ld\t%lf\t%lf\n" % (atom1,atom2,atom3,float(items[3]),float(items[4])))

    fo.close()

    # update dihedral
    bfile = '../../buildDna/in00_dihe.list'
    fin = open(bfile, 'r')
    rawdata = fin.readlines(); fin.close()

    fo = open('in00_dihe.list', 'w')
    for line in rawdata:
        items = line.split()
        atom1 = int(items[0]) + n_prot_atoms
        atom2 = int(items[1]) + n_prot_atoms
        atom3 = int(items[2]) + n_prot_atoms
        atom4 = int(items[3]) + n_prot_atoms
        fo.write("%ld\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%d\n" % (atom1,atom2,atom3,atom4,float(items[4]), float(items[5]), float(items[6]), float(items[7]), int(items[8])))

    fo.close()


if __name__ == '__main__':
    prot = lmp.Data()
    prot.read_from_file('../../sbm_scripts/protein/data.histone', peptideFlag=1)
    
    dna = lmp.Data()
    dna.read_from_file('../dnaCoord_from_xtal/dna_from_xtal.in', peptideFlag=1)

    update_dna_topology(len(prot.atoms), len(dna.atoms))
