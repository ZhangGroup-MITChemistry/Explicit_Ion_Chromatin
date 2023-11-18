#!/usr/bin/env python

# units in gromacs are in kilo jourles
# 
# sigma = 4 \AA; and epsilon = 1.0 for repulsive potential

import sys
sys.path.append('/Path_to_lmp_tools/lib/')
import lammps_tools as lmp

import numpy as np
import time

from gromacs_parser import *

kj2kcal = 0.239006

# http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
amino_acid_list = ['ALA','ARG', 'ASN',  'ASP',  
                   'CYS','GLU', 'GLN','GLY',
                   'HIS','ILE', 'LEU', 'LYS',    
                   'MET','PHE', 'PRO','SER',
                   'THR','TRP','TYR','VAL']

amino_acid_mass = [71.0788, 156.1875,114.1038,115.0886,
                  103.1388, 129.1155,128.1307,57.0519,
                  137.1411, 113.1594, 113.1594, 128.1741,
                  131.1926, 147.1766, 97.1167, 87.0782,
                  101.1051, 186.2132, 163.1760, 99.1326]

#   --------------------------------------------------  #
#   function starts here

def read_top(top_filename):
    gp = GromacsParser()
    gp.process_file(top_filename)

    return gp

def read_gro(gro_filename):
    '''
    read xyz from gro file
    '''

    fo = open(gro_filename, 'r')
    rawdata = fo.readlines()
    fo.close()

    n_atom = int(rawdata[1].strip())
    xyz = np.zeros((n_atom, 3), dtype=float)

    for line in rawdata[2:-1]:
        items = line.split()

        ia = int(items[2]) - 1

        xyz[ia,0] = float(items[3])*10
        xyz[ia,1] = float(items[4])*10
        xyz[ia,2] = float(items[5])*10

    return n_atom, xyz

def add_atom_types(topo):
    '''
    add the 20 amino acids that have different mass
    '''

    n_atom_types = len(amino_acid_mass)
    for iat in range(n_atom_types):
        topo.add_atom_type(amino_acid_mass[iat], None, None)

def add_bond(gmcs_topo, lmps_topo, energy_scale):
    '''
    all the bond between alpha-carbons will be modeled using list/ca

    bond length is nm in gromacs. need to change to A
    Also need to convert the K unit from kj/nm^2 to kcal/A^2
    '''

    comment='list/ca'
    coeffs = [10.0, 30.0]
    lmps_topo.add_bond_type(coeffs, comment)

    # write a bond list file
    fo = open('ca_bond_list.txt', 'w')

    for bond in gmcs_topo.current_molecule_type.bonds:
        atoml = (int(bond[0]), int(bond[1]))
        lmps_topo.add_bond(atoml, comment=None, atom_names=None, i=None, bond_type=lmps_topo.bond_types[0])

        # ;ai     aj      r0(nm)  Kb
        # gromacs has 1/2 in the code and lammps does not
        fo.write('%s\t%s\t%15.10e\t%15.10e\n'%(bond[0], bond[1], float(bond[3])*10, float(bond[4])*kj2kcal/100/2.0 * energy_scale))

    fo.close()

def add_angle(gmcs_topo, lmps_topo, energy_scale):
    '''
    all the angle between alpha-carbons will be modeled using list/ca
    '''

    # type 1 is for DNA
    comment='list/ca'
    coeffs = [10.0, 30.0]
    lmps_topo.add_angle_type(coeffs, comment)

    # write a angle list file
    fo = open('ca_angle_list.txt', 'w')

    for angle in gmcs_topo.current_molecule_type.angles:
        atoml = (int(angle[0]), int(angle[1]), int(angle[2]))
        lmps_topo.add_angle(atoml, comment=None, atom_names=None, i=None, angle_type=lmps_topo.angle_types[0])

        # ;ai     aj      r0(nm)  Kb
        # gromacs has 1/2 in the code and lammps does not
        fo.write('%s\t%s\t%s\t%s\t%15.10e\n'%(angle[0], angle[1], angle[2], angle[4], float(angle[5])*kj2kcal/2.0 * energy_scale))

    fo.close()

def add_dihedral(gmcs_topo, lmps_topo, energy_scale):
    '''
    all the dihedrals between alpha-carbons will be modeled using list/ca
    '''

    # type 1 is for DNA
    comment='list/ca'
    coeffs = [10.0, 30.0]
    lmps_topo.add_dihedral_type(coeffs, comment)

    # write a dihed list file
    fo = open('ca_dihed_list.txt', 'w')

    for idl, dihedral in enumerate(gmcs_topo.current_molecule_type.dihedrals):
        if np.mod(idl,2) == 0:
            # the two are redundant, only use one of them
            atoml = (int(dihedral[0]), int(dihedral[1]), int(dihedral[2]), int(dihedral[3]))
            lmps_topo.add_dihedral(atoml, comment=None, atom_names=None, i=None, dihedral_type=lmps_topo.dihedral_types[0])

            # ;ai     aj      r0(nm)  Kb
            fo.write('%s\t%s\t%s\t%s\t%s\t%15.10e\n'%(dihedral[0], dihedral[1], dihedral[2], dihedral[3], dihedral[5],float(dihedral[6])*kj2kcal * energy_scale))

    fo.close()

def write_pair_list(gmcs_topo, energy_scale):

    # write a angle list file
    fo = open('ca_pair_list.txt', 'w')
    fo2 = open('exclusion_list.txt', 'w')

    for pair in gmcs_topo.current_molecule_type.pairs:

        #; i j ftype A mu sigma a
        fo.write('%s\t%s\t%s\t%15.10e\t%15.10e\t%15.10e\t%15.10e\t%15.10e\n'%(pair[0], pair[1], 'ca', float(pair[3])*kj2kcal * energy_scale, float(pair[4])*10, float(pair[5])*10, float(pair[6])*1e12, float(pair[4])*10+6*float(pair[5])*10) )
        fo2.write('%s\t%s\n'%(pair[0], pair[1]))

    fo.close()
    fo2.close()


def build_sbm_4_lmps(top_filename, gro_filename, energy_scale=None):

    if not energy_scale:
        energy_scale = 1.0

    amino_acid_charge = {"GLY" : 0.0, "ALA" : 0.0, "LEU" : 0.0, "ILE" : 0.0,
                 "ARG" : 1.0, "LYS" : 1.0, "MET" : 0.0, "CYS" : 0.0,
                 "TYR" : 0.0, "THR" : 0.0, "PRO" : 0.0, "SER" : 0.0,
                 "TRP" : 0.0, "ASP" :-1.0, "GLU" :-1.0, "ASN" : 0.0,
                 "GLN" : 0.0, "PHE" : 0.0, "HIS" : 0.0, "VAL" : 0.0,
                 "M3L" : 0.0, "MSE" : 0.0, "CAS" : 0.0, "HSD" : 0.0 }

    lmps_topo = lmp.Data()
    lmps_topo.header = 'protein dna data file produced on: %s'%time.asctime()
    lmps_topo.box[0] = (-1000, 1000)
    lmps_topo.box[1] = (-1000, 1000)
    lmps_topo.box[2] = (-1000, 1000)

    # format: ;nr  type  resnr residue atom  cgnr
    residue_index = 3
    gmcs_topo   = read_top(top_filename)
    n_atom, xyz = read_gro(gro_filename)

    # self.add_atom(x, y, z, mol_i, charge, comment=comment, i=i, atom_type=self.atom_types[atom_type_i - 1], res_i=residue_i)
    #self.current_molecule_type.atoms.append(fields)

    resname_index = add_atom_types(lmps_topo)
    molecule_id = 1

    for ia in range(n_atom):
        resname  = gmcs_topo.current_molecule_type.atoms[ia][residue_index]
        charge = amino_acid_charge[resname]
        lmps_topo.add_atom(xyz[ia,0], xyz[ia,1], xyz[ia,2], molecule_id, charge, i=ia+1, atom_type=lmps_topo.atom_types[amino_acid_list.index(resname)], res_i=ia+1)

    # build topology
    add_bond(gmcs_topo, lmps_topo, energy_scale)

    add_angle(gmcs_topo, lmps_topo, energy_scale)

    add_dihedral(gmcs_topo, lmps_topo, energy_scale)

    write_pair_list(gmcs_topo, energy_scale)

    return lmps_topo

if __name__ == '__main__':

    build_sbm_4_lmps('smog.top', 'smog2.gro')

