#!/usr/bin/env python

import sys
sys.path.append('/Users/xl23/bin/lmp_tools/lib/')
import lammps_tools as lmp
import fileinput
import numpy as np
import time
from copy import *

def merge(prot, dna):
    '''
    merge the two data files
    '''

    topo = lmp.Data()
    topo.header = 'protein dna data file produced on: %s'%time.asctime()
    topo.box = prot.box

    n_prot_atoms = len(prot.atoms)
    n_prot_mols  = len(prot.mols)
    n_prot_res   = prot.atoms[-1]['res_i']
    n_dna_atom_types = len(dna.atom_types)

    # merge mols
    for mol in prot.mols:
        topo.mols.append(mol)

    for mol in dna.mols:
        mol_new = copy(mol)
        mol_new['i'] += n_prot_mols

        topo.mols.append(mol_new)

    # merge atom types; dna atom types must be listed first
    for atom_type in dna.atom_types:
        topo.atom_types.append(atom_type)

    for atom_type in prot.atom_types:
        atom_type_new = copy(atom_type)
        atom_type_new['i'] += n_dna_atom_types

        topo.atom_types.append(atom_type_new)

    # protein atoms must be listed first, but their types must be shifted
    for atom in prot.atoms:
        atom_new = copy(atom)
        atom_new['atom_type_i'] = topo.atom_types[atom['atom_type_i']+n_dna_atom_types-1]['i']
        #atom_new['atom_type'] = topo.atom_types[atom['atom_type_i']]

        topo.atoms.append(atom_new)

    for atom in dna.atoms:
        atom_new = copy(atom)
        atom_new['i'] += n_prot_atoms
        atom_new['mol_i'] += n_prot_mols
        atom_new['res_i'] += n_prot_res
        atom_new['charge'] /= 0.6
        topo.atoms.append(atom_new)

    #   ----    bonds   ----    #
    #   list DNA bonds first
    n_dna_bonds = len(dna.bonds)
    n_dna_bond_types = len(dna.bond_types)

    for bond_type in dna.bond_types:
        topo.bond_types.append(bond_type)

    for bond_type in prot.bond_types:
        bond_type_new = copy(bond_type)
        bond_type_new['i'] += n_dna_bond_types

        topo.bond_types.append(bond_type_new)

    for bond in dna.bonds:
        topo.add_bond((bond['atom1']['i'] + n_prot_atoms, bond['atom2']['i'] + n_prot_atoms) , comment=None, atom_names=None, i=bond['i'], bond_type=topo.bond_types[bond['bond_type_i']-1])

    for bond in prot.bonds:
        topo.add_bond((bond['atom1']['i'], bond['atom2']['i']) , comment=None, atom_names=None, i=bond['i']+n_dna_bonds, bond_type=topo.bond_types[bond['bond_type_i']+n_dna_bond_types-1])

    #   ----    angles   ----    #
    #   only DNA has angles 
    n_dna_angles = len(dna.angles)
    n_dna_angle_types = len(dna.angle_types)

    for angle_type in dna.angle_types:
        topo.angle_types.append(angle_type)

    for angle_type in prot.angle_types:
       angle_type_new = copy(angle_type)
       angle_type_new['i'] += n_dna_angle_types

       topo.angle_types.append(angle_type_new)

    for angle in dna.angles:
        topo.add_angle((angle['atom1']['i'] + n_prot_atoms,angle['atom2']['i'] + n_prot_atoms,angle['atom3']['i'] + n_prot_atoms), comment=None, atom_names=None, i=None, angle_type=topo.angle_types[angle['angle_type_i']-1])

    for angle in prot.angles:
        topo.add_angle((angle['atom1']['i'], angle['atom2']['i'], angle['atom3']['i']), comment=None, atom_names=None, i=None, angle_type=topo.angle_types[angle['angle_type_i']+n_dna_angle_types-1])

    #   ----    dihedrals   ----    #
    #   only DNA angles first
    n_dna_dihedals = len(dna.dihedrals)
    n_dna_dihedral_types = len(dna.dihedral_types)

    for dihedral_type in dna.dihedral_types:
        topo.dihedral_types.append(dihedral_type)

    for dihedral_type in prot.dihedral_types:
       dihedral_type_new = copy(dihedral_type)
       dihedral_type_new['i'] += n_dna_dihedral_types

       topo.dihedral_types.append(dihedral_type_new)

    for dihedral in dna.dihedrals:
        topo.add_dihedral((dihedral['atom1']['i'] + n_prot_atoms,dihedral['atom2']['i'] + n_prot_atoms,dihedral['atom3']['i'] + n_prot_atoms,dihedral['atom4']['i'] + n_prot_atoms ), comment=None, atom_names=None, i=None, dihedral_type=topo.dihedral_types[dihedral['dihedral_type_i']-1])

    for dihedral in prot.dihedrals:
        topo.add_dihedral((dihedral['atom1']['i'], dihedral['atom2']['i'], dihedral['atom3']['i'], dihedral['atom4']['i']), comment=None, atom_names=None, i=None, dihedral_type=topo.dihedral_types[dihedral['dihedral_type_i']+n_dna_dihedral_types-1])


    return topo

def replicate(topo, shiftz=True):
    '''
    shift the nucleosome in the z direction by 100 angstrom
    '''
    topo_2nucl = deepcopy(topo)

    n_atom_per_nucl = len(topo.atoms)
    n_mol_per_nucl  = len(topo.mols)
    n_res_per_nucl  = topo.atoms[-1]['res_i']
    n_bond_per_nucl = len(topo.bonds)
    n_angle_per_nucl = len(topo.angles)
    n_dihedral_per_nucl = len(topo.dihedrals)

    # add mols
    for mol in topo.mols:
        mol_new = copy(mol)
        mol_new['i'] += n_mol_per_nucl
        topo_2nucl.mols.append(mol_new)

    # add atoms
    for atom in topo.atoms:
        atom_new = copy(atom)
        
        if shiftz: 
            atom_new['z'] += 100.0
        else:
            atom_new['x'] += 100.0
            atom_new['y'] += 100.0

        atom_new['i'] += n_atom_per_nucl
        atom_new['mol_i'] += n_mol_per_nucl
        atom_new['res_i'] += n_res_per_nucl

        topo_2nucl.atoms.append(atom_new)

    # add bonds
    for bond in topo.bonds:
        bond_new = copy(bond)
        bond_new['atom1'] = topo_2nucl.atoms[bond['atom1']['i']+n_atom_per_nucl-1]
        bond_new['atom2'] = topo_2nucl.atoms[bond['atom2']['i']+n_atom_per_nucl-1]

        bond_new['i'] += n_bond_per_nucl

        topo_2nucl.bonds.append(bond_new)

    # add angles
    for angle in topo.angles:
        angle_new = copy(angle)
        angle_new['atom1'] = topo_2nucl.atoms[angle['atom1']['i']+n_atom_per_nucl-1]
        angle_new['atom2'] = topo_2nucl.atoms[angle['atom2']['i']+n_atom_per_nucl-1]
        angle_new['atom3'] = topo_2nucl.atoms[angle['atom3']['i']+n_atom_per_nucl-1]

        angle_new['i'] += n_angle_per_nucl

        topo_2nucl.angles.append(angle_new)

    # add dihedrals
    for dihedral in topo.dihedrals:
        dihedral_new = copy(dihedral)
        dihedral_new['atom1'] = topo_2nucl.atoms[dihedral['atom1']['i']+n_atom_per_nucl-1]
        dihedral_new['atom2'] = topo_2nucl.atoms[dihedral['atom2']['i']+n_atom_per_nucl-1]
        dihedral_new['atom3'] = topo_2nucl.atoms[dihedral['atom3']['i']+n_atom_per_nucl-1]
        dihedral_new['atom4'] = topo_2nucl.atoms[dihedral['atom4']['i']+n_atom_per_nucl-1]

        dihedral_new['i'] += n_dihedral_per_nucl

        topo_2nucl.dihedrals.append(dihedral_new)

    return topo_2nucl

if __name__ == '__main__':

    prot = lmp.Data()
    prot.read_from_file('../buildHistone/data.histoneCore', peptideFlag=1)

    dna = lmp.Data()
    dna.read_from_file('./dnaCoord_from_xtal/dna_from_xtal.in', peptideFlag=1)

    update_dna_topology(len(prot.atoms), len(prot.atoms)+len(dna.atoms))

    #topo = merge(prot, dna)
    #topo_2nucl = replicate(topo)

    #topo.write_to_file('data.prot_dna.one_nucl', peptideFlag=1)
    #topo_2nucl.write_to_file('data.prot_dna.forvmd', peptideFlag=0)
    #topo_2nucl.write_to_file('data.prot_dna', peptideFlag=1)

