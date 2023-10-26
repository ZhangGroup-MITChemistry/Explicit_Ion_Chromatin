#########################################################################
# Author: Charlse.Zhang
# Created Time: Thu 03 Apr 2014 05:47:21 PM CDT
# File Name: genConf.sh
# Description: 
#########################################################################
#!/bin/bash

spn2cdir="$HOME/Downloads/GitHub/LAMMPS-3SPN2/"
export X3DNA="$HOME/bin/x3dna-v2.4/"
export PATH="$HOME/bin/x3dna-v2.4/bin:$PATH"

pyPath="$HOME/bin/anaconda2/bin"

# modified the pdb2cg_dna.py from  ${spn2cdir}/utils to output residue id

PDBname=$1

$pyPath/python ../pdb2cg_dna.py $PDBname

${spn2cdir}/utils/replace_atoms.sh ../../buildDna/bdna_curv_conf_4list.in dna_conf.in dna_from_xtal_template.in

python replace_atom_type.py ../../buildDna/bdna_curv_conf_4list.in dna_from_xtal_template.in dna_from_xtal.in

