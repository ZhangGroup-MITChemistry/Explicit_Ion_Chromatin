#########################################################################
# Author: Xingcheng Lin
# Created Time: Thu 03 Apr 2014 05:47:21 PM CDT
# File Name: genConf.sh
# Description: 
#########################################################################
#!/bin/bash

spn2cdir="/Path_to_GitHub/USER-3SPN2/"
export X3DNA="/Path_to_x3dna-v2.4/"
export PATH="/Path_to_x3dna-v2.4/bin:$PATH"

# modified the pdb2cg_dna.py from  ${spn2cdir}/utils to output residue id
/Users/xl23/bin/anaconda2/bin/python ../pdb2cg_dna.py atomistic.pdb 

${spn2cdir}/utils/replace_atoms.sh ../../buildDna/bdna_curv_conf_4list.in dna_conf.in dna_from_xtal.in

