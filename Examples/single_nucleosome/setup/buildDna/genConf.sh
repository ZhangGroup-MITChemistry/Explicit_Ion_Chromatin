#########################################################################
# Author: Xingcheng Lin
# Created Time: Thu 03 Apr 2014 05:47:21 PM CDT
# File Name: genConf.sh
# Description: 
#########################################################################
#!/bin/bash

spn2cdir="Path_to_GitHub/LAMMPS-3SPN2/"

pyPath="Path_to_anaconda2/bin"

$pyPath/python ${spn2cdir}/utils/make_bp_params.py dnaSeq.txt 
x3dna_utils cp_std BDNA
rebuild -atomic bp_step.par atomistic.pdb

$pyPath/python ${spn2cdir}/utils/pdb2cg_dna.py atomistic.pdb 
${spn2cdir}/DSIM_ICNF/icnf.exe dnaSeq.txt 1 1 . 0 
${spn2cdir}/utils/replace_atoms.sh conf_lammps.in dna_conf.in bdna_curv_conf_4list.in

$pyPath/python ${spn2cdir}/utils/make_list_files.py bdna_curv_conf_4list.in 
