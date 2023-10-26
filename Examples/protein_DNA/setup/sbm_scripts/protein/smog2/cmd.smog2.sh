#########################################################################
# Author: Xingcheng Lin
# Created Time: Thu Feb 28 16:24:06 2019
# File Name: cmd.smog2.sh
# Description: 
#########################################################################
#!/bin/bash

export smog2dir=$HOME/bin/smog-2.2
export PDBname=$1


smog2 -i $PDBname -t $smog2dir/share/templates/SBM_AA -tCG $smog2dir/share/templates/SBM_calpha+gaussian -warn 1
