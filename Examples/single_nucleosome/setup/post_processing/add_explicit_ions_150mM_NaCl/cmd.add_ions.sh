#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Dec  1 20:27:23 2021
# File Name: cmd.add_ions.sh
# Description: Add ions based on the current data.prot_dna and 
# the desired ionic concentration/
#########################################################################
#!/bin/bash

num_NA=19649
num_MG=0
num_CL=19505
num_existing_atoms_orig=1854

python getSection.py data.prot_dna tmp.txt "Atoms" "Bonds"
# Remove the white space
gsed -i '/^$/d' tmp.txt

num_existing_atoms=$((num_existing_atoms_orig))
chainID_sta=10
resid_sta=1268


python place_ions.py $num_NA NA_ION.txt -299.0 -299.0 -299.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 35 1.0

# Update the number list
chainID_sta=$((chainID_sta+num_NA))
resid_sta=$((resid_sta+num_NA))
num_existing_atoms=$((num_existing_atoms+num_NA))


python place_ions.py $num_MG MG_ION.txt -295.0 -295.0 -295.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 36 2.0

# Update the number list
chainID_sta=$((chainID_sta+num_MG))
resid_sta=$((resid_sta+num_MG))
num_existing_atoms=$((num_existing_atoms+num_MG))

python place_ions.py $num_CL CL_ION.txt -292.0 -292.0 -292.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 37 -1.0

# Merge them together
cat tmp.txt NA_ION.txt MG_ION.txt CL_ION.txt > modified_ATOM.txt


python getSection.py data.prot_dna.4vmd tmp.txt "Atoms" "Bonds"
# Remove the white space
gsed -i '/^$/d' tmp.txt


num_existing_atoms=$((num_existing_atoms_orig))
chainID_sta=10
resid_sta=1268

python place_ions_4vmd.py $num_NA NA_ION_4vmd.txt -299.0 -299.0 -299.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 35 1.0

# Update the number list
chainID_sta=$((chainID_sta+num_NA))
resid_sta=$((resid_sta+num_NA))
num_existing_atoms=$((num_existing_atoms+num_NA))

python place_ions_4vmd.py $num_MG MG_ION_4vmd.txt -295.0 -295.0 -295.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 36 2.0

# Update the number list
chainID_sta=$((chainID_sta+num_MG))
resid_sta=$((resid_sta+num_MG))
num_existing_atoms=$((num_existing_atoms+num_MG))

python place_ions_4vmd.py $num_CL CL_ION_4vmd.txt -292.0 -292.0 -292.0 19.0 $chainID_sta $resid_sta $num_existing_atoms 37 -1.0

num_existing_atoms_updated=$((num_existing_atoms+num_CL))

# Merge them together
cat tmp.txt NA_ION_4vmd.txt MG_ION_4vmd.txt CL_ION_4vmd.txt > modified_ATOM_4vmd.txt

# Update the original data.prot_dna file to include the explicit ions
python updateData_Atoms.py data.prot_dna modified_ATOM.txt data.prot_dna_ions
python updateData_Atoms.py data.prot_dna.4vmd modified_ATOM_4vmd.txt data.prot_dna_ions.4vmd

# Replace the total number of atoms
gsed -i "s/$num_existing_atoms_orig atoms/$num_existing_atoms_updated atoms/g; s/34 atom types/37 atom types/g; s/-1000.000000 1000.000000 xlo xhi/-300.0 300.0 xlo xhi/g; s/-1000.000000 1000.000000 ylo yhi/-300.0 300.0 ylo yhi/g; s/-1000.000000 1000.000000 zlo zhi/-300.0 300.0 zlo zhi/g; s/34   99.1326/34   99.1326\n  35   22.9898\n  36   24.3050\n  37   35.4530/g" data.prot_dna_ions
gsed -i "s/$num_existing_atoms_orig atoms/$num_existing_atoms_updated atoms/g; s/34 atom types/37 atom types/g; s/-1000.000000 1000.000000 xlo xhi/-300.0 300.0 xlo xhi/g; s/-1000.000000 1000.000000 ylo yhi/-300.0 300.0 ylo yhi/g; s/-1000.000000 1000.000000 zlo zhi/-300.0 300.0 zlo zhi/g; s/34   99.1326/34   99.1326\n  35   22.9898\n  36   24.3050\n  37   35.4530/g" data.prot_dna_ions.4vmd

