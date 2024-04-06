# Explicit_Ion_Chromatin
Lammps Implementation of the Explicit-ion Residue-resolution Chromatin Model

## Update in 2024
We provide an example src folder of the Lammps 2016 version which include all the customized codes.

## Installation guidance
These codes have been tested to run correctly on the Lammps 2016 version.

Copy the codes in "Code/USER-SMOG" to the Lammps src/ folder, then compile for the executables.

## Example

* Change the path of the path of relevant environmental variables to the corresponding ones on your computer.

* Running an explicit-ion simulation with a protein-DNA complex structure (PDB code: 1APL).

* Step 1 Prepare the DNA sequence
* Go to the folder protein_DNA/setup/buildDna/ and use the command:
```
bash genConf.sh
```

* Step 2 Prepare the structure-based protein model based on the structure
* Go to the folder protein_DNA/setup/sbm_scripts/protein/smog2/ and use the command:
```
bash cmd.smog2.sh protein.pdb
cd ../
python gen_lmps_file.py
```

* Step 3 Prepare the protein-DNA complex input file for Lammps;
* Go to the folder protein_DNA/setup/post_processing and use the command:
```
bash build_one_nucl.sh
```

* Step 4 Prepare the explicit-ion input file
* Copy the data.prot_dna to the folder protein_DNA/setup/post_processing/add_explicit_ions_100mM_NaCl_200box and execute the command in the folder:
```
bash cmd.add_ions.sh 
```
The user needs to have a good estimate of the number of ions based on the targeted ionic concentration.


* Running an explicit-ion simulation with a nucleosome structure (PDB code: 1APL).

* The step is the same as the one used for running the protein-DNA complex structure, except for two additional steps
* Step 2.1 When preparing the protein input file, we need to delete the bonded interactions for the histone tails
* Execute this command in the folder single_nucleosome/setup/sbm_scripts/single_nucl/
```
bash cmd.remove_idr.sh
```

* Step 3.1 When preparing the Lammps input file for the nucleosome, we need to run this script in the folder single_nucleosome/setup/post_processing/groupsDefinition/ to generate the file for those domains applying the rigid-body approximation:
```
python prep_nucl_rigid.py
```
The user need to determine the number of rigid DNA base pairs by changing the parameters in the file prep_nucl_rigid.py

## Reference:
* Xingcheng Lin and Bin Zhang, Explicit Ion Modeling Predicts Physicochemical Interactions for Chromatin Organization, eLife (2023), https://doi.org/10.7554/eLife.90073
