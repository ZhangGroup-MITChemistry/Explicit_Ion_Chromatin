#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Apr 17 18:40:46 2019
# File Name: cmd.removeidr.sh
# Description: Remove the disordered region of PRC2 in smog.top
#########################################################################
#!/bin/bash

python getSection.py smog2/smog.top pairs.dat "[ pairs ]" "[ exclusions ]"
python removeIDR_for_pairs_section.py pairs.dat pairs_new.dat

python getSection.py smog2/smog.top exclusions.dat "[ exclusions ]" "[ system ]"
python removeIDR_for_exclusions_section.py exclusions.dat exclusions_new.dat

python getSection.py smog2/smog.top dihedrals.dat "[ dihedrals ]" "[ pairs ]"
python removeIDR_for_dihedrals_section.py dihedrals.dat dihedrals_new.dat

python updateTop_pairs.py smog2/smog.top pairs_new.dat smog_new.top
mv smog_new.top smog.top

python updateTop_exclusions.py smog.top exclusions_new.dat smog_new.top
mv smog_new.top smog.top

python updateTop_dihedrals.py smog.top dihedrals_new.dat smog_new.top
mv smog_new.top smog.top
