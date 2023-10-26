export startingRef1=50.0
export pyPath=$HOME/bin/anaconda2/bin

startID1=$1
endID1=$2

for ((i=$startID1; i<=$endID1; i++))
do
   ref1=$(bc<<<"$startingRef1+50.0*($i)" | awk '{printf "%1.1f\n", $0}')
   echo $ref1
   mkdir -p equil/umbreref_${ref1}
   cp starting_structure/start_job_temp300.dcd equil/umbreref_${ref1}/start_job_temp300.dcd
   sed "s/RefValue1/$ref1/g" template_restart1d.py > restart1d.py
   sed "s/RefValue1/$ref1/g" template_plumed_ini.dat > plumed_ini.dat
   sed "s/RefValue1/$ref1/g" template_plumed_again.dat > plumed_again.dat
   $pyPath/python restart1d.py 0
done

