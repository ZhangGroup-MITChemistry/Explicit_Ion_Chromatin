export startingRef=50.0
export pyPath=$HOME/bin/anaconda2/bin

startID=$1
endID=$2

for ((i=0; i<=11; i++));
do
	ref=$(bc<<<"$startingRef+50.0*($i)" | awk '{printf "%1.1f\n", $0}')
	echo $ref;
	sed "s/RefValue/$ref/g" template_restart1d.py > restart1d.py
        sed "s/RefValue/$ref/g" template_input.colvars_ini > input.colvars_ini.dat
        sed "s/RefValue/$ref/g" template_input.colvars_again > input.colvars_again.dat
	$pyPath/python restart1d.py $startID $endID
done

