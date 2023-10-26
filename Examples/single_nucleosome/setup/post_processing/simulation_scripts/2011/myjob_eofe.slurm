#!/bin/bash -l
##

#SBATCH --output=slurm-%j.out
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --partition=sched_mit_binz,newnodes,sched_mit_hill
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH --exclude=node[762-930,932,934]
#SBATCH --constraint=centos7

module purge
module load gcc
module load engaging/openmpi/1.8.8
export pyPath=$HOME/bin/anaconda2/bin

lammpsdir="/home/xclin/bin/awsem/bin"

curr=$SLURM_ARRAY_TASK_ID
prev=$(( $curr - 1 ))
curr=$(printf "%'02d" ${curr})
prev=$(printf "%'02d" ${prev})

cd run${curr}

if [ $curr -eq "00" ];
then
        cp ../../../start_job_temp*.dcd ./
	sed "s/RANDOMSEED/$RANDOM/g" ../../../template_template.lammps > ../../../template.lammps
	$pyPath/python  ../../../get_lammps_input.py ../run${prev} ./ ${curr}
else
        cp ../run${prev}/*.cpt ./        
	sed "s/RANDOMSEED/$RANDOM/g" ../../../template_template.lammps > ../../../template.lammps
        $pyPath/python  ../../../get_lammps_input.py ../run${prev} ./ ${curr}
fi

mpirun -np 12 $lammpsdir/lmp_openmpi_plumed246_3spn_2016_sbm_rigidtemper_tanh_expions -partition 12x1 -in proteinDna_sim.in

cd ../

