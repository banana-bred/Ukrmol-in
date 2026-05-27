#!/bin/bash -l

export exec=NICE-CI_ARCHER.x

export pwd=`pwd`

export name=$1

echo $pwd

echo $name

export nproc=$2



export wclim=$3


echo "Nproc=" $nproc

export num_nodes=$(((nproc-1)/24 + 1))

echo $num_nodes



echo "Nnodes=" $num_nodes, "Nproc=" $nproc, " Memory = "$MEM, "jobtype = " $jobtype, "wclimit = " $wclim
echo "Working dir is " $pwd



qsub -r n -N $name -j oe -e $name.e -A ecse0807 -l "select=$num_nodes,walltime=$wclim:00:00" \
     -v "name=$name,wd=$pwd,nproc=$nproc,exec=$exec" \
     $pwd/run_mpici.csh

