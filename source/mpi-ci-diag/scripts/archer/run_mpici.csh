#!/bin/bash
echo $name
echo $exec
echo $nproc
echo $wd
#Move to the working directory
cd $wd

module swap PrgEnv-cray PrgEnv-intel

module swap cray-mpich/7.2.6 cray-mpich/7.4.3

module load cray-tpsl/16.07.1

module load cray-petsc/3.7.2.1

export OMP_NUM_THREADS=1

export LD_LIBRARY_PATH=/home/ecse0807/ecse0807/alref/usr/lib:$LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH

aprun -n $nproc $wd/$exec $wd/$name > $wd/$name.out

