#!/bin/bash

#PBS -N TEST_LEE
#PBS -l nodes=gr01:ppn=72 amd:nodes=gr02:ppn=72
#PBS -q normal
#PBS -M hekun_lee@sjtu.edu.cn
#PBS -m abe
#PBS -o run.out
#PBS -e run.err
cd $PBS_O_WORKDIR
cat $PBS_O_WORKDIR
echo "------------This is the node file -------------"
cat $PBS_NODEFILE
echo "-----------------------------------------------"
cat $PBS_NODEFILE > host.mpd
np=`cat $PBS_NODEFILE | wc -l`
echo "The number of core is" $np
echo
echo
mpirun -machinefile host.mpd -np $np /home/hklee/work/code_test/mpi_test/test_on_cluster
