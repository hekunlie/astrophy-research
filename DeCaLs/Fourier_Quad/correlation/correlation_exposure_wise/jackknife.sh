#!/bin/bash

#PBS -N JKNF_LEE
#PBS -l nodes=1:ppn=72
#PBS -q normal
#PBS -M hekun_lee@sjtu.edu.cn
#PBS -m abe
#PBS -o run.out
#PBS -e run.err
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > host.mpd
np=`cat $PBS_NODEFILE | wc -l`
echo "The number of core is" $np
mpirun -machinefile host.mpd -np $np /home/hklee/work/CFHT/correlation/code/jackknife/get_corr /home/hklee/work/CFHT/correlation 200 360
rm host.mpd
