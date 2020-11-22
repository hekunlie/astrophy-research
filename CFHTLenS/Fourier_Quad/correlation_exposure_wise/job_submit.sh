#!/bin/bash

#PBS -N CORR_LEE
#PBS -l nodes=5:ppn=72
#PBS -q normal
#PBS -M hekun_lee@sjtu.edu.cn
#PBS -m abe
#PBS -o run.out
#PBS -e run.err
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > host.mpd
np=`cat $PBS_NODEFILE | wc -l`
echo "The number of core is" $np
mpirun -machinefile host.mpd -np $np /home/hklee/work/CFHT/correlation/code/tomo_cal/cal_cor /home/hklee/work/CFHT/correlation
rm host.mpd
