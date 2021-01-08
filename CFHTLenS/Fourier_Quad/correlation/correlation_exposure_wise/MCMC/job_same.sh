#!/bin/bash

#PBS -N MC1_LEE
#PBS -l nodes=1:ppn=72
#PBS -q normal
#PBS -M hekun_lee@sjtu.edu.cn
#PBS -m abe
#PBS -o run_same.out
#PBS -e run_same.err

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > host.mpd
np=`cat $PBS_NODEFILE | wc -l`
echo "The number of core is" $np
python correlation_mcmc.py 1 4234 10000 $np
rm host.mpd
