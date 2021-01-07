#!/bin/bash

#PBS -N MC0_LEE
#PBS -l nodes=1:ppn=72
#PBS -q normal
#PBS -M hekun_lee@sjtu.edu.cn
#PBS -m abe
#PBS -o run_diff.out
#PBS -e run_diff.err

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > host.mpd
np=`cat $PBS_NODEFILE | wc -l`
echo "The number of core is" $np
python correlation_mcmc.py 0 1234 10000 $np
rm host.mpd
