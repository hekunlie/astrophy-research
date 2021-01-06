#!/bin/bash

#PBS -N MCMC_LEE
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
python ./test.py $np
rm host.mpd
