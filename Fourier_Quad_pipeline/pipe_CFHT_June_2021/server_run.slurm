#!/bin/bash

#SBATCH --job-name=CFHT
#SBATCH --partition=cpu
#SBATCH -n 400
#SBATCH --exclusive
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=betajzhang@sjtu.edu.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err
##SBATCH --time=00:20:10

dir=/lustre/home/acct-phyzj/phyzj

srun --mpi=pmi2 -n 400 $dir/jzhang/pipe_CFHT_June_2021/main $dir/CFHT/i/source.list $dir/jzhang/pipe_CFHT_June_2021

