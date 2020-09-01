#!/bin/bash

#SBATCH --job-name=hkl_correlation
#SBATCH --partition=cpu
#SBATCH --mail-user=hekun_lee@sjtu.edu.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -n 160
#SBATCH --ntasks-per-node=40

ulimit -s unlimited
ulimit -l unlimited
#source /usr/share/Modules/init/bash

src_path=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/CFHT/correlation

srun --mpi=pmi2 $src_path/code/cal_cor $src_path 509


