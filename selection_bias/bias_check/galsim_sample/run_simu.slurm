#!/bin/bash

#SBATCH --job-name=hkl_egalsimu
#SBATCH --partition=cpu
#SBATCH --mail-user=hekun_lee@sjtu.edu.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -n 400
#SBATCH --ntasks-per-node=40

ulimit -s unlimited
ulimit -l unlimited

source /usr/share/Modules/init/bash
code_path=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/selection_bias/bias_check/galsim_sample/change_flux/epsf/imgs_5/code
data_path=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/selection_bias/bias_check/galsim_sample/change_flux/epsf/imgs_5

srun --mpi=pmi2 python $code_path/simulation.py $data_path 1


