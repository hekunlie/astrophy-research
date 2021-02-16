#!/bin/bash

#PBS -N accu_LEE
#PBS -l nodes=1:ppn=72
#PBS -q normal
#PBS -M hekun_lee@sjtu.edu.cn
#PBS -m abe
#PBS -o run.out
#PBS -e run.err

module load compiler/gcc-6.5.0
#module load lee_cfitsio/cfitsio-3420
#module load cfitsio/3.48
#module load lee_hdf5/hdf5-1.10.4
#module load lee_fftw/fftw-3.3.8
module load lee_mpich/mpich-3.2.1
module load lee_python/python-3.8.6

src_path=/home/hklee/work/new_PDF/accu_test

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > host.mpd
np=`cat $PBS_NODEFILE | wc -l`
echo "The number of core is" $np
mpirun -np $np python3 $src_path/galsim_accu_test.py c 500 1
mpirun -np $np python3 $src_path/galsim_accu_test.py e 500 1
rm -rf host.mpd
