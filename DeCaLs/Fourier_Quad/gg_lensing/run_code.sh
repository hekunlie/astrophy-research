data_path=/home/hklee/work/DECALS/DECALS_v210729/gg_lensing
core_num=55
jack_num=150
dz=0.2


for om in 0.20 0.25 0.30 0.35 0.4;
    do
        mpirun -np 5 python prepare_cata.py prepare_foreground $om

        mpirun -np 20 python prepare_cata.py prepare_background $om

        mpirun -np ${core_num} ./GGL_cal ${data_path} cata/background cata/foreground $jack_num ${dz}

        mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${om}.hdf5
    done
# for n in 200;
#     do
#         mpirun -np 5 python prepare_cata.py prepare_foreground 0.0 10 gal_jkf_dr9.hdf5 $n

#         mpirun -np ${core_num} ./GGL_cal ${data_path} cata/background cata/foreground $n ${dz}

#         mv ${data_path}/result/result.hdf5 ${data_path}/result/jkf_compare/dr9_${n}.hdf5
#     done


# for n in 200;
#     do
#         mpirun -np 5 python prepare_cata.py prepare_foreground 0.0 10 gal_jkf_dr8.hdf5 $n

#         mpirun -np ${core_num} ./GGL_cal ${data_path} cata/background cata/foreground $n ${dz}

#         mv ${data_path}/result/result.hdf5 ${data_path}/result/jkf_compare/dr8_${n}.hdf5
#     done

