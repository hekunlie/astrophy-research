data_path=/home/hklee/work/DECALS/gg_lensing
core_num=55
jack_num=100

z1=0.25
z2=0.45


mpirun -np 1 python prepare_cata.py prepare_foreground ${z1} ${z2} lowz_DECALS_overlap.hdf5
mpirun -np 1 python prepare_cata.py prepare_pdf 

dz=0.3
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_dz_${dz}_lowz.hdf5


dz=0.15
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_dz_${dz}_lowz.hdf5


dz=0.001
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_dz_${dz}_lowz.hdf5


# z1=0.4
# z2=0.55
# mpirun -np 1 python prepare_cata.py prepare_foreground ${z1} ${z2} cmass_DECALS_overlap.hdf5
# mpirun -np 1 python prepare_cata.py prepare_pdf 

# dz=0.3
# mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
# mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_dz_${dz}_cmass.hdf5

# dz=0.15
# mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
# mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_dz_${dz}_cmass.hdf5

# dz=0.001
# mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
# mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_dz_${dz}_cmass.hdf5
