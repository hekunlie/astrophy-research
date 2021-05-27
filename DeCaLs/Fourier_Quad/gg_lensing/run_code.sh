data_path=/home/hklee/work/DECALS/gg_lensing
core_num=40
jack_num=100

z1=0.2
z2=0.3
dz=0.002

m1=14
m2=14.5
mpirun -np 1 python prepare_cata.py prepare_foreground ${z1} ${z2} ${m1} ${m2}
mpirun -np 1 python prepare_cata.py prepare_pdf
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_${m1}_${m2}.hdf5


m1=14.5
m2=15
mpirun -np 1 python prepare_cata.py prepare_foreground ${z1} ${z2} ${m1} ${m2}
mpirun -np 1 python prepare_cata.py prepare_pdf
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_${m1}_${m2}.hdf5



z1=0.3
z2=0.4

m1=14
m2=14.5
mpirun -np 1 python prepare_cata.py prepare_foreground ${z1} ${z2} ${m1} ${m2}
mpirun -np 1 python prepare_cata.py prepare_pdf
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_${m1}_${m2}.hdf5


m1=14.5
m2=15
mpirun -np 1 python prepare_cata.py prepare_foreground ${z1} ${z2} ${m1} ${m2}
mpirun -np 1 python prepare_cata.py prepare_pdf
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_${m1}_${m2}.hdf5



z1=0.4
z2=0.5

m1=14
m2=14.5
mpirun -np 1 python prepare_cata.py prepare_foreground ${z1} ${z2} ${m1} ${m2}
mpirun -np 1 python prepare_cata.py prepare_pdf
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_${m1}_${m2}.hdf5


m1=14.5
m2=15
mpirun -np 1 python prepare_cata.py prepare_foreground ${z1} ${z2} ${m1} ${m2}
mpirun -np 1 python prepare_cata.py prepare_pdf
mpirun -np ${core_num} ./GGL_cal /home/hklee/work/DECALS/gg_lensing ${jack_num} ${dz}
mv ${data_path}/result/result.hdf5 ${data_path}/result/result_${z1}_${z2}_${m1}_${m2}.hdf5