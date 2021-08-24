
data_path=/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_1

for i in 0.00 0.05 0.10 0.15 0.20;
do
    for j in 0 1 2;
    do
    mpirun -np 11 python /home/hklee/work/Galaxy_Galaxy_lensing_test/code/ggl_test.py $j $i noise_free

    done

done
