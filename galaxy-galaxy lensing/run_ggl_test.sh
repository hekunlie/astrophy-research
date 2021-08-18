
data_path=/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z/result/dilution_test/dilution

for i in 0.00 0.05 0.10 0.15 0.20;
do
    for j in 0 1 2;
    do
    mpirun -np 15 python ggl_test.py $j $i 0

    cp -r ${data_path}/run/$j ${data_path}/$i
    done

done
