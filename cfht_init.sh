function prepare(){
    cd $1
    for dir in w?????
    do
        rm -r -f $dir/astrometry
        mkdir $dir/astrometry
        rm -r -f $dir/stamps
        mkdir $dir/stamps
        rm -r -f $dir/result
        mkdir $dir/result
    done
}

data_path=/mnt/ddnfs/data_users/hkli/CFHT/w1234/original
prepare $data_path
