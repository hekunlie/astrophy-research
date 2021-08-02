code_dir=/lustre/home/acct-phyzj/phyzj/jzhang/pipe_CFHT_June_2021
data_dir=/lustre/home/acct-phyzj/phyzj/CFHT/i
cat_dir=/lustre/home/acct-phyzj/phyzj/CFHT/galaxy
gaia_dir=/lustre/home/acct-phyzj/phyzj/CFHT/gaia_dr2
gen_list_dir=$code_dir/gen_list_CFHT

function prepare_fits(){
    cd $1
    for dir in w?????
    do
	cd $1/$dir  
	for file in *.fz
	do
	    funpack $file
	    echo $dir $file
	done
	mkdir science
	for file in *.fits
	do
	    python $2/python/separate_fits.py  $file science
	    echo $dir $file
	    rm $file
	done
    done
}

#prepare_fits $data_dir $code_dir

function prepare_list(){
    cd $1
    rm nname.dat

    for dir in w?????
    do
	echo $dir >> $1/nname.dat
	for name in $1/$dir/science/*.fits
	do
          echo ${name##*/} >> $1/nname.dat
	done
    done
    
    cd $2
    mpif77 *.f -o main
    ./main $1 $3 $4
}

#prepare_list $data_dir $gen_list_dir $cat_dir $gaia_dir

function prepare_run(){
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

#prepare_run $data_dir 

cd $code_dir
mpif77 *.f -o main -mcmodel=medium -lcfitsio

