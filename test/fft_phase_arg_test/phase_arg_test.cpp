#include<FQlib.h>
#include<hk_iolib.h>

#define MY_FLOAT double

int main(int argc, char*argv[])
{
    int size;
    int seed;

    size = atoi(argv[1]);
    seed = atoi(argv[2]);

    gsl_initialize(seed,0);
    
    MY_FLOAT *stamp = new MY_FLOAT[size*size]{};
    MY_FLOAT *stamp_pow = new MY_FLOAT[size*size]{};
    MY_FLOAT *stamp_pow_real = new MY_FLOAT[size*size]{};
    MY_FLOAT *stamp_pow_imag = new MY_FLOAT[size*size]{};

    MY_FLOAT *phase_arg = new MY_FLOAT[size*size]{};

    addnoise(stamp, size*size, 40, rng0);

    pow_spec(stamp, stamp_pow, stamp_pow_real, stamp_pow_imag, phase_arg, size, size);


    char data_path[300], set_name[40];

    sprintf(data_path, "test.hdf5");

    sprintf(set_name, "/img");
    write_h5(data_path, set_name, stamp, size, size, true);

    sprintf(set_name, "/img_pow");
    write_h5(data_path, set_name, stamp_pow, size, size, false);

    sprintf(set_name, "/img_pow_real");
    write_h5(data_path, set_name, stamp_pow_real, size, size, false);

    sprintf(set_name, "/img_pow_imag");
    write_h5(data_path, set_name, stamp_pow_imag, size, size, false);

    sprintf(set_name, "/img_arg");
    write_h5(data_path, set_name, phase_arg, size, size, false);

    return 0;
}