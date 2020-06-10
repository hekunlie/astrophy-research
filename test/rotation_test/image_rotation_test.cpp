#include<FQlib.h>
#include<hk_iolib.h>

#define MY_FLOAT float

int main(int argc, char **argv)
{
    int size;
    char data_path[300];
    
    size = 44;

    MY_FLOAT *image_in = new MY_FLOAT[size*size];
    MY_FLOAT *image_out = new MY_FLOAT[size*size];

    sprintf(data_path,"image_in.fits");
    read_fits(data_path, image_in);

    image_rotation(image_in, image_out, size);

    sprintf(data_path,"image_out.fits");
    write_fits(data_path, image_out, size, size);

    return 0;
}