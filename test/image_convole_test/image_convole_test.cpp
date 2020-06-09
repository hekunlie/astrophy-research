#include<FQlib.h>
#include<hk_iolib.h>

#define MY_FLOAT double

int main(int argc, char *argv[])
{   
    char img_path[200];

    int img_size = 44;
    int kernel_size = 5;
    
    MY_FLOAT *img = new MY_FLOAT[img_size*img_size];
    MY_FLOAT *img_con = new MY_FLOAT[img_size*img_size];
    MY_FLOAT *kernel = new MY_FLOAT[kernel_size*kernel_size];

    sprintf(img_path,"img.fits");
    read_fits(img_path, img);
    
    sprintf(img_path,"kernel.fits");
    read_fits(img_path, kernel);

    image_convole(img, img_con, img_size, kernel, kernel_size);

    sprintf(img_path, "img_con.fits");
    write_fits(img_path,img_con,img_size, img_size);

    return 0;
}