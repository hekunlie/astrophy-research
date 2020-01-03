#include<FQlib.h>
#include<hk_iolib.h>

int main(int argc, char**argv)
{
    int i,j,k;
    double psf_scale, theta, ellip;
    int stamp_size;
    double img_cent;
    char stamp_path[200];
    double *psf_img, *gal_img;
    double img_sum;

    double *pts = new double[8]{};

    stamp_size = atoi(argv[1]);
    psf_scale = atof(argv[2]);
    theta = atof(argv[3])*Pi;
    ellip = atof(argv[4]);
    
    img_cent = stamp_size*0.5 -0.5;
    psf_img = new double[stamp_size*stamp_size]{};
    gal_img = new double[stamp_size*stamp_size]{};


    std::cout<<"Size: "<<stamp_size<<" PSF scale: "<<psf_scale<<" theta: "<<theta/Pi<<"Pi ellipticity: "<<ellip<<std::endl;

    std::cout<<"Create PSF"<<std::endl;
    
    create_psf_e(psf_img, psf_scale, stamp_size, img_cent, ellip, theta, 2);
    
    std::cout<<"write to file"<<std::endl;
    sprintf(stamp_path, "!epsf.fits");
    write_fits(stamp_path, psf_img, stamp_size, stamp_size);
    
    img_sum = 0;
    for(i=0;i<stamp_size*stamp_size;i++){img_sum+=psf_img[i];}
    std::cout<<"Total flux: "<<img_sum<<std::endl;
    
    convolve_e(gal_img, pts, 1, stamp_size, img_cent, 4, 0, psf_scale, 0, 0, 2, ellip, theta);
    img_sum = 0;
    for(i=0;i<stamp_size*stamp_size;i++){img_sum+=gal_img[i];}
    std::cout<<"Total flux: "<<img_sum<<std::endl;

    for(i=0;i<stamp_size*stamp_size;i++){gal_img[i] -= psf_img[i];}
    show_arr(gal_img,stamp_size, stamp_size);
    sprintf(stamp_path, "!gal.fits");
    write_fits(stamp_path, gal_img, stamp_size, stamp_size);


    delete[] psf_img;

    return 0;
}