#include<FQlib.h>
#include<hk_iolib.h>

int main(int argc, char **argv)
{
    int i,j,k;
    int seed;
    int pts_num ;
    double max_radius;
    double flux_per_pt;
    int size;
    double g1, g2;

    int psf_type;
    double psf_scale;

    int rotation_times;
    int num_scale;
    int tag;

    para all_para;

    size = 64;
    psf_scale = 4;
    psf_type = 2;
    pts_num = 100;
    flux_per_pt = 10;
    max_radius = 10;
    all_para.stamp_size = size;

    rotation_times = 4;
    num_scale = 1;
    // input g1, g2
    g1 = atof(argv[1]);
    g2 = atof(argv[2]);
    seed = atoi(argv[3]);
    gsl_initialize(seed);

    double *pts = new double[2*pts_num]();

    double *check_img = new double[size*size*(rotation_times+1)]();
    double *gal_img = new double[size*size]();
    double *psf_img = new double[size*size]();
    double *gal_img_pow = new double[size*size]();
    double *psf_img_pow = new double[size*size]();

    // create random points(x,y)
    create_points(pts,pts_num, max_radius);
    // create psf
    create_psf(psf_img, psf_scale, size, psf_type);
    // power spectrum
    pow_spec(psf_img, psf_img_pow,size, size);
    // get half light radius
    get_psf_radius(psf_img_pow,&all_para, 2);
    std::cout<<"PSF HLR: "<<all_para.psf_hlr<<std::endl;

    stack(check_img, psf_img, 0, size, 1, rotation_times+1);

    // shear estimators
    double *mg1 = new double[rotation_times*num_scale]();
    double *mg2 = new double[rotation_times*num_scale]();
    double *mn = new double[rotation_times*num_scale]();
    
    for(j=0;j<num_scale;j++)
    {   
        tag = j*rotation_times;
        for(i=0;i<rotation_times;i++)
        {   
            initialize_arr(gal_img,size*size,0);
            initialize_arr(gal_img_pow,size*size,0);

            // draw the galaxy image        
            convolve(gal_img, pts, flux_per_pt, size, pts_num, i, psf_scale, g1, g2, psf_type, 1, &all_para);

            // power spectrum
            pow_spec(gal_img, gal_img_pow, size, size);
            shear_est(gal_img_pow, psf_img_pow, &all_para);

            mg1[tag+i] = all_para.n1;
            mg2[tag+i] = all_para.n2;
            mn[tag+i] = all_para.dn;

            if(j == 0)stack(check_img, gal_img, i+1,size, 1, rotation_times+1);
        }

    }

    double g1_s,g2_s, n_s;
    double est_g1, est_g2;
    g1_s = 0;
    g2_s = 0;
    n_s = 0;

    sum_arr(mg1, rotation_times*num_scale, 0, rotation_times*num_scale, g1_s);
    sum_arr(mg2, rotation_times*num_scale, 0, rotation_times*num_scale, g2_s);
    sum_arr(mn, rotation_times*num_scale, 0, rotation_times*num_scale, n_s);

    est_g1 = g1_s/n_s;
    est_g2 = g2_s/n_s;

    char inform[200];
    sprintf(inform,"Est g1: %8.5f, g2: %8.5f\nTrue g1: %8.5f, g2: %8.5f\ndiff g1: %8.5f, g2: %8.5f",
            est_g1,est_g2, g1, g2, g1-est_g1, g2-est_g2);
    std::cout<<inform<<std::endl;

    sprintf(inform, "!sample.fits");
    write_fits(inform, check_img, size, size*(rotation_times+1));

    delete[] mg1;
    delete[] mg2;
    delete[] mn;
    delete[] gal_img;
    delete[] gal_img_pow;
    delete[] psf_img;
    delete[] psf_img_pow;
    delete[] check_img;
    delete[] pts;
    gsl_free();
    return 0;
}