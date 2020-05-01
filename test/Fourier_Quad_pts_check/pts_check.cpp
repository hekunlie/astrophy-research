#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>

int main(int argc, char **argv)
{
    int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);


    int i,j,k;
    int seed1,seed2;
    int pts_num ;
    double max_radius;
    double flux_per_pt,total_flux, noise_sig;
    int size, img_cent;
    double g1, g2;

    int psf_type;
    double psf_scale;

    int sub_num, total_num;
    int rotation_times;
    double theta;
    int num_scale;
    int tag, check_label;
    char set_name[300], result_path[300];

    fq_paras all_para;

    size = 44;
    img_cent = size/2-0.5;
    psf_scale = 4;
    psf_type = 2;
    pts_num = 30;
    max_radius = 7;
    all_para.stamp_size = size;
    
    rotation_times = 4;

    // input g1, g2
    seed1 = atoi(argv[1]);
    seed2 = atoi(argv[2]);

    g1 = atof(argv[3]);
    g2 = atof(argv[4]);
    num_scale = atoi(argv[5]);
    total_flux = atof(argv[6]);
    flux_per_pt = total_flux/pts_num;

    noise_sig = atof(argv[7]);
    check_label = 30;

    char inform[300],inform1[300],inform2[300];
    sprintf(inform, "Seed: %d %d, g1: %.5f, g2: %.5f, Num_pts: %d, Flux_per_pt: %.5f, Noise sigma: %.5f",seed1,seed2, g1, g2, pts_num, flux_per_pt, noise_sig);
    std::cout<<inform<<std::endl;

    

    double *pts = new double[2*pts_num]();
    double *pts_r = new double[2*pts_num]();

    double *check_img = new double[size*size*rotation_times*check_label]();
    double *gal_img = new double[size*size]();
    double *psf_img = new double[size*size]();
    double *gal_img_pow = new double[size*size]();
    double *psf_img_pow = new double[size*size]();
    double *noise = new double[size*size]{};
    double *pnoise = new double[size*size]{};

    gsl_initialize(seed1,0);
    gsl_initialize(seed2+rank,1);

    // create random points(x,y)
    create_points(pts,pts_num, max_radius, 1, rng0);
    // if(rank ==0)
    // {
    //     sprintf(set_name,"/data");
    //     sprintf(result_path,"pts_point.hdf5");
    //     write_h5(result_path, set_name, pts, 2, pts_num, true);
    // }
    for(i=0;i<numprocs;i++)
    {
        if(rank == i)
        {std::cout<<"------------------"<<std::endl;show_arr(pts,2,pts_num);std::cout<<"------------------"<<std::endl;}
        MPI_Barrier(MPI_COMM_WORLD);
    }
        
    // create psf
    create_psf(psf_img, psf_scale, size, img_cent, psf_type);
    // power spectrum
    pow_spec(psf_img, psf_img_pow,size, size);
    // get half light radius
    get_psf_radius(psf_img_pow, &all_para, 2);
    std::cout<<"PSF HLR: "<<all_para.psf_hlr<<std::endl;

    // stack(check_img, psf_img, 0, size, 1, rotation_times+1);

    sub_num = rotation_times*num_scale;
    total_num = sub_num*numprocs;

    // shear estimators
    double *mg1 = new double[sub_num]();
    double *mg2 = new double[sub_num]();
    double *mn = new double[sub_num]();
    double *mnu1 = new double[sub_num]();
    double *mnu2 = new double[sub_num]();
    
    double *total_mg1;
    double *total_mg2;
    double *total_mn;
    double *total_mnu1;
    double *total_mnu2;
    int *gather_count = new int[numprocs];
    for(i=0;i<numprocs;i++){gather_count[i] = sub_num;}
    if(rank == 0)
    {
        total_mg1 = new double[total_num]();
        total_mg2 = new double[total_num]();
        total_mn = new double[total_num]();
        total_mnu1 = new double[total_num]();
        total_mnu2 = new double[total_num]();
    }

    for(j=0;j<num_scale;j++)
    {   
        tag = j*rotation_times;
        for(i=0;i<rotation_times;i++)
        {   
            initialize_arr(gal_img,size*size,0);
            initialize_arr(gal_img_pow,size*size,0);
            initialize_arr(noise, size*size, 0);
            initialize_arr(pnoise, size*size, 0);
            initialize_arr(pts_r, 2*pts_num, 0);
            // draw the galaxy image        
            
            theta = Pi/4*i;
            rand_uniform(0, Pi*2, theta, rng1);
            coord_rotation(pts, pts_num, theta, pts_r);
            convolve(gal_img, pts_r, flux_per_pt, size, img_cent, pts_num, 0, psf_scale, g1, g2, psf_type);
            
            // convolve(gal_img, pts, flux_per_pt, size, img_cent, pts_num, i, psf_scale, g1, g2, psf_type);

            // addnoise(gal_img, size*size, noise_sig, rng1);
            // // power spectrum
            pow_spec(gal_img, gal_img_pow, size, size);

            // addnoise(noise, size*size, noise_sig, rng1);
            // pow_spec(noise, pnoise, size,size);

            // noise_subtraction(gal_img_pow, pnoise, &all_para, 1, 0);


            shear_est(gal_img_pow, psf_img_pow, &all_para);

            mg1[tag+i] = all_para.n1;
            mg2[tag+i] = all_para.n2;
            mn[tag+i] = all_para.dn;
            mnu1[tag+i] = all_para.dn + all_para.du;
            mnu2[tag+i] = all_para.dn - all_para.du;

            // std::cout<<all_para.n1<<" "<<all_para.n2<<" "<<all_para.dn<<std::endl;
            if(tag+i < check_label*rotation_times)stack(check_img, gal_img, tag+i, size, check_label, rotation_times);
        }

    }
      
    sprintf(inform, "!sample_%d.fits", rank);
    write_fits(inform, check_img, check_label*size, size*rotation_times);
    sprintf(result_path,"data_%d.hdf5", rank);
    sprintf(set_name,"/mg1");
    write_h5(result_path,set_name, mg1, sub_num, 1, true);
    sprintf(set_name,"/mg2");
    write_h5(result_path,set_name, mg2, sub_num, 1, false);
    sprintf(set_name,"/mn");
    write_h5(result_path,set_name, mn, sub_num, 1, false);
    sprintf(set_name,"/mnu1");
    write_h5(result_path,set_name, mnu1, sub_num, 1, false);
    sprintf(set_name,"/mnu2");
    write_h5(result_path,set_name, mnu2, sub_num, 1, false);

    double gh1, gh1_sig;
    double gh2, gh2_sig;
    double chimin;
    double *check_chisq = new double[40]{};

    // find_shear_mean(mg1, mn, sub_num, gh1, gh1_sig, 100);
    // find_shear_mean(mg2, mn, sub_num, gh2, gh2_sig, 100);
    // sprintf(inform1,"Est g1: %8.5f(%8.5f), g2: %8.5f(%8.5f)\nTrue g1: %8.5f, g2: %8.5f\ndiff g1: %8.5f, g2: %8.5f",
    //             gh1,gh1_sig, gh2, gh2_sig, g1, g2, g1-gh1, g2-gh2);
    
    // find_shear(mg1, mnu1, sub_num, 10, gh1, gh1_sig, chimin, check_chisq,20,0,10,-0.05,0.05,40);
    // find_shear(mg2, mnu2, sub_num, 10, gh2, gh2_sig, chimin, check_chisq,20,0,10,-0.05,0.05,40);
    // sprintf(inform2,"Est g1: %8.5f(%8.5f), g2: %8.5f(%8.5f)\nTrue g1: %8.5f, g2: %8.5f\ndiff g1: %8.5f, g2: %8.5f",
    //             gh1,gh1_sig, gh2, gh2_sig, g1, g2, g1-gh1, g2-gh2);
    
    // for(i=0;i<numprocs;i++)
    // {
    //     if(rank==i)
    //     {
    //         std::cout<<"--------------------"<<rank<<"-------------------------"<<std::endl;
    //         std::cout<<inform1<<std::endl;
    //         std::cout<<inform2<<std::endl;
    //         std::cout<<"---------------------------------------------"<<std::endl;
            
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    MPI_Barrier(MPI_COMM_WORLD);
    my_Gatherv(mg1, gather_count, total_mg1, numprocs, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    my_Gatherv(mg2, gather_count, total_mg2, numprocs, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    my_Gatherv(mn, gather_count, total_mn, numprocs, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    my_Gatherv(mnu1, gather_count, total_mnu1, numprocs, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    my_Gatherv(mnu2, gather_count, total_mnu2, numprocs, rank);
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {

        sprintf(result_path,"data.hdf5");
        sprintf(set_name,"/mg1");
        write_h5(result_path,set_name, total_mg1, total_num, 1, true);
        sprintf(set_name,"/mg2");
        write_h5(result_path,set_name, total_mg2, total_num, 1, false);
        sprintf(set_name,"/mn");
        write_h5(result_path,set_name, total_mn, total_num, 1, false);
        sprintf(set_name,"/mnu1");
        write_h5(result_path,set_name, total_mnu1, total_num, 1, false);
        sprintf(set_name,"/mnu2");
        write_h5(result_path,set_name, total_mnu2, total_num, 1, false);


        find_shear_mean(total_mg1, total_mn, total_num, gh1, gh1_sig, 100);
        find_shear_mean(total_mg2, total_mn, total_num, gh2, gh2_sig, 100);

        sprintf(inform1,"Est g1: %8.5f(%8.5f), g2: %8.5f(%8.5f)\nTrue g1: %8.5f, g2: %8.5f\ndiff g1: %8.5f, g2: %8.5f",
                gh1,gh1_sig, gh2, gh2_sig, g1, g2, g1-gh1, g2-gh2);

        std::cout<<"--------------------total-------------------------"<<std::endl;
        std::cout<<inform1<<std::endl;

        find_shear(total_mg1, total_mnu1, total_num, 10, gh1, gh1_sig, chimin, check_chisq,20,0,100,-0.05,0.05,40);
        find_shear(total_mg2, total_mnu2, total_num, 10, gh2, gh2_sig, chimin, check_chisq,20,0,100,-0.05,0.05,40);

        sprintf(inform2,"Est g1: %8.5f(%8.5f), g2: %8.5f(%8.5f)\nTrue g1: %8.5f, g2: %8.5f\ndiff g1: %8.5f, g2: %8.5f",
                gh1,gh1_sig, gh2, gh2_sig, g1, g2, g1-gh1, g2-gh2);
        

        std::cout<<inform2<<std::endl;
        std::cout<<"---------------------------------------------"<<std::endl;
        

    }
    delete[] mg1;
    delete[] mg2;
    delete[] mn;
    delete[] gal_img;
    delete[] gal_img_pow;
    delete[] psf_img;
    delete[] psf_img_pow;
    delete[] check_img;
    delete[] pts;

    MPI_Finalize();
    return 0;
}