#include<hk_FQlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>
#define IMG_CHECK_LABEL 3
#define FLUX_PDF_UNI
#define CPSF

#define SAVE_MEM

#define DATA_SEP

#ifdef SAVE_MEM
#define MY_FLOAT float
#else
#define MY_FLOAT double
#endif

void arr_sqrt(MY_FLOAT *arr_in, MY_FLOAT *arr_out, const int length);

void arr_scale(MY_FLOAT *arr1, const MY_FLOAT scale,const int length);

void arr_add(MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length);

void arr_add(MY_FLOAT *arr1, const MY_FLOAT*arr2,const MY_FLOAT*arr3,const int length);

void arr_deduct(MY_FLOAT *result_buff, const MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length);
void arr_deduct(MY_FLOAT *result_buff, const MY_FLOAT *arr1, const MY_FLOAT*arr2, const MY_FLOAT *arr3, const int length);
void arr_copy(MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length);

void noise_subtraction_new(MY_FLOAT *pow_arr1, const MY_FLOAT*pow_arr2,const MY_FLOAT *noise_pow_arr, const int length);
void data_sep(const MY_FLOAT *data_arr, MY_FLOAT **sub_data, const int data_row, const int data_col);

void write_data(char *result_path, MY_FLOAT *total_data, char **mg_name, MY_FLOAT **mg_data, int total_data_row, int shear_data_cols);

void task_alloc_d(const int total_task_num, const int portion, const int my_part_id, int &my_start, int &my_end, int *task_count);


int main(int argc, char*argv[])
{	
	/* it is designed to run on the cluster (should save the memory) */
	/* the total chips should be divisible by the numprocs !!!! */
	/* loop the shear points, the chips will be distributed to threads */

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
	char parent_path[400], data_path[500], chip_path[500], para_path[500], shear_path[500], result_path[500], log_path[500];
	char buffer[500], log_inform[500], set_name[50];

	int i, j, k, ib, m,n;
	int seed_ini, seed_pts, seed_n1, seed_n2, seed_step;

	int num_p, size, source_tag, img_len;
    MY_FLOAT max_radius;
	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id;

    MY_FLOAT *radius, *flux;
    MY_FLOAT psf_scale, psf_type, psf_thresh_scale, psf_noise_sig;
    MY_FLOAT sig_level,gal_noise_sig, flux_i;
    MY_FLOAT g1, g2, ts, te, t1, t2;
    int gg_len, gg_count;
    MY_FLOAT *g1t, *g2t;
    MY_FLOAT psf_ellip, ellip_theta;
    MY_FLOAT img_cent;
    MY_FLOAT pts_step;
    int *rand_seed;
    MY_FLOAT *big_img[2];
    char src_type[50];

    strcpy(parent_path, argv[1]);
    source_tag = atoi(argv[2]);

	size = 56;//atoi(argv[4]);
    total_chips = atoi(argv[3]);
    strcpy(src_type, argv[4]);

	seed_step = 2;

    if(rank == 0)
    {
        std::cout<<parent_path<<std::endl;
        std::cout<<size<<" "<<source_tag<<std::endl;
    }
	pts_step = 3;

	num_p = 40;

	stamp_num = 10000;
	shear_data_cols = 5;

	psf_type = 2;
	psf_scale = 4;
	psf_ellip = 0.10;
	ellip_theta = 0;
	psf_thresh_scale = 2.;

	sig_level = 1.5;
	psf_noise_sig = 0;
    gal_noise_sig = 60;
	
    img_len = size*size;
	img_cent = size*0.5 - 0.5;
    stamp_nx = 100;


	
#ifdef SAVE_MEM
	fq_paras_float all_paras;
#else
	fq_paras all_paras;
#endif

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.gal_noise_sig = gal_noise_sig;
	
	sprintf(log_path, "%s/log/%02d.dat", parent_path, rank);    

    total_data_row = total_chips*stamp_num;


    if(rank == 0)
    {
        std::cout<<total_chips<<" chips"<<std::endl;
    }
	/////////////////// task distribution /////////////////////////////////////
	int *scatter_count,*gather_count;
	// for scatterring the flux to each thread
	scatter_count = new int[numprocs]{};
	// for the gatherv when finish in each shear point
	gather_count = new int[numprocs]{};
    // std::cout<<rank<<" "<<numprocs<<std::endl;
	task_alloc(total_chips, numprocs, rank, chip_st, chip_ed, scatter_count);
    // for(i=0;i<numprocs;i++){scatter_count[i] = 1;}
	sub_chip_num = scatter_count[rank];
	// the sub-data from the source processed by each thread
	sub_data_row = sub_chip_num * stamp_num;
	for(i=0;i<numprocs;i++)
	{
		// the real count of galaxies for each thread
		scatter_count[i] = scatter_count[i]*stamp_num;
		// the real amount of data of each thread
		gather_count[i] = scatter_count[i]*shear_data_cols;
	}

    if(rank == 0)
    {
        big_img[0] = new MY_FLOAT[stamp_nx*stamp_nx*size*size];
        big_img[1] = new MY_FLOAT[stamp_nx*stamp_nx*size*size];
    }

	MY_FLOAT *point = new MY_FLOAT[2 * num_p]{};

	MY_FLOAT *stamp_img[20];
	MY_FLOAT *stamp_pow_img[20];
	MY_FLOAT *psf_img[10];
	MY_FLOAT *noise_img[10];
	MY_FLOAT *noise_pow_img[10];

	
	// the shear estimators data matrix  
	MY_FLOAT *total_data;	
	MY_FLOAT *sub_noise_free_data;
    MY_FLOAT *sub_noisy_data;
    


	for(i=0;i<10;i++)
	{
		stamp_img[i] = new MY_FLOAT[img_len];
		stamp_pow_img[i] = new MY_FLOAT[img_len];
		noise_img[i] = new MY_FLOAT[img_len];
		noise_pow_img[i] = new MY_FLOAT[img_len];
		if(i<6){psf_img[i] = new MY_FLOAT[img_len]{};}
	}
    

	sub_noise_free_data = new MY_FLOAT[gather_count[rank]]{};
	sub_noisy_data = new MY_FLOAT[gather_count[rank]]{};
	
    g1t = new MY_FLOAT[total_chips*stamp_num];
	g2t = new MY_FLOAT[total_chips*stamp_num];
    radius = new MY_FLOAT[total_chips*stamp_num]{};
    flux = new MY_FLOAT[total_chips*stamp_num]{};
    rand_seed = new int[total_chips]{};


    sprintf(data_path,"%s/cata/background/continue_source_z", parent_path);

    sprintf(para_path,"%s/params/%s_para_%d.hdf5", data_path, src_type, source_tag);

    sprintf(set_name,"/gamma1");
    read_h5(para_path, set_name, g1t);
    sprintf(set_name,"/gamma2");
    read_h5(para_path, set_name, g2t);

    sprintf(set_name,"/flux");
    read_h5(para_path, set_name, flux);

    sprintf(set_name,"/radius");
    read_h5(para_path, set_name, radius);

    sprintf(set_name,"/seed");
    read_h5(para_path, set_name, rand_seed);


	create_psf(psf_img[0], psf_scale, size, img_cent, psf_type);

	pow_spec(psf_img[0], psf_img[1], size, size);
	get_psf_radius(psf_img[1], &all_paras, psf_thresh_scale);



	if (0 == rank)
	{	
		std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		std::cout << parent_path << std::endl;
		std::cout << "Total chip: " << total_chips<< ", Stamp size: " << size  << std::endl;
		std::cout << "Total cpus: " << numprocs << std::endl;
		std::cout <<"PSF Scale: "<<psf_scale<< " PSF THRESH: " << all_paras.psf_pow_thresh <<" PSF HLR: " << all_paras.psf_hlr << std::endl;
		std::cout <<"MAX RADIUS: "<< max_radius <<" , Step: "<<pts_step<< ", SIG_LEVEL: " << sig_level <<"sigma"<< std::endl;

        sprintf(chip_path,"%s/imgs/%s_psf.hdf5", data_path, src_type);
        sprintf(set_name,"/data");
        write_h5(chip_path, set_name, psf_img[0], size, size, true);

		std::cout<<"Gal Num of each thread: ";
		show_arr(scatter_count,1,numprocs);
		std::cout<<"---------------------------------------------------------------------------"<<std::endl<<std::endl;
	}
	for(i=0;i<10;i++)
	{
		if(i == rank)
		{
			std::cout<<rank<<" "<<chip_st<<" "<<chip_ed<<" "<<scatter_count[i]<<" "<<gather_count[i]<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
    MPI_Barrier(MPI_COMM_WORLD);
	
	// loop the shear points

    ts = clock();

    sprintf(log_inform, "size: %d, total chips: %d (%d cpus),  point num: %d , noise sigma: %.2f. Numprocs: %d", size, total_chips, numprocs, num_p, gal_noise_sig, numprocs);
    write_log(log_path, log_inform);
    sprintf(log_inform, "PSF scale: %.2f, max radius: %.2f, Step: %.4f", psf_scale, max_radius, pts_step);
    write_log(log_path, log_inform);
    sprintf(log_inform, "RANK: %03d, my chips: %d - %d", rank, chip_st, chip_ed);
    write_log(log_path, log_inform);
    
    MPI_Barrier(MPI_COMM_WORLD);
    



    // loop the chips
    for (i = chip_st; i < chip_ed; i++)
    {
        // initialize GSL
        seed_pts = rand_seed[i];
        seed_n1 = rand_seed[i] + i;
        gsl_initialize(seed_pts,0);
        gsl_initialize(seed_n1,1);

        t1 = clock();

        sprintf(log_inform, "RANK: %03d:, chip: %05d, start. seed: %d, %d", rank, i, seed_pts, seed_n1);
        write_log(log_path, log_inform);
        if (rank == 0)
        {
            std::cout << log_inform << std::endl;
        }

        row = (i-chip_st)*stamp_num*shear_data_cols;
        
        // loop the stamps
        for (j = 0; j < stamp_num; j++)
        {	

            g1 = g1t[i*stamp_num + j];
            g2 = g2t[i*stamp_num + j];

            max_radius = radius[i*stamp_num + j];
            flux_i = flux[i*stamp_num + j]/num_p;

            initialize_arr(point, num_p * 2, 0);				
            for(k=0;k<8;k++)
            {
                initialize_arr(stamp_img[k], size*size, 0);
                initialize_arr(noise_img[k], size*size, 0);
            }
            for(k=0;k<3;k++)
            {					
                addnoise(noise_img[k], img_len, gal_noise_sig, rng1);
                pow_spec(noise_img[k], noise_pow_img[k], size, size);
            }

            create_points(point, num_p, max_radius, pts_step, rng0);

            // convolve_e(point,num_p,flux_i, g1, g2, stamp_img[0], size, img_cent, psf_scale,psf_type,psf_ellip, ellip_theta);

            convolve(point, num_p, flux_i, g1, g2, stamp_img[0], size, img_cent, psf_scale, psf_type);


            // noise free
            pow_spec(stamp_img[0], stamp_pow_img[0], size, size);
            // noisy image
            arr_add(stamp_img[1], stamp_img[0], noise_img[0], img_len);
            pow_spec(stamp_img[1], stamp_pow_img[1], size, size);
            arr_deduct(stamp_pow_img[2], stamp_pow_img[1], noise_pow_img[1], img_len);
            

            if(rank == 0 and i < 2)
            {
                stack(big_img[0], stamp_img[0], j, size, stamp_nx, stamp_nx);
                stack(big_img[1], stamp_img[1], j, size, stamp_nx, stamp_nx);
            }


            /////////////////////// Noise free /////////////////////////
            shear_est(stamp_pow_img[0], psf_img[1], &all_paras);
            sub_noise_free_data[row + j * shear_data_cols] = all_paras.n1;
            sub_noise_free_data[row + j * shear_data_cols + 1] = all_paras.n2;
            sub_noise_free_data[row + j * shear_data_cols + 2] = all_paras.dn;
            sub_noise_free_data[row + j * shear_data_cols + 3] = all_paras.du;
            sub_noise_free_data[row + j * shear_data_cols + 4] = all_paras.dv;

            /////////////////////////// noisy image /////////////////////////////////
            shear_est(stamp_pow_img[2], psf_img[1], &all_paras);
            sub_noisy_data[row + j * shear_data_cols] = all_paras.n1;
            sub_noisy_data[row + j * shear_data_cols + 1] = all_paras.n2;
            sub_noisy_data[row + j * shear_data_cols + 2] = all_paras.dn;
            sub_noisy_data[row + j * shear_data_cols + 3] = all_paras.du;
            sub_noisy_data[row + j * shear_data_cols + 4] = all_paras.dv;

        }
        t2 = clock();
        sprintf(log_inform, "RANK: %03d, chip: %05d, done in %.2f s.", rank, i, (t2 - t1) / CLOCKS_PER_SEC);
        write_log(log_path, log_inform);
        if (rank == 0)
        {
            std::cout << log_inform << std::endl;
        }
        if(rank == 0 and i < 2)
        {
            sprintf(chip_path, "%s/imgs/%s_%d_gal_chip_%05d_noise_free.hdf5", data_path, src_type, source_tag, i);
            sprintf(set_name,"/data");
            write_h5(chip_path, set_name, big_img[0], stamp_nx*size, stamp_nx*size, true);

            sprintf(chip_path, "%s/imgs/%s_%d_gal_chip_%05d_noisy.hdf5", data_path, src_type, source_tag, i);
            write_h5(chip_path, set_name, big_img[1], stamp_nx*size, stamp_nx*size, true);

            // write_fits(chip_path, big_img_check[0], stamp_nx*size, stamp_nx*size);
        }

        gsl_free(0);
        gsl_free(1);
        // gsl_free(2);
    }


    sprintf(log_inform, "RANK: %03d, gather data", rank);
    write_log(log_path, log_inform);
    if (rank == 0)
    {
        std::cout << log_inform << std::endl;
    }

    // finish the chip loop
    // MPI_Barrier(MPI_COMM_WORLD);
    // for(i=0; i<numprocs;i++)
    // {
    //     if(rank == i)
    //     {
    //         sprintf(result_path, "%s/data/data_%d_noise_free.hdf5", parent_path, shear_tag);
    //         sprintf(set_name,"/data_%d", i);
    //         if(i == 0)
    //         {write_h5(result_path, set_name, sub_noise_free_data, sub_data_row, shear_data_cols,true);}
    //         else
    //         {write_h5(result_path, set_name, sub_noise_free_data, sub_data_row, shear_data_cols,false);}


    //         sprintf(result_path, "%s/data/data_%d_noisy_cpp.hdf5", parent_path, shear_tag);
    //         sprintf(set_name,"/data_%d", i);
    //         if(i == 0)
    //         {write_h5(result_path, set_name, sub_noisy_data, sub_data_row, shear_data_cols,true);}
    //         else
    //         {write_h5(result_path, set_name, sub_noisy_data, sub_data_row, shear_data_cols,false);}
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    if(rank == 0){total_data = new MY_FLOAT[total_chips*stamp_num*shear_data_cols];}
    MPI_Barrier(MPI_COMM_WORLD);

    my_Gatherv(sub_noise_free_data, gather_count, total_data, numprocs, rank);
    if (0 == rank)
    {
        sprintf(result_path, "%s/data/%s_data_%d_noise_free.hdf5", data_path, src_type, source_tag);
        sprintf(set_name,"/data");
        write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols,true);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    my_Gatherv(sub_noisy_data, gather_count, total_data, numprocs, rank);
    if (0 == rank)
    {
        sprintf(result_path, "%s/data/%s_data_%d_noisy_cpp.hdf5", data_path, src_type, source_tag);
        sprintf(set_name,"/data");
        write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols,true);
    }
    MPI_Barrier(MPI_COMM_WORLD);



    te = clock();
    if (rank == 0)
    {
        sprintf(buffer, "rank %d: done in %.2f sec ", rank, (te - ts) / CLOCKS_PER_SEC);
        std::cout << buffer<<std::endl;
        std::cout<<"---------------------------------------------------------------------------"<<std::endl;
    }
    // return 0;
	MPI_Finalize();
	return 0;
}



void arr_sqrt(MY_FLOAT *arr_in, MY_FLOAT *arr_out, const int length)
{
    for(int i=0;i<length;i++)
    {
        arr_out[i] = sqrt(arr_in[i]);
    }
}

void arr_scale(MY_FLOAT *arr1, const MY_FLOAT scale,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] += arr1[i]*scale;
    }
}


void arr_add(MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] += arr2[i];
    }
}

void arr_add(MY_FLOAT *arr1, const MY_FLOAT*arr2,const MY_FLOAT*arr3,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] = arr2[i] + arr3[i];
    }
}

void arr_deduct(MY_FLOAT *result_buff, const MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        result_buff[i] = arr1[i] - arr2[i];
    }
}
void arr_deduct(MY_FLOAT *result_buff, const MY_FLOAT *arr1, const MY_FLOAT*arr2, const MY_FLOAT *arr3, const int length)
{
    for(int i=0;i<length;i++)
    {
        result_buff[i] = arr1[i] - arr2[i] - arr3[i];
    }
}
void arr_copy(MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] = arr2[i];
    }
}

void noise_subtraction_new(MY_FLOAT *pow_arr1, const MY_FLOAT*pow_arr2,const MY_FLOAT *noise_pow_arr, const int length)
{
    for(int i=0;i<length;i++)
    {
        pow_arr1[i] = 2*pow_arr1[i] - pow_arr2[i];// + 2*noise_pow_arr[i];
    }
}
void data_sep(const MY_FLOAT *data_arr, MY_FLOAT **sub_data, const int data_row, const int data_col)
{
	int i, j;
	for(i=0;i<data_row;i++)
	{
		for(j=0;j<data_col;j++)
		{
			sub_data[j][i] = data_arr[i*data_col+j];
		}
	}
}

void write_data(char *result_path, MY_FLOAT *total_data, char **mg_name, MY_FLOAT **mg_data, int total_data_row, int shear_data_cols)
{	
	int k;
	data_sep(total_data, mg_data, total_data_row, shear_data_cols);
	for(k=0; k<shear_data_cols;k++)
	{	
		if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
		else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
	}		
}

void task_alloc_d(const int total_task_num, const int portion, const int my_part_id, int &my_start, int &my_end, int *task_count)
{	
	int i,j;
	int sub_num;
    int st,ed;
    int *temp_count = new int[portion]{};
	sub_num = total_task_num / portion;
	j = total_task_num%portion;
	for(i=0;i<portion;i++)
	{
		temp_count[i] = sub_num;
		// if(i<j){temp_count[i] = sub_num + 1;}
        // else{temp_count[i] = sub_num;}
	}
	st = 0;
	for(i=0;i<my_part_id;i++)
	{
		st += temp_count[i];
	}
	ed = st + temp_count[my_part_id];
    my_start = st;
    my_end = ed;
    for(i=0; i<portion; i++){task_count[i] = temp_count[i];}
    delete[] temp_count;
}
