#include <FQlib.h>
#include <mpi.h>
#include<hk_iolib.h>

#define SEX_NUM 6

int main(int argc, char**argv)
{
    /* calculate the gaussian weighted quadrupole                                   */
    /* it loops the shear points, each threads will calculate its own files         */
    /* argv[1]: the source name, one of "galsim_dimmer  pts_bright  pts_dimmer"     */
    /* it will read all "sex" data at once, then calculate the resolution factor for each      */
    /* "sex" catalog, say "sex2_1.5"..., because there is only one image catalog but there are */
    /* "sex" data.                                                                              */

    int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
	MPI_Status mpi_status;

    char total_path[200], img_path[300], sex_data_path[300], result_path[300],set_name[30];
    char log_path[200], log_inform[300];
    char time_now[50];


    char source_name[20];

    char *sex_folder_name[SEX_NUM];
    double *sex_data, *sex_data_shared;

    int total_chips;
    double *img, *stamp;
    double *chip_data, *final_data;

    int chip_st, chip_ed;
    int i,j,k,tag, gal_label;
    int nx, ny, size;
    int sex_row, sex_col,sex_block_size;

    int shear_num, shear_id;

    double st1, st2, st3, st4;

    double gal_quad;
    double eff_radius_sq;

    int quad_tag;
    char inform[40];

    size = 64;
    nx = 100;
    ny = 100;

    img = new double[nx*ny*size*size];
    stamp = new double[size*size];

    sex_row = 10000000;
    sex_col = 4;
    sex_block_size = sex_row*sex_col;
    sex_data = new double[sex_row*sex_col];

    
    strcpy(source_name,argv[1]);
    total_chips = 1000;
    
    sprintf(total_path,"/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s",source_name);
    sprintf(log_path,"logs/%d.dat", rank);
    
    for(i=0;i<SEX_NUM;i++)
    {
        sex_folder_name[i] = new char[20];
    }
    sprintf(sex_folder_name[0],"sex2_1.5");
    sprintf(sex_folder_name[1],"sex2_2");
    sprintf(sex_folder_name[2],"sex2_4");

    sprintf(sex_folder_name[3],"sex4_1.5");
    sprintf(sex_folder_name[4],"sex4_2");
    sprintf(sex_folder_name[5],"sex4_4");


    MPI_Win win_final_data,win_sex_data;
	MPI_Aint final_data_size, sex_data_size;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(SEX_NUM*nx*ny*total_chips* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data, &win_final_data);
		MPI_Win_allocate_shared(SEX_NUM*sex_row*sex_col*nx*ny*total_chips* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &sex_data_shared, &win_sex_data);

	}
	else
	{
		int dispu_total;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data, &win_final_data);
		MPI_Win_shared_query(win_final_data, 0, &final_data_size, &dispu_total, &final_data);

        MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data, &win_sex_data);
		MPI_Win_shared_query(win_sex_data, 0, &final_data_size, &dispu_total, &final_data);
	}

    final_data_size = nx*ny*total_chips;
    

    // task distribution
    task_alloc(total_chips, numprocs, rank, chip_st, chip_ed);
    if(rank == 0)
    {
        std::cout<<source_name<<" "<<sex_folder<<" "<<total_chips<<std::endl;
    }
    // for(i=0;i<numprocs;i++)
    // {
    //     if(rank == i)
    //     {
    //         std::cout<<rank<<" "<<chip_st<<" "<<chip_ed<<std::endl;
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    

    shear_num = 10;
    for(shear_id=0; shear_id<shear_num; shear_id++)
    {   
        st1 = clock();

        sprintf(sex_data_path, "%s/result/data/%s/sex_%d.hdf5",total_path, sex_folder, shear_id);
        sprintf(set_name,"/data");
        read_h5(sex_data_path, set_name, sex_data);

        if(rank == 0)
        {
            initialize_arr(final_data, final_data_size, 0);

            std::cout<<shear_id<<" "<<sex_data_path<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // loop the chip
        for(k=chip_st;k<chip_ed;k++)
        {      
            st3 = clock();

            get_time(time_now, 50);
            sprintf(log_inform,"%04d. %s",k, time_now);
            write_log(log_path, log_inform);

            sprintf(img_path,"%s/%d/gal_chip_%04d.fits",total_path, shear_id, k);
            read_fits(img_path, img);

            for(i=0;i<ny;i++)
            {
                for(j=0;j<nx;j++)
                {
                    tag = i*nx + j;
                    gal_label = k*nx*ny + tag;
                    segment(img, stamp, tag, size, nx, ny);
                    // \pi * r^2 = pixel number
                    eff_radius_sq = sex_data[gal_label*sex_col + 3]/Pi;           
            
                    if (eff_radius_sq > 0.00001)            
                    {
                        get_quad(stamp, size, eff_radius_sq, gal_quad);
                        final_data[gal_label] = gal_quad;
                    }
                }
            }
            st4 = clock();

            get_time(time_now, 50);
            sprintf(log_inform,"%04d. %s. %.2f",k, time_now, (st4-st3)/CLOCKS_PER_SEC);
            write_log(log_path, log_inform);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        st2 = clock();
        // write the data to disk
        if(rank == 0)
        {
            sprintf(result_path,"%s/result/data/%s/Rfactor_%d.hdf5",total_path,sex_folder,shear_id);
            write_h5(result_path, set_name,final_data,final_data_size,1,true);
            get_time(time_now,50);
            std::cout<<result_path<<std::endl;
            std::cout<<time_now<<" "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    delete[] sex_data;
    delete[] img;
    delete[] stamp;    
    
    MPI_Win_free(&win_final_data);
    MPI_Finalize();

    return 0;
}