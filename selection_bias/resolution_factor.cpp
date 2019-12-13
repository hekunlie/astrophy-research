#include <FQlib.h>
#include <hk_mpi.h>
#include<hk_iolib.h>

int main(int argc, char**argv)
{
    /* calculate the gaussian weighted quadrupole                                   */
    /* it loops the shear points, each threads will calculate its own files         */
    /* argv[1]: the source name, one of "galsim_dimmer  pts_bright  pts_dimmer"     */
    /* argv[2]: the folder name, "sex2_1.5", "sex4_2"...                            */
    /* argv[3]: total chip number in that folder                                    */

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
    int total_chips;
    char *sex_folder[6];
    char source_name[20];


    strcpy(total_path,argv[1]);
   
    sprintf(log_path,"%s/logs/%d.dat",total_path, rank);
    
    double *img, *stamp;
    double *sex_data[6];
    double *chip_data;
    double *final_data[6];

    int chip_st, chip_ed;
    int i, j, k, m, n, tag, gal_label;
    int nx, ny, size;
    int sex_row, sex_col;

    int shear_num, shear_id;

    double st1, st2, st3, st4;

    double gal_quad;
    double eff_radius_sq;

    int quad_tag;
    char inform[40];

    size = 64;
    nx = 100;
    ny = 100;
    total_chips = 1000;
    sex_row = nx*ny*total_chips;

    img = new double[nx*ny*size*size];
    stamp = new double[size*size];

    for(n=0;n<6;n++){sex_data[n] = new double[sex_row]; sex_folder[n] = new char[40];}
    sprintf(sex_folder[0],"sex2_1.5");
    sprintf(sex_folder[1],"sex2_2");
    sprintf(sex_folder[2],"sex2_4");
    sprintf(sex_folder[3],"sex4_1.5");
    sprintf(sex_folder[4],"sex4_2");
    sprintf(sex_folder[5],"sex4_4");

    MPI_Win win_final_data_1,win_final_data_2,win_final_data_3,win_final_data_4,win_final_data_5,win_final_data_6;
	MPI_Aint final_data_size;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(nx*ny*total_chips* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[0], &win_final_data_1);
		MPI_Win_allocate_shared(nx*ny*total_chips* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[1], &win_final_data_2);
		MPI_Win_allocate_shared(nx*ny*total_chips* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[2], &win_final_data_3);
		MPI_Win_allocate_shared(nx*ny*total_chips* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[3], &win_final_data_4);
		MPI_Win_allocate_shared(nx*ny*total_chips* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[4], &win_final_data_5);
		MPI_Win_allocate_shared(nx*ny*total_chips* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[5], &win_final_data_6);
	}
	else
	{
		int dispu_total;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[0], &win_final_data_1);
		MPI_Win_shared_query(win_final_data_1, 0, &final_data_size, &dispu_total, &final_data[0]);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[1], &win_final_data_2);
		MPI_Win_shared_query(win_final_data_2, 0, &final_data_size, &dispu_total, &final_data[1]);

        MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[2], &win_final_data_3);
		MPI_Win_shared_query(win_final_data_3, 0, &final_data_size, &dispu_total, &final_data[2]);

        MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[3], &win_final_data_4);
		MPI_Win_shared_query(win_final_data_4, 0, &final_data_size, &dispu_total, &final_data[3]);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[4], &win_final_data_5);
		MPI_Win_shared_query(win_final_data_5, 0, &final_data_size, &dispu_total, &final_data[4]);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &final_data[5], &win_final_data_6);
		MPI_Win_shared_query(win_final_data_6, 0, &final_data_size, &dispu_total, &final_data[5]);
	}

    final_data_size = nx*ny*total_chips;
    

    // task distribution
    task_alloc(total_chips, numprocs, rank, chip_st, chip_ed);
    if(rank == 0)
    {   
        for(n=0;n<6;n++){std::cout<<sex_folder[n]<<" ";}
        std::cout<<total_chips<<std::endl;
    }
    for(i=0;i<numprocs;i++)
    {
        if(rank == i)
        {
            std::cout<<rank<<" "<<chip_st<<" "<<chip_ed<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    

    shear_num = 10;
    for(shear_id=0; shear_id<shear_num; shear_id++)
    {   
        st1 = clock();
        for(n=0;n<6;n++)
        {
            sprintf(sex_data_path, "%s/result/data/%s/area_%d.hdf5", total_path, sex_folder[n], shear_id);
            sprintf(set_name,"/data");
            read_h5(sex_data_path, set_name, sex_data[n]);
            if(rank==0){std::cout<<shear_id<<" "<<sex_data_path<<std::endl;}
        }

        if(rank == 0)
        {
            for(n=0;n<6;n++){initialize_arr(final_data[n], final_data_size, 0);}
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // loop the chip
        for(k=chip_st;k<chip_ed;k++)
        {      
            st3 = clock();

            get_time(time_now, 50);
            sprintf(log_inform,"%04d. %s",k, time_now);
            write_log(log_path, log_inform);
            if(rank==0){std::cout<<shear_id<<" "<<log_inform<<std::endl;}

            sprintf(img_path,"%s/%d/gal_chip_%04d.fits",total_path, shear_id, k);
            read_fits(img_path, img);

            for(i=0;i<ny;i++)
            {
                for(j=0;j<nx;j++)
                {
                    tag = i*nx + j;
                    gal_label = k*nx*ny + tag;
                    segment(img, stamp, tag, size, nx, ny);

                    for(n =0; n<6; n++)
                    {   // \pi * r^2 = pixel number
                        eff_radius_sq = sex_data[n][gal_label]/Pi;           
                
                        if (eff_radius_sq > 0.000001)            
                        {
                            get_quad(stamp, size, eff_radius_sq, gal_quad);
                            final_data[n][gal_label] = gal_quad;
                        }
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
            for(n=0;n<6;n++)
            {   
                sprintf(result_path,"%s/result/data/%s/rfactor_%d.hdf5",total_path,sex_folder[n],shear_id);
                write_h5(result_path, set_name,final_data[n],final_data_size,1,true);
                get_time(time_now,50);
                std::cout<<result_path<<std::endl;
                std::cout<<time_now<<" "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    for(n=0;n<6;n++){delete[] sex_data[n];}
    delete[] img;
    delete[] stamp;    
    
    MPI_Win_free(&win_final_data_1);
    MPI_Win_free(&win_final_data_2);
    MPI_Win_free(&win_final_data_3);
    MPI_Win_free(&win_final_data_4);
    MPI_Win_free(&win_final_data_5);
    MPI_Win_free(&win_final_data_6);
    MPI_Finalize();

    return 0;
}