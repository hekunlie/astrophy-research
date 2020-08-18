// #include<FQlib.h>
#include<hk_mpi.h>
#include<functions.h>
#include<ctime>
#include<FQlib.h>

#define PRINT_INFO

int main(int argc, char *argv[])
{
    /*
        argv[1]: the parent path that includes "./data", "./result"
        argv[2] ~ argv[i]: the sky area name in the "catalog.hdf5", like "w1", "w2"...
    */

	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    char source_list[300];
    char data_path[300], result_path[300];
    char set_name[50], temp_char[50];

    data_info field_info;

    int i,j,k,m,n;
    int total_expos_num;
    int red_shift_bin_num, z1,z2;

    int ir, radius_bin_num;
    int my_fnm, target_fnm, total_field_num;
    float ra_z1, ra_z2, dec_z1, dec_z2, cos_dec_z1, cos_dec_z2;

    int ic, chi_num;
    float mg1_z1, mg2_z1, mn_z1, mu_z1;
    float mg1_z2, mg2_z2, mn_z2, mu_z2;

    float delta_ra, delta_dec, delta_radius;
    float sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    strcpy(source_list, argv[1]);
    strcpy(data_path,argv[2]);
    z1 = atoi(argv[3]);
    z2 = atoi(argv[4]);

    total_field_num = atoi(argv[5]);



    // read the information of each exposure file
    initialize(source_list, &field_info, total_field_num, numprocs, rank);
    if(rank == 0){std::cout<<"Initialization"<<std::endl;}
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef PRINT_INFO
    if(rank==0)
    {
        for(i=0;i<total_field_num;i++)
        {
            std::cout<<field_info.field_name[i]<<"  "<<field_info.exposure_num_of_field[i]
            <<"  "<<field_info.field_cen_ra[i]<<"  "<<field_info.field_cen_dec[i]
            <<"  "<<field_info.field_delta_ra[i]<<"  "<<field_info.field_delta_dec[i]
            <<"  "<<field_info.field_delta_len[i]<<"  "<<field_info.field_cen_cos_dec[i]<<std::endl;
        }
        for(i=0; i<numprocs; i++){std::cout<<field_info.field_num_each_rank[i]<<" ";}
        std::cout<<std::endl;
    }

    for(j=0; j<numprocs;j++)
    {
        if(j==rank)
        {
            std::cout<<rank<<" "<<field_info.my_field_st<<" "<<field_info.my_field_ed<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif


    if(rank == 0){std::cout<<"Redshifit: "<<z1<<" "<<z2<<std::endl;}

    // read the catalog of redshift bin z1 & z2
    read_field_data(&field_info, z1, z2);
    if(rank == 0)
    {
#ifdef PRINT_INFO
        for(i=0; i<total_field_num; i++)
        {
            std::cout<<"block num: "<<field_info.block_num[i]
            <<"  total gal nun z1: "<<field_info.total_gal_num_z1[i]
            <<"  total gal nun z1: "<<field_info.total_gal_num_z2[i]<<std::endl;
        }
#endif
        std::cout<<"Redshift bin:"<<std::endl;
        show_arr(field_info.zbin, 1, field_info.zbin_num+1);
        std::cout<<std::endl;

        std::cout<<"Radius bin:"<<std::endl;
        show_arr(field_info.theta_bin, 1, field_info.theta_bin_num+1);
        std::cout<<std::endl;

        std::cout<<"Chi guess: "<<field_info.chi_guess_num<<" points"<<std::endl;
        show_arr(field_info.chi_guess,1, field_info.chi_guess_num);

        std::cout<<"Bin num for PDF_SYM: "<<field_info.chi_bin_num<<std::endl;

    }


    // loop the fields
    for(my_fnm=field_info.my_field_st; my_fnm < field_info.my_field_ed; my_fnm++)
    {   ;
        ////////////////  search pairs in the current field  ////////////////////
        // for(i=0; i<field_info.block_num[my_fnm]; i++)
        // {   
        //     m = i*field_info.field_data_col;

        //     ra_z1 = field_info.field_data_z1[my_fnm][m+field_info.ra_idx];
        //     dec_z1 = field_info.field_data_z1[my_fnm][m+field_info.dec_idx];
        //     cos_dec_z1 = field_info.field_data_z1[my_fnm][m+field_info.cos_dec_idx];

        //     mg1_z1 = field_info.field_data_z1[my_fnm][n+field_info.mg1_idx];
        //     mg2_z1 = field_info.field_data_z1[my_fnm][n+field_info.mg1_idx];

        //     mn_z1 = field_info.field_data_z1[my_fnm][n+field_info.mn_idx];
        //     mu_z1 = field_info.field_data_z1[my_fnm][n+field_info.mu_idx];

        //     for(j=i; j<field_info.block_num[my_fnm]; j++)
        //     {   
        //         n = j*field_info.field_data_col;
        //         ra_z2 = field_info.field_data_z2[my_fnm][n+field_info.ra_idx];
        //         dec_z2 = field_info.field_data_z2[my_fnm][n+field_info.dec_idx];
        //         cos_dec_z2 = field_info.field_data_z2[my_fnm][n+field_info.cos_dec_idx];
                
        //         // the seperation angle (arc minute)
        //         delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
        //         delta_dec = dec_z2 - dec_z1;
        //         delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
 
        //         // shear estimators rotation (position angle defined as East of North)
        //         sin_theta = delta_ra/delta_radius;
        //         cos_theta = delta_dec/delta_radius;

        //         sin_2theta = 2*sin_theta*cos_theta;
        //         cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

        //         sin_4theta = 2*sin_2theta*cos_2theta;
        //         cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;

        //         mg1_z2 = field_info.field_data_z2[my_fnm][n+field_info.mg1_idx]*cos_2theta - 
        //             field_info.field_data_z2[my_fnm][n+field_info.mg2_idx]*sin_2theta;
        //         mg2_z2 = field_info.field_data_z2[my_fnm][n+field_info.mg1_idx]*sin_2theta + 
        //             field_info.field_data_z2[my_fnm][n+field_info.mg2_idx]*cos_2theta;

        //         mn_z2 = field_info.field_data_z2[my_fnm][n+field_info.mn_idx];

        //         mu_z2 = field_info.field_data_z2[my_fnm][n+field_info.mu_idx]*cos_4theta - 
        //             field_info.field_data_z2[my_fnm][n+field_info.mv_idx]*sin_4theta;
        //         for(ic=0; ic<chi_num; ic++)
        //         {
                    
        //         }
        //         // loop the radius bin
        //         for(ir=0;ir<radius_bin_num;ir++)
        //         {

        //         }
        //     }
        // }

        // ////////////////// search pairs in the target fields  //////////////////////
        // for(target_fnm=my_fnm; target_fnm<total_field_num; target_fnm++)
        // {
        //     // find the needed fields
        //     // dx^2 + dy^2 - sqrt(dra^2+ddec^2)*2 <= max radius
        //     delta_ra = field_info.field_cen_ra[my_fnm] - field_info.field_cen_ra[target_fnm];
        //     delta_dec = (field_info.field_cen_dec[my_fnm] - field_info.field_cen_dec[target_fnm])*field_info.field_cen_cos_dec[my_fnm];
        //     delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec) - field_info.field_delta_len[target_fnm] - field_info.field_delta_len[my_fnm];
            
        //     if(delta_radius <= field_info.theta_bin[field_info.theta_bin_num+1]);
        //     {
                
        //     }
        // }

    }

    
    MPI_Finalize();
    
    return 0;
}
