#include<hk_mpi.h>
#include<functions_expo_wise.h>
#include<ctime>


#define PRINT_INFO

int main(int argc, char *argv[])
{
    /*

    */

	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    char source_list[300];
    char data_path[300], result_path[300];
    char set_name[50], temp_char[50];
    char log_inform[500], log_path[600];

    double st1, st2, st3, st4, st5, tt;

    data_info expo_info;

    MY_FLOAT *ggcor_1, *ggcor_2, *gg_read;

    int i,j;
    int fnm_1, fnm_2, total_expo_num, label, num_label;


    strcpy(source_list, argv[1]);
    strcpy(expo_info.parent_path, argv[2]);
    total_expo_num = atoi(argv[3]);

    sprintf(log_path, "%s/log/%d_log.dat",expo_info.parent_path, rank);
    

    // read the information of each exposure file
    initialize(source_list, &expo_info, total_expo_num);

    // read the catalog of redshift bin z1 & z2
    read_data(&expo_info);

    // find all the potential expo pair for calculation (i, j),
    // does not include the expo itself i!= j
    task_prepare(numprocs, rank, &expo_info);

    

    ////////////////////////////////  PRINT INFO  ////////////////////////////////////////////
    if(rank == 0){std::cout<<"Initialization"<<std::endl;}
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0)
    {
#ifdef PRINT_INFO
        for(i=0;i<total_expo_num;i++)
        {
            std::cout<<expo_info.expo_name_path[i]<<"  "<<expo_info.expo_name[i]<<" "<<expo_info.expo_gal_num[i]
            <<"  "<<expo_info.expo_cen_ra[i]<<"  "<<expo_info.expo_cen_dec[i]
            <<"  "<<expo_info.expo_delta_ra[i]<<"  "<<expo_info.expo_delta_dec[i]
            <<"  "<<expo_info.expo_delta_len[i]<<"  "<<expo_info.expo_cen_cos_dec[i]<<std::endl;
        }
#endif
        std::cout<<std::endl;
        std::cout<<"Exposure pairs for each cpu:"<<std::endl;
        show_arr(expo_info.expo_pair_num_each_rank,1,numprocs);
        std::cout<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    sprintf(log_inform, "Rank %d. Total field pairs: %d. My field pairs: %d ~ %d.", rank, expo_info.expo_pair_num,
    expo_info.my_expo_pair_st,expo_info.my_expo_pair_ed);
    write_log(log_path, log_inform);
    for(j=0; j<numprocs;j++)
    {
        if(j==rank)
        {
            std::cout<<log_inform<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if(rank == 0)
    {
        std::cout<<"Redshift bin:"<<std::endl;
        show_arr(expo_info.zbin, 1, expo_info.zbin_num+1);
        std::cout<<log_inform<<std::endl;
        
        std::cout<<"Radius bin:"<<std::endl;
        show_arr(expo_info.theta_bin, 1, expo_info.theta_bin_num+1);
        std::cout<<std::endl;


        std::cout<<"Chi guess: "<<expo_info.chi_guess_num<<" points"<<std::endl;
        show_arr(expo_info.chi_guess, 1, expo_info.chi_guess_num);
        std::cout<<std::endl;
        std::cout<<"Bin num for PDF_SYM: "<<expo_info.mg_bin_num<<" "<<expo_info.mg_bin_num1<<" "<<expo_info.mg_bin_num2<<" "<<expo_info.mg_bin_num3<<std::endl;
        show_arr(expo_info.mg_bin, 1, expo_info.mg_bin_num+1);

        std::cout<<std::endl<<expo_info.mg_bin_num<<" "<<expo_info.chi_block_len<<" "<<expo_info.ir_chi_block_len<<" "<<
        expo_info.expo_chi_block_len<<std::endl<<std::endl;

        std::cout<<std::endl<<expo_info.parent_path<<std::endl;
        std::cout<<std::endl<<expo_info.gg_len<<std::endl<<std::endl;     
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ////////////////////////////////  PRINT INFO-end  ////////////////////////////////////////////



    initialize_total_chi_block(&expo_info);
    // initialize_expo_chi_block(&expo_info,my_fnm);

    // // loop the field pairs
    // st1 = clock();
    // for(i=expo_info.my_expo_pair_st; i < expo_info.my_expo_pair_ed; i++)
    // {   
    //     // field pair label
    //     fnm_1 = expo_info.expo_pair_label_1[i];
    //     fnm_2 = expo_info.expo_pair_label_2[i];

        
    //     ////////////////  search pairs in the current field  ////////////////////
    //     st2 = clock();
    //     sprintf(log_inform,"Start %d (%d ~ %d) field pair: %d-%s(%d) <-> %d-%s(%d)", i, expo_info.my_expo_pair_st,expo_info.my_expo_pair_ed, fnm_1, 
    //     expo_info.expo_name[fnm_1], expo_info.total_gal_num_z1[fnm_1], fnm_2, expo_info.expo_name[fnm_2],expo_info.total_gal_num_z2[fnm_2]);
    //     if(rank == 0){std::cout<<log_inform<<std::endl;}
    //     write_log(log_path, log_inform);
        
    //     if (fnm_1 == fnm_2)
    //     {
    //         find_pairs_same_field(&expo_info,fnm_1);
    //     }
    //     else
    //     {   
    //         find_pairs_diff_field(&expo_info, fnm_1, fnm_2);
    //     }

    //     // find_pairs(&expo_info, fnm_1, fnm_2);

    //     st3 = clock();
    //     tt =  (st3 - st2)/CLOCKS_PER_SEC;

    //     sprintf(log_inform,"Finish %d (%d ~ %d) field pair: %d-%s(%d) <-> %d-%s(%d) in %.2f sec.",i,expo_info.my_expo_pair_st,expo_info.my_expo_pair_ed,fnm_1, 
    //     expo_info.expo_name[fnm_1], expo_info.total_gal_num_z1[fnm_1], fnm_2, expo_info.expo_name[fnm_2],expo_info.total_gal_num_z2[fnm_2], tt);
    //     if(rank == 0)
    //     {
    //         std::cout<<log_inform<<std::endl;
    //         std::cout<<"========================================================================================="<<std::endl;
    //     }
    //     write_log(log_path, log_inform);
    //     if(num_label == 1)
    //     {
    //         ;// collect_chi_block(&expo_info, my_fnm);
    //         // save_expo_chi_block(&expo_info, my_fnm);
    //     }        
        
    // }
    // st4 = clock();
    // tt =  (st4 - st1)/CLOCKS_PER_SEC;
    // sprintf(log_inform,"All field pair finished in %.2f sec. (%d ~ %d)", tt,expo_info.my_expo_pair_st,expo_info.my_expo_pair_ed);
    // if(rank == 0)
    // {
    //     std::cout<<log_inform<<std::endl;
    //     std::cout<<"========================================================================================="<<std::endl;
    // }
    // write_log(log_path, log_inform);
    MPI_Finalize();
    
    return 0;
}
