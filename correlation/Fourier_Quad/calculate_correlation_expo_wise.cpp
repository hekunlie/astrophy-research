#include<hk_mpi.h>
#include<functions_expo_wise.h>
#include<ctime>


#define NOT_PRINT_INFO

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

    if(rank == 0)
    {
        std::cout<<"Redshift bin:"<<std::endl;
        show_arr(expo_info.zbin, 1, expo_info.zbin_num+1);
                
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



    // exit(0);

    // loop the expo pairs
    st1 = clock();
    for(i=0; i < expo_info.task_expo_num; i++)
    {   
        // expo pair label
        fnm_1 = expo_info.task_expo_label[i];

        initialize_expo_chi_block(&expo_info);

        for(fnm_2=fnm_1+1; fnm_2<expo_info.total_expo_num;fnm_2++)
        {    
            expo_distance(&expo_info,fnm_1, fnm_2, label);
            if(label == 1)
            {
                st2 = clock();
                sprintf(log_inform,"Start %d/%d expo pair: %d-%s(%d) <-> %d-%s(%d)", i+1, expo_info.task_expo_num, fnm_1, 
                expo_info.expo_name[fnm_1], expo_info.expo_gal_num[fnm_1], fnm_2, expo_info.expo_name[fnm_2],expo_info.expo_gal_num[fnm_2]);
                if(rank == 0){std::cout<<log_inform<<std::endl;}
                write_log(log_path, log_inform);

                ////////////////  search pairs ////////////////////
                find_pairs(&expo_info, fnm_1, fnm_2);

                st3 = clock();
                tt =  (st3 - st2)/CLOCKS_PER_SEC;
                sprintf(log_inform,"Finish in %.2f sec.", tt);
                
                if(rank == 0)
                {
                    std::cout<<log_inform<<std::endl;
                    std::cout<<"========================================================================================="<<std::endl;
                }
                write_log(log_path, log_inform);
            }
        
        }
        save_expo_chi_block(&expo_info, fnm_1);     
        
    }
    st4 = clock();
    tt =  (st4 - st1)/CLOCKS_PER_SEC;
    sprintf(log_inform,"All expo pairs finished in %.2f sec. (%d ~ %d)", tt,expo_info.my_expo_pair_st,expo_info.my_expo_pair_ed);
    if(rank == 0)
    {
        std::cout<<log_inform<<std::endl;
        std::cout<<"========================================================================================="<<std::endl;
    }
    write_log(log_path, log_inform);
    MPI_Finalize();
    
    return 0;
}
