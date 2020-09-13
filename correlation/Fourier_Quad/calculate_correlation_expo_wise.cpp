#include<hk_mpi.h>
#include<functions_expo_wise.h>
#include<ctime>


#define NOT_PRINT_INFO

int main(int argc, char *argv[])
{

	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    char source_list[300];
    char data_path[300], result_path[300];
    char set_name[50], temp_char[50];
    char log_inform[500], log_path[600];

    data_info expo_info;
    double st1, st2, st3, st4, st5, st6, tt;
    int i,j;
    int fnm_1, fnm_2, total_expo_num, label;
    double count_sum;


    strcpy(expo_info.parent_path, argv[1]);
    total_expo_num = atoi(argv[2]);

    sprintf(log_path, "%s/log/%d_log.dat",expo_info.parent_path, rank);
    
    // strcpy(result_path, argv[4]);

    // read the information of each exposure file
    initialize(&expo_info, total_expo_num);

    // // read the catalog of redshift bin z1 & z2
    // read_data(&expo_info);

    // find all the potential expo pair for calculation (i, j), i!= j
    // does not include the expo itself 
    task_prepare(numprocs, rank, &expo_info);

    if(rank == 0){initialize_thread_pool(&expo_info, numprocs);}

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
        // std::cout<<"Exposure pairs for each cpu:"<<std::endl;
        // show_arr(expo_info.expo_pair_num_each_rank,1,numprocs);
        std::cout<<"Total expo pairs: "<<expo_info.task_expo_num<<std::endl;
        std::cout<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        std::cout<<"Redshift bin: "<<expo_info.zbin_num<<std::endl;
        show_arr(expo_info.zbin, 1, expo_info.zbin_num+1);
                
        std::cout<<"Theta bin: "<<expo_info.theta_bin_num<<std::endl;
        show_arr(expo_info.theta_bin, 1, expo_info.theta_bin_num+1);
        std::cout<<std::endl;


        std::cout<<"Chi guess: "<<expo_info.chi_guess_num<<" points"<<std::endl;
        show_arr(expo_info.chi_guess, 1, expo_info.chi_guess_num);
        std::cout<<std::endl;
        std::cout<<"Bin num for PDF_SYM: "<<expo_info.mg_bin_num<<" "<<expo_info.mg_bin_num1<<" "<<expo_info.mg_bin_num2<<" "<<expo_info.mg_bin_num3<<std::endl;
        show_arr(expo_info.mg_bin, 1, expo_info.mg_bin_num+1);

        std::cout<<std::endl<<expo_info.mg_bin_num<<" "<<expo_info.chi_block_len<<" "
        <<expo_info.ir_chi_block_len<<" "<<expo_info.iz_chi_block_len<<" "<<
        expo_info.expo_chi_block_len<<std::endl<<std::endl;

        std::cout<<std::endl<<expo_info.parent_path<<std::endl;
        std::cout<<std::endl<<expo_info.gg_len<<std::endl<<std::endl;     
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ////////////////////////////////  PRINT INFO-end  ////////////////////////////////////////////


    ////////////////////////////////  Start  ////////////////////////////////////////////
    int task_end = 0;
    int thread_live = numprocs - 1;
    int task_labels[2];
    MPI_Status status;
    MPI_Request request;

    st1 = clock();

    if(rank>0)
    {   
        while(true)
        {   
            
            // then thread 0 will send the task, epxo_pair label, to it.
            MPI_Send(task_labels, 2, MPI_INT, 0, rank, MPI_COMM_WORLD);
            MPI_Recv(task_labels, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            // if it receives [-1,-1], no task left, it will break the loop.
            if(task_labels[0] >-1)
            {   
                st2 = clock();

                sprintf(log_inform,"expo pair: %d-%s(%d) <-> %d-%s(%d)", task_labels[0], expo_info.expo_name[task_labels[0]], expo_info.expo_gal_num[task_labels[0]], 
                        task_labels[1], expo_info.expo_name[task_labels[1]],expo_info.expo_gal_num[task_labels[1]]);

                if(rank == 0){std::cout<<log_inform<<std::endl;}
                write_log(log_path, log_inform);


                initialize_expo_chi_block(&expo_info);

                read_expo_data_1(&expo_info, task_labels[0]);
                read_expo_data_2(&expo_info, task_labels[1]);
                //////////////  search pairs ////////////////////
                find_pairs_new(&expo_info, task_labels[0], task_labels[1]);

                // if more 1 pair has been found, write into the result file
                if(expo_info.gg_pairs > 1){save_expo_data(&expo_info, task_labels[0], task_labels[1], rank);}
                

                st3 = clock();
                tt =  (st3 - st2)/CLOCKS_PER_SEC;
                sprintf(log_inform,"Finish in %.2f sec. %g pairs. Expo pairs got now: ", tt, expo_info.gg_pairs,expo_info.expo_pair_label_1.size());
                if(rank == 1)
                {                       
                    std::cout<<std::endl<<log_inform<<std::endl;
                }
                write_log(log_path, log_inform);
            }          
            else{break;}
        }
        
        st4 = clock();
        tt =  (st4 - st1)/CLOCKS_PER_SEC;
        sprintf(log_inform,"Finish in %.2f sec.", tt);
        write_log(log_path, log_inform);
        if(rank == 1)
        {   
            std::cout<<log_inform<<std::endl;
        }
    }
    else
    {   
        while(thread_live > 0)
        {    
            MPI_Irecv(task_labels, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);

            if(task_end<expo_info.task_expo_num)
            {
                task_labels[0] = expo_info.task_expo_pair_labels_1[task_end];
                task_labels[1] = expo_info.task_expo_pair_labels_1[task_end];
                task_end++;

                sprintf(log_inform,"Send %d-%s(%d) <-> %d-%s(%d) to CPU %d", task_labels[0], expo_info.expo_name[task_labels[0]], expo_info.expo_gal_num[task_labels[0]], 
                        task_labels[1], expo_info.expo_name[task_labels[1]],expo_info.expo_gal_num[task_labels[1]], status.MPI_SOURCE);
                write_log(log_path, log_inform);
            }
            else
            {
                task_labels[0] = -1;
                task_labels[1] = -1;
                thread_live--;

                sprintf(log_inform,"No task left for CPU %d", status.MPI_SOURCE);
                write_log(log_path, log_inform);
            }                
            MPI_Send(task_labels, 2, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);    
        }

        st5 = clock(); 
        tt =  (st5 - st1)/CLOCKS_PER_SEC;
        sprintf(log_inform,"Finish in %.2f sec.", tt);
        write_log(log_path, log_inform);
        if(rank == 0)
        {   
            std::cout<<log_inform<<std::endl;
        }      
    }
        
    ////////////////////////////////  End  ////////////////////////////////////////////



    // ////////////////////////////////// loop the expo pairs  ////////////////////////////////
    // st1 = clock();
    // for(i=0; i < expo_info.task_expo_num; i++)
    // {   
    //     st2 = clock();
    //     // expo pair label
    //     fnm_1 = expo_info.task_expo_label[i];

    //     read_expo_data_1(&expo_info,fnm_1);

        
    //     sprintf(log_inform,"Start %d/%d. expo  %d-%s(%d)", i+1, expo_info.task_expo_num, 
    //             fnm_1, expo_info.expo_name[fnm_1], expo_info.expo_gal_num[fnm_1]);
    //     if(rank == 0){std::cout<<log_inform<<std::endl;}
    //     write_log(log_path, log_inform);

    //     // only search pairs from the exposures have a bigger exposure label
    //     // to avoid double counting
    //     for(fnm_2=fnm_1+1; fnm_2<expo_info.total_expo_num;fnm_2++)
    //     {    
    //         expo_distance(&expo_info, fnm_1, fnm_2, label);

    //         if(label == 1)
    //         {   
    //             st3 = clock();
                
    //             sprintf(log_inform,"expo pair: %d-%s(%d) <-> %d-%s(%d)", fnm_1, expo_info.expo_name[fnm_1], expo_info.expo_gal_num[fnm_1], 
    //                     fnm_2, expo_info.expo_name[fnm_2],expo_info.expo_gal_num[fnm_2]);
    //             if(rank == 0){std::cout<<log_inform<<std::endl;}
    //             write_log(log_path, log_inform);

    //             initialize_expo_chi_block(&expo_info);

    //             read_expo_data_2(&expo_info, fnm_2);
    //             //////////////  search pairs ////////////////////
    //             find_pairs_new(&expo_info, fnm_1, fnm_2);

    //             if(expo_info.gg_pairs > 1){save_expo_data(&expo_info, fnm_1, fnm_2, rank);}

    //             st4 = clock();
    //             tt =  (st4 - st3)/CLOCKS_PER_SEC;
    //             sprintf(log_inform,"Finish in %.2f sec. %g pairs. Expo pairs now: ", tt, expo_info.gg_pairs,expo_info.expo_pair_label_1.size());
    //             if(rank == 0)
    //             {                       
    //                 std::cout<<std::endl<<log_inform<<std::endl;
    //             }
    //             write_log(log_path, log_inform);
                
    //         }
    //         if(fnm_2 == fnm_1+5){exit(0);}
    //     }

    //     st5 = clock();
    //     tt =  (st5 - st2)/CLOCKS_PER_SEC;
    //     sprintf(log_inform,"Finish %d/%d. expo %d-%s(%d) in %.2f sec.", i+1, expo_info.task_expo_num,fnm_1,
    //             expo_info.expo_name[fnm_1], expo_info.expo_gal_num[fnm_1], tt);
    //     write_log(log_path, log_inform);
    //     if(rank == 0)
    //     {   
    //         std::cout<<log_inform<<std::endl;
    //         std::cout<<"========================================================================================="<<std::endl<<std::endl;
    //     }
    // }

    // st6 = clock();
    // tt =  (st6 - st1)/CLOCKS_PER_SEC;
    // sprintf(log_inform,"All expo pairs finished in %.2f sec.", tt);
    // if(rank == 0)
    // {
    //     std::cout<<log_inform<<std::endl;
    //     std::cout<<"========================================================================================="<<std::endl;
    // }
    // write_log(log_path, log_inform);
    // //////////////////////////////// loop the expo pairs-end ////////////////////////////////


    MPI_Finalize();
    
    return 0;
}
