#include<hk_mpi.h>
#include<hk_correlation_functions_expo_wise.h>
#include<ctime>


#define NOT_PRINT_INFO

int main(int argc, char *argv[])
{

	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    char cata_sub_path[80],result_sub_path[80];
    char log_inform[500], log_path[600];

    data_info expo_info;
    double st1, st2, st3, st4, st5, st6, tt;
    int i,j;
    int pair_num, label;
    double count_sum;

    int task_end = 0;
    int thread_live;
    int task_labels[2];
    MPI_Status status;
    MPI_Request request;

    
    strcpy(expo_info.parent_path, argv[1]);
    strcpy(cata_sub_path, argv[2]);
    strcpy(result_sub_path, argv[3]);

    expo_info.my_rank = rank;

    sprintf(expo_info.cata_path,"%s/%s", expo_info.parent_path, cata_sub_path);
    sprintf(expo_info.result_path,"%s/%s", expo_info.parent_path, result_sub_path);

    sprintf(log_path, "%s/log/%d_log.dat",expo_info.parent_path, rank);
    
    // read the information of each exposure file
    initialize(&expo_info);

    // find all the potential expo pair for calculation, 
    // (i, j), i!= j, does not include the expo itself 
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
        std::cout<<"Totally "<<expo_info.total_expo_num<<" exposures, "<<expo_info.task_expo_num<<" exposure pairs"<<std::endl;
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

        std::cout<<std::endl<<"G bins: "<<expo_info.mg_bin_num<<". Minimum Chi block len: "<<expo_info.chi_block_len<<std::endl;
        std::cout<<"Chi block len of one point: "<<expo_info.ir_chi_block_len<<std::endl;
        std::cout<<"Chi block len of each Z bin: "<<expo_info.iz_chi_block_len<<std::endl;
        std::cout<<"Total Chi block len: "<<expo_info.expo_chi_block_len<<std::endl<<std::endl;

        std::cout<<"Buffer size: "<<expo_info.max_buffer_size<<" elements. Actual size: "<<expo_info.actual_buffer_size<<" elements"<<std::endl;
        std::cout<<"Max block num in buffer: "<<expo_info.max_block_in_buffer <<std::endl;
        std::cout<<"Block size: "<<expo_info.block_size_in_buffer<<" elements"<<std::endl;
        std::cout<<"Now "<<expo_info.total_buffer_num<<" buffers & "<<expo_info.block_count<<" blocks"<<std::endl;

        std::cout<<std::endl<<expo_info.parent_path<<std::endl;
        std::cout<<std::endl<<expo_info.cata_path<<std::endl;
        std::cout<<std::endl<<expo_info.result_path<<std::endl;

        std::cout<<std::endl<<expo_info.gg_len<<std::endl<<std::endl;     
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ////////////////////////////////  PRINT INFO-end  ////////////////////////////////////////////



    ////////////////////////////////  Start  ////////////////////////////////////////////
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

                if(rank == 1){std::cout<<log_inform<<std::endl;}
                write_log(log_path, log_inform);

                initialize_expo_chi_block(&expo_info);

                read_expo_data_1(&expo_info, task_labels[0]);
                read_expo_data_2(&expo_info, task_labels[1]);
                
                //////////////  search pairs ////////////////////
                find_pairs_diff_expo_dev(&expo_info, task_labels[0], task_labels[1]);
                //find_pairs_same_expo(&expo_info, task_labels[0], task_labels[1]);
                //find_pairs_stack_expo(&expo_info, task_labels[0], task_labels[1]);


                // if more 1 pair has been found, write into the result file
                // if(expo_info.gg_pairs > 1){save_expo_data(&expo_info, task_labels[0], task_labels[1], rank);}

                if(expo_info.gg_pairs > 1000)
                {
                    expo_info.expo_pair_label_1.push_back(task_labels[0]);
                    expo_info.expo_pair_label_2.push_back(task_labels[1]);
                    save_expo_data_new(&expo_info,  rank, 0);
                }
                

                st3 = clock();
                tt =  (st3 - st2)/CLOCKS_PER_SEC;
                pair_num = expo_info.expo_pair_label_1.size();
                sprintf(log_inform,"Finish in %.2f sec. %g pairs. Expo pairs got now: %d\nNow %d buffers, %d block in current buffer", 
                        tt, expo_info.gg_pairs,pair_num, expo_info.total_buffer_num,expo_info.block_count);
                if(rank == 1)
                {                       
                    std::cout<<log_inform<<std::endl<<std::endl;
                }
                write_log(log_path, log_inform);
            }          
            else
            {   
                save_expo_data_new(&expo_info,  rank, 1);
                
                sprintf(log_inform,"End, Write & Exit. Expo pairs got now: %d", pair_num);
                if(rank == 1)
                {                       
                    std::cout<<std::endl<<log_inform<<std::endl;
                }
                write_log(log_path, log_inform);

                break;
            }
        }
        save_expo_pair_label(&expo_info, rank);

        st4 = clock();
        tt =  (st4 - st1)/CLOCKS_PER_SEC;
        sprintf(log_inform,"CPU %d. Finish in %.2f sec.", rank, tt);
        write_log(log_path, log_inform);
        // if(rank == 1)
        // {   
            std::cout<<log_inform<<std::endl;
        // }
    }
    else
    {   
        // CPU 0 is the master for task distribution
        thread_live = numprocs - 1;

        while(thread_live > 0)
        {    
            MPI_Irecv(task_labels, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);

            if(task_end<expo_info.task_expo_num)
            {
                task_labels[0] = expo_info.task_expo_pair_labels_1[task_end];
                task_labels[1] = expo_info.task_expo_pair_labels_2[task_end];
                task_end++;

                sprintf(log_inform,"Send %d-%s(%d) <-> %d-%s(%d) to CPU %d. %d/%d", task_labels[0], expo_info.expo_name[task_labels[0]], expo_info.expo_gal_num[task_labels[0]], 
                        task_labels[1], expo_info.expo_name[task_labels[1]],expo_info.expo_gal_num[task_labels[1]], status.MPI_SOURCE,task_end, expo_info.task_expo_num);
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
        sprintf(log_inform,"CPU %d: Finish in %.2f sec.", rank, tt);
        write_log(log_path, log_inform);
        // if(rank == 0)
        // {   
            std::cout<<log_inform<<std::endl;
        // }      
    }
        
    ////////////////////////////////  End  ////////////////////////////////////////////

    MPI_Finalize();
    
    return 0;
}
