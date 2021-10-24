#include<hk_science_cal_lib.h>
#include<iostream>


int main(int argc, char *argv[])
{
    int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int i, j, hist_total_len;
    ggl_data_info data_info;
    int foreground_expo_label = -1;
    int thread_live;
    int task_label;
    MPI_Status status;
    MPI_Request request;
    double time_st, time_ed;
    char back_sub_path[100], fore_sub_path[100], result_sub_path[100];
    char set_name[50];
    char result_file_nm[100];

    time_st = clock();

    // initialization
    strcpy(data_info.ggl_total_path, argv[1]);
    strcpy(back_sub_path, argv[2]);
    strcpy(fore_sub_path, argv[3]);

    data_info.jack_num = atoi(argv[4]);
    strcpy(result_file_nm, argv[5]);

    data_info.back_dz = 0.01;//atof(argv[5]);
    data_info.rank = rank;
    data_info.numprocs = numprocs;

    sprintf(data_info.ggl_foreground_inform_path, "%s/%s/foreground_source_list.dat", data_info.ggl_total_path,fore_sub_path);
    sprintf(data_info.ggl_background_inform_path, "%s/%s/background_source_list.dat", data_info.ggl_total_path,back_sub_path);


    // read exposure informs from back/foreground list file
    // and somethings for PDF_SYM
    ggl_initialize(&data_info);

    data_info.dz_hist_bin_num = 400;
    hist_total_len = data_info.dz_hist_bin_num*data_info.sep_bin_num;

    data_info.dz_bin = new MY_FLOAT[data_info.dz_hist_bin_num+1]{};

    data_info.dz_hist = new double[hist_total_len]{};
    double *dz_hist_buff = new double[hist_total_len]{};

    initialize_arr(data_info.dz_hist, hist_total_len, 0);
    initialize_arr(dz_hist_buff, hist_total_len, 0);

    linspace(0, 3, data_info.dz_hist_bin_num+1, data_info.dz_bin); 

    // if(rank == 0)
    // {
    //     show_arr(data_info.dz_bin, 1, data_info.dz_hist_bin_num+1);
    //     std::cout<<std::endl;
    //     show_arr(data_info.dz_hist, 1, data_info.dz_hist_bin_num);
    //     std::cout<<std::endl;
    //     show_arr(dz_hist_buff, 1, data_info.dz_hist_bin_num);
    //     std::cout<<std::endl;
    // }

    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank > 0)
    {
       while(true)
       {    
            foreground_expo_label = -1;
            // then thread 0 will send the task, epxo_pair label, to it.
            MPI_Send(&foreground_expo_label, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
            MPI_Recv(&foreground_expo_label, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

            sprintf(data_info.ggl_log_inform,"%d. receive %d th foreground exposure\n", rank, foreground_expo_label);
            // std::cout<<data_info.ggl_log_inform;

            if(foreground_expo_label > -1)
            {
                ggl_dz_hist(&data_info, foreground_expo_label);
            }
            else
            {
                sprintf(data_info.ggl_log_inform,"%d Break. receive %d th foreground exposure\n", rank, foreground_expo_label);
                std::cout<<data_info.ggl_log_inform;
                break;
            }
       } 
    }
    else
    {   
        // CPU 0 is the master for task distribution
        thread_live = numprocs - 1;
        task_label = 0;

        sprintf(data_info.ggl_log_inform,"Start\n");
        std::cout<<data_info.ggl_log_inform;

        while (thread_live>0)
        {
            MPI_Irecv(&foreground_expo_label, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);

            if(task_label < data_info.len_expo_num)
            {
                foreground_expo_label = task_label;
                task_label ++;
            }
            else
            {
                foreground_expo_label = -1;
                thread_live --;
            }
            MPI_Send(&foreground_expo_label, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD); 
            sprintf(data_info.ggl_log_inform,"send %d th foreground exposure to %d worker %d\n", foreground_expo_label,status.MPI_SOURCE, thread_live);
            // std::cout<<data_info.ggl_log_inform; 
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank > 0)
    {MPI_Send(data_info.dz_hist, hist_total_len, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);}
    else
    {
        for(i=1;i<numprocs;i++)
        {
            MPI_Recv(dz_hist_buff, hist_total_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);

            for(j=0; j<hist_total_len; j++)
            {data_info.dz_hist[j] += dz_hist_buff[j];}
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // sprintf(data_info.ggl_result_path,"%s/result/%s.hdf5", data_info.ggl_total_path, result_file_nm);

    // for(i=0; i< numprocs; i++)
    // {
    //     if(rank == i)
    //     {            
    //         if(i == 0)
    //         {
    //             sprintf(set_name,"/dz_hist_bin");
    //             write_h5(data_info.ggl_result_path, set_name, data_info.dz_bin, 1, data_info.dz_hist_bin_num+1, true);
    //         }
    //         sprintf(set_name,"/dz_hist_%d", i);
    //         write_h5(data_info.ggl_result_path, set_name, data_info.dz_hist, 1, hist_total_len, false);
    //     }

    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    // MPI_Barrier(MPI_COMM_WORLD);
    
    time_ed = clock();

    if(rank == 0)
    {   
        int total_num = 0;
        for(i=0; i<data_info.len_expo_num; i++){total_num += data_info.len_data_row[i];}

        std::cout<<total_num<<"Lens"<<std::endl;


        sprintf(data_info.ggl_result_path,"%s/result/%s.hdf5", data_info.ggl_total_path, result_file_nm);
        sprintf(set_name,"/dz_hist_bin");
        write_h5(data_info.ggl_result_path, set_name, data_info.dz_bin, 1, data_info.dz_hist_bin_num+1, true);
        sprintf(set_name,"/dz_hist");
        write_h5(data_info.ggl_result_path, set_name, data_info.dz_hist, 1, hist_total_len, false);
        sprintf(set_name,"/dz_hist_re");
        write_h5(data_info.ggl_result_path, set_name, data_info.dz_hist, data_info.sep_bin_num, data_info.dz_hist_bin_num, false);


        sprintf(set_name,"/total_len_num");
        write_h5(data_info.ggl_result_path, set_name, &total_num, 1, 1, false);

        sprintf(data_info.ggl_log_inform,"worker %d. Finish in %.2f sec\n", data_info.rank, (time_ed-time_st)/CLOCKS_PER_SEC);
        std::cout<<data_info.ggl_log_inform; 
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
    return 0;
}