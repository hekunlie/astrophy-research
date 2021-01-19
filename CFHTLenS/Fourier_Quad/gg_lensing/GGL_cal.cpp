#include<science_cal_lib.h>
#include<iostream>


int main(int argc, char *argv[])
{
    int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    ggl_data_info data_info;
    int foreground_expo_label = -1;
    int thread_live;
    int task_label;
    MPI_Status status;
    MPI_Request request;

    // initialization
    strcpy(data_info.ggl_total_path, argv[1]);
    data_info.jack_num = atoi(argv[2]);
    data_info.rank = rank;
    data_info.numprocs = numprocs;
    // read exposure informs from back/foreground list file
    // and somethings for PDF_SYM
    ggl_initialize(&data_info);


    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank > 0)
    {
       while(true)
       {    
            foreground_expo_label = -1;
            // then thread 0 will send the task, epxo_pair label, to it.
            MPI_Send(&foreground_expo_label, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
            MPI_Recv(&foreground_expo_label, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

            sprintf(data_info.log_infrom,"receive %d th foreground exposure\n", foreground_expo_label);
            std::cout<<data_info.log_infrom;

            if(foreground_expo_label > -1)
            {
                ggl_find_pair(&data_info, foreground_expo_label);
            }
            else
            {
                sprintf(data_info.log_infrom,"%d Break. receive %d th foreground exposure\n", rank, foreground_expo_label);
                std::cout<<data_info.log_infrom;
                break;
            }
       } 
    }
    else
    {   
        // CPU 0 is the master for task distribution
        thread_live = numprocs - 1;
        task_label = 0;

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
            sprintf(data_info.log_infrom,"send %d th foreground exposure to %d worker %d\n", foreground_expo_label,status.MPI_SOURCE, thread_live);
            std::cout<<data_info.log_infrom; 
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout<<rank<<std::endl;
    ggl_collect_chi(&data_info);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
    return 0;
}