#include<mpi.h>
#include<iostream>
#include<cmath>

int main(int argc, char *argv[])
{
    int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int i, num;
    int task_labels_s[2],task_labels_r[2];
    MPI_Status status;
    MPI_Request request;
    int *labels = new int[numprocs-1]{};
    int sum = 0;
    if(rank > 0 )
    {   
        task_labels_s[0] = rank;
        task_labels_s[1] = rank;

        if(rank == 1){num=1000;}
        else{num=10000000000*rank*rank*rank;}

        for(i=0;i<num;i++){sin(i)*cos(i)*tan(i)*tan(i)*tan(i);}

        // MPI_Isend(task_labels, 2, MPI_INT, 0, rank, MPI_COMM_WORLD, &request);
        
        MPI_Send(task_labels_s, 2, MPI_INT, 0, rank, MPI_COMM_WORLD);

        MPI_Recv(task_labels_r, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        // std::cout<<rank<<" I sent "<<task_labels_s[0]<<" "<<task_labels_s[1]<<" and I receive "<<task_labels_r[0]<<" "<<task_labels_r[1]<<std::endl;
    }
    if(rank == 0)
    {   
        while(sum < numprocs-1)
        {


            MPI_Irecv(task_labels_r, 2, MPI_INT, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);

            // std::cout<<MPI_ANY_SOURCE<<" "<<MPI_ANY_TAG<<" "<<request<<" "<<status.MPI_SOURCE<<" "<<status.MPI_TAG<<std::endl;
            std::cout<<"Receive from "<<status.MPI_SOURCE<<": "<<task_labels_r[0]<<" "<<task_labels_r[1]<<std::endl; 
            task_labels_s[0] = status.MPI_SOURCE;
            task_labels_s[1] = status.MPI_SOURCE;
            MPI_Send(task_labels_s, 2, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
            labels[status.MPI_SOURCE-1] = 1;
            sum =0;
            for(i=0;i<numprocs-1;i++){sum+=labels[i];}
            // for(i=0;i<numprocs-1; i++){std::cout<<labels[i]<<" ";}
            
            // std::cout<<sum<<std::endl;
        }
        std::cout<<"Finish"<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}