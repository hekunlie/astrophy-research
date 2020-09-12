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
    int task_labels[2];
    MPI_Status status;

    if(rank > 0 )
    {   
        task_labels[0] = rank;
        task_labels[1] = rank;

        if(rank == 1){num=1000;}
        else{num=1000000000*rank*rank*rank;}

        for(i=0;i<num;i++){sin(i)*cos(i)*tan(i)*tan(i)*tan(i);}

        MPI_Send(task_labels, 2, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Recv(task_labels, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    if(rank == 0)
    {
        for(i=numprocs-1; i>0; i--)
        {   
            MPI_Recv(task_labels, 2, MPI_INT, i,i, MPI_COMM_WORLD, &status);
            MPI_Send(task_labels, 2, MPI_INT, i, 0, MPI_COMM_WORLD);  
            std::cout<<i<<" "<<task_labels[0]<<" "<<task_labels[1]<<std::endl;  
        }
    }
    MPI_Finalize();
    return 0;
}