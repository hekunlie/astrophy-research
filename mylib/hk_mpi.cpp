#include"hk_mpi.h"

void my_Scatterv(const double *send_buf, const int*send_count, double *rec_buf, const int num_pros, const int rank,const int send_rank)
{
    int *displ = new int[num_pros]();
    int i,j,k;

    for(i=0; i<num_pros; i++)
    {
        k=0;  
        for(j=0; j<i; j++)
        {
            k+=send_count[j];
        }
        displ[i] = k;
    }
    MPI_Scatterv(send_buf, send_count, displ, MPI_DOUBLE, rec_buf, send_count[rank], MPI_DOUBLE, send_rank, MPI_COMM_WORLD);

    delete[] displ;
}

void my_Gatherv(const double *send_buf, const int*send_count, double *rec_buf, const int num_pros, const int rank,const int rev_rank)
{
    int *displ = new int[num_pros]();
    int i,j,k;

    for(i=0; i<num_pros; i++)
    {
        k=0;  
        for(j=0; j<i; j++)
        {
            k+=send_count[j];
        }
        displ[i] = k;
    }
    MPI_Gatherv(send_buf, send_count[rank], MPI_DOUBLE, rec_buf, send_count, displ, MPI_DOUBLE, rev_rank, MPI_COMM_WORLD);

    delete[] displ;
}

void my_Send_Recv(const double send_val, double *rev_buf, const int num_pros, const int rank, const int rev_rank)
{
    int i;
    MPI_Status status;

    if (rank != rev_rank)
    {
        MPI_Send(&send_val,1, MPI_DOUBLE, rev_rank, rank, MPI_COMM_WORLD);
    }
    else
    {
        rev_buf[rank] = send_val;
        for(i=0;i<num_pros;i++)
        {
            if(i != rev_rank)
            {
                MPI_Recv(&rev_buf[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            }
        }
    }
}

void my_Send_Recv(const int send_val, int *rev_buf, const int num_pros, const int rank, const int rev_rank)
{
    int i;
    MPI_Status status;

    if (rank != rev_rank)
    {
        MPI_Send(&send_val,1, MPI_INT, rev_rank, rank, MPI_COMM_WORLD);
    }
    else
    {
        rev_buf[rank] = send_val;
        for(i=0;i<num_pros;i++)
        {
            if(i != rev_rank)
            {
                MPI_Recv(&rev_buf[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            }
        }
    }
}