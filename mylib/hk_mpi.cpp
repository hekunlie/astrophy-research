#include"hk_mpi.h"

void my_Scatterv(const double *send_buf, const int*send_count, const int num_pros, const int rank, double *rec_buf, const int send_rank)
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