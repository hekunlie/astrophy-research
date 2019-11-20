#ifndef __HK_MPI__
#define __HK_MPI__
#pragma once
#include<mpi.h>

void my_Scatterv(const double *send_buf, const int*send_count, const int num_pros, const int rank, double *rec_buf, const int send_rank=0);
/* scatter the parts of a array to threads */
#endif