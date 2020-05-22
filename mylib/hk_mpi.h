#ifndef __HK_MPI__
#define __HK_MPI__
#pragma once
#include<mpi.h>

void my_Scatterv(const double *send_buf, const int*send_count, double *rec_buf, const int num_pros, const int rank,const int send_rank=0);//checked 2019-11-21
/* scatter the parts of a array to threads */
void my_Scatterv_all(const double *send_buf, const int total_len, double *rec_buf, const int num_pros, const int send_rank=0);
/* scatter the whole array to each thread, each thread gets a copy of the 'send_buf' */

void my_Gatherv(const double *send_buf, const int*send_count, double *rec_buf, const int num_pros, const int rank,const int rev_rank=0);//checked 2019-11-21
void my_Gatherv(const float *send_buf, const int*send_count, float *rec_buf, const int num_pros, const int rank,const int rev_rank=0);
/* gather arrays into a large one */

void my_Send_all(const double *send_buf, double *recv_buf, const int buf_len, const int num_pros, const int send_rank, const int my_rank);
/* the send_rank send the "send_buf" to each thread */

void my_Send_Recv(const double send_val, double *rev_buf, const int num_pros, const int rank, const int rev_rank=0);//checked 2019-11-21
void my_Send_Recv(const int send_val, int *rev_buf, const int num_pros, const int rank, const int rev_rank=0);//checked 2019-11-21
/* send data to rev_rank */


#endif
