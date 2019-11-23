#include<hk_mpi.h>
#include<FQlib.h>

int main(int argc, char **argv)
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	MPI_Status status;

	double t1, t2, t3, s;
	int i, j, num,sub_num;
	num = atoi(argv[1]);
	sub_num = num/numprocs;

	double *send_buf, *re_gather_buf;
	double *recv_buf = new double[sub_num] {};
	int *disp = new int[num];
	int *send_num = new int[numprocs];

	if (rank == 0)
	{			
		std::cout<<"Total: "<<num<<" Each rank: "<<sub_num<<std::endl;
		send_buf = new double[num]();
		for (i = 0; i < numprocs; i++)
		{	
			for(j=0;j<sub_num;j++)
			{
				send_buf[i*sub_num+j] = i+1;
			}
			std::cout<<"Sum of each sub-array: "<<(i+1)*sub_num<<std::endl;
		}
		
	}

	for (i = 0; i < numprocs; i++)
	{
		send_num[i] = sub_num;
 	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	t1 = clock();


	// scatter
	if(rank == 0)
	{
		std::cout<<"Scatter"<<std::endl;
	}
	my_Scatterv(send_buf, send_num, recv_buf, numprocs, rank);
	sum_arr(recv_buf, sub_num, 0, sub_num, s);
	for (i = 0; i < numprocs; i++)
	{
		if (rank == i)
		{		
			std::cout <<rank<<" "<< s << std::endl;
			//show_arr(recv_buf,1,sub_num);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);


	// gather the data scattered before from each threads
	if(rank == numprocs-1)
	{
		re_gather_buf = new double[num]();
	}
	my_Gatherv(recv_buf, send_num, re_gather_buf, numprocs, rank, numprocs-1);

	if(rank == numprocs -1)
	{
		for(i=0;i<numprocs;i++)
		{	
			s=0;
			sum_arr(&re_gather_buf[i*sub_num],sub_num, 0, sub_num, s);
			std::cout<<"Sum of each sub-array: "<<s<<std::endl;
		}
	}


	// send & receive
	if(rank == 0)
	{
		std::cout<<"Sned & Recv"<<std::endl;
	}
	double send_val_d, *rev_vals_d;
	int send_val_i, *rev_vals_i;

	send_val_d = rank+1;
	send_val_i = rank+1;

	if(rank == numprocs-1)
	{
		rev_vals_d = new double[numprocs];
		rev_vals_i = new int[numprocs];
	}
	my_Send_Recv(send_val_d, rev_vals_d,numprocs, rank, numprocs-1);
	my_Send_Recv(send_val_i, rev_vals_i,numprocs, rank, numprocs-1);
	if(rank == numprocs-1)
	{
		for(i=0;i<numprocs;i++)
		{
			std::cout<<"Rev Double: "<<rev_vals_d[i]<<" Real: "<<i+1<<std::endl;
			std::cout<<"Rev Int: "<<rev_vals_i[i]<<" Real: "<<i+1<<std::endl;

		}
	}

	t2 = clock();


	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
	{
		std::cout <<"Finish "<< (t2 - t1) / CLOCKS_PER_SEC<<" Sec"<<std::endl;
	}
	return 0;
}