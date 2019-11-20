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

	double *send_buf;
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

	// scatter
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = clock();
	//MPI_Scatterv(sendbuf, send_num, disp, MPI_DOUBLE, recvbuf, num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	my_Scatterv(send_buf, send_num, numprocs, rank, recv_buf);
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

	// send & receive
	t2 = clock();
	// if (rank == 0)
	// {
	// 	for (i = 1; i < numprocs; i++)
	// 	{
	// 		MPI_Send(sendbuf, num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
	// 	}
	// }
	// else
	// {
	// 	MPI_Recv(sendbuf, num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
	// }

	// for (i = 0; i < numprocs; i++)
	// {
	// 	if (rank == i)
	// 	{
	// 		sum_arr(sendbuf, num, 0, num, s);
	// 		std::cout << s << std::endl;
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);

	t3 = clock();

	std::cout << (t2 - t1) / CLOCKS_PER_SEC<<" "<<(t3 - t2) / CLOCKS_PER_SEC << std::endl;
	return 0;
}