#include <mpi.h>
#include<FQlib.h>

int main(int argc, char **argv)
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	int i, j, num;
	num = 50000000;
	double *sendbuf = new double[num] {};
	double *recvbuf = new double[num] {};
	if (rank == 0)
	{
		for (i = 0; i < num; i++)
		{
			sendbuf[i] = 1;
		}
	}
	int *disp = new int[num];
	int *send_num = new int[num];
	for (i = 0; i < num; i++)
	{
		disp[i] = 0;
		send_num[i] = num;
 	}

	MPI_Status status;

	double t1, t2, t3;
	t1 = clock();
	MPI_Scatterv(sendbuf, send_num, disp, MPI_DOUBLE, recvbuf, num, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// stactter
	double s;
	for (i = 0; i < numprocs; i++)
	{
		if (rank == i)
		{
			sum_arr(recvbuf, num, 0, num, s);
			std::cout << s << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// send & receive
	t2 = clock();
	if (rank == 0)
	{
		for (i = 1; i < numprocs; i++)
		{
			MPI_Send(sendbuf, num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(sendbuf, num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
	}

	for (i = 0; i < numprocs; i++)
	{
		if (rank == i)
		{
			sum_arr(sendbuf, num, 0, num, s);
			std::cout << s << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	t3 = clock();

	std::cout << (t2 - t1) / CLOCKS_PER_SEC<<" "<<(t3 - t2) / CLOCKS_PER_SEC << std::endl;
	return 0;
}