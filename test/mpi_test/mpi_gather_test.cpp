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
	num = rank;
	double *sendbuf = new double[num] {};
	for (i = 0; i < num; i++)
	{
		sendbuf[i] = rank;
	}
	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{
			if (num > 0)
			{
				show_arr(sendbuf, 1, num);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	int *count = new int[numprocs];
	MPI_Allgather(&num, 1, MPI_INT, count, 1, MPI_INT, MPI_COMM_WORLD);
	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{
			show_arr(count, 1, numprocs);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	

	int total_num;
	int *disp = new int[numprocs] {};

	sum_arr(count, numprocs, 0, numprocs, total_num);

	double *recvbuf;
	if (rank == 0)
	{
		recvbuf = new double[total_num] {};
	}

	for (i = 0; i < numprocs; i++)
	{
		for (j = 0; j < i; j++)
		{
			disp[i] += count[j];
		}
	}
	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{
			if (num > 0)
			{
				show_arr(disp, 1, numprocs);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Gatherv(sendbuf, num, MPI_DOUBLE, recvbuf, count, disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		show_arr(recvbuf, 1, total_num);
	}
	MPI_Finalize();
	return 0;
}