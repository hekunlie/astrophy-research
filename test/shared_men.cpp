#include<FQlib.h>
#include"mpi.h"
#include<thread>
#include<chrono>


int main(int argc, char*argv[])
{
	int myid, numprocs, namelen;	

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int arr_len = 50000000;
	//double *check = new double[arr_len];
	MPI_Win win;
	MPI_Aint size = arr_len *sizeof(double);
	double *baseptr;
	int i, j, k;
	if (myid == 0)
	{		
		MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &baseptr, &win);
	}
	else
	{

		int disp_unit;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &baseptr, &win);
		MPI_Win_shared_query(win, 0, &size, &disp_unit, &baseptr);
		
	}
	std::cout << size << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);

	// myid 0 initialize the array

	if (myid == 0)
	{
		for (i = 0; i < arr_len; i++)
		{
			//std::cout << i << std::endl;
			baseptr[i] = -1;
			//check[i] = -1;
		}
		std::cout << "Begin: " << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// each thread show the array in turns
	//for (i = 0; i < numprocs; i++)
	//{
	//	if (i == myid)
	//	{
	//		std::cout << myid << std::endl;
	//		for (j = 0; j < arr_len; j++)
	//		{					
	//			std::cout << baseptr[j] << ", ";
	//		}
	//		std::cout << std::endl;
	//	}
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}

	// change the part of the array

	//for (i = 0; i < numprocs; i++)
	//{
	//	//MPI_Win_lock_all(0, win);
	//	if (i == myid)
	//	{
	//		for (j = 0; j < 5; j++)
	//		{
	//			baseptr[i * 5 + j] = myid;
	//		}
	//	}
	//	//MPI_Win_unlock_all(0, win);
	//}

	int step = arr_len / numprocs;
	std::cout << step << std::endl;
	for (i = 0; i < 30; i++)
	{
		for (j = 0; j < step; j++)
		{
			baseptr[myid * step + j] = pow(myid + 1.1, 10)*sin(myid);
			//check[myid * step + j] = pow(myid + 1.1, 10)*sin(myid);

		}
	}
	std::this_thread::sleep_for(std::chrono::seconds(5));
	if (myid == 0)
	{
		std::cout << "After: " << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// each thread show the array in turns
	//for (i = 0; i < numprocs; i++)
	//{
	//	if (i == myid)
	//	{	

	//		std::cout << myid << std::endl;
	//		for (j = 0; j < arr_len; j++)			{
	//			
	//			std::cout << baseptr[j] << ", ";
	//		}
	//		std::cout << std::endl;
	//	}
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}
	//delete[] check;
	MPI_Win_free(&win);
	MPI_Finalize();
	return 0;
}