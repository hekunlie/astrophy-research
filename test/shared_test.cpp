#include<FQlib.h>
#include<mpi.h>

int main(int argc, char* argv[])
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
	int disp_unit;

	long i, j;
	double *foregal_data_shared, *backgal_data_shared;
	long foregal_num = 8000, radius_num = 12;
	long backgal_num = 1000000000;
	double count_1=0, count_2 = 0;
	double count = 0;
	long num1 = foregal_num, num2 = backgal_num;
	for (j = 0; j < 4; j++)
	{
		count_1 = 0;
		count_2 = 0;
		MPI_Win win_1, win_2;
		MPI_Aint size_1, size_2;
		size_1 = num1 * sizeof(double);
		size_2 = num2 * sizeof(double);
		if (0 == rank)
		{
			MPI_Win_allocate_shared(size_1, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &foregal_data_shared, &win_1);
			MPI_Win_allocate_shared(size_2, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &backgal_data_shared, &win_2);
		}
		else
		{
			
			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &foregal_data_shared, &win_1);
			MPI_Win_shared_query(win_1, 0, &size_1, &disp_unit, &foregal_data_shared);

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &backgal_data_shared, &win_2);
			MPI_Win_shared_query(win_2, 0, &size_2, &disp_unit, &backgal_data_shared);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == rank)
		{
			for (i = 0; i < num1; i++)
			{
				foregal_data_shared[i] = 1+j;
			}

			for (i = 0; i < num2; i++)
			{
				backgal_data_shared[i] = 1 + j;
				
			}

			std::cout << "Initialize the buffer 1 & 2" << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (2 == rank)
		{
			for (i = 0; i < num1; i++)
			{
				count_1 += foregal_data_shared[i];
			}
			count = 0;
			
			for (i = 0; i < num2; i++)
			{
				count_2 += backgal_data_shared[i];
				count += 1;
			}
			backgal_data_shared[0] = count;
			std::cout << count_1<<" "<<count_2 <<" "<<count<<" "<< backgal_data_shared[0] << std::endl;
		}

		MPI_Win_free(&win_1);
		MPI_Win_free(&win_2);
		if (2 == rank)
		{
			std::cout << "Free the buffer 1 & 2" << std::endl;
		}
	}
	MPI_Finalize();
	return 0;
}