#include<FQlib.h>
#include<mpi.h>

double chisq_2d(const double *hist_arr, const int size)
{	
	int h = size / 2, i, j, s1, s2=size*size;
	double chi=0, n,m;
	for (i = 0; i < h; i++)
	{	
		s1 = i * size;
		for (j = 0; j < h; j++)
		{
			m = hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] - (hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j]);
			n = hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] + hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j];
			chi += m * m / n;
		}
	}
	return chi * 0.5;
}

int main(int argc, char *argv[])
{
	int myid, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int size = 5, psize = 30;
	int i, j,  m,n,count = 0, turns=1000000;
	double *pool = new double[psize]{};
	double *buffer = new double[size]{};



	for (i = 0; i < turns; i++)
	{	
		for (j = 0; j < size; j++)
		{
			buffer[j] ++;
			count++;
			check_buffer(pool, buffer, size, size, count, 50000);
			if (i == turns - 1)
			{
				count = 50001;
				check_buffer(pool, buffer, size, size, count, 50000);
			}
			if (count == 0)
			{
				std::cout << "Re-count: " << i << ", " << j <<", "<<count<< std::endl;
				show_arr(pool, 6, 5);
				std::cout << std::endl;
				show_arr(buffer, 1, 5);
			}

		}
	}
	


	delete[] pool;
	delete[] buffer;
	MPI_Finalize();
	return 0;
}