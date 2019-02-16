#include<FQlib.h>
#include<mpi.h>

int main(int argc, char *argv[])
{
	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	// initialize GSL
	int seed;
	gsl_initialize(123);

	// read the parameters of the data structure	
	int max_area[1]{}; //how many sky areas
	int max_buffer[1]{}; // the maximum data number of all sky areas (for convenience)
	int max_block_length[1]{}; // the maximum data number of all blocks
	int max_grid_length[1]{}; // the maximum grid number of all sky areas
	int grid_info[2]{}; // the grid shape to be read of each area
	int grid_ny, grid_nx, grid_num; // number of blocks of each area along each axis

	int data_col = 7; // columns of data 
	// the column of each quantity in the data array
	int ra_idx = 0, dec_idx = 1;
	int mg1_idx = 2, mg2_idx = 3, mn_idx = 4, mu_idx = 5, mv_idx = 6;

	char h5f_path[150], set_name[50], attrs_name[50];
	char log_inform[200];

	sprintf(h5f_path, "......hdf5");
	// of the grid
	sprintf(set_name, "/grid");
	sprintf(attrs_name, "max_area");
	read_h5_attrs(h5f_path, set_name, attrs_name, max_area);

	// loop the areas
	int i, j, k;
	int ix, iy;

	for (i = 0; i < max_area[0]; i++)
	{
		sprintf(attrs_name, "max_blcok_length");
		read_h5_attrs(h5f_path, set_name, attrs_name, max_block_length);
		sprintf(attrs_name, "max_grid_length");
		read_h5_attrs(h5f_path, set_name, attrs_name, max_grid_length);
		sprintf(attrs_name, "max_buffer");
		read_h5_attrs(h5f_path, set_name, attrs_name, max_buffer);


		// the shared buffer in the memory
		int arr_len = data_col * max_buffer[0];
		MPI_Win win, win_err;
		MPI_Aint size, size_err;
		double *buffer_ptr;
		int *err_ptr;
		if (rank == 0)
		{
			size = arr_len * sizeof(double);
			size_err = sizeof(int);
			MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &buffer_ptr, &win);
			MPI_Win_allocate_shared(size_err, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &err_ptr, &win_err);
		}
		else
		{
			int disp_unit, disp_init_err;

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &buffer_ptr, &win);
			MPI_Win_shared_query(win, 0, &size, &disp_unit, &buffer_ptr);

			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &err_ptr, &win_err);
			MPI_Win_shared_query(win_err, 0, &size_err, &disp_init_err, &err_ptr);		
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == rank)
		{
			err_ptr[0] = 1;
		}

		// read the shape of the grid "w_i"
		sprintf(set_name, "/grid/w_%d",i);
		sprintf(attrs_name, "grid_info");
		read_h5_attrs(h5f_path, set_name, attrs_name, grid_info);
		grid_ny = grid_info[0];
		grid_nx = grid_info[1];
		grid_num = grid_ny * grid_nx;		

		// initialize the data in the buffer for all threads
		if (0 == rank)
		{	
			int temp_size = max_block_length[0] * data_col;
			int block_shape[2]{}, num_in_block;
			int block_id;
			double *temp = new double[temp_size];

			sprintf(log_inform, "INITIALIZE THE BIG BUFFER FOR W_%d: ", i);
			std::cout << log_inform << std::endl;
			sprintf(log_inform, "the biggest blcok: %d x %d", max_block_length[0] , data_col);
			std::cout << log_inform << std::endl;
			for (iy = 0; iy < grid_ny; iy++)
			{
				if (0 == err_ptr[0])
				{
					break;
				}
				for (ix = 0; ix < grid_nx; ix++)
				{
					// initialize temp
					initialize_arr(temp, temp_size, -1.);

					sprintf(set_name, "/grid/w_%d/%d/%d", i, iy, ix);
					sprintf(attrs_name, "shape");
					read_h5_attrs(h5f_path, set_name, attrs_name, block_shape);
					read_h5(h5f_path, set_name, temp);
					num_in_block = block_shape[0] * block_shape[1];

					sprintf(log_inform, "Initialzie the block [%d, %d], The true shape (%d, %d)", iy, ix, block_shape[0], block_shape[1]);
					std::cout << log_inform << std::endl;
					// check
					if (grid_info[1] != data_col)
					{
						sprintf(log_inform, "The data shape, (%d, %d), does not match the data_col, %d ", iy, ix, block_shape[0], block_shape[1], data_col);
						std::cout << log_inform << std::endl;
						err_ptr[0] = 0;
						break;
					}
					if (-1. != temp[num_in_block])
					{
						sprintf(log_inform, "Something may be wrong!!! The first elemet, %g, after the target elements is not zero!!", temp[num_in_block]);
						std::cout << log_inform << std::endl;
						err_ptr[0] = 0;
						break;
					}
					
					block_id = iy * grid_nx + ix;
					for (k = 0; k < temp_size; k++)
					{
						buffer_ptr[block_id*temp_size + k] = temp[k];
					}
				}
			}
			delete[] temp;
			sprintf(log_inform, "FINISH BUFFER INITIALZIATION");
			std::cout << log_inform << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == err_ptr[0])
		{
			exit(0);
		}

		// mask for non-repeating calculation
		int *mask = new int[max_block_length[0]]{};
		initialize_arr(mask, max_block_length[0] , 1);

		// the bounds of each block
		int *boundy = new int[grid_num * 4]{};
		int *boundx = new int[grid_num * 4]{};

		sprintf(set_name, "/grid/w_%d/boundy", i);
		read_h5(h5f_path, set_name, boundy);
		sprintf(set_name, "/grid/w_%d/boundx", i);
		read_h5(h5f_path, set_name, boundx);

		// each thread gets its own targets blocks
		int *task_list = new int[grid_num] {};
		int *my_tasks = new int[grid_num] {};
		for (k = 0; k < grid_num; k++)
		{
			task_list[k] = k;
		}
		// -1 labels the end
		task_alloc(task_list, grid_num, numprocs, rank, my_tasks);

		// loop the target blocks 
		for (k = 0; k < grid_num; k++)
		{
			if (my_tasks[k] > -1)
			{
				// the block (iy, ix) of area "w_i"
				iy = my_tasks[k] / grid_nx;
				ix = my_tasks[k] % grid_ny;
				!!!

			}
		}

		delete[] mask;
		delete[] boundy;
		delete[] boundx;
		delete[] task_list;
		delete[] my_tasks;

	}

	MPI_Win_free(&win);
	MPI_Win_free(&win_err);
	// free the GSL
	gsl_free();
	MPI_Finalize();
	return 0;
}