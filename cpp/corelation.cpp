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

	pts_info pts_data;

	int i, j, k;
	int ix, iy, ik, ir, ib, ib_;
	int label, label_, tag;

	// read the parameters of the data structure	

	int max_area[1]{}; // sky areas
	int radius_num[1]; // bin number 

	int grid_info[2]{}; // the grid shape to be read of each area
	int grid_ny, grid_nx, grid_num; // number of blocks of each area along each axis
	int max_block_length[2]{}; // [row, col] the maximum data number of all blocks in one sky area
	double block_scale[1]{}; // the length of the block side

	int data_col, block_size; // columns of data & size of array of the biggest block

	// the column of each quantity in the data array
	int ra_idx = 0, dec_idx = 1;
	int mg1_idx = 2, mg2_idx = 3, mn_idx = 4, mu_idx = 5, mv_idx = 6;

	double ra, dec, distance_sq;

	char h5f_path[150], set_name[50], attrs_name[50];
	char log_inform[200];

	sprintf(h5f_path, "......hdf5");
	// read the number of the sky areas
	sprintf(set_name, "/grid");
	sprintf(attrs_name, "max_area");
	read_h5_attrs(h5f_path, set_name, attrs_name, max_area);
	// the bin number
	sprintf(set_name, "/radius");
	sprintf(attrs_name, "bin_num");
	read_h5_attrs(h5f_path, set_name, attrs_name, radius_num);
	// the bins of radian for correlation function
	double *radius = new double[radius_num[0] + 1]{};
	read_h5(h5f_path, set_name, radius);


	// loop the areas
	for (i = 0; i < max_area[0]; i++)
	{		
		sprintf(set_name, "/grid/w_%d", i);

		// read the shape of the grid "w_i"
		sprintf(attrs_name, "grid_info");
		read_h5_attrs(h5f_path, set_name, attrs_name, grid_info);
		grid_ny = grid_info[0];
		grid_nx = grid_info[1];
		grid_num = grid_ny * grid_nx;

		// of the single block
		sprintf(attrs_name, "max_block_length");
		read_h5_attrs(h5f_path, set_name, attrs_name, max_block_length);
		data_col = max_block_length[1];
		// the maximum block size. row * col
		block_size = max_block_length[0] * max_block_length[1];

		sprintf(attrs_name, "block_scale");
		read_h5_attrs(h5f_path, set_name, attrs_name, block_scale);

	
		// mask for non-repeating calculation
		int *mask = new int[max_block_length[0]]{};

		// the bounds of each block
		double *boundy = new double[grid_num * 4]{};
		double *boundx = new double[grid_num * 4]{};
		// the origin of the coordinates is the first cross of the grid (array)
		block_bound(block_scale[0], grid_ny, grid_nx, boundy, boundx);

		// each thread gets its own targets blocks in my_tasks
		int *task_list = new int[grid_num] {};
		int *my_tasks = new int[grid_num] {};
		int *search_blocks = new int[grid_num] {};

		// initialize the information of pts_info
		pts_data.scale = block_scale[0];
		pts_data.ny = grid_ny;
		pts_data.nx = grid_nx;
		pts_data.blocks_num = grid_num;

		// the shared buffer in the memory
		int arr_len = block_size *grid_num;
		MPI_Win win, win_num, win_err;
		MPI_Aint size, size_num, size_err;
		double *buffer_ptr;
		int *num_ptr;
		int *err_ptr;
		if (rank == 0)
		{
			size = arr_len * sizeof(double);	
			MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &buffer_ptr, &win);
			size_num = grid_num * sizeof(int);
			MPI_Win_allocate_shared(size_num, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &num_ptr, &win_num);
			size_err = sizeof(int);
			MPI_Win_allocate_shared(size_err, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &err_ptr, &win_err);
		}
		else
		{
			int disp_unit, disp_unit_num, disp_unit_err;

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &buffer_ptr, &win);
			MPI_Win_shared_query(win, 0, &size, &disp_unit, &buffer_ptr);

			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &num_ptr, &win_num);
			MPI_Win_shared_query(win_num, 0, &size_num, &disp_unit_num, &num_ptr);

			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &err_ptr, &win_err);
			MPI_Win_shared_query(win_err, 0, &size_err, &disp_unit_err, &err_ptr);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == rank)
		{
			err_ptr[0] = 1;
			sprintf(log_inform, "Create the shared buffer for program. %d", err_ptr[0]);
			std::cout << log_inform << std::endl;
		}


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

			// read the number of points in each block
			sprintf(set_name, "/grid/w_%d/num_in_block", i);
			read_h5(h5f_path, set_name, num_ptr);

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
					if (fabs( temp[num_in_block] + 1)>0.00001)
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
			
			   
		for (k = 0; k < grid_num; k++)
		{
			task_list[k] = k;
		}
		// -1 labels the end
		task_alloc(task_list, grid_num, numprocs, rank, my_tasks);

		// loop the target blocks 
		for (k = 0; k < grid_num; k++)
		{
			// initialize the mask for non-repeating calculation
			initialize_arr(mask, max_block_length[0], 1);
			if (my_tasks[k] > -1)
			{
				// the block (iy, ix) of area "w_i"
				iy = my_tasks[k] / grid_nx;
				ix = my_tasks[k] % grid_ny;
				label_ = my_tasks[k] * block_size;

				// loop the point in the block
				for (ik = 0; ik < block_size; ik++)
				{
					// remember 8 parameters in pts_data must be initialized
					pts_data.idy = iy;
					pts_data.idx = ix;
					label = label_ + ik * data_col;
					pts_data.y = buffer_ptr[label + dec_idx];
					pts_data.x = buffer_ptr[label + ra_idx];

					// loop the radius scales 
					for (ir = 0; ir < radius_num[0]; ir++)
					{
						find_block(&pts_data, radius[ir], radius[ir + 1], boundy, boundx, search_blocks);
						
						// find the pairs in the searched blocks
						for (ib = 0; ib < grid_num; ib++)
						{
							if (search_blocks[ib] > -1)
							{
								for (ib_ = 0; ib_ < num_ptr[search_blocks[ib]]; ib_++)
								{
									tag = search_blocks[ib] * block_size + ib_ * data_col;
									distance = 
								}
							}
						}

					}

				}

			}
		}

		delete[] mask;
		delete[] boundy;
		delete[] boundx;
		delete[] task_list;
		delete[] my_tasks;
		delete[] search_blocks;

		MPI_Win_free(&win);
		MPI_Win_free(&win_num);
		MPI_Win_free(&win_err);

	}

	delete[] radius;
	// free the GSL
	gsl_free();
	MPI_Finalize();
	return 0;
}