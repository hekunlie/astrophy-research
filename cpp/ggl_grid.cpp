#include<FQlib.h>
#include<mpi.h>

int main(int argc, char *argv[])
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	int i, j, k;
	int area_num = 4;
	int area_id;

	int z_id = 0, ra_id = 1, dec_id = 2, mg1_id = 3, mg2_id = 4, mn_id = 5, mu_id = 6, mv_id = 7;
	int nib_id = 8, bs_id = 9, be_id = 10, bdy_id = 11, bdx_id = 12;
	// for the initial data
	int num_ini;
	int data_col = 8;	
	double *data_ini[13];
	// for others
	double *data_temp[13];
	int *mask, count;

	int *num_in_block, *block_start, *block_end;
	double *block_boundy, *block_boundx;

	// redshift range for foreground galaxy
	double low_z, high_z, block_scale;
	low_z = atof(argv[1]);
	high_z = atof(argv[2]);
	block_scale = atof(argv[3]);

	char h5f_path_1[150], h5f_path_2[150], log_inform[200];
	char data_path[150];

	char  attrs_name[50];
	int shape[2];
	char set_name[50],*names[13];
	for (i = 0; i < 13; i++)
	{
		names[i] = new char[20];
	}
	sprintf(names[0], "Z");
	sprintf(names[1], "RA");
	sprintf(names[2], "DEC");
	sprintf(names[3], "G1");
	sprintf(names[4], "G2");
	sprintf(names[5], "N");
	sprintf(names[6], "U");
	sprintf(names[7], "V");
	sprintf(names[8], "num_in_block");
	sprintf(names[9], "block_start");
	sprintf(names[10], "block_end");
	sprintf(names[11], "block_boundy");
	sprintf(names[12], "block_boundx");

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/");
	sprintf(h5f_path_1, "%scata_result_ext_cut.hdf5", data_path);
	sprintf(h5f_path_2, "%scata_result_ext_grid.hdf5", data_path);

	if (0 == rank)
	{
		sprintf(set_name, "/foreground/");
		create_h5_group(h5f_path_2, set_name, TRUE);
		sprintf(set_name, "/background/");
		create_h5_group(h5f_path_2, set_name, FALSE);
	}
	//sprintf(log_inform, "Foreground: %d galaxies in W_%d", area_id, count + 1);
	//std::cout << log_inform << std::endl;
	
	for (area_id = 1; area_id < area_num+1; area_id++)
	{
		sprintf(attrs_name, "shape");
		sprintf(set_name, "/w_%d/Z", area_id);
		read_h5_attrs(h5f_path_1, set_name, attrs_name, shape, "d");
		num_ini = shape[0];
		mask = new int[num_ini] {};
		for (i = 0; i < data_col; i++)
		{
			data_ini[i] = new double[num_ini];
			sprintf(set_name, "/w_%d/%s", area_id, names[i]);
			read_h5(h5f_path_1, set_name, data_ini[i]);
		}
	
		// the foreground galaxy
		if (0 == rank)
		{	
			count = 0;
			for (i = 0; i < num_ini; i++)
			{
				// find the foreground galaxy
				if (low_z <= data_ini[z_id][i] and data_ini[z_id][i] <= high_z)
				{
					mask[i] = 1;
					count += 1;
				}
			}
			shape[0] = count;
			shape[1] = 1;
			
			for (i = 0; i < data_col; i++)
			{	
				// alloc the array for foreground galaxy
				data_temp[i] = new double[count];
			}
			count = 0;
			for (i = 0; i < num_ini; i++)
			{
				
				if (1 == mask[i])
				{
					for (j = 0; j < data_col; j++)
					{
						data_temp[j][count] = data_ini[j][i];
					}
					count += 1;
				}				
			}
			delete[] mask;

			// create the file
			sprintf(attrs_name, "shape");
			for (j = 0; j < data_col; j++)
			{
				// write to file
				sprintf(set_name, "/foreground/w_%d/%s", area_id, names[j]);
				write_h5(h5f_path_2, set_name, data_temp[j], count, 1, FALSE);
				write_h5_attrs(h5f_path_2, set_name, attrs_name, shape, 2, "d");
				// free the memory
				delete[] data_temp[j];
			}

			sprintf(log_inform, "Foreground: %d galaxies in W_%d", count, area_id);
			std::cout << log_inform << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// build the grid


	}



	return 0;
}