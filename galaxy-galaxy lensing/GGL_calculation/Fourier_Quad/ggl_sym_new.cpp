#include<FQlib.h>
#include<mpi.h>

#define MAX_DATA_COL 20
#define MAX_AREA 20
#define MG_BIN_NUM 8
#define ALLOT_THRESH 200000

//#define SMALL_CATA

#ifdef SMALL_CATA
#define MY_INT int
#else
#define MY_INT long
#endif 

int main(int argc, char ** argv)
{
	/* crit_coeff is the coefficient of the critical density, see the note for detials */
	
	/* argv[1]: the file name of foreground source												*/											
	/* argv[2] -- argv[n]: area labels, like 1,3, 4														*/

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	char total_path[300], data_path[300], set_name[30], result_path[300];
	char fore_source[50];
	char logs[300];

	int i, j, k;
	double st1, st2;

	int ia, areas[MAX_AREA]{}, area_num, area_id;

	int radius_num, radius_id;
	double *radius_bin;

	int *pair_num, *pair_each_radius, dataset_num, set_id, pair_count;
	int data_num, data_col, data_row;
	int num_label, small_num;
	double *data[MAX_DATA_COL], *total_data;

	
	// the foreground name
	strcpy(fore_source, argv[1]);
	// the radius label
	area_num = argc - 2;
	std::cout << fore_source << ": ";
	for (i = 2; i < argc; i++)
	{
		areas[i - 2] = atoi(argv[i]);
		std::cout << areas[i - 2] << " ";
	}
	std::cout<<std::endl;

	sprintf(total_path, "/mnt/perc/hklee/CFHT/gg_lensing/result/%s/fourier_cata_old/", fore_source);	

	// read the radius bin
	sprintf(data_path, "%sw_1/radius_0.hdf5", total_path);
	sprintf(set_name, "/radius_bin");
	read_h5_datasize(data_path, set_name, radius_num);
	radius_bin = new double[radius_num];
	read_h5(data_path, set_name, radius_bin);
	radius_num = radius_num - 1;

	pair_each_radius = new int[radius_num];
	initialize_arr(pair_each_radius, radius_num, 0);

	// count the total pair num in each radius bin
	for(radius_id =0; radius_id <radius_num; radius_id++)
	{
		for (ia = 0; ia < area_num; ia++)
		{
			sprintf(data_path, "%sw_%d/radius_%d.hdf5", total_path, areas[ia], radius_id);
			sprintf(set_name, "/pair_count");
			read_h5_datasize(data_path, set_name, dataset_num);
			pair_num = new int[dataset_num];
			read_h5(data_path, set_name, pair_num);
			sum_arr(pair_num, dataset_num, 0, dataset_num, pair_count);
			pair_each_radius[radius_id] += pair_count;
		}
	}
	std::cout << "Radius bin:" << std::endl;
	show_arr(radius_bin, 1, radius_num+1);
	std::cout << "Pair count:" << std::endl;
	show_arr(pair_each_radius, 1, radius_num);

	//////////////////////////////////////////////// calculate ////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// first, calculate the signal of the radius bin in which pair is too little
	// it will allot those to each threads
	small_num = 0;
	for (radius_id = 0; radius_id < radius_num; radius_id++)
	{
		if (pair_each_radius[radius_id] <= ALLOT_THRESH)
		{
			small_num++;
		}
	}
	if (0 < small_num <= numprocs)
	{
		if (rank < small_num)
		{
			for (ia = 0; ia < area_num; ia++)
			{
				sprintf(data_path, "%sw_%d/radius_%d.hdf5", total_path, areas[ia], rank);
				sprintf(set_name, "/pair_count");
				read_h5_datasize(data_path, set_name, dataset_num);
				pair_num = new int[dataset_num];
				read_h5(data_path, set_name, pair_num);
				sum_arr(pair_num, dataset_num, 0, dataset_num, pair_count);
				pair_each_radius[radius_id] += pair_count;
			}
			for (set_id = 0; set_id < dataset_num; set_id++)
			{
				data_row = pair_num[set_id];
				sprintf(set_name, "/pair_data_%d", set_id);
				read_h5_datasize(data_path, set_name, data_num);
				data_col = data_num / data_row;

			}

		}
		
	}
	else
	{

	}

	for (radius_id = small_num; radius_id < radius_num; radius_id++)
	{
		for (ia = 0; ia < area_num; ia++)
		{
			sprintf(data_path, "%sw_%d/radius_%d.hdf5", total_path, areas[ia], radius_id);
			sprintf(set_name, "/pair_count");
			read_h5_datasize(data_path, set_name, dataset_num);
			pair_num = new int[dataset_num];
			read_h5(data_path, set_name, pair_num);
			sum_arr(pair_num, dataset_num, 0, dataset_num, pair_count);
			pair_each_radius[radius_id] += pair_count;
		}
		for (set_id = 0; set_id < dataset_num; set_id++)
		{
			data_row = pair_num[set_id];
			sprintf(set_name, "/pair_data_%d", set_id);
			read_h5_datasize(data_path, set_name, data_num);
			data_col = data_num / data_row;

		}


	}

	

	MPI_Finalize();
	return 0;
}