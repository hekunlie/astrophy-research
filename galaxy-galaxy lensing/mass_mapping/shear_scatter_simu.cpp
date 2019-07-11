#include<FQlib.h>
#include<mpi.h>

int main(int argc, char**argv)
{
	/* generate the galaxies on the shear suface, i.e. g = f(x,y)										*/
	/* to see the scatter of the estimated shears.																*/			
	/* because the shear itself has the scatter on the surface, which simulates				*/
	/* the real case that the shear measured in a area has more than one value.			*/


	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	int i, j, k, temp;
	int total_chip_num, my_chip_s, my_chip_e, my_chip_num;
	double st1, st2, st1c, st2c;

	int data_row, data_col = 5;

	int total_data_num, my_data_num, data_start;

	int stamp_nx;
	int gal_num, size, point_num = 80;
	int psf_type=2;
	double psf_scale = 4, max_radius =8;
	double flux_i;

	int shear_num;
	double *g1, *g2;
	
	double *points;
	double *psf_stamp, *gal_stamp, *big_img;
	float *big_img_f;
	double *psf_pow, *gal_pow;
	double *shear_data_shared;
	int *data_num_shared;

	char shear_path[200], chip_path[200], parent_path[200], data_path[200];
	char set_name[30], time_curr[50];
	char inform[200];

	// intput the total chip number of each shear point
	// each chip contains 10000 galaxies
	size = 48;
	stamp_nx = 100;
	gal_num = stamp_nx* stamp_nx;
	total_chip_num = atoi(argv[1]);
	total_data_num = total_chip_num * gal_num * data_col;

	sprintf(parent_path, "/mnt/perc/hklee/CFHT/multi_shear_/");
	sprintf(shear_path, "/mnt/perc/hklee/CFHT/multi_shear_/data/shear.hdf5");

	// read the shear
	sprintf(set_name, "/g1");
	read_h5_datasize(shear_path, set_name, shear_num);

	g1 = new double[shear_num];
	g2 = new double[shear_num];

	sprintf(set_name, "/g1");
	read_h5(shear_path, set_name, g1);
	sprintf(set_name, "/g2");
	read_h5(shear_path, set_name, g2);

	MPI_Win win_shear_data, win_data_num;
	MPI_Aint  data_size;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(total_data_num * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_data_shared, &win_shear_data);
		MPI_Win_allocate_shared(numprocs * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &data_num_shared, &win_data_num);
	}
	else
	{
		int dispu_total;

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_data_shared, &win_shear_data);
		MPI_Win_shared_query(win_shear_data, 0, &data_size, &dispu_total, &shear_data_shared);

		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &data_num_shared, &win_data_num);
		MPI_Win_shared_query(win_data_num, 0, &data_size, &dispu_total, &data_num_shared);
	}

	MPI_Barrier(MPI_COMM_WORLD);


	my_chip_s = 0;
	my_chip_e = 0;
	if (total_chip_num <= numprocs)
	{	
		// if cpus >= chip number, some cpus have nothing to do
		if (rank < total_chip_num)
		{
			my_chip_s = rank;
			my_chip_e = rank + 1;
		}
	}
	else
	{
		i = total_chip_num / numprocs;
		j = total_chip_num % numprocs;

		my_chip_s = rank * i;
		my_chip_e = (rank + 1)*i;

		if (rank == numprocs - 1)
		{
			my_chip_e = (rank + 1)*i + j;
		}
	}
	my_chip_num = my_chip_e - my_chip_s;
	my_data_num = my_chip_num * gal_num * data_col;
	data_num_shared[rank] = my_data_num;

	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(inform, "Total: %d shears. Total: %d chips. Rank %d gets [%d, %d], %d data_point",shear_num, total_chip_num, rank, my_chip_s, my_chip_e, my_data_num);
	std::cout << inform << std::endl;
	
	para all_paras;
	all_paras.gal_noise_sig = 0;
	all_paras.psf_noise_sig = 0;
	all_paras.stamp_size = size;
	all_paras.max_source = 30;
	all_paras.area_thres = 5;
	all_paras.detect_thres = 0;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = max_radius; 

	psf_stamp = new double[size*size];
	psf_pow = new double[size*size];

	points = new double[point_num * 2];

	gal_stamp = new double[size*size];
	gal_pow = new double[size*size];
	
	big_img = new double[gal_num*size*size];
	big_img_f = new float[gal_num*size*size];

	create_psf(psf_stamp, psf_scale, size, 2);
	pow_spec(psf_stamp, psf_pow, size, size);
	get_psf_radius(psf_pow, &all_paras, 2);

	if (rank == 0)
	{
		sprintf(chip_path, "!%sdata/psf.fits", parent_path);
		write_fits(chip_path, psf_stamp, size, size);
		sprintf(chip_path, "!%sdata/psf_pow.fits", parent_path);
		write_fits(chip_path, psf_pow, size, size);
	}

	for (i = 0; i < shear_num; i++)
	{	
		st1 = clock();
		if (rank == 0)
		{
			initialize_arr(shear_data_shared, total_data_num, 0);
			
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		for (j = my_chip_s; j < my_chip_e; j++)
		{	
			st1c = clock();
			gsl_initialize(rank*i+j + i + 100000);

			data_row = gal_num * data_col * j;

			sprintf(chip_path, "!%sdata/shear_%d_chip_%d.fits", parent_path, i, j);
			initialize_arr(big_img, size*size*gal_num, 0);

			for (k = 0; k < gal_num; k++)
			{
				initialize_arr(gal_stamp, size*size, 0);
				initialize_arr(gal_pow, size*size, 0);
				initialize_arr(points, point_num * 2, 0);

				initialize_para(&all_paras);

				create_points(points, point_num, max_radius);
				flux_i = 100;

				convolve(gal_stamp, points, flux_i, size, point_num, 0, psf_scale, g1[i], g2[i], psf_type, 0, &all_paras);

				pow_spec(gal_stamp, gal_pow, size, size);

				shear_est(gal_pow, psf_pow, &all_paras);

				stack(big_img, gal_stamp, k, size, stamp_nx, stamp_nx);

				shear_data_shared[data_row + k * data_col] = all_paras.n1;
				shear_data_shared[data_row + k * data_col + 1] = all_paras.n2;
				shear_data_shared[data_row + k * data_col + 2] = all_paras.dn;
				shear_data_shared[data_row + k * data_col + 3] = all_paras.du;
				shear_data_shared[data_row + k * data_col + 4] = all_paras.dv;
			}
			for (k = 0; k < stamp_nx*size* stamp_nx*size; k++)
			{
				big_img_f[k] = big_img[k];
			}
			write_fits(chip_path, big_img, stamp_nx*size, stamp_nx*size);
			gsl_free();
			st2c = clock();
			if (rank == 0)
			{	
				get_time(time_curr, 50);
				std::cout << "Chip " << j << ": " << (st2c - st1c) / CLOCKS_PER_SEC<<" "<<time_curr << std::endl;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		st2 = clock();
		// write down the shear data of each shear point
		if (rank == 0)
		{
			sprintf(data_path, "%sresult/shear_%d.hdf5", parent_path, i);
			sprintf(set_name, "/data");
			write_h5(data_path, set_name, shear_data_shared, total_data_num / data_col, data_col, TRUE);
			get_time(time_curr, 50);
			std::cout << "Shear " << i << ": " << (st2 - st1) / CLOCKS_PER_SEC << " " << time_curr<< std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
	MPI_Barrier(MPI_COMM_WORLD);

	delete[] big_img;
	delete[] gal_stamp;
	delete[] psf_stamp;
	delete[] gal_pow;
	delete[] psf_pow;
	delete[] points;
	delete[] g1;
	delete[] g2;

	MPI_Finalize();
	return 0;
}