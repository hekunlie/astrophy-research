#include<FQlib.h>
#include<mpi.h>

int main(int argc, char**argv)
{
	/* generate the galaxies on the shear suface, i.e. g = f(x), and measure their shape								*/

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	int i, j, k, temp, seed;
	int chip_num, expo_num;
	double st1, st2, st1c, st2c;

	int total_data_num, my_num_st, my_num_ed;

	double *shear, g1, g2;
	int shear_num;
	double *flux, flux_i;

	int stamp_nx, stamp_ny, stamp_num;
	int source_num, size, point_num;
	int psf_type;
	double psf_scale , max_radius;
	double gal_noise_sigma;

	double *points;
	double *psf_stamp, *gal_stamp, *big_img;
	double *psf_pow, *gal_pow;
	double *noise, *noise_pow;
	double *shear_data_shared, *result_data;
	int *data_num_shared;
	int data_row, data_col;

	char param_path[200], chip_path[200], parent_path[200], data_path[200];
	char set_name[30], time_curr[50];
	char inform[200];

	sprintf(parent_path, "/mnt/perc/hklee/multi_shear_dectect");

	sprintf(param_path, "%s/param_0.hdf5", parent_path);

	read_h5_datasize(param_path, "/flux", source_num);

	   
	size = 50;
	point_num = 80;
	max_radius = 8;
	
	gal_noise_sigma = 70;

	psf_type = 2;
	psf_scale = 4;

	stamp_nx = 10;
	stamp_ny = 10;
	stamp_num = stamp_nx * stamp_ny;

	data_row = source_num;
	data_col = 7; // 5 shear estimators + g, x
	total_data_num = data_row * data_col;

	i = source_num / numprocs;
	j = source_num % numprocs;
	my_num_st = i * rank;
	my_num_ed = i * (rank + 1);
	if (rank == numprocs - 1)
	{
		my_num_ed += j;
	}
	shear_num = 3;
	shear = new double[2 * shear_num]{-0.04,0.01,0.03,0.04,0.03,-0.02};
	// for (i = 0; i < shear_num; i++)
	// {
	// 	shear[i] = -0.05 + i * 0.01;
	// 	shear[i + shear_num] = 0.05 - i * 0.01;
	// }



	para all_paras;
	all_paras.gal_noise_sig = gal_noise_sigma;
	all_paras.psf_noise_sig = 0;
	all_paras.stamp_size = size;
	all_paras.max_source = 30;
	all_paras.area_thres = 5;
	all_paras.detect_thres = 0;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = max_radius;

	flux = new double[source_num];

	psf_stamp = new double[size*size];
	psf_pow = new double[size*size];

	points = new double[point_num * 2];

	gal_stamp = new double[size*size];
	gal_pow = new double[size*size];

	noise = new double[size*size];
	noise_pow = new double[size*size];

	big_img = new double[stamp_nx*stamp_ny*size*size];

	create_psf(psf_stamp, psf_scale, size, psf_type);
	pow_spec(psf_stamp, psf_pow, size, size);
	get_psf_radius(psf_pow, &all_paras, 2);

	MPI_Win win_data;
	MPI_Aint  data_size;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(total_data_num * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &result_data, &win_data);
	}
	else
	{
		int dispu_total;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &result_data, &win_data);
		MPI_Win_shared_query(win_data, 0, &data_size, &dispu_total, &result_data);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0)
	{

		std::cout << "Total: " << source_num << std::endl;
		sprintf(chip_path, "!%s/psf.fits", parent_path);
		write_fits(chip_path, psf_stamp, size, size);
		sprintf(chip_path, "!%s/psf_pow.fits", parent_path);
		write_fits(chip_path, psf_pow, size, size);
	}

	//for (i = 0; i < numprocs; i++)
	//{
	//	if (i == rank)
	//	{
	//		std::cout << rank << " " << my_num_st << " " << my_num_ed << std::endl;
	//	}
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}

	
	for (j = 0; j < shear_num; j++)
	{			
		if (rank == 0)
		{
			initialize_arr(result_data, total_data_num, 0);
		}
		seed = 1233 + rank*50 + rank + j;
		gsl_initialize(seed);
		k = 0;

		g1 = shear[j];
		g2 = shear[j + shear_num];

		sprintf(param_path, "%s/param_%d.hdf5", parent_path, j);
		sprintf(set_name, "/flux");
		read_h5(param_path, set_name, flux);

		st1c = clock();
		for (i = my_num_st; i < my_num_ed; i++)
		{
			initialize_arr(gal_stamp, size*size, 0);
			initialize_arr(gal_pow, size*size, 0);
			initialize_arr(points, point_num * 2, 0);
			initialize_arr(noise, size*size, 0);
			initialize_arr(noise_pow, size*size, 0);
			initialize_para(&all_paras);

			create_points(points, point_num, max_radius);
			flux_i = flux[i] / point_num;

			convolve(gal_stamp, points, flux_i, size, point_num, 0, psf_scale, g1, g2, psf_type, 1, &all_paras);

			addnoise(gal_stamp, size*size, gal_noise_sigma);

			addnoise(noise, size*size, gal_noise_sigma);
			pow_spec(noise, noise_pow, size, size);

			pow_spec(gal_stamp, gal_pow, size, size);
			noise_subtraction(gal_pow, noise_pow, &all_paras, 1, 1);

			shear_est(gal_pow, psf_pow, &all_paras);

			result_data[i * data_col] = all_paras.n1;
			result_data[i * data_col + 1] = all_paras.n2;
			result_data[i * data_col + 2] = all_paras.dn;
			result_data[i * data_col + 3] = all_paras.du;
			result_data[i * data_col + 4] = all_paras.dv;

			if (rank == 0 and k < 100)
			{
				stack(big_img, gal_stamp, k, size, stamp_ny, stamp_nx);
			}
			k++;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		st2c = clock();

		if (rank == 0)
		{
			sprintf(chip_path, "!%s/check_%d.fits", parent_path,j);
			write_fits(chip_path, big_img, stamp_ny*size, stamp_nx*size);
			// write down the shear data of each shear point
			sprintf(data_path, "%s/result_%d.hdf5", parent_path, j);
			sprintf(set_name, "/data");
			write_h5(data_path, set_name, result_data, total_data_num / data_col, data_col, TRUE);
			get_time(time_curr, 50);
			std::cout << (st2c - st1c) / CLOCKS_PER_SEC << " " << time_curr << std::endl;
		}
		gsl_free();
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Win_free(&win_data);

	delete[] big_img;
	delete[] gal_stamp;
	delete[] psf_stamp;
	delete[] gal_pow;
	delete[] psf_pow;
	delete[] points;
	delete[] shear;
	delete[] flux;

	MPI_Finalize();
	return 0;
}