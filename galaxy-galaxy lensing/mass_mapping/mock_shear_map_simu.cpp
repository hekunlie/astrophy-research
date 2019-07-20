#include<FQlib.h>
#include<mpi.h>

int main(int argc, char**argv)
{
	/* generate the galaxies on the shear suface, i.e. g = f(x,y)										*/
	/* to see the scatter of the estimated shears.																*/
	/* because the shear itself has the scatter on the surface, which simulates				*/
	/* the real case that the shear measured in a area has more than one value.			*/
	/* each rank creates an exposure																				*/			

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	int i, j, k, temp, seed;
	int chip_num, expo_num;
	double st1, st2, st1c, st2c;

	int total_data_num, my_data_num, data_start;

	int my_expo_s, my_expo_e;

	int stamp_nx, stamp_ny, stamp_num;
	int source_num, size, point_num;
	int psf_type;
	double psf_scale , max_radius;
	double sigma;

	double *g1, *g2, *g;
	double *ra, *dec;
	double *flux, flux_i;

	double *points;
	double *psf_stamp, *gal_stamp, *big_img;
	float *big_img_f;
	double *psf_pow, *gal_pow;
	double *noise, *noise_pow;
	double *shear_data_shared, *result_data;
	int *data_num_shared;
	int data_row, data_col;

	char shear_path[200], chip_path[200], parent_path[200], data_path[200];
	char set_name[30], time_curr[50];
	char inform[200];

	// intput the total chip number of each shear point
	// each chip contains 10000 stamps
	size = 48;
	point_num = 80;
	max_radius = 8;
	
	sigma = 30;

	psf_type = 2;
	psf_scale = 4;

	stamp_nx = 10;
	stamp_ny = 10;
	stamp_num = stamp_nx * stamp_nx;
	source_num = 30000;
	expo_num = 50;
	chip_num = 1;
	data_row = source_num;
	data_col = 10; // 5 shear estimators + RA, DEC , g, g1 , g2
	total_data_num = chip_num * data_row * data_col;

	seed = 12335812+ rank*10 + rank;

	i = expo_num / numprocs;
	j = expo_num % numprocs;
	my_expo_s = i * rank;
	my_expo_e = i * (rank + 1);
	if (rank == numprocs - 1)
	{
		my_expo_e += j;
	}

	sprintf(parent_path, "/mnt/perc/hklee/CFHT/multi_shear/cluster_field/");
	sprintf(shear_path, "%sparam.hdf5", parent_path);

	g1 = new double[source_num];
	g2 = new double[source_num];
	g = new double[source_num];

	ra = new double[source_num];
	dec = new double[source_num];
	
	flux = new double[source_num];
	result_data = new double[total_data_num];

	// read the shear
	sprintf(set_name, "/g");
	read_h5(shear_path, set_name, g);
	sprintf(set_name, "/g1");
	read_h5(shear_path, set_name, g1);
	sprintf(set_name, "/g2");
	read_h5(shear_path, set_name, g2);

	// read  the  RA and DEC
	sprintf(set_name, "/ra");
	read_h5(shear_path, set_name, ra);
	sprintf(set_name, "/dec");
	read_h5(shear_path, set_name, dec);
	// read flux
	sprintf(set_name, "/flux");
	read_h5(shear_path, set_name, flux);

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

	noise = new double[size*size];
	noise_pow = new double[size*size];

	big_img = new double[stamp_nx*stamp_ny*size*size];
	big_img_f = new float[stamp_nx*stamp_ny*size*size];

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

	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{
			std::cout << rank << " " << my_expo_s << " " << my_expo_e << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	for (j = my_expo_s; j < my_expo_e; j++)
	{
		st1c = clock();
		gsl_initialize(seed);

		sprintf(chip_path, "!%sdata/expo_%d.fits", parent_path, j);
		initialize_arr(big_img, size*size*stamp_nx*stamp_ny, 0);

		for (k = 0; k < source_num; k++)
		{
			initialize_arr(gal_stamp, size*size, 0);
			initialize_arr(gal_pow, size*size, 0);
			initialize_arr(points, point_num * 2, 0);
			initialize_arr(noise, size*size, 0);
			initialize_arr(noise_pow, size*size, 0);
			initialize_para(&all_paras);

			create_points(points, point_num, max_radius);
			flux_i = flux[k]/point_num;

			convolve(gal_stamp, points, flux_i, size, point_num, 0, psf_scale, g1[k], g2[k], psf_type, 0, &all_paras);

			addnoise(gal_stamp, size*size, sigma);

			addnoise(noise, size*size, sigma);
			pow_spec(noise, noise_pow, size, size);

			pow_spec(gal_stamp, gal_pow, size, size);
			noise_subtraction(gal_pow, noise_pow, &all_paras, 1, 1);

			shear_est(gal_pow, psf_pow, &all_paras);

			if (k < 100)
			{
				stack(big_img, gal_stamp, k, size, stamp_nx, stamp_nx);
			}

			result_data[k * data_col] = all_paras.n1;
			result_data[k * data_col + 1] = all_paras.n2;
			result_data[k * data_col + 2] = all_paras.dn;
			result_data[k * data_col + 3] = all_paras.du;
			result_data[k * data_col + 4] = all_paras.dv;
			result_data[k * data_col + 5] = g[k];
			result_data[k * data_col + 6] = g1[k];
			result_data[k * data_col + 7] = g2[k];
			result_data[k * data_col + 8] = ra[k];
			result_data[k * data_col + 9] = dec[k];
		}
		for (k = 0; k < stamp_nx*size* stamp_nx*size; k++)
		{
			big_img_f[k] = big_img[k];
		}
		write_fits(chip_path, big_img_f, stamp_nx*size, stamp_nx*size);
		gsl_free();
		
		// write down the shear data of each shear point
		sprintf(data_path, "%sresult/expo_%d.hdf5", parent_path, j);
		sprintf(set_name, "/data");
		write_h5(data_path, set_name, result_data, total_data_num / data_col, data_col, TRUE);
		st2c = clock();

		get_time(time_curr, 50);
		std::cout << "EPXO: " << j <<" "<< (st2c - st1c) / CLOCKS_PER_SEC << " " << time_curr << std::endl;

	}



	delete[] big_img;
	delete[] gal_stamp;
	delete[] psf_stamp;
	delete[] gal_pow;
	delete[] psf_pow;
	delete[] points;
	delete[] g1;
	delete[] g2;
	delete[] flux;
	delete[] ra;
	delete[] dec;

	MPI_Finalize();
	return 0;
}