#include"FQlib.h"
#include "mpi.h"
#include<hdf5.h>

int main(int argc, char*argv[])
{
	int myid, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
	para all_paras;
	int size = 48, data_col = 20, chip_num_total =10, nx=100, stamp_num=10000, chip_num, detect_label;
	chip_num = chip_num_total / numprocs;
	int data_row = chip_num*stamp_num, row_s;
	int data_row_t = chip_num_total*stamp_num;
	int chip_id_s, chip_id_e;
	int i, j, k,m,n;
	std::cout << sizeof(int) << std::endl;

	double noise_sig = 60, sig_level=1.5;

	double *gal = new double[size*size]{};
	double *pgal = new double[size*size]{};
	int *mask = new int[size*size]{};
	double *img = new double[size*size * stamp_num]{};
	int *mask_img = new int[size*size*stamp_num]{};
	double *data = new double[data_row*data_col]{};
	double *data_t = new double[data_row_t*data_col]{};
	std::cout << sizeof(mask[0]) << std::endl;
	char img_path[150], data_path[150], set_name[20];

	all_paras.gal_noise_sig = noise_sig;
	all_paras.psf_noise_sig = 0;
	all_paras.stamp_size = size;
	all_paras.max_source = 30;
	all_paras.area_thres = 5;
	all_paras.detect_thres = noise_sig * sig_level;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = 5.5; // because the max half light radius of the galsim source is 5.5 pixels

	chip_id_s = chip_num * myid;
	chip_id_e = chip_num * (myid + 1);
	
	for (i = chip_id_s; i < chip_id_e; i++)
	{	
		initialize_arr(img, size*size*stamp_num);
		initialize_arr(mask_img, size*size*stamp_num);
		sprintf(img_path, "/mnt/ddnfs/data_users/hkli/simu_test/scatter/%d.fits", i);
		read_fits(img, img_path);

		row_s = (i - chip_id_s)*stamp_num*data_col;

		for (j = 0; j < stamp_num; j++)
		{
			initialize_arr(gal, size*size);
			initialize_arr(pgal, size*size);
			initialize_arr(mask, size*size);

			segment(img, gal, j, size, nx, nx);
			pow_spec(gal, pgal, size, size);
			detect_label = galaxy_finder(gal, mask, &all_paras, false);
			if (0 == myid && i == 0 & j == 10)
			{
				sprintf(img_path, "!/mnt/ddnfs/data_users/hkli/simu_test/scatter/mask.fits");
				write_fits(mask, size, size, img_path);
				for (m = 0; m < size; m++)
				{
					for (n = 0; n < size; n++)
					{

						std::cout << mask[m*size + n] << " ";
						if (n == size - 1)
						{
							std::cout << std::endl;
						}						
					}					
				}				
			}
			stack(mask_img, mask, j, size, nx, nx);
			snr_est(pgal, &all_paras, 2);

			data[row_s + j * data_col + 0] = all_paras.gal_flux2;
			data[row_s + j * data_col + 1] = all_paras.gal_flux_alt;
			data[row_s + j * data_col + 2] = all_paras.gal_flux;
			data[row_s + j * data_col + 3] = all_paras.gal_osnr;

			data[row_s + j * data_col + 4] = all_paras.gal_flux2_ext[0];
			data[row_s + j * data_col + 5] = all_paras.gal_flux2_ext[1];
			data[row_s + j * data_col + 6] = all_paras.gal_flux2_ext[2];
			data[row_s + j * data_col + 7] = all_paras.gal_flux2_ext[3];
			data[row_s + j * data_col + 8] = all_paras.gal_flux2_ext[4];

			data[row_s + j * data_col + 9] = all_paras.gal_flux_ext[0];
			data[row_s + j * data_col + 10] = all_paras.gal_flux_ext[1];
			data[row_s + j * data_col + 11] = all_paras.gal_flux_ext[2];
			data[row_s + j * data_col + 12] = all_paras.gal_flux_ext[3];
			data[row_s + j * data_col + 13] = all_paras.gal_flux_ext[4];

			data[row_s + j * data_col + 14] = all_paras.gal_size_ext[0];
			data[row_s + j * data_col + 15] = all_paras.gal_size_ext[1];
			data[row_s + j * data_col + 16] = all_paras.gal_size_ext[2];
			data[row_s + j * data_col + 17] = all_paras.gal_size_ext[3];
			data[row_s + j * data_col + 18] = all_paras.gal_size_ext[4];

			data[row_s + j * data_col + 19] = detect_label;
		
		}

		sprintf(img_path, "!/mnt/ddnfs/data_users/hkli/simu_test/scatter/%d_mask.fits", i);
		write_fits(mask_img, nx*size, nx*size, img_path);

	}
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(set_name, "/data");
	MPI_Gather(data, data_row*data_col, MPI_DOUBLE, data_t, data_row*data_col, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (0 == myid)
	{
		sprintf(data_path, "/mnt/ddnfs/data_users/hkli/simu_test/scatter/result.hdf5");
		write_h5(data_path, set_name, data_row_t, data_col, data_t, NULL);
	}


	delete[] gal;
	delete[] pgal;
	delete[] mask;
	delete[] img;
	delete[] mask_img;
	delete[] data;
	delete[] data_t;
	return 0;
}