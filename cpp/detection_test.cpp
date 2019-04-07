#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"
#include<hdf5.h>
#include<stdio.h>
#include<string>


using namespace std;

int main(int argc, char*argv[])
{
	int myid, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	para all_paras;
	ifstream fin;
	string s, detect_info;

	/* 14 (g1,g2) points and each pairs contain 500 chips which cotians 10000 gals */
	int total_chip_num = 1000, chip_num, stamp_num = 10000, shear_pairs = 10;
	/* remember to change the data_cols when you change the number of estimators recorded */
	int i, j, seed, data_rows, shear_esti_data_cols = 7, snr_para_data_cols = 7, chip_id, shear_id, detect_label;
	int size = 64, num_p = 40, stamp_nx = 100, psf_type = 2;
	double psf_scale = 4., max_radius = 9, st, ed, s1, s2;
	double g1 = 0., g2 = 0.;
	double gal_noise_sig = 60, psf_noise_sig = 0., scale = 2.;
	int total_num = total_chip_num * stamp_num;
	
	chip_num = 1000;
	data_rows = chip_num*stamp_num;

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.detect_thres = gal_noise_sig * 1.5;
	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = psf_noise_sig;
	all_paras.max_distance = max_radius;
	initialize_para(&all_paras);

	double *big_img = new double[stamp_nx*stamp_nx*size*size]();
	int *mask_img = new int[stamp_nx*stamp_nx*size*size]();
	double *check_img = new double[stamp_nx*stamp_nx*size*size]();

	double *point = new double[2 * num_p]();
	double *gal = new double[size*size]();
	double *gpow = new double[size*size]();
	double *psf = new double[size*size]();
	double *ppow = new double[size*size]();
	double *noise = new double[size*size]();
	double *pnoise = new double[size*size]();
	int *mask = new int[size*size];
	char para_inform[500], data_path[60], shear_path[120], para_path[120], log_path[120], log_inform[120], chip_path[120], buffer[200];
	char h5_path[120], snr_h5_path[120];
	char set_name1[20], set_name2[20];
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer_m2/");// the total path of all the data



	///* the shear estimators data matrix data[i][j] */
	//double *data_m = new double[data_rows*shear_esti_data_cols]();
	//double **data = new double*[data_rows];
	//for (i = 0; i < data_rows; i++)
	//{
	//	data[i] = data_m + i*shear_esti_data_cols;
	//}

	///* the snr parameters data matrix data_snr[i][j] */
	//double *data_s = new double[data_rows*snr_para_data_cols]();
	//double **data_snr = new double*[data_rows];
	//for (i = 0; i < data_rows; i++)
	//{
	//	data_snr[i] = data_s + i*snr_para_data_cols;
	//}

	// initialize gsl
	int sss1, sss2;
	sss1 = 2430;
	sss2 = 130;
	seed = myid *sss1 + sss2;

	gsl_initialize(seed);


	st = clock();
	for (shear_id = 0; shear_id < 1; shear_id++)
	{
		for (i = 0; i < 1; i++)
		{
			s1 = clock();
			sprintf(chip_path, "%s%d/gal_chip_%04d.fits", data_path, shear_id, i);
			
			read_fits(chip_path, big_img);
			std::cout << chip_path << std::endl;
			for (j = 0; j < stamp_num; j++)
			{
				initialize_arr(gal, size*size, 0);
				initialize_para(&all_paras);
				initialize_arr(mask, size*size, 0);

				segment(big_img, gal, j, size, stamp_nx, stamp_nx);

				galaxy_finder(gal, mask, &all_paras, false, detect_label, detect_info);
				stack(check_img, gal, j, size, stamp_nx, stamp_nx);
				stack(mask_img, mask, j, size, stamp_nx, stamp_nx);
				//std::cout << detect_label << std::endl;
				if (detect_info == "Too many source!")
				{
					sprintf(h5_path, "!/home/hkli/work/mask.fits");
					write_fits(h5_path, mask, size, size);
					sprintf(h5_path, "!/home/hkli/work/gal.fits");
					write_fits(h5_path, gal, size, size);
					//show_arr(mask, size, size);
					std::cout << shear_id << " " << i << " " << j <<" "<< detect_info<< std::endl;
					break;
				}
			}
			sprintf(h5_path, "!/home/hkli/work/mask.fits");
			write_fits(h5_path, mask_img, stamp_nx*size, stamp_nx*size);
			sprintf(h5_path, "!/home/hkli/work/gal.fits");
			write_fits(h5_path, check_img, stamp_nx*size, stamp_nx*size);
			s2 = clock();
			sprintf(log_inform, "Thread: %d, chip: %d, done in %.2f s.", myid, i, (s2 - s1) / CLOCKS_PER_SEC);
			cout << log_inform << endl;

		}
	}
	//sprintf(set_name1, "/data");

	//sprintf(snr_h5_path, "%sdata_%d.hdf5", data_path, myid);
	//write_h5(snr_h5_path, set_name1, data_rows, snr_para_data_cols, data_snr[0], NULL);


	ed = clock();


	delete[] big_img;
	delete[] point;
	delete[] gal;
	delete[] gpow;
	delete[] psf;
	delete[] ppow;
	delete[] noise;
	delete[] pnoise;
	gsl_free();
	MPI_Finalize();
	return 0;
}
