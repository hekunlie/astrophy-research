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
#include<string>

int main(int argc, char*argv[])
{
		int myid, numprocs, namelen;
		char processor_name[MPI_MAX_PROCESSOR_NAME];

		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Get_processor_name(processor_name, &namelen);

		para all_paras;

		std::ifstream fin;
		std::string s, str_stampsize = "stamp_size", str_total_num = "total_num", str_noise = "noise_sig", str_shear_num = "shear_num", str_nx = "stamp_col";
		char data_path[100], chip_path[150], snr_h5_path[150], para_path[150], shear_path[150], buffer[200], h5_path[150], set_1[50],set_2[50], log_path[150], log_inform[250];
		sprintf(data_path, "/mnt/ddnfs/data_users/hkli/selection_bias_64/");
		std::string str_data_path = "/mnt/ddnfs/data_users/hkli/selection_bias_64/";
		std::string str_paraf_path = str_data_path + "parameters/para.ini";
		sprintf(log_path, "%slogs/m_%02d.dat", data_path, myid);

		int num_p = 100, size, total_chips, chip_num, shear_pairs, data_row, total_data_row;
		int stamp_num = 10000, stamp_nx, shear_esti_data_cols = 7, snr_para_data_cols = 7;		
		int row, row_s, seed, chip_id_s, chip_id_e, shear_id, psf_type = 2;
		double max_radius=8, psf_scale=4., psf_thres_scale = 2., sig_level = 1.5, psf_noise_sig = 0, gal_noise_sig, psf_peak = 0;
		int i, j, k;
		double g1, g2, ts, te, t1, t2;

		read_para(str_paraf_path, str_stampsize, size);
		read_para(str_paraf_path, str_total_num, total_chips);
		read_para(str_paraf_path, str_noise, gal_noise_sig);
		read_para(str_paraf_path, str_shear_num, shear_pairs);
		read_para(str_paraf_path, str_nx, stamp_nx);

		chip_num = total_chips / numprocs;
		total_data_row = total_chips * stamp_num;
		data_row = chip_num * stamp_num;

		chip_id_s = chip_num * myid;
		chip_id_e = chip_num * (myid + 1);

		all_paras.gal_noise_sig = gal_noise_sig;
		all_paras.psf_noise_sig = psf_noise_sig;
		all_paras.stamp_size = size;
		all_paras.max_source = 30;
		all_paras.area_thres = 5;
		all_paras.detect_thres = gal_noise_sig * sig_level;
		all_paras.img_x = size;
		all_paras.img_y = size;
		all_paras.max_distance = 8; // because the max half light radius of the galsim source is 5.5 pixels

		initialize_para(&all_paras);

		double *big_img = new double[stamp_nx*stamp_nx*size*size]();
		double *point = new double[2 * num_p]();
		double *gal = new double[size*size]();
		double *pgal = new double[size*size]();
		double *psf = new double[size*size]();
		double *ppsf = new double[size*size]();
		double *noise = new double[size*size]();
		double *pnoise = new double[size*size]();
		double *shear = new double[2 * shear_pairs](); // [...g1,...,..g2,...]
		/* the shear estimators data matrix  */
		double *data = new double[data_row*shear_esti_data_cols]();
		/* the snr parameters data matrix */
		double *data_s = new double[data_row*snr_para_data_cols]();
		double *flux = new double[total_data_row];
		double *mag = new double[total_data_row];

		const char *str;
		sprintf(shear_path, "%sparameters/shear.dat", data_path);
		fin.open(shear_path);
		i = 0;
		while (!fin.eof())
		{
			getline(fin, s);
			str = s.c_str();
			shear[i] = atof(str);
			i++;
		}
		fin.close();

		#ifndef PRECISION
		DATA_TYPE *cp = new DATA_TYPE[stamp_nx*stamp_nx*size*size]();
		#endif

		// initialize gsl
		int sss1, sss2;
		sss1 = 2430;
		sss2 = 130;
		//sss1 = 5811430;
		//sss2 = 7161130;
		seed = myid *sss1 +sss2;
		//seed = myid *380 + 1401;// no bias
		gsl_rng_initialize(seed);


		// read parameters

		//PSF
		create_psf(psf, psf_scale, size, psf_type);

		pow_spec(psf, ppsf, size, size);

		get_psf_radius(ppsf, &all_paras, psf_thres_scale);

		for (shear_id = 0; shear_id < shear_pairs; shear_id++)
		{
			sprintf(para_path, "%sparameters/para_%d.hdf5", data_path, shear_id);
			sprintf(set_1, "/flux");
			sprintf(set_2, "/mag");
			read_h5(para_path, set_1, flux, set_2, mag, NULL, NULL);			

			ts = clock();

			g1 = shear[shear_id];
			g2 = shear[shear_id + shear_pairs];

			sprintf(log_path, "%slogs/%d_log.dat", data_path, myid );

			sprintf(para_inform, "myid: %03d, chip_id: %d, shear_id: %d", myid, chip_id, shear_id);
			write_log(log_path, para_inform);

			sprintf(para_inform, "size: %d, chip_num: %d,  point num: %d , noise sigma: %.2f, seed: myid*%d + %d, ", size, chip_num, num_p, gal_noise_sig, sss1, sss2);
			write_log(log_path, para_inform);

			sprintf(para_inform, "PSF scale: %.2f, max radius: %.2f", psf_scale, max_radius);
			write_log(log_path, para_inform);

			for (i = chip_id_s; i < chip_id_e; i++)
			{
				t1 = clock();

				sprintf(chip_path, "!%s%d/gal_chip_%04d.fits", data_path, shear_id, i);
				initialize_arr(big_img, stamp_nx*stamp_nx*size*size);

				sprintf(log_inform, "Thread: %d, chip: %d, start.", myid, i);
				write_log(log_path, log_inform);

				row = (i - chip_id_s) *stamp_num*shear_esti_data_cols;
				row_s = (i - chip_id_s) *stamp_num*snr_para_data_cols;

				for (j = 0; j < stamp_num; j++)
				{
					initialize_arr(gal, size*size);
					initialize_arr(pgal, size*size);
					initialize_arr(point, num_p * 2);
					initialize_arr(noise, size*size);
					initialize_arr(pnoise, size*size);
					initialize_para(&all_paras);

					create_points(point, num_p, max_radius);

					convolve(gal, point, flux[(chip_id + i)*stamp_num + j] / num_p, size, num_p, 0, psf_scale, g1, g2, psf_type);

					addnoise(gal, size*size, gal_noise_sig);

					stack(big_img, gal, j, size, stamp_nx, stamp_nx);

					galaxy_finder(gal, &all_paras, false);

					pow_spec(gal, pgal, size, size);

					snr_est(pgal, &all_paras, 1);

					addnoise(noise, size*size, gal_noise_sig);
					pow_spec(noise, pnoise, size, size);

					noise_subtraction(pgal, pnoise, &all_paras, 1, 1);
					shear_est(pgal, ppsf, &all_paras);

					data[row + j * shear_esti_data_cols + 0] = g1;
					data[row + j * shear_esti_data_cols + 1] = g2;
					data[row + j * shear_esti_data_cols + 2] = all_paras.n1;
					data[row + j * shear_esti_data_cols + 3] = all_paras.n2;
					data[row + j * shear_esti_data_cols + 4] = all_paras.dn;
					data[row + j * shear_esti_data_cols + 5] = all_paras.du;
					data[row + j * shear_esti_data_cols + 6] = all_paras.dv;

					data_s[row + j * snr_para_data_cols + 0] = all_paras.gal_osnr;
					data_s[row + j * snr_para_data_cols + 1] = all_paras.gal_flux;
					data_s[row + j * snr_para_data_cols + 2] = all_paras.gal_flux_alt;
					data_s[row + j * snr_para_data_cols + 3] = all_paras.gal_flux2;
					data_s[row + j * snr_para_data_cols + 4] = all_paras.gal_snr;
					data_s[row + j * snr_para_data_cols + 5] = all_paras.gal_size;
					data_s[row + j * snr_para_data_cols + 6] = mag[i*stamp_num + j];

				}

#ifdef PRECISION
				write_img(big_img, size*stamp_nx, size*stamp_nx, chip_path);
#else
				copy(big_img, big_img + stamp_nx * stamp_nx*size*size, cp);
				write_img(cp, size*stamp_nx, size*stamp_nx, chip_path);
#endif	

				t2 = clock();
				sprintf(log_inform, "Thread: %d, chip: %d, done in %.2f s.", myid, i, (t2 - t1) / CLOCKS_PER_SEC);
				write_log(log_path, log_inform);
				if (myid == 0)
				{
					std::cout << log_inform << std::endl;
				}
			}

			sprintf(set_1, "/data");

			sprintf(h5_path, "%sresult/data/data_%d_%d.hdf5", data_path, shear_id, myid / shear_pairs);
			write_h5(h5_path, set_1, data_row, shear_esti_data_cols, data, NULL);

			sprintf(snr_h5_path, "%sresult/data/data_2sig/data_%d_%d.hdf5", data_path, shear_id, myid / shear_pairs);
			write_h5(snr_h5_path, set_1, data_row_s, snr_para_data_cols, data_s, NULL);


			te = clock();
			if (myid == 0)
			{
				sprintf(buffer, "myid %d:  done in %g \n", myid, (te - ts) / CLOCKS_PER_SEC);
				std::cout << buffer<<std::endl;
			}
		}		

		delete[] big_img;
		delete[] point;
		delete[] gal;
		delete[] pgal;
		delete[] psf;
		delete[] ppsf;
		delete[] noise;
		delete[] pnoise;
		delete[] data;
		delete[] data_s;
		delete[] shear;
		delete[] flux;
		delete[] mag;
		gsl_rng_free();
		MPI_Finalize();
		return 0;
}
