#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include<hdf5.h>
#include<stdio.h>
#include<string>
#include "FQlib.h"

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
	char data_path[100], chip_path[150], snr_h5_path[150], para_path[150], buffer[200], h5_path[150], set_name[50], log_path[150], log_inform[250],coeff_path[50];
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/simu_test_64/");
	std::string str_data_path = "/mnt/ddnfs/data_users/hkli/simu_test_64/";
	std::string str_paraf_path = str_data_path + "parameters/para.ini";
	sprintf(log_path, "%slogs/m_%02d.dat", data_path, myid);

	int size, total_chips, chip_num, shear_pairs, data_row, total_data_row;
	int stamp_num = 10000, stamp_nx, shear_esti_data_cols = 7, snr_para_data_cols = 10;
	int i, j, k=0, row, row_s, seed, chip_id_s, chip_id_e, shear_id, temp_s=myid, detect_label, h;
	double psf_thres_scale = 2., sig_level = 2.0, psf_noise_sig = 0, gal_noise_sig, ts, te, t1, t2, psf_peak = 0, temp_flux = 0;;

	int cmd = 1;

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

	if (0 == myid)
	{	
		std::cout<<data_path<<std::endl;
		if (0 == cmd)
		{	
			std::cout << "OPERATION: detect & measure, SIG_LEVEL: " << sig_level << " sigma" << std::endl;
			
		}
		if (1 == cmd)
		{	
			std::cout << "OPERATION: detect , SIG_LEVEL: " << sig_level << " sigma" << std::endl;
		}
		std::cout << "Total chip: " << total_chips<<"Total cpus: "<<numprocs <<" Stamp size: "<<size <<std::endl;
		sprintf(log_inform, "RANK: %03d,  thread: %d, total cips: %d, individual chip: %d , size：%d, stamp_col: %d", myid, numprocs, total_chips, chip_num, size, stamp_nx);
		write_log(log_path, log_inform);
	}

	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = psf_noise_sig;
	all_paras.stamp_size = size;
	all_paras.max_source = 30;
	all_paras.area_thres = 5;
	all_paras.detect_thres = gal_noise_sig * sig_level;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = 5.5; // because the max half light radius of the galsim source is 5.5 pixels

	double *psf = new double[size*size]();
	double *ppsf = new double[size*size]();
	double *ppsf_cp = new double[size*size]();
	double *big_img = new double[size*size*stamp_num]();
	int *check_img = new int[size*size*stamp_num]();
	int *mask = new int[size*size]{};
	double *gal = new double[size*size]();
	double *pgal = new double[size*size]();
	double *noise = new double[size*size]();
	double *pnoise = new double[size*size]();
	double *mag = new double[total_data_row]();
	double *recvbuf, *recvbuf_s;
	double *data = new double[data_row*shear_esti_data_cols]();
	// the snr parameters data matrix data_snr[i][j]
	double *data_s = new double[data_row*snr_para_data_cols]();
	if (0 == myid)
	{
		recvbuf = new double[total_chips*stamp_num*shear_esti_data_cols];
		recvbuf_s = new double[total_chips*stamp_num*snr_para_data_cols];
	}

	//double *coeff = new double[3750]();
	//sprintf(coeff_path, "coeffs.hdf5");
	//sprintf(set_name, "/data");
	//read_h5(coeff_path, set_name, coeff, NULL, NULL, NULL, NULL);

	sprintf(chip_path, "%spsf.fits", data_path);
	read_fits(chip_path, psf);
	pow_spec(psf, ppsf, size, size);

	seed = 12300;
	gsl_rng_initialize(seed);	

	//addnoise(noise, size*size, gal_noise_sig);
	//pow_spec(noise, pnoise, size, size);	
	//smooth(ppsf, coeff, &all_paras);
	//smooth(pnoise, coeff, &all_paras);		
	//noise_subtraction(ppsf, pnoise, &all_paras, 1, 1);	
	//normalize_arr(ppsf, size);

	get_psf_radius(ppsf, &all_paras, psf_thres_scale);

	gsl_rng_free();

	if (0 == myid)
	{
		std::cout << "PSF THRES: " << all_paras.psf_pow_thres << std::endl << all_paras.psf_hlr << std::endl;
	}

	for (shear_id = 0; shear_id < shear_pairs; shear_id++)
	{
		ts = clock();
		sprintf(log_inform, "RANK: %03d, SHEAR %02d: my chips: %d - %d, total chips: %d (%d cpus)", myid, shear_id, chip_id_s, chip_id_e, total_chips, numprocs);
		write_log(log_path, log_inform);
		
		sprintf(para_path, "%sparameters/para_%d.hdf5", data_path, shear_id);
		sprintf(set_name, "/mag");
		read_h5(para_path, set_name, mag);
		
		if (0 == myid)
		{
			
			initialize_arr(recvbuf, total_chips*stamp_num*shear_esti_data_cols);
			initialize_arr(recvbuf_s, total_chips*stamp_num*snr_para_data_cols);
		}
		initialize_arr(data, data_row*shear_esti_data_cols);
		initialize_arr(data_s, data_row*snr_para_data_cols);


		for (i = chip_id_s; i < chip_id_e; i++)
		{
			seed = myid * i + shear_id + 1+ i + temp_s;
			temp_s++;
			gsl_rng_initialize(seed);
			t1 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: %04d 's chip start...", myid, shear_id, i);
			write_log(log_path, log_inform);
			if (0 == myid)
			{
				std::cout << log_inform << std::endl;
			}

			sprintf(chip_path, "%s%d/gal_chip_%04d.fits", data_path, shear_id, i);
			initialize_arr(big_img, stamp_nx*stamp_nx*size*size);
			initialize_arr(check_img, stamp_nx*stamp_nx*size*size);

			read_fits(chip_path, big_img);

			row = (i - chip_id_s) *stamp_num*shear_esti_data_cols;
			row_s = (i - chip_id_s) *stamp_num*snr_para_data_cols;

			for (j = 0; j < stamp_num; j++)			
			{				
				initialize_arr(noise, size*size);
				initialize_arr(pnoise, size*size);
				initialize_arr(gal, size*size);
				initialize_arr(pgal, size*size);
				initialize_para(&all_paras);
				initialize_arr(mask, size*size);

				segment(big_img, gal, j, size, stamp_nx, stamp_nx);
				pow_spec(gal, pgal, size, size);

				detect_label = galaxy_finder(gal, mask, &all_paras, false);
				if (i<chip_id_s+2 && 0 == myid && 0 == shear_id)
				{
					stack(check_img, mask, j, size, stamp_nx, stamp_nx);
				}
				snr_est(pgal, &all_paras, 2);

				if (cmd == 0)
				{
					addnoise(noise, size*size, gal_noise_sig);
					pow_spec(noise, pnoise, size, size);
					//smooth(pnoise, ppsf, coeff, &all_paras);
					//smooth(pgal, ppsf, coeff, &all_paras);
					
					noise_subtraction(pgal, pnoise, &all_paras, 1, 1);
					shear_est(pgal, ppsf, &all_paras);

					data[row + j * shear_esti_data_cols + 0] = i;
					data[row + j * shear_esti_data_cols + 1] = j;
					data[row + j * shear_esti_data_cols + 2] = all_paras.n1;
					data[row + j * shear_esti_data_cols + 3] = all_paras.n2;
					data[row + j * shear_esti_data_cols + 4] = all_paras.dn;
					data[row + j * shear_esti_data_cols + 5] = all_paras.du;
					data[row + j * shear_esti_data_cols + 6] = all_paras.dv;
				}

				data_s[row_s + j * snr_para_data_cols + 0] = all_paras.gal_flux2;
				data_s[row_s + j * snr_para_data_cols + 1] = all_paras.gal_flux_alt;
				data_s[row_s + j * snr_para_data_cols + 2] = all_paras.gal_flux;
				data_s[row_s + j * snr_para_data_cols + 3] = all_paras.gal_osnr;

				data_s[row_s + j * snr_para_data_cols + 4] = all_paras.gal_flux2_ext[0];
				data_s[row_s + j * snr_para_data_cols + 5] = all_paras.gal_flux2_ext[1];
				data_s[row_s + j * snr_para_data_cols + 6] = all_paras.gal_flux2_ext[2];
				data_s[row_s + j * snr_para_data_cols + 7] = all_paras.gal_flux2_ext[3];
				data_s[row_s + j * snr_para_data_cols + 8] = -mag[i*stamp_num + j];
				data_s[row_s + j * snr_para_data_cols + 9] = detect_label;

			 }		

			t2 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: %04d 's chip finish in %.2f sec", myid, shear_id, i, (t2 - t1) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (0 == myid)
			{
				std::cout << log_inform << std::endl;
			}
			gsl_rng_free();
		}
		MPI_Barrier(MPI_COMM_WORLD);
		sprintf(set_name, "/data");
		// the shear estimators
		if (0 == cmd)
		{
			MPI_Gather(data, data_row*shear_esti_data_cols, MPI_DOUBLE, recvbuf, data_row*shear_esti_data_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (0 == myid)
			{
				sprintf(h5_path, "%sresult/data/data_%d.hdf5", data_path, shear_id);
				write_h5(h5_path, set_name, recvbuf, total_data_row, shear_esti_data_cols);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// the SNR.. parameters
		MPI_Gather(data_s, data_row*snr_para_data_cols, MPI_DOUBLE, recvbuf_s, data_row*snr_para_data_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (0 == myid)
		{
			sprintf(snr_h5_path, "%sresult/data/data_%.1fsig/data_%d.hdf5", data_path, sig_level, shear_id);
			write_h5(snr_h5_path, set_name, recvbuf_s, total_data_row, snr_para_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		te = clock();
		sprintf(log_inform, "RANK: %03d, SHEAR %02d: write data to file and finish jobs in %.2f sec", myid, shear_id, (te - ts) / CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == myid)
		{
			std::cout << log_inform << std::endl;
		}
	}

	if (0 == myid)
	{
		delete[] recvbuf;
		delete[] recvbuf_s;
	}		
	delete[] psf;
	delete[] ppsf;
	delete[] ppsf_cp;
	delete[] big_img;
	delete[] mask;
	delete[] gal;
	delete[] pgal;
	delete[] noise;
	delete[] pnoise;
	delete[] data;
	delete[] data_s;
	delete[] mag;
	//delete[] coeff;
	delete[] check_img;
	MPI_Finalize();
	return 0;
}