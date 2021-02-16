#include<FQlib.h>
#include<hk_iolib.h>

#define MY_FLOAT float
int main(int argc, char **argv)
{

	// fq_paras all_paras;
	fq_paras_float all_paras;


	int size, i, j, k;
	int seed = 123;
	int num_p, psf_scale;
	MY_FLOAT theta,ellip, max_radius;
	char data_path[200], set_name[50];
	MY_FLOAT amplitude;
	MY_FLOAT img_cent;

	size = atoi(argv[1]);
	num_p = atoi(argv[2]);
	psf_scale = atof(argv[3]);
	ellip = atof(argv[4]);
	theta = atof(argv[5])*Pi;
	amplitude = atof(argv[6]);
	seed = atoi(argv[7]);

	img_cent = size/2;
	max_radius = 6;

	std::cout << "STAMP: " << size << std::endl;
	std::cout << "Number of Points: " << num_p << std::endl;
	std::cout << "PSF scale: " << psf_scale << std::endl;
	std::cout << "Ellipticity of PSF: " << ellip << std::endl;
	std::cout << "Position angle: " << theta/Pi*180 << std::endl;
	std::cout << "Amplitude: " << amplitude << std::endl;

	MY_FLOAT *img = new MY_FLOAT[size*size * 4]{};
	MY_FLOAT *img_ = new MY_FLOAT[size*size]{};
	MY_FLOAT *pts = new MY_FLOAT[2 * num_p]{};
	MY_FLOAT flux = 0;

	gsl_initialize(seed,0);

	sprintf(data_path,"./img.hdf5");

	create_points(pts, num_p, max_radius, 1, rng0);

	create_psf(img_, psf_scale, size, img_cent, 2);
	flux = 0;
	for (i = 0; i < size*size; i++)
	{
		flux += img_[i];
	}	
	std::cout << "PSF: " <<flux<< std::endl;
	sprintf(set_name, "/cpsf");
	write_h5(data_path, set_name, img_, size, size, true);
	
	convolve(pts, num_p, amplitude, 0,0, img_, size, img_cent, psf_scale, 2);
	sprintf(set_name, "/gal_cpsf");
	write_h5(data_path, set_name, img_, size, size, false);

	// elliptical psf
	// create_psf(img_, psf_scale, size, ellip, theta, amplitude,2);
	initialize_arr(img, size*size*4,0);
	for(i=0; i<4; i++)
	{	
		create_psf_e(img_, psf_scale, size, img_cent, ellip, Pi/4*i+theta, 2);

		stack(img, img_, i, size, 1, 4);

		flux = 0;
		for (j = 0; j < size*size; j++)
		{
			flux += img_[j];
		}	
		std::cout << "PSF: " <<flux<< std::endl;
		sprintf(set_name, "/epsf_%d",i);
		write_h5(data_path, set_name, img_, size, size, false);
	}
	
	
	sprintf(set_name, "/epsf");
	write_h5(data_path, set_name, img, size, size*4, false);

	for(i=0; i<4; i++)
	{
		create_psf_e(img_, psf_scale, size, img_cent, ellip, Pi/4*i, 2);
		convolve_e(pts, num_p, amplitude, 0,0, img_, size, img_cent, psf_scale, 2, ellip, Pi/4*i+theta);
		stack(img, img_, i, size, 1, 4);
	}
	sprintf(set_name, "/gal_epsf");
	write_h5(data_path, set_name, img, size, size*4, false);

	sprintf(set_name, "/pts");
	write_h5(data_path, set_name, pts, 2, num_p, false);






	// // convolve with elliptical psf
	// convolve(img_, pts, 1., size, num_p, 0, psf_scale, 0, 0, 2, 1, ellip, theta, amplitude, &all_paras);
	// stack(img, img_, 1, size, 2, 2);
	// for (i = 0; i < size*size; i++)
	// {
	// 	flux += img_[i];
	// }
	// std::cout << "Galaxy 1: " << flux<<std::endl;
	// flux = 0;
	// initialize_arr(img_, size*size, 0);

	// // circle psf
	// create_psf(img_, psf_scale, size, 2);
	// stack(img, img_, 2, size, 2, 2);
	// for (i = 0; i < size*size; i++)
	// {
	// 	flux += img_[i];
	// }
	// std::cout << "PSF: " << flux << std::endl;
	// flux = 0;
	// initialize_arr(img_, size*size, 0);
	// //convolve with circle psf
	// convolve(img_, pts, 1., size, num_p, 0, psf_scale, 0, 0, 2, 1, &all_paras);
	// stack(img, img_, 3, size, 2, 2);
	// for (i = 0; i < size*size; i++)
	// {
	// 	flux += img_[i];
	// }
	// std::cout << "Galaxy 2: " << flux << std::endl;

	// sprintf(data_path, "!sample.fits");
	// write_fits(data_path, img, size*2, size * 2);

	gsl_free(0);

	return 0;
}