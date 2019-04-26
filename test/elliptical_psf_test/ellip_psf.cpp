#include<FQlib.h>

int main(int argc, char **argv)
{
	para all_paras;

	int size, i, j, k;
	int seed = 123;
	int num_p, psf_scale;
	double theta,ellip;
	char data_path[200];
	double amplitude;

	size = atoi(argv[1]);
	num_p = atoi(argv[2]);
	psf_scale = atof(argv[3]);
	ellip = atof(argv[4]);
	theta = atof(argv[5]);
	amplitude = atof(argv[6]);
	seed = atoi(argv[7]);

	std::cout << "STAMP: " << size << std::endl;
	std::cout << "Number of Points: " << num_p << std::endl;
	std::cout << "PSF scale: " << psf_scale << std::endl;
	std::cout << "Ellipticity of PSF: " << ellip << std::endl;
	std::cout << "Position angle: " << theta/Pi*180 << std::endl;
	std::cout << "Amplitude: " << amplitude << std::endl;

	double *img = new double[size*size * 4]{};
	double *img_ = new double[size*size]{};
	double *pts = new double[2 * num_p]{};
	double flux = 0;

	gsl_initialize(seed);

	create_points(pts, num_p, 9);

	// elliptical psf
	create_psf(img_, psf_scale, size, ellip, theta, amplitude,2);
	stack(img, img_, 0, size, 2, 2);
	for (i = 0; i < size*size; i++)
	{
		flux += img_[i];
	}	
	std::cout << "PSF: " <<flux<< std::endl;
	flux = 0;
	initialize_arr(img_, size*size, 0);	
	// convolve with elliptical psf
	convolve(img_, pts, 1., size, num_p, 0, psf_scale, 0, 0, 2, 1, ellip, theta, amplitude, &all_paras);
	stack(img, img_, 1, size, 2, 2);
	for (i = 0; i < size*size; i++)
	{
		flux += img_[i];
	}
	std::cout << "Galaxy 1: " << flux<<std::endl;
	flux = 0;
	initialize_arr(img_, size*size, 0);

	// circle psf
	create_psf(img_, psf_scale, size, 2);
	stack(img, img_, 2, size, 2, 2);
	for (i = 0; i < size*size; i++)
	{
		flux += img_[i];
	}
	std::cout << "PSF: " << flux << std::endl;
	flux = 0;
	initialize_arr(img_, size*size, 0);
	//convolve with circle psf
	convolve(img_, pts, 1., size, num_p, 0, psf_scale, 0, 0, 2, 1, &all_paras);
	stack(img, img_, 3, size, 2, 2);
	for (i = 0; i < size*size; i++)
	{
		flux += img_[i];
	}
	std::cout << "Galaxy 2: " << flux << std::endl;

	sprintf(data_path, "!sample.fits");
	write_fits(data_path, img, size*2, size * 2);

	gsl_free();

	return 0;
}