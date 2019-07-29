#include<FQlib.h>

void com_dist_old(const double low_z, const double high_z, const double omg_m, const double omg_lam, double &result, const double precision_thresh)
{
	int i, j, bin_num;
	int max_run = 500;
	double result_1, result_2, d_res;
	double dz, z1, z2;

	bin_num = 20;
	result_1 = 0;
	for (i = 0; i < max_run; i++)
	{
		dz = (high_z - low_z) / (bin_num-1);
		result_2 = 0;
		for (j = 0; j < bin_num; j++)
		{
			z1 = low_z + j * dz;
			z2 = z1 + dz;
			result_2 = result_2 + 1. / sqrt(pow(1 + z1, 3)*omg_m + omg_lam) + 1. / sqrt(pow(1 + z2, 3)*omg_m + omg_lam);
			//result_2 = result_2 + 1. / sqrt(z1*omg_m + pow(z1, 4)*omg_lam) + 1. / sqrt(z2*omg_m + pow(z2, 4)*omg_lam);
		}
		result_2 = result_2 * dz*0.5;

		if (fabs(result_2 - result_1) <= precision_thresh)
		{
			// comoving distance  [ Mpc/h ]
			result = (result_2 + result_1) *0.5 * 1000 * C_0_hat;
			break;
		}
		else
		{
			result_1 = result_2;
			bin_num *= 2;
		}
	}
	if (i == max_run)
	{
		std::cout << "Max iteration, doesn't converge, " << result_2 - result_1 << "." << std::endl;
		exit(0);
	}
}
void com_dist_new(const double low_z, const double high_z, const double omg_m, const double omg_lam, double &result, const double precision_thresh = 1.e-8, const bool integ_only = FALSE)
{
	// Int f(x) = (f_1+f_2)/2*dx + (f_2+f_3)/2*dx + (f_3+f_4)/2*dx...
	//				= (f_1 + 2*f_2 + 2*f_3 + .. 2*f_{n-1} + f_n)/2*dx
	//				= ((f_1+f_n)/2 + f_2 + f_3 + .. + f_{n-1})*dx
	int i, j, bin_num, bin_num_0;
	int max_iters = 50, iters=0;
	double res_1, res_2, d_res;
	double dz, z;
	double scale = 1000 * C_0_hat;

	// run first time, if it converges, it returns the result
	// else, do the while block
	bin_num = 20;
	dz = (high_z - low_z) / (bin_num - 1);
	res_1 = 0;
	// (f_1 + f_n)/2
	res_2 = (1. / sqrt(pow(1 + low_z, 3)*omg_m + omg_lam) + 1. / sqrt(pow(1 + high_z, 3)*omg_m + omg_lam))*0.5;
	// add the median value
	for (j = 1; j < bin_num - 1; j++)
	{
		z = low_z + j * dz;
		res_2 = res_2 + 1. / sqrt(pow(1 + z, 3)*omg_m + omg_lam);
	}	
	while (TRUE)
	{
		d_res = fabs((res_2* dz - res_1)*scale);
		if (d_res < precision_thresh)
		{
			break;
		}
		if (iters > max_iters)
		{
			std::cout << "com_dist() Max iteration, doesn't converge, " << d_res << "( " << precision_thresh << " )." << std::endl;
			exit(0);
		}
		res_1 = res_2*dz;

		dz = dz *0.5;

		for (j = 0; j < bin_num - 1; j++)
		{
			z = low_z + (2 * j + 1) * dz;
			res_2 = res_2 + 1. / sqrt(pow(1 + z, 3)*omg_m + omg_lam);
		}			
		bin_num = bin_num + bin_num - 1;
		iters++;
	}
	if (integ_only)
	{
		result = res_2 * dz;
	}
	else
	{
		result = res_2 * dz*scale;
	}
}


int main(int argc, char**argv)
{

	double low_z, high_z;
	double omg_m0, omg_lam0;
	double dist1, dist2;
	double dist_precision;
	double t1, t2, t3, dt1, dt2;
	int i, times=1;

	char inform[200];


	low_z = atof(argv[1]);
	high_z = atof(argv[2]);
	omg_m0 = atof(argv[3]);
	omg_lam0 = 1 - omg_m0;

	dist_precision = atof(argv[4]);

	sprintf(inform, "Z: [%.4f, %.4f]. Omega_M0: %.4f, Omega_Lam0: %.4f.", low_z, high_z, omg_m0, omg_lam0);
	std::cout << inform << std::endl;

	com_distance(low_z, high_z, omg_m0, omg_lam0, dist1, dist_precision);

	sprintf(inform, "Dist_out: %.5f Mpc/h. Time: %.4f", dist1, dt1);
	std::cout << inform << std::endl;

	int num = 100;
	double *distance, *redshift;
	distance = new double[num];
	redshift = new double[num];

	for (i = 0; i < num; i++)
	{
		redshift[i] = (10. - 0) / (num - 1)*i;
		com_distance(0, redshift[i], omg_m0, omg_lam0, distance[i], dist_precision);
		sprintf(inform, "Z: %.6f. Dist: %.6f Mpc/h", redshift[i], distance[i]);
		std::cout << inform << std::endl;
	}
	char filename[100], set_name[40];
	sprintf(filename, "test.hdf5");
	sprintf(set_name, "/Z");
	write_h5(filename, set_name, redshift, num, 1, TRUE);
	sprintf(set_name, "/DISTANCE");
	write_h5(filename, set_name, distance, num, 1, FALSE);
	/*t1 = clock();
	for (i = 0; i < times; i++)
	{
		com_distance(low_z, high_z, omg_m0, omg_lam0, dist1, dist_precision);
	}
	t2 = clock();
	dt1 = (t2 - t1) / CLOCKS_PER_SEC;
	sprintf(inform, "Dist_out: %.5f Mpc/h. Time: %.4f",  dist1, dt1);
	std::cout << inform << std::endl;

	for (i = 0; i < times; i++)
	{
		com_dist_new(low_z, high_z, omg_m0, omg_lam0, dist2, dist_precision);
	}
	t3 = clock();
	dt2 = (t3 - t2) / CLOCKS_PER_SEC;
	sprintf(inform, "Dist_in: %.5f Mpc/h. Time: %.4f",  dist2, dt2);
	std::cout << inform << std::endl;
*/

	return 0;
}