#include<FQlib.h>

int main()
{
	/* calculate the comoving distance with omega_m0 = 0.31, omega_lambda0 = 0.69 */
	/* only calculate the integrate without the factor c/H0                                                */
	/* redshift: 0 ~ 10                                                                                                          */
	double omega_m = 0.31, omeg_lambda = 0.69;
	double precision = 0.00000001;
	int i, step = 1000;
	double z_step = 0.01;
	double distance;
	double *redshift = new double[step + 1]{};
	double *co_dist = new double[step + 1]{};
	for (i = 1; i < step + 1; i++)
	{
		redshift[i] = z_step * i;
		com_distance(0, redshift[i], precision, omega_m, omeg_lambda, distance);
		co_dist[i] = distance;
		std::cout << "Redshift: " << redshift[i] << "  Comoving distance: " << distance << std::endl;
	}
	
	char data_path[150];
	char set_name[50], attrs_name[50];
	int shape[2];
	double parameter[1];
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5");

	shape[0] = step + 1;
	shape[1] = 1;
	sprintf(attrs_name, "shape");

	sprintf(set_name, "/redshift");
	write_h5(data_path, set_name, redshift, step + 1, 1, TRUE);
	write_h5_attrs(data_path, set_name, attrs_name, shape, 2, "d");

	sprintf(attrs_name, "omega_m0");
	parameter[0] = omega_m;
	write_h5_attrs(data_path, set_name, attrs_name, parameter, 1, "d");
	sprintf(attrs_name, "omega_lam0");
	parameter[0] = omeg_lambda;
	write_h5_attrs(data_path, set_name, attrs_name, parameter, 1, "d");

	sprintf(set_name, "/distance");
	write_h5(data_path, set_name, co_dist, step + 1, 1, FALSE);
	write_h5_attrs(data_path, set_name, attrs_name, shape, 2, "d");
	
	return 0;
}