#include<FQlib.h>

int main(int argc, char **argv)
{
	int i, j, k;
	double z_max, dist, hubble_h, c0, precision;
	double parameters[4];
	double *redshift, *com_dist_py, *com_dist_c;
	char data_path[100], set_name[30], inform[150];
	int num;

	sprintf(data_path, "dist.hdf5");
	sprintf(set_name, "/Z");
	read_h5_datasize(data_path, set_name, num);

	redshift = new double[num];
	com_dist_c = new double[num];
	com_dist_py = new double[num];

	read_h5(data_path, set_name, redshift);
	sprintf(set_name, "/DISTANCE");
	read_h5(data_path, set_name, com_dist_py);

	sprintf(set_name, "/PARAMETER");
	read_h5(data_path, set_name, parameters);

	hubble_h = parameters[2];
	c0 = parameters[3];
	precision = 1.e-8;
	for (i = 0; i < num; i++)
	{
		com_distance(0, redshift[i], parameters[0], parameters[1], dist, precision);
		com_dist_c[i] = dist * 1000 * c0 / hubble_h;
		sprintf(inform, "Z: %.3f. Dist_py (Dist_cpp): %.4f (%.4f) Mpc. Diff: %.4f", redshift[i], com_dist_py[i], com_dist_c[i], com_dist_py[i] - com_dist_c[i]);
		std::cout << inform << std::endl;
	}
	sprintf(set_name, "/DISTANCE_C");
	write_h5(data_path, set_name, com_dist_c, num, 1, FALSE);


	return 0;
}