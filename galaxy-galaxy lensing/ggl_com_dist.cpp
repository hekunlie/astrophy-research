#include<FQlib.h>
#include<mpi.h>

// assign the comoving distance to the data according to their redshifts

int main(int argc, char *argv[])
{
	int i, j, k;
	int area_id;
	char data_path[200], set_name[30], attrs_name[30], inform[200];
	
	int refer_num = 100001;
	double *redshift_refer = new double[refer_num];
	double *dist_refer = new double[refer_num];

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5");
	sprintf(set_name, "redshift");
	read_h5(data_path, set_name, redshift_refer);

	sprintf(set_name, "distance");
	read_h5(data_path, set_name, dist_refer);
	std::cout << redshift_refer[refer_num - 1] << " " << dist_refer[refer_num - 1] << std::endl;

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_cut.hdf5");
	sprintf(attrs_name, "shape");

	int shape[2];

	for (area_id = 1; area_id < 5; area_id++)
	{		
		sprintf(set_name, "/w_%d/Z", area_id);
		
		read_h5_attrs(data_path, set_name, attrs_name, shape, "d");
		double *redshift = new double[shape[0]];
		double *dist = new double[shape[0]];

		read_h5(data_path, set_name, redshift);

		std::cout << shape[0] << std::endl;
		for (i = 0; i < shape[0]; i++)
		{
			find_near(redshift_refer, redshift[i], refer_num, k);
			dist[i] = dist_refer[k];
			if (0 == i % 500000)
			{
				sprintf(inform, "Z: %7.5f, (%7.5f). Distance: %7.5f, %d", redshift[i], redshift_refer[k], dist_refer[k], k);
				std::cout << inform << std::endl;
			}
		}

		sprintf(set_name, "/w_%d/DISTANCE", area_id);
		write_h5(data_path, set_name, dist, shape[0], 1, FALSE);
		write_h5_attrs(data_path, set_name, attrs_name, shape, 2, "d");

		delete[] redshift;
		delete[] dist;

	}

	delete[] redshift_refer;
	delete[] dist_refer;

	return 0;
}
