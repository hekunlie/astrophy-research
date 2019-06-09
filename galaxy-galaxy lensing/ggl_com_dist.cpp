#include<FQlib.h>

// assign the comoving distance to the data according to their redshifts
// argv[1]: directory to the data file

int main(int argc, char *argv[])
{
	int i, j, k;
	int data_num;
	char data_path[200], set_name[30], attrs_name[30], inform[200];
	char target_set_name[20];

	int refer_num = 100001;
	double *redshift_refer = new double[refer_num];
	double *dist_refer = new double[refer_num];

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5");
	sprintf(set_name, "redshift");
	read_h5(data_path, set_name, redshift_refer);

	sprintf(set_name, "distance");
	read_h5(data_path, set_name, dist_refer);
	std::cout << redshift_refer[refer_num - 1] << " " << dist_refer[refer_num - 1] << std::endl;

	strcpy(data_path, argv[1]);
	strcpy(target_set_name, argv[2]);

	if (file_exist(data_path))
	{
		sprintf(set_name, "%s/Z", target_set_name);
		read_h5_datasize(data_path, set_name, data_num);
		if (data_num > 0)
		{
			std::cout << "Data number: " << data_num << std::endl;
			double *redshift = new double[data_num];
			double *dist = new double[data_num];

			read_h5(data_path, set_name, redshift);

			for (i = 0; i < data_num; i++)
			{
				find_near(redshift_refer, redshift[i], refer_num, k);
				dist[i] = dist_refer[k];

				if (data_num > 100)
				{
					if (0 == i % (int)(data_num / 10))
					{
						sprintf(inform, "Z: %7.5f, (%7.5f). Distance: %7.5f, %d", redshift[i], redshift_refer[k], dist_refer[k], k);
						std::cout << inform << std::endl;
					}
				}
				else
				{
					sprintf(inform, "Z: %7.5f, (%7.5f). Distance: %7.5f, %d", redshift[i], redshift_refer[k], dist_refer[k], k);
					std::cout << inform << std::endl;
				}
			}

			sprintf(set_name, "%s/DISTANCE", target_set_name);
			write_h5(data_path, set_name, dist, data_num, 1, FALSE);
			
			delete[] redshift;
			delete[] dist;
		}
		else
		{
			std::cout << "Wrong data number: " << data_num << std::endl;
			exit(0);
		}
	}
	else
	{
		std::cout << "Can't find the file. " << data_path << std::endl;
	}


	delete[] redshift_refer;
	delete[] dist_refer;

	return 0;
}
