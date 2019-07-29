#include<FQlib.h>

// assign the comoving distance to the data according to their redshifts
// argv[1]: directory to the data file

int main(int argc, char *argv[])
{
	int i, j, k;
	int data_num;
	char data_path[200], set_name[30], attrs_name[30], inform[200];
	char target_set_name[20];

	int refer_num;


	//sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5");
	sprintf(data_path, "/mnt/perc/hklee/CFHT/gg_lensing/data/redshift.hdf5");
	sprintf(set_name, "/Z");
	read_h5_datasize(data_path, set_name, refer_num);

	double *redshift_refer = new double[refer_num];
	double *dist_refer = new double[refer_num];
	double *dist_integ_refer = new double[refer_num];

	read_h5(data_path, set_name, redshift_refer);

	sprintf(set_name, "/DISTANCE");
	read_h5(data_path, set_name, dist_refer);

	sprintf(set_name, "/DISTANCE_INTEG");
	read_h5(data_path, set_name, dist_integ_refer);

	std::cout << redshift_refer[refer_num-1] << " " << dist_refer[refer_num-1] << " " << dist_integ_refer[refer_num-1]<<std::endl;

	sprintf(data_path, "%s",argv[1]);
	sprintf(target_set_name, "%s", argv[2]);

	if (file_exist(data_path))
	{
		sprintf(set_name, "%s/Z", target_set_name);
		read_h5_datasize(data_path, set_name, data_num);
		if (data_num > 0)
		{
			std::cout << "Data number: " << data_num << std::endl;
			double *redshift = new double[data_num];
			double *dist = new double[data_num];
			double *dist_integ = new double[data_num];

			read_h5(data_path, set_name, redshift);

			for (i = 0; i < data_num; i++)
			{
				find_near(redshift_refer, redshift[i], refer_num, k);
				dist[i] = dist_refer[k];
				dist_integ[i] = dist_integ_refer[k];

				if (data_num > 100)
				{
					if (0 == i % (int)(data_num / 20))
					{
						sprintf(inform, "Z: %.7f, (%.7f). Distance: %.7f, %.7f. %d", redshift[i], redshift_refer[k], dist_refer[k], dist_integ_refer[k], k);
						std::cout << inform << std::endl;
					}
				}
				else
				{
					sprintf(inform, "Z: %.7f, (%.7f). Distance: %.7f, %.7f. %d", redshift[i], redshift_refer[k], dist_refer[k], dist_integ_refer[k],k);
					std::cout << inform << std::endl;
				}
			}

			sprintf(set_name, "%s/DISTANCE", target_set_name);
			write_h5(data_path, set_name, dist, data_num, 1, FALSE);

			sprintf(set_name, "%s/DISTANCE_INTEG", target_set_name);
			write_h5(data_path, set_name, dist_integ, data_num, 1, FALSE);
			
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
