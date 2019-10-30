#include<FQlib.h>

// assign the comoving/physical distance to the data according to their redshifts
// argv[1]: directory to the data file

int main(int argc, char *argv[])
{
	int i, j, k;
	int data_num;
	char data_path[200], set_name[30], inform[300];
	char target_set_name[20];

	int refer_num;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// read reference data ///////////////////////////////////////////////////////////////

	//sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5");
	sprintf(data_path, "/mnt/perc/hklee/CFHT/gg_lensing/data/redshift.hdf5");
	sprintf(set_name, "/Z");
	read_h5_datasize(data_path, set_name, refer_num);

	double *redshift_refer = new double[refer_num];
	// comoving distance
	double *co_dist_refer = new double[refer_num];
	// angular diameter distance
	double *phy_dist_refer = new double[refer_num];
	double *dist_integ_refer = new double[refer_num];

	read_h5(data_path, set_name, redshift_refer);

	sprintf(set_name, "/COM_DISTANCE");
	read_h5(data_path, set_name, co_dist_refer);

	sprintf(set_name, "/PHY_DISTANCE");
	read_h5(data_path, set_name, phy_dist_refer);

	sprintf(set_name, "/DISTANCE_INTEG");
	read_h5(data_path, set_name, dist_integ_refer);

	std::cout << redshift_refer[refer_num-1] << " " << co_dist_refer[refer_num-1] <<
	 phy_dist_refer[refer_num-1] << " " << dist_integ_refer[refer_num-1]<<std::endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// read the target catalog
	sprintf(data_path, "%s",argv[1]);
	// the distance will be signed to the dataset like /target_set_name/COM_DISTANCE ...
	sprintf(target_set_name, "%s", argv[2]);

	if (file_exist(data_path))
	{
		// read the target catalog
		sprintf(set_name, "%s/Z", target_set_name);
		read_h5_datasize(data_path, set_name, data_num);
		if (data_num > 0)
		{
			std::cout << "Data number: " << data_num << std::endl;
			double *redshift = new double[data_num];

			double *co_dist = new double[data_num];
			double *phy_dist = new double[data_num];
			double *dist_integ = new double[data_num];

			read_h5(data_path, set_name, redshift);

			for (i = 0; i < data_num; i++)
			{
				find_near(redshift_refer, redshift[i], refer_num, k);
				co_dist[i] = co_dist_refer[k];
				phy_dist[i] = phy_dist_refer[k];
				dist_integ[i] = dist_integ_refer[k];

				if (data_num > 100)
				{
					if (0 == i % (int)(data_num / 20))
					{
						sprintf(inform, "Z: %.5f, (%.5f). Distance: %.6f (%.6f) Mpc/h, %.5f. %d", redshift[i], redshift_refer[k], 
						co_dist[i], phy_dist[i], dist_integ[i], k);
						std::cout << inform << std::endl;
					}
				}
				else
				{
					sprintf(inform, "Z: %.5f, (%.5f). Distance: %.6f (%.6f) Mpc/h, %.5f. %d", redshift[i], redshift_refer[k], 
					co_dist[i], phy_dist[i], dist_integ[i], k);
					std::cout << inform << std::endl;
				}
			}

			sprintf(set_name, "%s/COM_DISTANCE", target_set_name);
			write_h5(data_path, set_name, co_dist, data_num, 1, FALSE);
			
			sprintf(set_name, "%s/PHY_DISTANCE", target_set_name);
			write_h5(data_path, set_name, phy_dist, data_num, 1, FALSE);

			sprintf(set_name, "%s/DISTANCE_INTEG", target_set_name);
			write_h5(data_path, set_name, dist_integ, data_num, 1, FALSE);

			delete[] redshift;
			delete[] co_dist;
			delete[] phy_dist;
			delete[] dist_integ;

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
	delete[] co_dist_refer;
	delete[] phy_dist_refer;
	delete[] dist_integ_refer;

	return 0;
}
