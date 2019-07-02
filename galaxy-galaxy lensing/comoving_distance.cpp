#include<FQlib.h>
#include<mpi.h>

int main(int argc, char ** argv)
{
	/* calculate the comoving distance with omega_m0 = 0.31, omega_lambda0 = 0.69 */
	/* only calculate the integrate without the factor c/H0                                                */
	/* redshift: 0 ~ 10                                                                                                          */

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);


	double omega_m, omeg_lambda, z_max;
	double distance, precision, z_step;
	double *redshift, *com_dist, *my_com_dist;
	int i, j, step;	
	int my_num_i, my_num_e, my_num, total_num;
	double st1, st2;

	st1 = clock();

	// \Omega_m0 & \Omega_lanbda0
	omega_m = 0.268;
	omeg_lambda = 0.732;

	// 0 ~ z_max
	z_max = 10;
	// precision for integrate
	precision = 1.e-7;
	// point num [0, 1_max]
	step = 500000;
	z_step = z_max / step;

	total_num = step + 1;
	
	if (rank == 0)
	{
		std::cout <<"Z:  0 ~ "<<z_max<<" Z_step: "<<z_step<<".  Total num: " << total_num <<".  Precision:  "<<precision<<std::endl;
	}

	if (numprocs - 1 == rank)
	{
		my_num_i = total_num / numprocs * rank;
		my_num_e = total_num / numprocs * (rank + 1) + total_num % numprocs;
	}
	else
	{
		my_num_i = total_num / numprocs * rank;
		my_num_e = total_num / numprocs * (rank + 1);
	}
	for (i = 0; i < numprocs; i++)
	{
		if(rank == i)
		{
			std::cout << "RANK: " << rank << " " << my_num_i << " ~ " << my_num_e << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	


	MPI_Win win_dist_share, win_z_share;
	MPI_Aint  point_num;

	point_num = total_num;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(point_num * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &redshift, &win_z_share);
		MPI_Win_allocate_shared(point_num * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &com_dist, &win_dist_share);
	}
	else
	{
		int dispu_total;

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &redshift, &win_z_share);
		MPI_Win_shared_query(win_z_share, 0, &point_num, &dispu_total, &redshift);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &com_dist, &win_dist_share);
		MPI_Win_shared_query(win_dist_share, 0, &point_num, &dispu_total, &com_dist);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
	{
		for (i = 0; i < total_num; i++)
		{
			redshift[i] = i * z_step;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	my_num = my_num_e - my_num_i;
	
	my_com_dist = new double[my_num];

	for (i = my_num_i; i < my_num_e; i++)
	{
		com_distance(0, redshift[i], precision, omega_m, omeg_lambda, distance);
		my_com_dist[i - my_num_i] = distance;
	}
	
	for (i = 0; i < numprocs; i++)
	{
		if (rank == i)
		{
			for (j = 0; j < my_num; j++)
			{
				com_dist[j + my_num_i] = my_com_dist[j];
				std::cout << "Redshift: " << redshift[j + my_num_i] << "  Comoving distance: " << my_com_dist[j] << std::endl;
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	st2 = clock();

	if (rank == 0)
	{
		char data_path[150];
		char set_name[50], end_time[50];
		double parameter[1];
		//sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5");
		sprintf(data_path, "/mnt/perc/hklee/CFHT/gg_lensing/data/redshift.hdf5");


		sprintf(set_name, "/Z");
		write_h5(data_path, set_name, redshift, total_num, 1, TRUE);

		sprintf(set_name, "/OMEGA_M0");
		parameter[0] = omega_m;
		write_h5(data_path, set_name, parameter, 1, 1, FALSE);

		sprintf(set_name, "/OMEGA_LAMBDA0");
		parameter[0] = omeg_lambda;
		write_h5(data_path, set_name, parameter, 1, 1, FALSE);

		sprintf(set_name, "/DISTANCE");
		write_h5(data_path, set_name, com_dist, total_num, 1, FALSE);

		get_time(end_time, 50);
		std::cout << (st2 - st1) / CLOCKS_PER_SEC << " " << end_time << std::endl;
	}

	MPI_Win_free(&win_z_share);
	MPI_Win_free(&win_dist_share);

	delete[] my_com_dist;

	MPI_Finalize();
	
	return 0;
}