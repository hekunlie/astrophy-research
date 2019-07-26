#include<FQlib.h>
#include<mpi.h>

int main(int argc, char ** argv)
{
	/* calculate the comoving distance [Mpc/h]																*/
	/*	argv[1]: Omega_m0, (Omega_lambda0 = 1 - Omega_m0)	 								  */
	/*  argv[2]: Z_MAX, 0~ Z_MAX, the distance to be calculated                                   */
	
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
	int num_step, m,n, my_num_i, my_num_e, my_num, total_num;
	double st1, st2;

	st1 = clock();

	// \Omega_m0 & \Omega_lanbda0
	omega_m = atof(argv[1]);
	omeg_lambda = 1 - omega_m;

	// 0 ~ z_max
	z_max = atof(argv[2]);
	// precision for integrate
	precision = 1.e-9;
	// point num [0, Z_max]
	total_num = 500000;
	z_step = z_max / (total_num -1);
	
	if (rank == 0)
	{
		std::cout <<"Z:  0 ~ "<<z_max<<" Z_step: "<<z_step<<".  Total num: " << total_num <<".  Precision:  "<<precision<<std::endl;
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
	
	num_step = total_num / numprocs;
	m = int(num_step / 10);

	if (rank == 0)
	{
		std::cout << "|";
	}
	for (i = 0; i < num_step; i++)
	{
		com_distance(0, redshift[rank + i* numprocs], omega_m, omeg_lambda, distance, precision);
		com_dist[rank + i* numprocs] = distance;
		if (i % m == 0 and rank ==0)
		{
			std::cout << "|#";			
		}
	}

	if (rank == 0)
	{
		std::cout << "||"<<std::endl;
		for (i=num_step*numprocs; i<total_num;i++)
		{
			com_distance(0, redshift[i], omega_m, omeg_lambda, distance, precision);
			com_dist[i] = distance;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

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
	
		for (j = 0; j < total_num; j++)
		{
			std::cout << "Redshift: " << redshift[j] << "  Comoving distance: " << com_dist[j] <<" Mpc/h"<< std::endl;
		}		
		get_time(end_time, 50);
		std::cout << (st2 - st1) / CLOCKS_PER_SEC << " " << end_time << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Win_free(&win_z_share);
	MPI_Win_free(&win_dist_share);

	MPI_Finalize();
	
	return 0;
}