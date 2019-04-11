#include<FQlib.h>
#include<mpi.h>


int main(int argc, char *argv[])
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	char data_path[200], set_name[30], attrs_name[30];
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_cut.hdf5");
	sprintf(set_name, "Z");
	sprintf(attrs_name, "shape");
	int shape[2];
	read_h5_attrs(data_path, set_name, attrs_name, shape, "d");
	double *redshift = new double[shape[0]];
	double *dist = new double[shape[0]];

	read_h5(data_path, set_name, redshift);

	double *redshift_refer = new double[1001];
	double *dist_refer = new double[1001];

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5");
	sprintf(set_name, "redshift");
	read_h5(data_path, set_name, redshift_refer);

	sprintf(set_name, "distance");
	read_h5(data_path, set_name, dist_refer);

	int i, j, k;
	for (i = 0; i < shape[0]; i++)
	{
		find_near(redshift_refer, redshift[i], shape[0], k);
		dist[i] = dist_refer[k];
	}


	delete[] redshift;
	delete[] dist;
	delete[] redshift_refer;
	delete[] dist_refer;
	return 0;
}
