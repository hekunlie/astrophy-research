#include<FQlib.h>
#include<mpi.h>

int main(int argc, char *argv[])
{
	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	gsl_initialize(123);


	int i, j, k;
	int sky_areas;

	char data_path[100], log_inform[250], set_name[50], attrs_name[50];
	char data_set_path[100];

	sprintf(data_path, "/home/hkli/work/cpp/check.hdf5");

	int attrs_len = 2;
	double dshape[2]{};
	float fshape[2]{};
	int ishape[2]{};
	//sprintf(set_name, "/data");
	//read_h5(data_path, set_name, dshape);
	//std::cout << dshape[0] << ", " << dshape[1] << std::endl;
	sprintf(attrs_name, "double");
	printf(set_name, "/");
	read_h5_attrs(data_path, set_name, attrs_name, dshape);
	for (i = 0; i < attrs_len; i++)
	{
		std::cout << dshape[i] << std::endl;
	}

	/*sprintf(attrs_name, "float");
	read_h5_attrs(data_path, set_name, attrs_name, fshape);
	for (i = 0; i < attrs_len; i++)
	{
		std::cout << fshape[i] << std::endl;
	}

	sprintf(attrs_name, "int");
	read_h5_attrs(data_path, set_name, attrs_name, ishape);
	for (i = 0; i < attrs_len; i++)
	{
		std::cout << ishape[i] << std::endl;
	}*/
	// loop the sky areas
	for (i = 0; i < 1; i++)
	{
		;


	}



	gsl_free();
	MPI_Finalize();
	return 0;
}