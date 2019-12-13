#include<FQlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>

int main(int argc, char **argv)
{
    int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
    
    char data_path[200], set_name[30];
    char inform[200];
    int i,j,k;
    int seed;
    int data_num;
    double sigma, signal;
    double gh, gh_sig;
    int chi_fit_num;
    double *chi_fit_check;
    double *result_data;
    int bin_num[12]{2,4,8, 12, 14, 16,20,32,64, 96,128,192};
    // data number, sigma of Gauss, mu (signal) of Gauss, seed for rand()
    data_num = atoi(argv[1]);
    sigma = atof(argv[2]);
    signal = atof(argv[3]);
    seed = atoi(argv[4]);
    chi_fit_num = 20;
    chi_fit_check = new double[chi_fit_num*2];

    gsl_initialize(seed, 0);
    gsl_initialize(seed, 1);

    MPI_Win win_result_data;
	MPI_Aint result_data_size;

	if (0 == rank)
	{   
        //[[true g], [measured g], [sigma], [Num*sigma^2], [bin_num],  [Num]], 6 rows
		MPI_Win_allocate_shared(6*numprocs* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &result_data, &win_result_data);
	}
	else
	{
		int dispu_total;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &result_data, &win_result_data);
		MPI_Win_shared_query(win_result_data, 0, &result_data_size, &dispu_total, &result_data);
	}

    double *mg = new double[data_num];
    double *mn = new double[data_num];
    double noise;
    for(i=0;i<data_num;i++)
    {
        rand_gauss(sigma, signal, mg[i], rng0);
        //rand_gauss(mg[i]/50, 0, noise, rng1);
        //mg[i] +=noise;
        mn[i] = 1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    try
    {
        find_shear(mg, mn, data_num, bin_num[rank], gh, gh_sig, chi_fit_check, chi_fit_num,0,100, -0.05,0.05);
    }
    catch(const char*img)
    {
        std::cout<<rank<<" "<<bin_num[rank]<<" PDF is going wrong "<<img<<std::endl;
    }
    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {
            sprintf(inform,"Rank %d. Data num: %d, Signal: %7.5f, Sigma: %.4f, Est: %10.8f(%10.8f), bin_num: %d",rank, data_num, signal, sigma, gh, gh_sig, bin_num[rank]);
            std::cout<<inform<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // sprintf(inform,"Rank %d. Data num: %d, Signal: %7.5f, Sigma: %.4f, Est: %10.8f(%10.8f), bin_num: %d",rank, data_num, signal, sigma, gh, gh_sig, bin_num[rank]);
    // std::cout<<inform<<std::endl;
    // exit(0);
    result_data[rank             ] = signal;
    result_data[rank +   numprocs] = gh;
    result_data[rank + 2*numprocs] = gh_sig;
    result_data[rank + 3*numprocs] = data_num * gh_sig* gh_sig;
    result_data[rank + 4*numprocs] = bin_num[rank];
    result_data[rank + 5*numprocs] = data_num;

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
    {
        sprintf(set_name,"/data");
        sprintf(data_path,"data.hdf5");
        write_h5(data_path, set_name, result_data, 6,numprocs,true);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    delete[] mg;
    delete[] mn;
    delete[] chi_fit_check;

    gsl_free(0);
    gsl_free(1);
    MPI_Win_free(&win_result_data);

    MPI_Finalize();
    return 0;
}