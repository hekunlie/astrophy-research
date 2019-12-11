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
    
    char inform[200];
    int i,j,k;
    int seed;
    int data_num;
    double sigma, signal;
    double gh, gh_sig;
    int chi_fit_num;
    double *chi_fit_check;

    int bin_num[16]{2,4,8,12,14, 16,20,30,40,60, 80,100,120,140,160,200};

    data_num = atoi(argv[1]);
    sigma = atof(argv[2]);
    signal = atof(argv[3]);
    seed = atoi(argv[4]);
    chi_fit_num = 20;
    chi_fit_check = new double[chi_fit_num*2];

    gsl_initialize(seed, 0);
    gsl_initialize(seed, 1);



    double *mg = new double[data_num];
    double *mn = new double[data_num];
    double noise;
    for(i=0;i<data_num;i++)
    {
        rand_gauss(sigma, signal, mg[i], rng0);
        rand_gauss(mg[i]/50, 0, noise, rng1);
        mg[i] +=noise;
        mn[i] = 1;
    }

    find_shear(mg, mn, data_num, bin_num[rank],gh, gh_sig, chi_fit_check, chi_fit_num);
    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {
            sprintf(inform,"Data num: %d, Signal: %7.5f, Sigma: %.4f, Est: %10.8f(%10.8f), bin_num: %d",data_num, signal, sigma, gh, gh_sig, bin_num[rank]);
            std::cout<<inform<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    delete[] mg;
    delete[] mn;
    delete[] chi_fit_check;

    gsl_free(0);
    gsl_free(1);

    MPI_Finalize();
    return 0;
}