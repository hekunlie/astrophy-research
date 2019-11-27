#include<FQlib.h>
#include<hk_mpi.h>

int main(int argc, char **argv)
{
    int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

    int i,j,k;
    int seed;
    int test_num = 30;
    double *rand_val = new double[test_num];
    seed = rank;

    // before initialization
    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {   
            std::cout<<"Rank "<<rank<<" Before setup, GSL_SETUP_LABEL="<<GSL_SETUP_LABEL<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    // first initialization
    gsl_initialize(seed*100+100,0);

    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {   
            std::cout<<"Rank "<<rank<<" After setup, GSL_SETUP_LABEL="<<GSL_SETUP_LABEL<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    for(i=0;i<test_num;i++)
    {
        rand_gauss(10, 0, rand_val[i], rng0);
    }
    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {
            show_arr(rand_val,1, test_num);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    gsl_free(0);

    MPI_Barrier(MPI_COMM_WORLD);


    // initialize the rng's and free them
    for(k=0;k<numprocs;k++)
    {
        if(k==rank)
        {   
            std::cout<<"==================== "<<rank<<" ======================="<<std::endl;
            for(i=0;i<3;i++)
            {   
                // rng0
                gsl_initialize(seed+1,0);
                if(rank == 0)
                {
                    std::cout<<"GSL_SETUP_LABEL="<<GSL_SETUP_LABEL<<std::endl;
                }
                for(i=0;i<test_num;i++)
                {
                    rand_gauss(10, 0, rand_val[i], rng0);
                }
                show_arr(rand_val, 1, test_num);
                gsl_free(0);

                // rng1
                gsl_initialize(seed+1,1);
                for(i=0;i<test_num;i++)
                {
                    rand_gauss(10, 0, rand_val[i], rng1);
                }
                show_arr(rand_val, 1, test_num);
                gsl_free(1);

                // rng2
                gsl_initialize(seed+1,2);
                for(i=0;i<test_num;i++)
                {
                    rand_gauss(10, 0, rand_val[i], rng2);
                }
                show_arr(rand_val, 1, test_num);
                gsl_free(2);

                // rng3
                gsl_initialize(seed+1,3);
                for(i=0;i<test_num;i++)
                {
                    rand_gauss(10, 0, rand_val[i], rng3);
                }
                show_arr(rand_val, 1, test_num);
                gsl_free(3);
            }
            std::cout<<"==========================================="<<std::endl<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {   
            std::cout<<"Rank "<<rank<<" At the end, GSL_SETUP_LABEL="<<GSL_SETUP_LABEL<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}