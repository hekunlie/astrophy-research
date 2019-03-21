#include<FQlib.h>

# define NUM 10

def chisq_1d(const double *hist_num, const int size, double &chi_sq)
{   
    // the size must be an even number
    int i,j;
    int mid = size/2;
    double chi_count = 0;
    for(i=mid;i<size;i++)
    {
        chi_count += (hist_num[i] - hist[size - i - 1])*(hist_num[i] - hist[size - i - 1])/(hist_num[i] + hist[size - i - 1]);
    }
    chi_sq = chi_count*0.5;
}


int mian(int argc, char *argv[])
{
    double g1[NUM], g2[NUM];
    int i, j, k;
    
    int bin_num = 6, gh_num = 100;
    int file_num = argc-1;
    int data_num = 10000000, data_col=7;
    int data_stack_num = data_num*(argc-1); 
    int size = data_num*data_col;
    int size_s = data_stack_num*data_col;

    double chi_sq;

    double *chi_1 = new double[gh_num];
    double *chi_2 = new double[gh_num];

    double *mg_bin = new double[bin_num+1];
    int *num_in_bin = new int[bin_num];

    double *gh = new double[gh_num];
    for(i=0;i<gh_num;i++)
    {
        gh[i] = -0.08 + 0.16/100*i;
    }

    int *file_label = new int[file_num];
    for(i=1;i<argc;i++)
    {
        file_label[i] = atoi(argv[i]);
    }

    double *data = new double[size];
    double *data_stack = new double[size_s]

    double *mg1 = new double[data_stack_num];
    double *mg2 = new double[data_stack_num];
    double *mnu1 = new double[data_stack_num];
    double *mnu2 = new double[data_stack_num];
    double *dg1 = new double[data_stack_num];
    double *dg2 = new double[data_stack_num];

    char data_path[200], set_name[50];

    sprintf(set_name,"/data");
    for(i=0;i<argc-1;i++)
    {
        sprintf(data_path,"/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/result/data_%d.hdf5",file_label[i]);
        read_h5(data_path, set_name, data);
        
        for(j=i*size;j<(i+1)*size;j++)   
        {
            data_stack[j] = data[j - i*size];
        }
    }

    for(i=0;i<data_stack_num;i++)
    {
        mg1[i] = data_stack[i*data_col + 2];
        mg2[i] = data_stack[i*data_col + 3];
        mnu1[i] = data_stack[i*data_col + 4] + data_stack[i*data_col + 5];
        mnu2[i] = data_stack[i*data_col + 4] - data_stack[i*data_col + 5];
    }
    
    set_bin(mg1, data_stack_num, mg_bin, bin_num, 100);

    for(i=0;i<gh_num;i++)
    {
        for(j=0;j<data_stack_num;j++)
        {
            dg1[j] = mg1[j] - gh[i]*mnu1[j];
            dg2[j] = mg1[j] - gh[i]*mnu2[j];
        }
        
        histogram(dg1, mg_bin, num_in_bin, data_stack_num, bin_num);
        chisq_1d(num_in_bin, bin_num, chi_sq);
        chi_1[i] = chi_sq;

        histogram(dg2, mg_bin, num_in_bin, data_stack_num, bin_num);
        chisq_1d(num_in_bin, bin_num, chi_sq);
        chi_2[i] = chi_sq;
    }

    sprintf(set_name, "/chi_1");
    sprintf(data_paht, "/home/hklee/work/test/chisq.hdf5");
    write_h5(data_path, set_name, chi_1, gh_num, 1, TRUE);

    sprintf(set_name, "/chi_2");
    write_h5(data_path, set_name, chi_2, gh_num, 1, FALSE);
    
    delete[] dg2;
    delete[] dg1;
    delete[] mnu2;
    delete[] mnu1;
    delete[] mg2;
    delete[] mg1;
    delete[] data_stack;
    delete[] data;
    delete[] file_lable;
    delete[] gh;
    delete[] num_in_bin;
    delete[] mg_bin;
    delete[] chi_1;
    delete[] chi_2;

    return 0;
}