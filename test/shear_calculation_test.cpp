#include<FQlib.h>


void histogram_(double *data, const double *bins, int *num, const int data_num, const int bin_num)
{	
	// slower than the old version
	// initialize
	int i, j;
	int sl, sm, sr, ds, tag;

	initialize_arr(num, bin_num, 0);

	gsl_histogram * h = gsl_histogram_alloc(bin_num);
	gsl_histogram_set_ranges(h, bins, bin_num + 1);

	for (i = 0; i < data_num; i++)
	{
		gsl_histogram_increment(h, data[i]);	
	}
	for (i = 0; i < bin_num; i++)
	{
		num[i] = gsl_histogram_get(h, i);
	}
	gsl_histogram_free(h);
}

void histogram_new(double *data, const double *bins, int *num, const int data_num, const int bin_num)
{	// run slowly than the old version...
	// initialize
	int i, j;
	int sl, sm, sr, ds, tag;

	sort_arr(data, data_num, 1);
	initialize_arr(num, bin_num, 0);
	for (i = 0; i < data_num; i++)
	{
		//////////////////////
		//  run slowly due to this part
		sl = 0;
		sr = bin_num - 1;
		sm = int((sr - sl)*0.5);
		ds = bin_num;
		///////////////////////
		//while (sr > sl+2)
		while (ds > 3)
		{
			if (bins[sm] <= data[i])
			{
				sl = sm;
			}
			else
			{
				sr = sm;
			}
			sm = int((sr + sl)*0.5);
			ds = sr - sl + 1;
		}
		
		for (j = sl; j < sr+1; j++)
		{
			if (bins[j] <= data[i] and data[i] < bins[j + 1])
			{
				num[j] += 1;
				break;
			}
		}
		//std::cout << sl << " " << sr << " " << bins[j] <<" "<< data[i] <<" "<< bins[j + 1] << std::endl;
	}
}

void chisq_1d_(const int *hist_num, const int size, double &chi_sq)
{   
    // the size must be an even number
    int i,j;
    int mid = size/2;
    double chi_count = 0;
	double dn, sn;
    for(i=mid;i<size;i++)
    {
		dn = (hist_num[i] - hist_num[size - i - 1]);
		sn = (hist_num[i] + hist_num[size - i - 1]);
		chi_count += dn * dn / sn;
    }
    chi_sq = chi_count*0.5;
}


int main(int argc, char **argv)
{    
    int i, j, k;
	double st1, st2, st3, st4, st5, st6, sts[2];
    int bin_num = 8, gh_num = 50, shear_num=10;
    int file_num = argc-1;
    int data_num = 10000000, data_col=7;
    int data_stack_num = data_num*(argc-1); 
    int size = data_num*data_col;
    int size_s = data_stack_num*data_col;

    double chi_sq_1, chi_sq_2;
	double gh1, gh2, gh1_sig, gh2_sig;

	double *shear = new double[shear_num * 2];

	double *data = new double[size];
	double *data_stack = new double[size_s];

	double *mg1 = new double[data_stack_num];
	double *mg2 = new double[data_stack_num];
	double *mnu1 = new double[data_stack_num];
	double *mnu2 = new double[data_stack_num];
	double *dg1 = new double[data_stack_num];
	double *dg2 = new double[data_stack_num];

    double *chi_1 = new double[gh_num];
    double *chi_2 = new double[gh_num];

    double *mg_bin = new double[bin_num+1];

    int *num_in_bin_1 = new int[bin_num];
	int *num_in_bin_2 = new int[bin_num];

    double *gh = new double[gh_num];
    for(i=0;i<gh_num;i++)
    {
        gh[i] = -0.03 + 0.06/ gh_num *i;
    }

    int *file_label = new int[file_num];
    for(i=1;i<argc;i++)
    {
        file_label[i-1] = atoi(argv[i]);
    }

    char data_path[200], set_name[50];
	char log_inform[250];

	std::string data_path_s;
	std::stringstream str_1, str_2;

	data_path_s = "/mnt/ddnfs/data_users/hkli/selection_bias_64_1/parameters/shear.dat";
	read_text(data_path_s, shear, 2* shear_num);

    sprintf(set_name,"/data");
    for(i=0;i<argc-1;i++)
    {		
        sprintf(data_path,"/mnt/ddnfs/data_users/hkli/selection_bias_64_1/result/data/data_%d.hdf5",file_label[i]);
        read_h5(data_path, set_name, data);
		std::cout << file_label[i] << std::endl;
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
		st1 = clock();

        for(j=0;j<data_stack_num;j++)
        {
            dg1[j] = mg1[j] - gh[i]*mnu1[j];
            dg2[j] = mg2[j] - gh[i]*mnu2[j];
        }

		st2 = clock();


        histogram(dg1, mg_bin, num_in_bin_1, data_stack_num, bin_num);
		st3 = clock();
		//show_arr(num_in_bin, 1, bin_num);
		

        histogram(dg2, mg_bin, num_in_bin_2, data_stack_num, bin_num);
		//histogram_(dg2, mg_bin, num_in_bin_2, data_stack_num, bin_num);
		st4 = clock();

		chisq_1d(num_in_bin_1, bin_num, chi_sq_1);
		chi_1[i] = chi_sq_1;	
		st5 = clock();

        chisq_1d(num_in_bin_2, bin_num, chi_sq_2);
        chi_2[i] = chi_sq_2;	

		st6 = clock();

		sts[1] += (st4 - st3) / CLOCKS_PER_SEC;
		sts[0] += (st3 - st2) / CLOCKS_PER_SEC;
		sprintf(log_inform, "%d: chi_1: %.4f,  chi_2: %.4f.  t1: %.3f. t2: %.3f. t3: %.3f. t4: %.3f. t5: %.3f. ", i, chi_1[i], chi_2[i], 
			(st2 - st1) / CLOCKS_PER_SEC, (st3 - st2) / CLOCKS_PER_SEC,
			(st4 - st3) / CLOCKS_PER_SEC, (st5 - st4) / CLOCKS_PER_SEC, (st6 - st5) / CLOCKS_PER_SEC);
		std::cout <<log_inform << std::endl;
    }
	std::cout << sts[0] / gh_num << " " << sts[1]/gh_num<<" " <<(sts[0]   - sts[1]  )/sts[0] << std::endl;
	str_1 << "g1: ";
	str_2 << "g2: ";
	for (i = 0; i < argc - 1; i++)
	{
		str_1 << file_label[i] << " " << shear[file_label[i]] << ", ";
		str_2 << file_label[i] + shear_num << " " << shear[file_label[i] + shear_num] << ", ";
	}
	std::cout << str_1.str() << str_2.str() << std::endl;

	st1 = clock();
	fit_shear(gh, chi_1, gh_num, gh1, gh1_sig);
	st2 = clock();
	fit_shear(gh, chi_2, gh_num, gh2, gh2_sig);
	st3 = clock();
	sprintf(log_inform, "g1: %.5f (%.5f).  g2:  %.5f (%.5f).", gh1, gh1_sig, gh2, gh2_sig,(st2 - st1) / CLOCKS_PER_SEC, (st3 - st2) / CLOCKS_PER_SEC);
	std::cout << log_inform << std::endl;

    sprintf(set_name, "/chi_1");
    sprintf(data_path, "/home/hkli/work/test/chisq.hdf5");
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
    delete[] file_label;
    delete[] gh;
    delete[] num_in_bin_1;
	delete[] num_in_bin_2;
    delete[] mg_bin;
    delete[] chi_1;
    delete[] chi_2;

    return 0;
}