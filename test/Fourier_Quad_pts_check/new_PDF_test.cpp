#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>

void task_alloc(const int total_task_num, const int division_num, const int my_part_id, int &my_st_id, int &my_ed_id, int *task_count)
{
    int i,j,m,n;
    m = total_task_num/division_num;
    n = total_task_num%division_num;

    for(i=0;i<division_num;i++)
    {
        task_count[i] = m;
        if(i<n){task_count[i] +=1;}
    }
    m=0;
    n=0;
    for(i=0;i<my_part_id;i++)
    {
        m+=task_count[i];
    }
    n = m+task_count[my_part_id];
    my_st_id = m;
    my_ed_id = n;
}



void set_radius_bin(double *radius_arr, const int data_num, double * bins, const int bin_num, const double max_scale)
{
    
	int i;
	int step, num;
	double data_max;
	
    sort_arr(radius_arr, data_num, 1);
	step = data_num / bin_num;

	// make the boundary big enough to enclose all the data
	bins[0] = 0;
	bins[bin_num] = radius_arr[data_num-1] * max_scale;
	for (i = 1; i < bin_num; i++)
	{
		bins[i] = radius_arr[step*i];
	}
}


void cal_radial_chisq(const double *tri_ratio, const int data_num, const double *radius_bin, const int bin_num, double &chisq)
{
    int i,j,k;
    double temp=0;

    double *tri_ratio_cumulate = new double[bin_num]{};
    int *num_count = new int[bin_num]{};

    for(i=0;i<data_num;i++)
    {   
        for(j=0;j<bin_num;j++)
        {
            if(radius_bin[j] <= tri_ratio[i] <radius_bin[j+1])
            {
                num_count[j] += 1;
                tri_ratio_cumulate[j] += tri_ratio[i];
            }
        }
    }
    show_arr(num_count,1,bin_num);
    show_arr(tri_ratio_cumulate,1,bin_num);

    for(j=0;j<bin_num;j++)
    {
        if(num_count[j]==0){std::cout<<"Radial chi squared divided by ZERO !!!"<<std::endl;throw "Radial chi squared divided by ZERO !!!";}
        else{temp += (tri_ratio_cumulate[j]*tri_ratio_cumulate[j])/num_count[j];}        
    }
    chisq = temp;
    delete[] tri_ratio_cumulate;
    delete[] num_count;
}

void find_shear_dipole(const double *mg1, const double *mg2, const double *mn, const int data_num, const int bin_num, 
                        double &gh1, double &gh1_sig,double &gh2, double &gh2_sig, double *chi_check_g1,double *chi_check_g2,const int chi_fit_num=20, 
                        const double left_guess=-0.1, const double right_guess=0.1, const double chi_gap=40)
{
    int i, j, k;
    int change, iters;
    double g1, g1_sig, g2, g2_sig;

    double left, right, step;
	double chi_left, chi_right, chi_mid, chi_left_mid, chi_right_mid;
	double gh_left, gh_right, gh_mid, gh_left_mid, gh_right_mid;

    double *gh_fit = new double[chi_fit_num];
    double *chisq_fit = new double[chi_fit_num];

    double mg1n, mg2n, temp;
    double *ratio = new double[data_num];
    double *radius = new double[data_num];
    double *radius_bin = new double[bin_num+1]{};

    for(i=0;i<data_num;i++)
    {
        radius[i] = sqrt(mg1[i]*mg1[i] + mg2[i]*mg2[i]);
    }

    set_radius_bin(radius, data_num, radius_bin, bin_num, 100);
    show_arr(radius_bin,1, bin_num+1);
    // find g1
    iters = 0;
    change = 1;
    left = left_guess;
    right = right_guess;
    while (change == 1)
	{		
		change = 0;
		gh_mid = (left + right) *0.5;
		gh_left = left;
		gh_right = right;

		try
		{
            for(i=0;i<data_num;i++)
            {
                temp = mn[i]*gh_left;
                mg1n = mg1[i] - temp;
                mg2n = mg2[i] - temp;
                temp = mg2n/mg1n;
                radius = sqrt()
                ratio[i] = 1./sqrt(1 + temp*temp);
            }
			cal_radial_chisq(radius, data_num, radius_bin, bin_num, chi_left);
            
            for(i=0;i<data_num;i++)
            {
                temp = mn[i]*gh_mid;
                mg1n = mg1[i] - temp;
                mg2n = mg2[i] - temp;
                ratio = mg2n/mg1n;
                radius[i] = 1./sqrt(1 + ratio*ratio);
            }
			cal_radial_chisq(radius, data_num, radius_bin, bin_num, chi_mid);

            for(i=0;i<data_num;i++)
            {
                temp = mn[i]*gh_right;
                mg1n = mg1[i] - temp;
                mg2n = mg2[i] - temp;
                ratio = mg2n/mg1n;
                radius[i] = 1./sqrt(1 + ratio*ratio);
            }
			cal_radial_chisq(radius, data_num, radius_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}

		//std::cout << left << " "<< gh_left<<" "<< gh_mid<<" "<< gh_right <<" "<< right << std::endl;

		if (chi_left > chi_mid + chi_gap)
		{
			left = (gh_mid + gh_left) *0.5;
			change = 1;
		}
		if (chi_right > chi_mid + chi_gap)
		{
			right = (gh_mid + gh_right)*0.5;
			change = 1;
		}

		iters += 1;
		if (iters > 13)
		{
			break;
		}
	}
	step = (right - left) / (chi_fit_num - 1);

	for (i = 0; i < chi_fit_num; i++)
	{	
        gh_fit[i] = left + step * i;
        for(j=0;j<data_num;j++)
        {
            temp = mn[j]*gh_fit[i];
            mg1n = mg1[i] - temp;
            mg2n = mg2[i] - temp;
            ratio = mg2n/mg1n;
            radius[i] = 1./sqrt(1 + ratio*ratio);
        }
		try
		{
			cal_radial_chisq(radius, data_num, radius_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
		if (chi_check_g1)
		{	// for checking
			chi_check_g1[i] = chi_right;
			chi_check_g1[chi_fit_num + i] = gh_fit[i];
		}
	}
    fit_shear(gh_fit, chisq_fit, chi_fit_num, gh1, gh1_sig, -1);

    // find g2
    iters = 0;
    change = 1;
    left = left_guess;
    right = right_guess;
    while (change == 1)
	{		
		change = 0;
		gh_mid = (left + right) *0.5;
		gh_left = left;
		gh_right = right;

		try
		{
            for(i=0;i<data_num;i++)
            {   
                temp = mn[i]*gh_left;
                mg1n = mg1[i] - temp;
                mg2n = mg2[i] - temp;
                ratio = mg1n/mg2n;
                radius[i] = 1./sqrt(1 + ratio*ratio);
            }
			cal_radial_chisq(radius, data_num, radius_bin, bin_num, chi_left);
            
            for(i=0;i<data_num;i++)
            {   
                temp = mn[i]*gh_mid;
                mg1n = mg1[i] - temp;
                mg2n = mg2[i] - temp;
                ratio = mg1n/mg2n;
                radius[i] = 1./sqrt(1 + ratio*ratio);
            }
			cal_radial_chisq(radius, data_num, radius_bin, bin_num, chi_mid);

            for(i=0;i<data_num;i++)
            {
                temp = mn[i]*gh_right;
                mg1n = mg1[i] - temp;
                mg2n = mg2[i] - temp;
                ratio = mg1n/mg2n;
                radius[i] = 1./sqrt(1 + ratio*ratio);
            }
			cal_radial_chisq(radius, data_num, radius_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}

		//std::cout << left << " "<< gh_left<<" "<< gh_mid<<" "<< gh_right <<" "<< right << std::endl;

		if (chi_left > chi_mid + chi_gap)
		{
			left = (gh_mid + gh_left) *0.5;
			change = 1;
		}
		if (chi_right > chi_mid + chi_gap)
		{
			right = (gh_mid + gh_right)*0.5;
			change = 1;
		}

		iters += 1;
		if (iters > 13)
		{
			break;
		}
	}

	step = (right - left) / (chi_fit_num - 1);

	for (i = 0; i < chi_fit_num; i++)
	{	
        gh_fit[i] = left + step * i;

        for(j=0;j<data_num;j++)
        {   
            temp = mn[j]*gh_fit[i];
            mg1n = mg1[i] - temp;
            mg2n = mg2[i] - temp;
            ratio = mg1n/mg2n;
            radius[i] = 1./sqrt(1 + ratio*ratio);
        }
		try
		{
			cal_radial_chisq(radius, data_num, radius_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
		if (chi_check_g2)
		{	// for checking
			chi_check_g2[i] = chi_right;
			chi_check_g2[chi_fit_num + i] = gh_fit[i];
		}
	}
    fit_shear(gh_fit, chisq_fit, chi_fit_num, gh2, gh2_sig, -1);

    delete[] radius;
    delete[] radius_bin;

       // int grid_num = 30;
    // double chisq_ij;
    // double *g1_grid = new double[grid_num*grid_num];
    // double *g2_grid = new double[grid_num*grid_num];
    // double *chisq = new double[grid_num*grid_num];
    // find_shear_mean(mg1, mn, data_num, g1, g1_sig, 1000, 100);
    // find_shear_mean(mg2, mn, data_num, g2, g2_sig, 1000, 100);

    // // 10 sigma width
    // for(i=0;i<grid_num;i++) // y G1
    // {   
    //     m = g1 - g1_sig*5 + g1_sig*10/(grid_num-1)*i;        

    //     for(k=0;k<data_num;k++){mg1n[k] = mg1[k] - m*mn[i];}

    //     for(j=0;j<grid_num;j++) // x G2
    //     {
    //         n = g2 - g2_sig*5 + g2_sig*10/(grid_num-1)*i;
    //         g1_grid[i*grid_num+j] = m;
    //         g2_grid[i*grid_num+j] = n;
            
    //         for(k=0;k<data_num;k++){mg2n[k] = mg2[k] - n*mn[i];}

    //         cal_radial_chisq(mg1n, mg2n, data_num, radius_bin, bin_num, chisq_ij);
    //         chisq[i*grid_num + j] = chisq_ij;
    //     }
    // }
    // delete[] g1_grid;
    // delete[] g2_grid;

}

int main(int argc, char**argv)
{   
    int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

    char parent_path[200], data_path[200], set_name[30];
    char inform[200];
    char data_nm[30];

    double *data, *mg1, *mg2, *mn, *mnu1, *mnu2;
    double *check = new double[40];

    double *check1 = new double[40];
    double *check2 = new double[40];

    strcpy(parent_path, argv[1]);
    strcpy(data_nm, argv[2]);

    int ratio;
    ratio = atoi(argv[3]);

    int i,j,k;
    int data_row = 20000000, sub_row;
    double g1, g1_sig, g2, g2_sig;

    sub_row = data_row/ratio;

    data = new double[data_row*4];
    mg1 = new double[sub_row*4];
    mg2 = new double[sub_row*4];
    mn = new double[sub_row*4];
    mnu1 = new double[sub_row*4];
    mnu2 = new double[sub_row*4];

    sprintf(data_path, "%s/shear.hdf5", parent_path);
    
    sprintf(set_name,"/g1");
    read_h5_datasize(data_path, set_name,k);
    double *shear = new double[k];
    read_h5(data_path, set_name, shear);
    show_arr(shear, 1, k);
    std::cout<<std::endl;

    sprintf(set_name,"/g2");
    read_h5_datasize(data_path, set_name,k);
    read_h5(data_path, set_name, shear);
    show_arr(shear, 1, k);
    std::cout<<std::endl;

    sprintf(data_path, "%s/%s",parent_path, data_nm);
    sprintf(set_name,"/data");
    read_h5(data_path, set_name, data);

    for(i=0;i<sub_row;i++)
    {
        mg1[i] = data[i*4];
        mg2[i] = data[i*4+1];
        mn[i] = data[i*4+2];
        mnu1[i] = data[i*4+2] + data[i*4+3];
        mnu2[i] = data[i*4+2] - data[i*4+3];
    }

    find_shear_mean(mg1, mn, sub_row, g1, g1_sig, 10, 1);
    find_shear_mean(mg2, mn, sub_row, g2, g2_sig, 10, 1);

    for(i=0;i<numprocs;i++)
    {
        if(i ==  rank)
        {   
            sprintf(inform, "Ave g1: %8.6f (%8.6f), g2: %8.6f (%8.6f)", g1, g1_sig, g2, g2_sig);
            std::cout<<inform<<std::endl;
        }
    }


    find_shear(mg1, mnu1, sub_row, 10, g1, g1_sig, check);
    find_shear(mg2, mnu2, sub_row, 10, g2, g2_sig, check);

    for(i=0;i<numprocs;i++)
    {
        if(i ==  rank)
        {   
            sprintf(inform, "PDF g1: %8.6f (%8.6f), g2: %8.6f (%8.6f)", g1, g1_sig, g2, g2_sig);
            std::cout<<inform<<std::endl;
        }
    }



    find_shear_dipole(mg1, mg2, mn, sub_row, 10, g1, g1_sig, g2, g2_sig, check1, check2);
    for(i=0;i<numprocs;i++)
    {
        if(i ==  rank)
        {   
            sprintf(inform, "PDF-dipole g1: %8.6f (%8.6f), g2: %8.6f (%8.6f)", g1, g1_sig, g2, g2_sig);
            std::cout<<inform<<std::endl;
        }
    }

    

    MPI_Finalize();
    return 0;

}