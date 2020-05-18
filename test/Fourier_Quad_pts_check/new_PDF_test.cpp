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

//,const int chi_grid_num, double *g1_grid, double* g2_grid, double* chisq_grid
void find_shear_2d(const double *mg1, const double* mg2, const double *mn, const double *mu, const int data_num, const int bin_num, 
                        double &gh1, double &gh1_sig,double &gh2, double &gh2_sig, const int chi_grid_num,double *g1_grid, double *g2_grid,
                        double *chisq1_grid, double *chisq2_grid, double  g1_true, double  g2_true, double delta_g,double* coeff_check)

{
    int i, j, k;
    int change, iters, tag;
    double g1, g1_sig, g2, g2_sig;
    double m, n;
    double dgn, dgu;

    int *num_in_bin = new int[bin_num];

    double *gh1_fit = new double[chi_grid_num*chi_grid_num];
    double *gh2_fit = new double[chi_grid_num*chi_grid_num];
    double *chisq1_fit = new double[chi_grid_num*chi_grid_num];
    double *chisq2_fit = new double[chi_grid_num*chi_grid_num];
    double *fit_x, *fit_y, *fit_z;
    int fit_num1, fit_num2;
    double chisq_ij, chisq_thresh;

    double *coeff = new double[6];

    double *mg = new double[data_num];
    double *temp = new double[data_num];
    double *mg_bin = new double[bin_num+1]{};  
    
    char inform[200];

    //show_arr(radius_bin,1, bin_num+1);
    //std::cout<<std::endl;
    
       
    // dg = 0.005
    dgn = delta_g;
    dgu = delta_g*5;
    chisq_thresh = 30;
    // for g1
    set_bin(mg1, data_num, mg_bin, bin_num, 100);
    find_shear_mean(mg1, mn, data_num, g1, g1_sig, 100, 1);
    for(i=0;i<chi_grid_num;i++) // y 
    {   
        m = g1 - dgn + dgn*2/(chi_grid_num-1)*i;        

        for(k=0;k<data_num;k++){temp[k]=mg1[k]-mn[k]*m;}

        for(j=0;j<chi_grid_num;j++) // x 
        {   
            tag = i*chi_grid_num+j;

            n = g1 - dgu + dgu*2/(chi_grid_num-1)*j;

            gh1_fit[tag] = m;
            gh2_fit[tag] = n;
            
            for(k=0;k<data_num;k++){mg[k] = temp[k]- n*mu[k];}
            
            histogram(mg, mg_bin, num_in_bin, data_num, bin_num);
            cal_chisq_1d(num_in_bin, bin_num, chisq_ij);        
            chisq1_fit[tag] = chisq_ij;
        }
    }
    fit_num1 = 0;
    for(i=0; i<chi_grid_num*chi_grid_num; i++)
    {
        g1_grid[i] = gh1_fit[i];
        g1_grid[i+chi_grid_num*chi_grid_num] = gh2_fit[i];
        chisq1_grid[i] = chisq1_fit[i];
        if(chisq1_fit[i] <= chisq_thresh)
        {
            fit_num1 ++;
        }
    }
    fit_x = new double[fit_num1];
    fit_y = new double[fit_num1];
    fit_z = new double[fit_num1];
    fit_num1 = 0;
    for(i=0; i<chi_grid_num*chi_grid_num; i++)
    {
        if(chisq1_fit[i] <= chisq_thresh)
        {
            fit_x[fit_num1] = gh1_fit[i];
            fit_y[fit_num1] = gh2_fit[i];
            fit_z[fit_num1] = chisq1_fit[i];
            fit_num1 ++;
        }
    }

    poly_fit_2d(fit_x, fit_y, fit_z, fit_num1, 2, coeff);
    gh1 = -coeff[1]/2/coeff[3];
    gh1_sig = sqrt(0.5/coeff[3]);
    for(i=0;i<6;i++){coeff_check[i] = coeff[i];}
    // gh2 = -coeff[2]/2/coeff[5];
    // gh2_sig = sqrt(0.5/coeff[5]);
    delete[] fit_x;
    delete[] fit_y;
    delete[] fit_z;

    // for g2
    set_bin(mg2, data_num, mg_bin, bin_num, 100);
    find_shear_mean(mg2, mn, data_num, g2, g2_sig, 100, 1);
    for(i=0;i<chi_grid_num;i++) // y 
    {   
        m = g2 - dgn + dgn*2/(chi_grid_num-1)*i; 

        for(k=0;k<data_num;k++){temp[k]=mg2[k]-mn[k]*m;}

        for(j=0;j<chi_grid_num;j++) // x 
        {   
            tag = i*chi_grid_num+j;

            n = g2 - dgu + dgu*2/(chi_grid_num-1)*j;

            gh1_fit[tag] = m;
            gh2_fit[tag] = n;
            
            for(k=0;k<data_num;k++){mg[k] = temp[k] - n*mu[k];}
            
            histogram(mg, mg_bin, num_in_bin, data_num, bin_num);
            cal_chisq_1d(num_in_bin, bin_num, chisq_ij);        
            chisq2_fit[tag] = chisq_ij;
        }
    }

    fit_num2 = 0;
    for(i=0; i<chi_grid_num*chi_grid_num; i++)
    {
        g2_grid[i] = gh1_fit[i];
        g2_grid[i+chi_grid_num*chi_grid_num] = gh2_fit[i];
        chisq2_grid[i] = chisq2_fit[i];
        if(chisq2_fit[i] <= chisq_thresh)
        {
            fit_num2 ++;
        }
    }
    fit_x = new double[fit_num2];
    fit_y = new double[fit_num2];
    fit_z = new double[fit_num2];
    fit_num2 = 0;
    for(i=0; i<chi_grid_num*chi_grid_num; i++)
    {
        if(chisq2_fit[i] <= chisq_thresh)
        {
            fit_x[fit_num2] = gh1_fit[i];
            fit_y[fit_num2] = gh2_fit[i];
            fit_z[fit_num2] = chisq2_fit[i];
            fit_num2 ++;
        }
    }

    poly_fit_2d(fit_x, fit_y, fit_z, fit_num2, 2, coeff);
    gh2 = -coeff[1]/2/coeff[3];
    gh2_sig = sqrt(0.5/coeff[3]);
    for(i=0;i<6;i++){coeff_check[i+6] = coeff[i];}
    sprintf(inform,"%9.6f %9.6f (%9.6f) %9.6f %9.6f (%9.6f) %d %d",g1_true, g1, g1_sig, g2_true, g2, g2_sig,fit_num1,fit_num2);
    std::cout<<inform<<std::endl;
    // gh2 = -coeff[2]/2/coeff[5];
    // gh2_sig = sqrt(0.5/coeff[5]);

    delete[] fit_x;
    delete[] fit_y;
    delete[] fit_z;

    delete[] num_in_bin;
    delete[] gh1_fit;
    delete[] gh2_fit;
    delete[] chisq1_fit;
    delete[] chisq2_fit;
    delete[] coeff;
    delete[] temp;
    delete[] mg;
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


void cal_radial_chisq(const double *radius, const double *tri_ratio, const int data_num, const double *radius_bin, const int bin_num, double &chisq)
{
    int i,j,k;
    double temp=0;

    double *tri_ratio_cumulate = new double[bin_num]{};
    int *num_count = new int[bin_num]{};
    // show_arr(num_count,1,bin_num);
    // show_arr(tri_ratio_cumulate,1,bin_num);
    // show_arr(radius_bin, 1, bin_num+1);
    // show_arr(radius, 1, data_num);
    // show_arr(tri_ratio,1, data_num);
    for(i=0;i<data_num;i++)
    {   
        for(j=0;j<bin_num;j++)
        {
            if(radius_bin[j] <= radius[i] and radius[i] <radius_bin[j+1])
            {   
                //std::cout<<radius_bin[j] <<" "<< radius[i] <<" "<<radius_bin[j+1]<<std::endl;
                num_count[j] += 1;
                tri_ratio_cumulate[j] += tri_ratio[i];
            }
        }
    }
    // show_arr(num_count,1,bin_num);
    // show_arr(tri_ratio_cumulate,1,bin_num);

    for(j=0;j<bin_num;j++)
    {
        if(num_count[j]==0){std::cout<<"Radial chi squared divided by ZERO !!!"<<std::endl;throw "Radial chi squared divided by ZERO !!!";}
        else{temp += (tri_ratio_cumulate[j]*tri_ratio_cumulate[j])/num_count[j];}        
    }
    chisq = temp;
    delete[] tri_ratio_cumulate;
    delete[] num_count;
}

void cal_radial_chisq(const double *mg1n, const double *mg2n, const double *radius, const int data_num, const double *radius_bin, const int bin_num, double &chisq)
{
    int i,j,k;
    double temp=0;

    double *mg1_ratio = new double[bin_num]{};
    double *mg2_ratio = new double[bin_num]{};

    int *num_count = new int[bin_num]{};
    // show_arr(num_count,1,bin_num);
    // show_arr(tri_ratio_cumulate,1,bin_num);
    // show_arr(radius_bin, 1, bin_num+1);
    // show_arr(radius, 1, data_num);
    // show_arr(tri_ratio,1, data_num);

    for(i=0;i<data_num;i++)
    {   
        for(j=0;j<bin_num;j++)
        {
            if(radius_bin[j] <= radius[i] and radius[i] <radius_bin[j+1])
            {   
                //std::cout<<radius_bin[j] <<" "<< radius[i] <<" "<<radius_bin[j+1]<<std::endl;
                num_count[j] += 1;
                mg1_ratio[j] += mg1n[i]/radius[i];
                mg2_ratio[j] += mg2n[i]/radius[i];

            }
        }
    }
    // show_arr(num_count,1,bin_num);
    // show_arr(tri_ratio_cumulate,1,bin_num);

    for(j=0;j<bin_num;j++)
    {
        if(num_count[j]==0){std::cout<<"Radial chi squared divided by ZERO !!!"<<std::endl;throw "Radial chi squared divided by ZERO !!!";}
        else{temp += (mg1_ratio[j]*mg1_ratio[j] + mg2_ratio[j]*mg2_ratio[j])/num_count[j];}        
    }
    chisq = temp;
    delete[] mg1_ratio;
    delete[] mg2_ratio;
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
    //show_arr(radius_bin,1, bin_num+1);
    //std::cout<<std::endl;

    //////////////////// find g1 ////////////////////////////////////
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
                mg2n = mg2[i];// - temp;
                //temp = mg2n/mg1n;
                radius[i] = sqrt(mg1n*mg1n + mg2n*mg2n);
                ratio[i] = mg1n/radius[i];//1./sqrt(1 + temp*temp);
            }
			cal_radial_chisq(radius, ratio, data_num, radius_bin, bin_num, chi_left);
            // exit(0);
            for(i=0;i<data_num;i++)
            {
                temp = mn[i]*gh_mid;
                mg1n = mg1[i] - temp;
                mg2n = mg2[i];// - temp;
                //temp = mg2n/mg1n;
                radius[i] = sqrt(mg1n*mg1n + mg2n*mg2n);
                ratio[i] = mg1n/radius[i];//1./sqrt(1 + temp*temp);
            }
			cal_radial_chisq(radius, ratio, data_num, radius_bin, bin_num, chi_mid);

            for(i=0;i<data_num;i++)
            {
                temp = mn[i]*gh_right;
                mg1n = mg1[i] - temp;
                mg2n = mg2[i];// - temp;
                //temp = mg2n/mg1n;
                radius[i] = sqrt(mg1n*mg1n + mg2n*mg2n);
                ratio[i] = mg1n/radius[i];//1./sqrt(1 + temp*temp);
            }
			cal_radial_chisq(radius, ratio, data_num, radius_bin, bin_num, chi_right);
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
            mg1n = mg1[j] - temp;
            mg2n = mg2[j];// - temp;
            //temp = mg2n/mg1n;
            radius[j] = sqrt(mg1n*mg1n + mg2n*mg2n);
            ratio[j] = mg1n/radius[j];//1./sqrt(1 + temp*temp);
        }
		try
		{
			cal_radial_chisq(radius, ratio, data_num, radius_bin, bin_num, chi_right);
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
    /////////////////////////  find g1 ////////////////////////////////////
    
    ////////////////////////  find g2 ////////////////////////////////////
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
                mg1n = mg1[i];// - temp;
                mg2n = mg2[i] - temp;
                //temp = mg2n/mg1n;
                radius[i] = sqrt(mg1n*mg1n + mg2n*mg2n);
                ratio[i] = mg2n/radius[i];//1./sqrt(1 + temp*temp);
            }
			cal_radial_chisq(radius, ratio, data_num, radius_bin, bin_num, chi_left);
            
            for(i=0;i<data_num;i++)
            {   
                temp = mn[i]*gh_mid;
                mg1n = mg1[i];// - temp;
                mg2n = mg2[i] - temp;
                //temp = mg2n/mg1n;
                radius[i] = sqrt(mg1n*mg1n + mg2n*mg2n);
                ratio[i] = mg2n/radius[i];//1./sqrt(1 + temp*temp);
            }
			cal_radial_chisq(radius, ratio, data_num, radius_bin, bin_num, chi_mid);

            for(i=0;i<data_num;i++)
            {
                temp = mn[i]*gh_right;
                mg1n = mg1[i];// - temp;
                mg2n = mg2[i] - temp;
                //temp = mg2n/mg1n;
                radius[i] = sqrt(mg1n*mg1n + mg2n*mg2n);
                ratio[i] = mg2n/radius[i];//1./sqrt(1 + temp*temp);
            }
			cal_radial_chisq(radius, ratio, data_num, radius_bin, bin_num, chi_right);
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
            mg1n = mg1[j];// - temp;
            mg2n = mg2[j] - temp;
            //temp = mg2n/mg1n;
            radius[j] = sqrt(mg1n*mg1n + mg2n*mg2n);
            ratio[j] = mg2n/radius[j];//1./sqrt(1 + temp*temp);
        }
		try
		{
			cal_radial_chisq(radius, ratio, data_num, radius_bin, bin_num, chi_right);
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
    ////////////////////////  find g2 ////////////////////////////////////

    
    delete[] radius;
    delete[] radius_bin;
}

void find_shear_dipole_2d(const double *mg1, const double *mg2, const double *mn, const int data_num, const int bin_num, 
                        double &gh1, double &gh1_sig,double &gh2, double &gh2_sig,const int chi_grid_num, double *g1_grid, double* g2_grid, double* chisq_grid)
{
    int i, j, k;
    int change, iters, tag;
    double g1, g1_sig, g2, g2_sig;
    double m, n;
    double dg;

    double *gh1_fit = new double[chi_grid_num*chi_grid_num];
    double *gh2_fit = new double[chi_grid_num*chi_grid_num];
    double *chisq_fit = new double[chi_grid_num*chi_grid_num];
    double chisq_ij;

    double *coeff = new double[6];

    double *mg1n = new double[data_num];
    double *mg2n = new double[data_num];
    double *radius = new double[data_num];

    double *radius_bin = new double[bin_num+1]{};


    for(i=0;i<data_num;i++)
    {
        radius[i] = sqrt(mg1[i]*mg1[i] + mg2[i]*mg2[i]);
    }

    set_radius_bin(radius, data_num, radius_bin, bin_num, 100);
    //show_arr(radius_bin,1, bin_num+1);
    //std::cout<<std::endl;
    
    find_shear_mean(mg1, mn, data_num, g1, g1_sig, 1000, 100);
    find_shear_mean(mg2, mn, data_num, g2, g2_sig, 1000, 100);

    // dg = 0.005
    dg = 0.002;
    for(i=0;i<chi_grid_num;i++) // y G1
    {   
        m = g1 - dg + dg*2/(chi_grid_num-1)*i;        

        for(k=0;k<data_num;k++){mg1n[k] = mg1[k] - m*mn[k];}

        for(j=0;j<chi_grid_num;j++) // x G2
        {   
            tag = i*chi_grid_num+j;

            n = g2 - dg + dg*2/(chi_grid_num-1)*j;

            gh1_fit[tag] = m;
            gh2_fit[tag] = n;
            
            for(k=0;k<data_num;k++)
            {
                mg2n[k] = mg2[k] - n*mn[k]; 
                radius[k] = sqrt(mg1n[k]*mg1n[k] + mg2n[k]*mg2n[k]);
            }

            cal_radial_chisq(mg1n, mg2n, radius, data_num, radius_bin, bin_num, chisq_ij);
            chisq_fit[tag] = chisq_ij;
        }
    }

    for(i=0; i<chi_grid_num*chi_grid_num; i++)
    {
        g1_grid[i] = gh1_fit[i];
        g2_grid[i] = gh2_fit[i];
        chisq_grid[i] = chisq_fit[i];
    }
    poly_fit_2d(gh1_fit, gh2_fit, chisq_fit, chi_grid_num*chi_grid_num, 2, coeff);
    gh1 = -coeff[1]/2/coeff[3];
    gh1_sig = sqrt(0.5/coeff[3]);

    gh2 = -coeff[2]/2/coeff[5];
    gh2_sig = sqrt(0.5/coeff[5]);
    
    delete[] mg1n;
    delete[] mg2n;
    delete[] gh1_fit;
    delete[] gh2_fit;
    delete[] chisq_fit;
    delete[] radius;
    delete[] radius_bin;
    delete[] coeff;
}



int main(int argc, char**argv)
{   
    int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

    char parent_path[200], data_path[200], set_name[30], data_name[40], result_path[200];
    char inform[300];
    char data_nm[30];
    char time_now[50],time_now_1[50];

    double *data, *mg1, *mg2, *mn, *mu,  *mnu1, *mnu2;
    double *check = new double[40];
    double *g1t, *g2t, *shear_result, *chisq_check;
    int gh_grid_num;
    double *check1 = new double[40];
    double *check2 = new double[40];
    double *gh1_grid;
    double *gh2_grid;
    double *chisq1_grid;
    double *chisq2_grid;


    strcpy(parent_path, argv[1]);
    strcpy(data_nm, argv[2]);
    sprintf(data_name, data_nm, rank);
    int ratio;
    ratio = atoi(argv[3]);
    gh_grid_num = atoi(argv[4]);

    int add_time, add_scale;
    add_time = 100;//atoi(argv[5]);
    add_scale = 100;//atoi(argv[6]);
    double delta_g;
    delta_g = atof(argv[5]);

    int i,j,k,shear_num;
    int data_row = 20000000, sub_row;
    double g1, g1_sig, g2, g2_sig;


    gh1_grid = new double[gh_grid_num*gh_grid_num*2];
    gh2_grid = new double[gh_grid_num*gh_grid_num*2];
    chisq1_grid = new double[gh_grid_num*gh_grid_num];
    chisq2_grid = new double[gh_grid_num*gh_grid_num];
    double *coeff_check = new double[12];

    sub_row = data_row/ratio;

    data = new double[data_row*4];
    mg1 = new double[sub_row];
    mg2 = new double[sub_row];
    mn = new double[sub_row];
    mu = new double[sub_row];
    mnu1 = new double[sub_row];
    mnu2 = new double[sub_row];

    sprintf(result_path, "%s/new_pdf_%.3f.hdf5",parent_path, delta_g);

    sprintf(data_path, "%s/shear.hdf5", parent_path);
    
    sprintf(set_name,"/g1");
    read_h5_datasize(data_path, set_name,shear_num);
    g1t = new double[shear_num];
    g2t = new double[shear_num];
    double *mc = new double[4];
    double *mc_all = new double[8];

    double *fit_val = new double[shear_num];
    double *fit_err = new double[shear_num];

    read_h5(data_path, set_name, g1t);

    sprintf(set_name,"/g2");
    read_h5_datasize(data_path, set_name,shear_num);
    read_h5(data_path, set_name, g2t);
    if(rank == 0)
    {   
        show_arr(g1t, 1, shear_num);
        show_arr(g2t, 1, shear_num);
        std::cout<<rank<<" "<<sub_row<<" "<<gh_grid_num<<std::endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Win win_shear, win_chisq_check;
	MPI_Aint shear_size, chisq_size;

	shear_size = shear_num *4*sizeof(double);
    chisq_size = 2*shear_num*20*2*sizeof(double);
	if (0 == rank)
	{	
		MPI_Win_allocate_shared(shear_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_result, &win_shear);
		MPI_Win_allocate_shared(chisq_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chisq_check, &win_chisq_check);

	}
	else
	{
		int disp_unit_s, disp_unit_c, disp_unit_num;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_result, &win_shear);
		MPI_Win_shared_query(win_shear, 0, &shear_size, &disp_unit_s, &shear_result);
        MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chisq_check, &win_chisq_check);
		MPI_Win_shared_query(win_chisq_check, 0, &chisq_size, &disp_unit_s, &chisq_check);
	}

    sprintf(data_path, "%s/%s",parent_path, data_name);
    sprintf(set_name,"/data");
    read_h5(data_path, set_name, data);

    for(i=0;i<sub_row;i++)
    {
        mg1[i] = data[i*4];
        mg2[i] = data[i*4+1];
        mn[i] = data[i*4+2];
        mu[i] = data[i*4+3];

        mnu1[i] = mn[i] + mu[i];
        mnu2[i] = mn[i] - mu[i];
    }
    if(rank==0){std::cout<<"Read data"<<std::endl;}
    MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////// average ////////////////////////
    find_shear_mean(mg1, mn, sub_row, g1, g1_sig, add_time, add_scale);
    find_shear_mean(mg2, mn, sub_row, g2, g2_sig, add_time, add_scale);
    shear_result[rank] = g1;
    shear_result[rank + shear_num] = g1_sig;
    shear_result[rank + shear_num*2] = g2;
    shear_result[rank + shear_num*3] = g2_sig;
    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0;i<numprocs;i++)
    {
        if(i ==  rank)
        {   
            sprintf(inform, "Ave g1: (true %9.6f) %9.6f (%9.6f), g2: (true %9.6f) %9.6f (%9.6f)", g1t[rank],g1, g1_sig, g2t[rank], g2, g2_sig);
            std::cout<<inform<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == numprocs-1)
    {
        for(k=0;k<2;k++)
        {
            for(j=0;j<shear_num;j++)
            {
                fit_val[j] = shear_result[j + k*shear_num*2];
                fit_err[j] = shear_result[j + shear_num + k*shear_num*2];
            }
            if(k ==0){poly_fit_1d(g1t, fit_val, fit_err, shear_num, mc, 1);}
            else{poly_fit_1d(g2t, fit_val, fit_err, shear_num, mc, 1);}

            mc_all[k*4] = mc[2]-1;// m
            mc_all[k*4 + 1] = mc[3];// m_sig
            mc_all[k*4 + 2] = mc[0];// c
            mc_all[k*4 + 3] = mc[1];// c_sig
        }
        std::cout<<"Ave m & c:"<<std::endl;
        show_arr(mc_all,2,4);
        std::cout<<std::endl;
        initialize_arr(shear_result, 4*shear_num, 0);
        initialize_arr(chisq_check, 2*shear_num*40,0);
    }
    MPI_Barrier(MPI_COMM_WORLD);



    // /////////////////////  PDF //////////////////////////
    // find_shear(mg1, mnu1, sub_row, 10, g1, g1_sig, check);
    // // for(i=0;i<40;i++){chisq_check[rank*40 + i] = check[i];}
    // find_shear(mg2, mnu2, sub_row, 10, g2, g2_sig, check);
    // // for(i=0;i<40;i++){chisq_check[shear_num * 40 + rank*40 + i] = check[i];}
    
    // shear_result[rank] = g1;
    // shear_result[rank + shear_num] = g1_sig;
    // shear_result[rank + shear_num*2] = g2;
    // shear_result[rank + shear_num*3] = g2_sig;
    // MPI_Barrier(MPI_COMM_WORLD);
    // for(i=0;i<numprocs;i++)
    // {
    //     if(i ==  rank)
    //     {   
    //         sprintf(inform, "PDF g1: (true %9.6f) %9.6f (%9.6f), g2: (true %9.6f) %9.6f (%9.6f)", g1t[rank],g1, g1_sig, g2t[rank], g2, g2_sig);
    //         std::cout<<inform<<std::endl;
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // get_time(time_now,50);
    
    // if(rank == numprocs-1)
    // {   
    //     // get_time(time_now_1,50);
    //     // std::cout<<rank<<" "<<time_now_1<<std::endl;
    //     for(k=0;k<2;k++)
    //     {
    //         for(j=0;j<shear_num;j++)
    //         {
    //             fit_val[j] = shear_result[j + k*shear_num*2];
    //             fit_err[j] = shear_result[j + shear_num + k*shear_num*2];
    //         }
    //         if(k ==0){poly_fit_1d(g1t, fit_val, fit_err, shear_num, mc, 1);}
    //         else{poly_fit_1d(g2t, fit_val, fit_err, shear_num, mc, 1);}

    //         mc_all[k*4] = mc[2]-1;// m
    //         mc_all[k*4 + 1] = mc[3];// m_sig
    //         mc_all[k*4 + 2] = mc[0];// c
    //         mc_all[k*4 + 3] = mc[1];// c_sig
    //     }
    //     std::cout<<"PDF m & c:"<<std::endl;
    //     show_arr(mc_all,2,4);
    //     std::cout<<std::endl;
        
        
    //     sprintf(set_name, "/PDF/shear");
    //     write_h5(result_path, set_name, shear_result, 4, shear_num, true);
    //     sprintf(set_name, "/PDF/chisq");
    //     write_h5(result_path, set_name, chisq_check, 2*shear_num, 40, false);
    //     sprintf(set_name, "/PDF/mc");
    //     write_h5(result_path, set_name, mc_all, 2, 4, false);
        
    //     initialize_arr(shear_result, 4*shear_num, 0);
    //     initialize_arr(chisq_check, 2*shear_num*40,0);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // std::cout<<rank<<time_now<<std::endl;


    // //////////////////  PDF-dipole ///////////////////////////////
    // // find_shear_dipole(mg1, mg2, mn, sub_row, 10, g1, g1_sig, g2, g2_sig, check1, check2);
    // // find_shear_dipole_2d(mg1, mg2, mn, sub_row, 10, g1, g1_sig, g2, g2_sig, gh_grid_num, gh1_grid, gh2_grid, chisq_grid);
    find_shear_2d(mg1, mg2, mn, mu, sub_row, 10, g1, g1_sig, g2, g2_sig, gh_grid_num,gh1_grid, gh2_grid, chisq1_grid, chisq2_grid,g1t[rank],g2t[rank], delta_g, coeff_check);
    
    shear_result[rank] = g1;
    shear_result[rank + shear_num] = g1_sig;
    shear_result[rank + shear_num*2] = g2;
    shear_result[rank + shear_num*3] = g2_sig;
    
    for(i=0;i<40;i++){chisq_check[rank*40 + i] = check1[i];}
    for(i=0;i<40;i++){chisq_check[shear_num * 40 + rank*40 + i] = check2[i];}
    
    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0;i<numprocs;i++)
    {
        if(i ==  rank)
        {   
            sprintf(inform, "PDF_new g1: (true %9.6f) %9.6f (%9.6f), g2: (true %9.6f) %9.6f (%9.6f)", g1t[rank],g1, g1_sig, g2t[rank], g2, g2_sig);
            std::cout<<inform<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == numprocs-1)
    {
        for(k=0;k<2;k++)
        {
            for(j=0;j<shear_num;j++)
            {
                fit_val[j] = shear_result[j + k*shear_num*2];
                fit_err[j] = shear_result[j + shear_num + k*shear_num*2];
            }
            if(k ==0){poly_fit_1d(g1t, fit_val, fit_err, shear_num, mc, 1);}
            else{poly_fit_1d(g2t, fit_val, fit_err, shear_num, mc, 1);}

            mc_all[k*4] = mc[2]-1;// m
            mc_all[k*4 + 1] = mc[3];// m_sig
            mc_all[k*4 + 2] = mc[0];// c
            mc_all[k*4 + 3] = mc[1];// c_sig
        }
        std::cout<<"PDF_new m & c:"<<std::endl;
        show_arr(mc_all,2,4);
        std::cout<<std::endl;
        
        sprintf(set_name, "/PDF_new/shear");
        write_h5(result_path, set_name, shear_result, 4, shear_num, true);
        sprintf(set_name, "/PDF_new/chisq");
        write_h5(result_path, set_name, chisq_check, 2*shear_num, 40, false);
        sprintf(set_name, "/PDF_new/mc");
        write_h5(result_path, set_name, mc_all, 2, 4, false);

        initialize_arr(shear_result, 4*shear_num, 0);
        initialize_arr(chisq_check, 2*shear_num*40,0);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {
            sprintf(set_name, "/PDF_new/%d/gh1_grid", rank);
            write_h5(result_path, set_name, gh1_grid, gh_grid_num*2, gh_grid_num, false);
            sprintf(set_name, "/PDF_new/%d/gh2_grid", rank);
            write_h5(result_path, set_name, gh2_grid, gh_grid_num*2, gh_grid_num, false);
            sprintf(set_name, "/PDF_new/%d/chisq1_grid", rank);
            write_h5(result_path, set_name, chisq1_grid, gh_grid_num, gh_grid_num, false);
            sprintf(set_name, "/PDF_new/%d/chisq2_grid", rank);
            write_h5(result_path, set_name, chisq2_grid, gh_grid_num, gh_grid_num, false);
            sprintf(set_name, "/PDF_new/%d/coeff", rank);
            write_h5(result_path, set_name, coeff_check, 2, 6, false);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;

}