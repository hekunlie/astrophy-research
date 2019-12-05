#include<FQlib.h>
#include<hk_iolib.h>

int main(int argc, char* argv[])
{	
	int i,j,k,num;

	num = 15;
	double *x = new double[num] {};
	double *y = new double[num] {};
	double *yerr = new double[num]{};
	double *data = new double[6*num]{};
	double *mc_data = new double[2*4]{};
	double *mc_fit = new double[4]{};
	double *coeff = new double[4]{};

	char data_path[100], set_name[30];

	sprintf(data_path, "fit_shear_test.hdf5");
	sprintf(set_name, "/mean_result");
	read_h5(data_path, set_name, data);
	sprintf(set_name, "/mean_mc");
	read_h5(data_path, set_name, mc_data);

	initialize_arr(mc_fit,4,0);
	for(i=0;i<num;i++)
	{
		x[i] = data[i];
		y[i] = data[i + num];
		yerr[i] = data[i + 2*num];
	}
	poly_fit_1d(x, y, yerr, num, coeff, 1);
	mc_fit[0] = coeff[2]-1;
	mc_fit[1] = coeff[3];
	mc_fit[2] = coeff[0];
	mc_fit[3] = coeff[1];
	
	std::cout<<"m1 & c1:"<<std::endl;
	show_arr(mc_fit, 1, 4);

	initialize_arr(mc_fit,4,0);
	for(i=0;i<num;i++)
	{
		x[i] = data[i + 3*num];
		y[i] = data[i + 4*num];
		yerr[i] = data[i + 5*num];
	}
	poly_fit_1d(x, y, yerr, num, coeff, 1);
	mc_fit[0] = coeff[2]-1;
	mc_fit[1] = coeff[3];
	mc_fit[2] = coeff[0];
	mc_fit[3] = coeff[1];
	std::cout<<"m2 & c2:"<<std::endl;
	show_arr(mc_fit, 1, 4);
	
	std::cout<<"Result of Python:"<<std::endl;
	show_arr(mc_data, 2, 4);
	return 0;
}