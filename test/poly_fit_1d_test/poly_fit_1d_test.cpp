#include<FQlib.h>

int main(int argc, char* argv[])
{
	int order, num;
	order = atoi(argv[1]);
	num = atoi(argv[2]);
	double *x = new double[num] {};
	double *fx = new double[num] {};
	double *coeff = new double[order + 1]{};
	char data_path[150];
	char set_name[50];
	sprintf(data_path, "test.hdf5");
	sprintf(set_name, "/x");
	read_h5(data_path, set_name, x);
	sprintf(set_name, "/fx");
	read_h5(data_path, set_name, fx);
	//show_arr(x, 1, num);
	//show_arr(fx, 1, num);
	std::cout << "Order: " << order << ",  Num: " << num << std::endl;
	poly_fit_1d(x, fx, num, order, coeff);
	
	std::cout << "The coeff: ";
	show_arr(coeff, 1, order + 1);

	return 0;
}