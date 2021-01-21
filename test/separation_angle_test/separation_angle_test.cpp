#include<FQlib.h>
#include<science_cal_lib.h>
// #include<hk_iolib.h>

int main(int argc, char**argv)
{
	double RA1, RA2, DEC1, DEC2;
	float RA1f, RA2f, DEC1f, DEC2f;

	double sep_angle1;
	float sep_anglef1;
	double sep_angle2;
	float sep_anglef2;
	
	double st1, st2, tt1, tt2;
	int i;
	int num = 10000000;

	RA1 = atof(argv[1]);
	DEC1 = atof(argv[2]);
	RA2 = atof(argv[3]);
	DEC2 = atof(argv[4]);
	
	RA1f = RA1;
	DEC1f = DEC1;
	RA2f = RA2;
	DEC2f = DEC2;

	st1 = clock();
	for(i=0; i<num; i++)
	{
		separation_angle_1(RA1, DEC1, RA2, DEC2, sep_angle1);
		separation_angle_1(RA1f, DEC1f, RA2f, DEC2f, sep_anglef1);
	}
	st2 = clock();
	tt1 = (st2-st1)/CLOCKS_PER_SEC;

	std::cout<<sep_angle1<<" "<<sep_anglef1<<"  "<<tt1<<std::endl;

	st1 = clock();
	for(i=0; i<num; i++)
	{
		separation_angle_2(RA1, DEC1, RA2, DEC2, sep_angle2);
		separation_angle_2(RA1f, DEC1f, RA2f, DEC2f, sep_anglef2);
	}
	st2 = clock();
	tt2 = (st2-st1)/CLOCKS_PER_SEC;

	std::cout<<sep_angle2<<" "<<sep_anglef2<<"  "<<tt2<<std::endl;

	std::cout<<sep_angle1-sep_angle2<<" "<<sep_anglef1-sep_anglef2<<std::endl;
	std::cout<<(tt2-tt1)/tt2<<std::endl;
	
	// char inform[200];
	// sprintf(inform, "(%.4f, %.4f) <--  %.6f (%.4f deg) --> (%.4f, %.4f)", RA1, DEC1,  sep_angle, sep_angle / Pi * 180, RA2, DEC2);
	// std::cout << inform << std::endl;

	// if (argc > 5)
	// {
	// 	char path[200], set_name[20];
	// 	int i, j, k, arr_size, num;
	// 	int data_col = 6;

	// 	sprintf(path, "./data/sep_test.hdf5");
	// 	sprintf(set_name, "/data");

	// 	read_h5_datasize(path, set_name, arr_size);
	// 	num = arr_size / data_col;

	// 	double *data = new double[arr_size];
	// 	double *diff = new double[num * 2];
	// 	double diff_theta;

	// 	read_h5(path, set_name, data);

	// 	for (i = 0; i < num; i++)
	// 	{
	// 		separation(data[i*data_col], data[i*data_col + 2], data[i*data_col + 1], data[i*data_col + 3], diff_theta);
	// 		diff[i * 2] = diff_theta;
	// 		diff[i * 2 + 1] = diff_theta / Pi * 180;
	// 		if (fabs(diff_theta - data[i*data_col + 4]) > 0.001)
	// 		{
	// 			sprintf(inform, "P1: (%.4f, %.4f). P2: (%.4f, %.4f). Separation: %.6f (%.4f deg) %.6f",
	// 				data[i*data_col], data[i*data_col + 2], data[i*data_col + 1], data[i*data_col + 3], diff[i * 2], diff[i * 2 + 1] / Pi * 180, data[i*data_col + 4]);

	// 			std::cout << inform << std::endl;
	// 		}
	// 	}

	// 	sprintf(path, "./data/sep_test_cpp.hdf5");
	// 	sprintf(set_name, "/data");
	// 	write_h5(path, set_name, diff, num, 2, TRUE);
	// }
	return 0;
}