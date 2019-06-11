#include<FQlib.h>
#include<iostream>

//#define SMALL_CATA

int main(int argc, char **argv)
{
	double start = atof(argv[1]);
	double end = atof(argv[2]);
	int num = atoi(argv[3]);

	double *test_bin = new double[num];
	linspace(start, end, num, test_bin);
	show_arr(test_bin, 1, num);

#if defined(SMALL_CATA)
	std::cout << "123" << std::endl;
#endif //SMALL_CATA

	return 0;
}