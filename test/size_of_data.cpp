#include<FQlib.h>
#include<climits>
int main(int argc, char **argv)
{
	std::cout << "char    :" << CHAR_MIN << "  ~  " << CHAR_MAX <<". "<<sizeof(char)<<" Bytes"<< '\n';
	std::cout << "short   :" << SHRT_MIN << "  ~  " << SHRT_MAX << ". " << sizeof(short) << " Bytes" << '\n';
	std::cout << "int     :" << INT_MIN << "  ~  " << INT_MAX << ". " << sizeof(int) << " Bytes" << '\n';
	std::cout << "long    :" << LONG_MIN <<"  ~  " << LONG_MAX << ". " << sizeof(long) << " Bytes" << '\n';
	std::cout << "float    :" << LONG_MIN << "  ~  " << LONG_MAX << ". " << sizeof(float) << " Bytes" << '\n';
	std::cout << "double    :" << LONG_MIN << "  ~  " << LONG_MAX << ". " << sizeof(double) << " Bytes" << '\n';

	long m, n, k;
	int i, j;
	double *data = new double[100];
	for (i = 0; i < 100; i++)
	{
		data[i] = i;
	}

	m = 10;
	show_arr(data, 1, 100);
	show_arr(data + m, 1, 100 - m);
	return 0;
}