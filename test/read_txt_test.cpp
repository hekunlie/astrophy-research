#include<FQlib.h>

void read_text(const char *filename, double *data_buf, const int data_col, int &line_count, const int skipline=0)
{
	std::ifstream infile;
	std::string file_path, str;
	std::stringstream strs;
	int i, j;

	//char_to_str(filename, file_path);

	infile.open(filename);
	line_count = 0;
	if (skipline > 0)
	{
		while (!infile.eof())
		{
			str.clear();
			strs.clear();
			getline(infile, str);

			if (line_count > 0)
			{
				strs << str;
				for (i = 0; i < data_col; i++)
				{
					strs >> data_buf[(line_count - 1)*data_col + i];
				}
			}
			line_count += 1;
		}
	}
	else
	{
		while (!infile.eof())
		{
			str.clear();
			strs.clear();
			getline(infile, str);
			std::cout << "READ: " << str <<" The first"<<str[0]<< std::endl;
			strs << str;
			for (i = 0; i < data_col; i++)
			{
				strs >> data_buf[line_count*data_col + i];
			}

			line_count += 1;
		}
	}
	infile.close();
}

int main(int argc, char **argv)
{
	double st1, st2;
	char *filename = new char[50];

	int lines;
	int i, j, k;
	i = atoi(argv[2]);
	j = atoi(argv[3]);
	k = atoi(argv[4]);

	double *data = new double[i*j];

	st1 = clock();
	read_text(argv[1], data, j, lines, k);
	st2 = clock();
	if (lines < 10)
	{
		show_arr(data, i, j);
	}
	
	std::cout << lines << std::endl;
	int m, n;
	for (m =0; m < j; m++)
	{
		std::cout << data[(i - 1)*j + m] << " ";
	}
	std::cout << std::endl;
	std::cout << (st2 - st1) / CLOCKS_PER_SEC << std::endl;
	return 0;
}