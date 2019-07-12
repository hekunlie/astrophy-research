#include<FQlib.h>

void read_text_(const char *filename, double *data_buf, const int data_col,  const int skipline=0)
{
	std::ifstream infile;
	std::string file_path, str, temp;
	std::stringstream strs;

	int i, j, line_count;

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
			std::cout << str << std::endl;
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
			temp.clear();
			getline(infile, str);
			std::cout << str << std::endl;
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

void write_text_(const char*filename, const double *data_buf, const int data_row, const int data_col, const int mode)
{
	int i, j, k;
	std::ofstream fw;

	if (mode == 0)
	{
		fw.open(filename, std::ios::out | std::ios::trunc);
	}
	else
	{
		fw.open(filename, std::ios::out | std::ios::app);
	}

	for (i = 0; i < data_row; i++)
	{
		for (j = 0; j < data_col - 1; j++)
		{
			fw << data_buf[i*data_col + j] << "\t";
		}
		fw << data_buf[i*data_col + j]<< std::endl;
	}
	fw.close();
}

void write_text_(const char*filename, double *data_buf, const int data_row, const int data_col, const char * comment, const int mode)
{
	int i, j, k;
	std::ofstream fw;

	if (mode == 0)
	{
		fw.open(filename, std::ios::out | std::ios::trunc);
	}
	else
	{
		fw.open(filename, std::ios::out | std::ios::app);
	}

	fw << comment << std::endl;
	for (i = 0; i < data_row; i++)
	{
		for (j = 0; j < data_col - 1; j++)
		{
			fw << data_buf[i*data_col + j] << "\t";
		}
		fw << data_buf[i*data_col + j] << std::endl;
	}
	fw.close();
}

int main(int argc, char **argv)
{
	double st1, st2;
	char *filename = new char[50];
	char *commet = new char[50];

	int lines;
	int i, j, k, data_row, data_col, skip_line;

	data_row = atoi(argv[2]);
	data_col = atoi(argv[3]);
	skip_line = atoi(argv[4]);

	double *data = new double[data_row*data_col];

	st1 = clock();
	read_text(argv[1], data, data_col, skip_line);
	st2 = clock();

	sprintf(filename, "%s_w", argv[1]);
	write_text(filename, data, data_row, data_col,0);

	sprintf(filename, "%s_w_c", argv[1]);
	sprintf(commet, "# this is the test");
	write_text(filename, data, data_row, data_col, commet,0);

	if (lines < 10)
	{
		show_arr(data, data_row, data_col);
	}
	
	std::cout << lines<<" lines"<< std::endl;
	std::cout << (st2 - st1) / CLOCKS_PER_SEC << std::endl;
	return 0;
}