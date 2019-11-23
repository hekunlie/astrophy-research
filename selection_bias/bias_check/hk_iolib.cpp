#include "hk_iolib.h"

std::ofstream default_loggers;

/********************************************************************************************************************************************/
/* file reading and writting*/
/********************************************************************************************************************************************/

void char_to_str( const char *char_in, std::string &string_out)
{
	std::stringstream media;
	media << char_in;
	string_out = media.str();
}

void write_log(char*filename, char *inform)
{
	char time_now[40];
	time_t timer;
	time(&timer);
	strftime(time_now, 40, "%Y-%m-%d %H:%M:%S", localtime(&timer));

	default_loggers.open(filename, std::ios::out | std::ios::app);
	default_loggers << time_now << " ---- " << inform << std::endl;
	default_loggers.close();
}

bool file_exist(const char *filename)
{
	struct stat my_stat;
	return (stat(filename, &my_stat) == 0);
}

void read_config(const std::string path, const std::string target_section, const std::string target_para, std::string &para)
{
	// find the content or value of a parameter, "name", in section,"section".
	// the config file muse be written as follow:
	// [section a]			 # section name must be enclosed by "[ ]"
	// para_a = aaaaa    #  the content in each section must consist of parameter name, "=" and the content of the parameter
	//								# and there is a space in each side of "=" to separate the the name and the content

	std::ifstream infile;
	std::string str, str1, str2, str3, sect_name;
	std::stringstream strs;
	int sect_found = 0, para_found = 0, line_label = 0;
	int pos1, pos2;
	infile.open(path);
	while (!infile.eof())
	{
		str.clear();
		str1.clear();
		str2.clear();
		str3.clear();
		strs.clear();

		getline(infile, str);

		// find the section name which in a "[ ]"
		pos1 = str.find("[");
		pos2 = str.find("]");
		if (pos1 > -1 and pos2 > pos1)
		{
			sect_name = str.substr(pos1 + 1, pos2 - pos1 - 1);
			// find the section
			if (target_section == sect_name)
			{
				// label the found of section
				sect_found = 1;
				// label this line as the line of section name not the parameter lines
				line_label = 1;
				// the section has been found and skip this line in file
				continue;
			}
		}
		else
		{
			line_label = 0;
		}
		if (1 == sect_found and 0 == line_label)
		{
			strs << str;
			strs >> str1 >> str2 >> str3;			
			if (target_para == str1)
			{
				para = str3;
				para_found = 1;
				break;
			}
		}
	 }
	infile.close();
	if (0 == sect_found)
	{
		str.clear();
		str = "Section '" + target_section + "' can not be found!!";
		std::cout << str << std::endl;
		exit(0);
	}
	if (0 == para_found)
	{
		str.clear();
		str = "Parameter '" + target_para + "' can not be found!!";
		std::cout << str << std::endl;
		exit(0);
	}
}

void read_para(const std::string path, const std::string name, double &para)
{
	std::ifstream infile;
	std::string str, str1, str2, str3;
	std::stringstream strs;
	int f = 0;
	infile.open(path);
	while (!infile.eof())
	{
		str.clear();
		str1.clear();
		str2.clear();
		str3.clear();
		strs.clear();

		getline(infile, str);
		strs << str;
		strs >> str1 >> str2 >> str3;

		if (str1 == name)
		{
			para = std::stod(str3);
			f = 1;
			break;
		}
	}
	infile.close();
	if (f == 0)
	{
		str.clear();
		str = name + " can not be found!!";
		std::cout << str << std::endl;
		exit(0);
	}
}

void read_para(const std::string path, const std::string name, float &para)
{
	std::ifstream infile;
	std::string str, str1, str2, str3;
	std::stringstream strs;
	int f = 0;
	infile.open(path);
	while (!infile.eof())
	{
		str.clear();
		str1.clear();
		str2.clear();
		str3.clear();
		strs.clear();

		getline(infile, str);
		strs << str;
		strs >> str1 >> str2 >> str3;

		if (str1 == name)
		{
			para = std::stof(str3);
			f = 1;
			break;
		}
	}
	infile.close();
	if (f == 0)
	{
		str.clear();
		str = name + " can not be found!!";
		std::cout << str << std::endl;
		exit(0);
	}
}

void read_para(const std::string path, const std::string name, int &para)
{
	std::ifstream infile;
	std::string str, str1, str2, str3;
	std::stringstream strs;
	infile.open(path);
	int f = 0;
	while (!infile.eof())
	{
		str.clear();
		str1.clear();
		str2.clear();
		str3.clear();
		strs.clear();

		getline(infile, str);
		strs << str;
		strs >> str1 >> str2 >> str3;

		if (str1 == name)
		{
			para = std::stoi(str3);
			f = 1;
			break;
		}
	}
	infile.close();
	if (f == 0)
	{
		str.clear();
		str = name + " can not be found!!";
		std::cout << str << std::endl;
		exit(0);
	}
}


void read_text(const char *filename, double *data_buf, const int data_col, const int skipline)
{
	std::ifstream infile;
	std::string file_path, str, temp;
	std::stringstream strs;

	int i, j;

	//char_to_str(filename, file_path);
	int line_count = 0;

	infile.open(filename);

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

void read_text(const std::string path, double *arr, const int start_line, const int read_lines, const int read_cols)
{
	// start_line counts from 0 not 1 !!!
	int i = start_line, ie = start_line+read_lines, j;
	std::ifstream infile;
	std::string str;
	std::string *strs= new std::string[read_cols];
	std::stringstream strs_;

	infile.open(path);

	while (i<ie)
	{	
		// clear the string before assignment
		str.clear();	

		// read the line
		getline(infile, str);

		if (i >= start_line)
		{
			strs_.clear();
			// assign the contents to stringstream
			strs_ << str;
			for (j = 0; j < read_cols; j++)
			{
				strs[j].clear();
				strs_ >> strs[j];
				// to the array
				arr[(i - start_line)*read_cols + j] = std::stod(strs[j]);
			}
		}
		i++;
	}

	infile.close();

	delete[] strs;
}

void read_text(const std::string path, double *arr, const int read_lines)
{
	std::ifstream infile;
	int i = 0;
	infile.open(path);
	while (i<read_lines && infile>>arr[i])
	{
		i++;
	}	
	infile.close();
}

void read_text(const std::string path, float *arr, const int read_lines)
{
	std::ifstream infile;
	int i = 0;
	infile.open(path);
	while (i<read_lines && infile >> arr[i])
	{
		i++;
	}
	infile.close();
}

void read_text(const std::string path, int *arr, const int read_lines)
{
	std::ifstream infile;
	int i = 0;
	infile.open(path);
	while (i<read_lines && infile >> arr[i])
	{
		i++;
	}
	infile.close();
}


void write_text(const char*filename, const double *data_buf, const int data_row, const int data_col, const int mode)
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
		fw << data_buf[i*data_col + j] << std::endl;
	}
	fw.close();
}

void write_text(const char*filename, double *data_buf, const int data_row, const int data_col, const char * comment, const int mode)
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


void read_h5_datasize(const char *filename, const char *set_name, int &elem_num)
{
	int num;
	herr_t status;
	hid_t file_id, dataset_id, space_id;

	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	space_id = H5Dget_space(dataset_id);
	num = H5Sget_simple_extent_npoints(space_id);
	if (num > 0)
	{
		elem_num = num;
	}
	else
	{
		elem_num = -1;
		std::cout << "Failed in reading the size of " << set_name << "." << std::endl;
	}
	status = H5Sclose(space_id);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void read_h5_datasize(const char *filename, const char *set_name, long &elem_num)
{
	long num;
	herr_t status;
	hid_t file_id, dataset_id, space_id;

	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	space_id = H5Dget_space(dataset_id);
	num = H5Sget_simple_extent_npoints(space_id);
	if (num > 0)
	{
		elem_num = num;
	}
	else
	{
		elem_num = -1;
		std::cout << "Failed in reading the size of " << set_name << "." << std::endl;
	}
	status = H5Sclose(space_id);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void read_h5(const char *filename, const char *set_name, double *arr)
{
	hid_t file_id;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t dataset_id;
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void read_h5(const char *filename, const char *set_name, float *arr)
{
	hid_t file_id;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t dataset_id;
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void read_h5(const char *filename, const char *set_name, int *arr)
{
	hid_t file_id;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t dataset_id;
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void read_h5(const char *filename, const char *set_name, long *arr)
{
	hid_t file_id;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t dataset_id;
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void read_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, double *buff, std::string flag)
{
	hid_t file_id, dataset_id, attr;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (flag == "d")
	{
		dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
		attr = H5Aopen_name(dataset_id, attrs_name);
		status = H5Aread(attr, H5T_NATIVE_DOUBLE, buff);
		status = H5Aclose(attr);
		status = H5Dclose(dataset_id);
	}
	else
	{
		dataset_id = H5Gopen(file_id, set_name, H5P_DEFAULT);
		attr = H5Aopen_name(dataset_id, attrs_name);
		status = H5Aread(attr, H5T_NATIVE_DOUBLE, buff);
		status = H5Aclose(attr);
		status = H5Gclose(dataset_id);
	}


	status = H5Fclose(file_id);
}

void read_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, float *buff, std::string flag)
{
	hid_t file_id, dataset_id, attr;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (flag == "d")
	{
		dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
		attr = H5Aopen_name(dataset_id, attrs_name);
		status = H5Aread(attr, H5T_NATIVE_FLOAT, buff);
		status = H5Aclose(attr);
		status = H5Dclose(dataset_id);
	}
	else
	{
		dataset_id = H5Gopen(file_id, set_name, H5P_DEFAULT);
		attr = H5Aopen_name(dataset_id, attrs_name);
		status = H5Aread(attr, H5T_NATIVE_FLOAT, buff);
		status = H5Aclose(attr);
		status = H5Gclose(dataset_id);
	}


	status = H5Fclose(file_id);
}

void read_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, int *buff, std::string flag)
{
	hid_t file_id, dataset_id, attr;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (flag == "d")
	{
		dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
		attr = H5Aopen_name(dataset_id, attrs_name);
		status = H5Aread(attr, H5T_NATIVE_INT, buff);
		status = H5Aclose(attr);
		status = H5Dclose(dataset_id);
	}
	else
	{
		dataset_id = H5Gopen(file_id, set_name, H5P_DEFAULT);
		attr = H5Aopen_name(dataset_id, attrs_name);
		status = H5Aread(attr, H5T_NATIVE_INT, buff);
		status = H5Aclose(attr);
		status = H5Gclose(dataset_id);
	}
	
	
	status = H5Fclose(file_id);
}

void write_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, const double *attrs_buffer, const int buffer_len, std::string flag)
{
	hid_t file_id, dataset_id, dataspace_id, attrs_id;
	herr_t status;
	hsize_t num = buffer_len;
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

	if (flag == "d")
	{
		dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
		dataspace_id = H5Screate_simple(1, &num, NULL);
		attrs_id = H5Acreate2(dataset_id, attrs_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrs_id, H5T_NATIVE_DOUBLE, attrs_buffer);
		status = H5Aclose(attrs_id);
		status = H5Sclose(dataspace_id);
		status = H5Dclose(dataset_id);
	}
	else
	{
		dataset_id = H5Gopen(file_id, set_name, H5P_DEFAULT);
		dataspace_id = H5Screate_simple(1, &num, NULL);
		attrs_id = H5Acreate2(dataset_id, attrs_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrs_id, H5T_NATIVE_DOUBLE, attrs_buffer);
		status = H5Aclose(attrs_id);
		status = H5Sclose(dataspace_id);
		status = H5Gclose(dataset_id);
	}
	status = H5Fclose(file_id);
}

void write_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, const int *attrs_buffer, const int buffer_len, std::string flag)
{
	hid_t file_id, dataset_id, dataspace_id, attrs_id;
	herr_t status;
	hsize_t num = buffer_len;
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

	if (flag == "d")
	{
		dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
		dataspace_id = H5Screate_simple(1, &num, NULL);
		attrs_id = H5Acreate2(dataset_id, attrs_name, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrs_id, H5T_NATIVE_INT, attrs_buffer);
		status = H5Aclose(attrs_id);
		status = H5Sclose(dataspace_id);
		status = H5Dclose(dataset_id);
	}
	else
	{
		dataset_id = H5Gopen(file_id, set_name, H5P_DEFAULT);
		dataspace_id = H5Screate_simple(1, &num, NULL);
		attrs_id = H5Acreate2(dataset_id, attrs_name, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrs_id, H5T_NATIVE_INT, attrs_buffer);
		status = H5Aclose(attrs_id);
		status = H5Sclose(dataspace_id);
		status = H5Gclose(dataset_id);
	}
	status = H5Fclose(file_id);
}

void create_h5_group(const char *filename, const char *set_name, const bool trunc)
{
	int i, j, m, count = 0, s_count = 0;

	int *slash;
	char *name;

	hid_t file_id, group_id;
	herr_t status;
	unsigned rank = 2;
	hsize_t dims[2];
	hid_t dataspace_id;
	hid_t dataset_id;

	// the length of the "set_name"
	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}
	slash = new int[count] {};
	name = new char[count + 1]{};

	// find the position of "/"
	for (i = 0; i < count; i++)
	{
		if (set_name[i] == '/')
		{
			slash[s_count] = i;
			s_count++;
		}
	}
	// the position of '\0', the end
	slash[s_count] = count;
	if (trunc)// truncate it
	{
		file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	}
	else
	{
		if (file_exist(filename))
		{
			file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
		}
		else
		{
			file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		}
	}

	for (i = 0; i < s_count; i++)
	{
		for (j = 0; j < slash[i + 1]; j++)
		{
			name[j] = set_name[j];
		}
		name[j] = '\0';// label the end
		status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
		status = H5Lget_info(file_id, name, NULL, H5P_DEFAULT);
		if (status != 0)// doesn't exist
		{
			group_id = H5Gcreate(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Gclose(group_id);
		}
	}
	status = H5Eset_auto(H5E_DEFAULT, (H5E_auto2_t)H5Eprint2, stderr);
	status = H5Fclose(file_id);
	delete[] slash;
	delete[] name;

}

void write_h5(const char *filename, const char *set_name, const double*data, const int row, const int column, const bool trunc)
{
	int i, count, s_count;
	int *slash;
	char *name, *new_name;
	// the length of the set_name
	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}
	slash = new int[count] {};
	name = new char[count + 1];
	new_name = new char[count + 1];
	// find the slashes in the set_name
	s_count = 0;
	for (i = 0; i < count; i++)
	{
		if (set_name[i] == '/')
		{
			slash[s_count] = i;
			s_count++;
		}
	}
	// label the end of the set_name
	slash[s_count] = count;
	// the set_name is something like /a/b/c
	// the dataset will be under c,while /a/b must be created before
	// the "c" is in the "new_names",  "/a/b" in "names"
	if (s_count > 1)
	{
		// the case like /a/b, /a/b/c...
		for (i = 0; i < slash[s_count - 1]; i++)
		{
			name[i] = set_name[i];
		}
		name[slash[s_count - 1]] = '\0';
	}
	else
	{
		// the case like /a
		name[0] = '/';
		name[1] = '\0';
	}

	for (i = slash[s_count - 1] + 1; i < slash[s_count]; i++)
	{
		new_name[i - slash[s_count - 1] - 1] = set_name[i];
	}
	new_name[slash[s_count] - slash[s_count - 1] - 1] = '\0';
	//std::cout << set_name << std::endl;
	//std::cout << name << std::endl;
	//std::cout << new_name << std::endl;
	//show_arr(slash, 1, count);
	//std::cout << s_count << std::endl;
	// try to create /a/b
	create_h5_group(filename, name, trunc);

	hid_t file_id, group_id, dataset_id, dataspace_id;
	herr_t status;
	unsigned rank;

	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	group_id = H5Gopen1(file_id, name);

    if(row == 1 or column == 1)
    {
        hsize_t dims[1];
	    rank = 1;
        dims[0] = column*row;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }
    else
	{
        hsize_t dims[2];
	    rank = 2;
	    dims[0] = row;
	    dims[1] = column;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }
	
	dataset_id = H5Dcreate(group_id, new_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Gclose(group_id);
	status = H5Fclose(file_id);

	delete[] slash;
	delete[] name;
	delete[] new_name;
}

void write_h5(const char *filename, const char *set_name, const float*data, const int row, const int column, const bool trunc)
{
	int i, count, s_count;
	int *slash;
	char *name, *new_name;
	// the length of the set_name
	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}
	slash = new int[count] {};
	name = new char[count + 1];
	new_name = new char[count + 1];
	// find the slashes in the set_name
	s_count = 0;
	for (i = 0; i < count; i++)
	{
		if (set_name[i] == '/')
		{
			slash[s_count] = i;
			s_count++;
		}
	}
	// label the end of the set_name
	slash[s_count] = count;
	// the set_name is something like /a/b/c
	// the dataset will be under c,while /a/b must be created before
	// the "c" is in the "new_names",  "/a/b" in "names"
	if (s_count > 1)
	{
		for (i = 0; i < slash[s_count - 1]; i++)
		{
			name[i] = set_name[i];
		}
		name[slash[s_count - 1]] = '\0';
	}
	else
	{
		name[0] = '/';
		name[1] = '\0';
	}

	for (i = slash[s_count - 1] + 1; i < slash[s_count]; i++)
	{
		new_name[i - slash[s_count - 1] - 1] = set_name[i];
	}
	new_name[slash[s_count] - slash[s_count - 1] - 1] = '\0';
	//std::cout << set_name << std::endl;
	//std::cout << name << std::endl;
	//std::cout << new_name << std::endl;
	//show_arr(slash, 1, count);
	//std::cout << std::endl;
	// try to create /a/b if it doesn't exist 
	create_h5_group(filename, name, trunc);

	hid_t file_id, group_id, dataset_id, dataspace_id;
	herr_t status;
	hsize_t dims[2];
	unsigned rank;

	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	group_id = H5Gopen1(file_id, name);

    if(row == 1 or column == 1)
    {
        hsize_t dims[1];
	    rank = 1;
        dims[0] = column*row;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }
    else
	{
        hsize_t dims[2];
	    rank = 2;
	    dims[0] = row;
	    dims[1] = column;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }
	dataset_id = H5Dcreate(group_id, new_name, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Gclose(group_id);
	status = H5Fclose(file_id);

	delete[] slash;
	delete[] name;
	delete[] new_name;
}

void write_h5(const char *filename, const char *set_name, const int*data, const int row, const int column, const bool trunc)
{
	int i, count, s_count;
	int *slash;
	char *name, *new_name;
	// the length of the set_name
	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}
	slash = new int[count] {};
	name = new char[count + 1];
	new_name = new char[count + 1];
	// find the slashes in the set_name
	s_count = 0;
	for (i = 0; i < count; i++)
	{
		if (set_name[i] == '/')
		{
			slash[s_count] = i;
			s_count++;
		}
	}
	// label the end of the set_name
	slash[s_count] = count;
	// the set_name is something like /a/b/c
	// the dataset will be under c,while /a/b must be created before
	// the "c" is in the "new_names",  "/a/b" in "names"
	if (s_count > 1)
	{
		for (i = 0; i < slash[s_count - 1]; i++)
		{
			name[i] = set_name[i];
		}
		name[slash[s_count - 1]] = '\0';
	}
	else
	{
		name[0] = '/';
		name[1] = '\0';
	}

	for (i = slash[s_count - 1] + 1; i < slash[s_count]; i++)
	{
		new_name[i - slash[s_count - 1] - 1] = set_name[i];
	}
	new_name[slash[s_count] - slash[s_count - 1] - 1] = '\0';
	//std::cout << set_name << std::endl;
	//std::cout << name << std::endl;
	//std::cout << new_name << std::endl;
	//show_arr(slash, 1, count);
	//std::cout << std::endl;
	// try to create /a/b if it doesn't exist 
	create_h5_group(filename, name, trunc);

	hid_t file_id, group_id, dataset_id, dataspace_id;
	herr_t status;
	hsize_t dims[2];
	unsigned rank;
	
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	group_id = H5Gopen1(file_id, name);

    if(row == 1 or column == 1)
    {
        hsize_t dims[1];
	    rank = 1;
        dims[0] = column*row;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }
    else
	{
        hsize_t dims[2];
	    rank = 2;
	    dims[0] = row;
	    dims[1] = column;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }

	dataset_id = H5Dcreate(group_id, new_name, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Gclose(group_id);
	status = H5Fclose(file_id);

	delete[] slash;
	delete[] name;
	delete[] new_name;
}

void write_h5(const char *filename, const char *set_name, const long *data, const int row, const int column, const bool trunc)
{
	int i, count, s_count;
	int *slash;
	char *name, *new_name;
	// the length of the set_name
	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}
	slash = new int[count] {};
	name = new char[count + 1];
	new_name = new char[count + 1];
	// find the slashes in the set_name
	s_count = 0;
	for (i = 0; i < count; i++)
	{
		if (set_name[i] == '/')
		{
			slash[s_count] = i;
			s_count++;
		}
	}
	// label the end of the set_name
	slash[s_count] = count;
	
	// the set_name is something like /a/b/c
	// the dataset will be under c,while /a/b must be created before
	// the "c" is in the "new_names",  "/a/b" in "names"
	if (s_count > 1)
	{
		for (i = 0; i < slash[s_count - 1]; i++)
		{
			name[i] = set_name[i];
		}
		name[slash[s_count - 1]] = '\0';
	}
	else
	{
		name[0] = '/';
		name[1] = '\0';
	}

	for (i = slash[s_count - 1] + 1; i < slash[s_count]; i++)
	{
		new_name[i - slash[s_count - 1] - 1] = set_name[i];
	}
	new_name[slash[s_count] - slash[s_count - 1] - 1] = '\0';
	//std::cout << set_name << std::endl;
	//std::cout << name << std::endl;
	//std::cout << new_name << std::endl;
	//show_arr(slash, 1, count);
	//std::cout << std::endl;
	// try to create /a/b if it doesn't exist 
	create_h5_group(filename, name, trunc);

	hid_t file_id, group_id, dataset_id, dataspace_id;
	herr_t status;
	hsize_t dims[2];
	unsigned rank;

	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	group_id = H5Gopen1(file_id, name);

    if(row == 1 or column == 1)
    {
        hsize_t dims[1];
	    rank = 1;
        dims[0] = column*row;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }
    else
	{
        hsize_t dims[2];
	    rank = 2;
	    dims[0] = row;
	    dims[1] = column;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }
	
	dataset_id = H5Dcreate(group_id, new_name, H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Gclose(group_id);
	status = H5Fclose(file_id);

	delete[] slash;
	delete[] name;
	delete[] new_name;
}


// void read_fits(const char *filename, double *arr)
// {
// 	fitsfile *fptr;															/* pointer to the FITS file, defined in fitsio.h */
// 	int status = 0, nfound, anynull;
// 	long naxes[2], fpixel = 1, nbuffer, npixels, ii;
// 	double datamin, datamax, nullval = 0;

// 	fits_open_file(&fptr, filename, READONLY, &status);
// 	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);		/* read the NAXIS1 and NAXIS2 keyword to get image size */
// 	npixels = naxes[0] * naxes[1];	/* number of pixels in the image, python_arr[naxes[1], naxes[0]] */
// 	double *buffer = new double[npixels];										/* create a new array */
// 	fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval, buffer, &anynull, &status);
// 	for (ii = 0; ii < npixels; ii++) arr[ii] = buffer[ii];
// 	delete[] buffer;														/* (have to) delete the array */
// 	fits_close_file(fptr, &status);
// }

// void read_fits(const char *filename, float *arr)
// {
// 	fitsfile *fptr;			/* pointer to the FITS file, defined in fitsio.h */
// 	int status = 0, nfound, anynull;
// 	long naxes[2], fpixel = 1, nbuffer, npixels, ii;
// 	double datamin, datamax, nullval = 0;

// 	fits_open_file(&fptr, filename, READONLY, &status);
// 	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);		/* read the NAXIS1 and NAXIS2 keyword to get image size */
// 	npixels = naxes[0] * naxes[1];	/* number of pixels in the image, python_arr[naxes[1], naxes[0]] */
// 	float *buffer = new float[npixels];	/* create a new array */
// 	fits_read_img(fptr, TFLOAT, fpixel, npixels, &nullval, buffer, &anynull, &status);
// 	for (ii = 0; ii < npixels; ii++) arr[ii] = buffer[ii];
// 	delete[] buffer;		/* (have to) delete the array */
// 	fits_close_file(fptr, &status);
// }

// void read_fits(const char *filename, int *arr)
// {
// 	fitsfile *fptr;			/* pointer to the FITS file, defined in fitsio.h */
// 	int status = 0, nfound, anynull;
// 	long naxes[2], fpixel = 1, nbuffer, npixels, ii;
// 	double datamin, datamax, nullval = 0;

// 	fits_open_file(&fptr, filename, READONLY, &status);
// 	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);		/* read the NAXIS1 and NAXIS2 keyword to get image size */
// 	npixels = naxes[0] * naxes[1];	/* number of pixels in the image, python_arr[naxes[1], naxes[0]] */
// 	int *buffer = new int[npixels];	/* create a new array */
// 	fits_read_img(fptr, TINT, fpixel, npixels, &nullval, buffer, &anynull, &status);
// 	for (ii = 0; ii < npixels; ii++) arr[ii] = buffer[ii];
// 	delete[] buffer;		/* (have to) delete the array */
// 	fits_close_file(fptr, &status);
// }

// void write_fits(const char *filename, double *img, const int ysize, const int xsize)
// {
// 	fitsfile *fptr;		/* pointer to the FITS file; defined in fitsio.h */
// 	int status, ii, jj;
// 	long fpixel = 1, naxis = 2, nelements, exposure;
// 	long naxes[2] ;	    /* x, y */
// 	naxes[0] = xsize;
// 	naxes[1] = ysize;	 /* the same as numpy array */	
// 	status = 0;			/* initialize status before calling fitsio routines */
// 	fits_create_file(&fptr, filename, &status);		/* create new file */
// 	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);		 /* Create the primary array image (16-bit short integer pixels */
// 	//fits_update_key(fptr, TDOUBLE, "EXPOSURE", &exposure, "Total Exposure Time", &status);  /* Write a keyword; must pass the ADDRESS of the value */	
// 	nelements = xsize*ysize;         /* number of pixels to write */
// 	fits_write_img(fptr, TDOUBLE, fpixel, nelements, img, &status);     /* Write the array of integers to the image */
// 	fits_close_file(fptr, &status);              /* close the file */
// 	fits_report_error(stderr, status);      /* print out any error messages */	
// }

// void write_fits(const char *filename, float *img, const int ysize, const int xsize)
// {
// 	fitsfile *fptr;				/* pointer to the FITS file; defined in fitsio.h */
// 	int status, ii, jj;
// 	long fpixel = 1, naxis = 2, nelements, exposure;
// 	long naxes[2];			/* x, y */
// 	naxes[0] = xsize;
// 	naxes[1] = ysize;   	/* the same as numpy array */
// 	status = 0;		/* initialize status before calling fitsio routines */
// 	fits_create_file(&fptr, filename, &status);		/* create new file */
// 	fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);	 /* Create the primary array image (16-bit short integer pixels */
// 	nelements = xsize * ysize;         /* number of pixels to write */
// 	fits_write_img(fptr, TFLOAT, fpixel, nelements, img, &status);     /* Write the array of integers to the image */
// 	fits_close_file(fptr, &status);              /* close the file */
// 	fits_report_error(stderr, status);      /* print out any error messages */
// }

// void write_fits(const char *filename, int *img, const int ysize, const int xsize)
// {
// 	fitsfile *fptr;			/* pointer to the FITS file; defined in fitsio.h */
// 	int status, ii, jj;
// 	long fpixel = 1, naxis = 2, nelements;
// 	long naxes[2];		/* x, y */
// 	naxes[0] = xsize;
// 	naxes[1] = ysize; 	/* the same as numpy array */
// 	status = 0;	/* initialize status before calling fitsio routines */
// 	fits_create_file(&fptr, filename, &status);		/* create new file */
// 	fits_create_img(fptr, LONG_IMG, naxis, naxes, &status);	 /* Create the primary array image (16-bit short integer pixels */
// 	nelements = xsize * ysize;         /* number of pixels to write */
// 	fits_write_img(fptr, TINT, fpixel, nelements, img, &status);     /* Write the array of integers to the image */
// 	fits_close_file(fptr, &status);              /* close the file */
// 	fits_report_error(stderr, status);      /* print out any error messages */
// }
