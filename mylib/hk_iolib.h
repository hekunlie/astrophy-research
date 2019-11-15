#ifndef __HK_IOLIB__
#define __HK_IOLIB__

//#pragma once

#include<iostream>
#include<fitsio.h>
#include<fstream> // ifstream,ofstream
#include<sstream> //std::stringstream
#include<hdf5.h>
#include<string>
#include<string.h>
#include<ctime>
#include<ciso646> // for "and, not, or, ..."
#include<sys/stat.h> // for stat()


/********************************************************************************************************************************************/
/* file reading and writting*/
/********************************************************************************************************************************************/

void char_to_str(const char *char_in, std::string &string_out);//checked
/* convert a char to string */

bool file_exist(const char *filename);//checked
/* check the existence of a file , will called by create_h5_group()*/

void write_log(char *filename, char *inform); //checked
/* write char to log file */

void read_config(const std::string path, const std::string target_section, const std::string target_para, std::string &para);//checked
/*  find the content or value of a parameter, "name", in section,"section".
	 the config file muse be written as follow:
	 [section a]			 # section name must be enclosed by "[ ]"
	 para_a = aaaaa    #  the content in each section must consist of parameter name, "=" and the content of the parameter
								# and there is a space in each side of "=" to separate the the name and the content

	 path: the config file path
	 target_section: the section name of which to be read
	 target_para: the target parameter to read
	 para: the value of the target_para
*/

void read_para(const std::string path, const std::string name, int &para);//checked
void read_para(const std::string path, const std::string name, double &para);
void read_para(const std::string path, const std::string name, float &para);
/* read the parameters ("name") value from parameter file */

void read_text(const char *filename, double *data_buf, const int data_col, const int skipline = 0);
/* read the data in the text file into the data_buf.														*/
/* skipline: 1, skip the first line of annonation,   else 0												*/
/* the function itself does not check the type of the content in the file,					*/
/* so nothing except the data and the first line of annotation should be in the file	*/

void read_text(const std::string path, double *arr, const int start_line, const int read_lines, const int read_cols);
void read_text(const std::string path, double *arr, const int read_lines);//checked
void read_text(const std::string path, float *arr, const int read_lines);
void read_text(const std::string path, int *arr, const int read_lines);
/* read data from txt file which should be just one column.
	the read_lines limits the maximum lines to read
*/
void write_text(const char*filename, double *data_buf, const int data_row, const int data_col, const int mode=0);
void write_text(const char*filename, double *data_buf, const int data_row, const int data_col, const char * comment, const int mode=0);
/* write array into file															*/
/* mode: 0, truncate the file before writting							*/
/*				1, append to the end of the file							*/
/* comment: something like "#....", write to the first line	*/	

void read_h5_datasize(const char *filename, const char *set_name, int &elem_num);//checked
void read_h5_datasize(const char *filename, const char *set_name, long &elem_num);//checked
/* read the element number in the dataset 
	if it fails, elem_num=-1 
*/
void read_h5(const char *filename, const char *set_name, double *data);//checked
void read_h5(const char *filename, const char *set_name, float *data);
void read_h5(const char *filename, const char *set_name, int *data);//checked
void read_h5(const char *filename, const char *set_name, long *data);//checked
/* if the length of arr is longer than the data, the rest element  of "arr" will not be changed 
	initializing the arr before each reading is highly recommended.
*/
void read_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, double *buff, std::string flag);//checked
void read_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, float *buff, std::string flag);//checked
void read_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, int *buff, std::string flag);//checked

void write_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, const double *attrs_buffer, const int buffer_len, std::string flag);
void write_h5_attrs(const char *filename, const char *set_name, const char *attrs_name, const int *attrs_buffer, const int buffer_len, std::string flag);
/* the attributes must be attached to the non-root directory, or it will raise the error "/".
	flag: "d" for data set, "g" for data group
*/
void create_h5_group(const char *filename, const char *set_name, const bool trunc);//checked
/* create the data group 
*/
void write_h5(const char *filename, const char *set_name, const double *data, const int row, const int column, const bool trunc);//checked
void write_h5(const char *filename, const char *set_name, const float *data, const int row, const int column, const bool trunc);//checked
void write_h5(const char *filename, const char *set_name, const int *data, const int row, const int column, const bool trunc);//checked
void write_h5(const char *filename, const char *set_name,  const long *data, const int row, const int column, const bool trunc);//checked
/* write the hdf5 file
	if trunc == true, the file will be truncated before data are writted into the file
	else, write directly
*/

void read_fits(const char *filename, double *arr);//checked
void read_fits(const char *filename, float *arr);
void read_fits(const char *filename, int *arr);//checked
void write_fits(const char *filename, double *img, const int ysize, const int xsize);//checked
void write_fits(const char *filename, float *img, const int ysize, const int xsize);
void write_fits(const char *filename, int *img, const int ysize, const int xsize);//checked
/* read and write the array to  fits file, 
	be careful with the datetype "TINT" and "LONG_IMG"!!! 
	the length of INT may be different in different platform,
	however, the "TINT32BIT" doesn't work with "LONG_IMG".

	the size of each axises should be inputted，ARRAY(y, x)
*/

#endif