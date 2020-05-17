#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<FQlib.h>


int main(int argc, char** argv)
{
	std::string path = "/home/hkli/work/envs/envs.dat";
	std::string section, para, para_con;
	
	char_to_str(argv[1], section);
	char_to_str(argv[2], para);

	std::cout << "To find the " << para << " in " << section << std::endl;
	read_config(path, section, para, para_con);
	std::cout <<"Found: "<<section<<" "<<para<<" "<< para_con << std::endl;

	return 0;
}