#include<iostream>
#include<fstream>
#include<sstream>

void read_file(const char *filename, char **buffer, int &lines)
{
    std::ifstream infile;
	std::string str;
	std::stringstream strs;

	int i, j;
    int line_count;

	infile.open(filename);
	line_count = 0;

    while (!infile.eof())
    {
        str.clear();
        strs.clear();
        getline(infile, str);
                
        strs << str;
        strs >> buffer[line_count];

        // std::cout << str << std::endl;
        line_count += 1;
    }
    lines = line_count;
}

int main()
{
    int lines;
    char *buffer[100];
    char filename[200];
    int i;
    for(i=0;i<100;i++){buffer[i] = new char[20];}

    sprintf(filename,"D:/shear.dat");
    read_file(filename, buffer, lines);
    std::cout<<"Lines: "<<lines<<std::endl;

    for(i=0; i<lines; i++)
    {
        std::cout<<"line: "<<i<<" "<<buffer[i]<<std::endl;
    }

    return 0;
}