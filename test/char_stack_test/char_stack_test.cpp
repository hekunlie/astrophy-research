#include<iostream>
#include<sstream>
#include<string>
#include<hk_iolib.h>

void char_stack_(char **char_in, const int chars_num, char *char_out,std::string linker)
{
    std::string ss, temp;
    char_to_str(char_in[0],ss);
    for(int i=1; i<chars_num; i++)
    {   
        char_to_str(char_in[i],temp);
	    ss += linker + temp;
    }
    strcpy(char_out, ss.c_str());
}

int main(int argc, char *argv[])
{
    char *nm[10];
    int i,j,k;

    for(i=1;i<argc;i++)
    {   
        nm[i-1] = new char[50];
        strcpy(nm[i-1], argv[i]);
    }
    char final_char[500];
    std::string linker;
    linker = "=";

    char_stack(nm, argc-1, final_char, linker);
    std::cout<<"After "<<final_char<<std::endl;

    return 0;
}

