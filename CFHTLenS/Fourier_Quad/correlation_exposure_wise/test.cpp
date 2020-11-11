#include<functions_expo_wise.h>
#include<hk_iolib.h>

int main(int argc, char **argv)
{
    char file_path[400], hdf5_path[400];
    char set_name[20];
    sprintf(set_name,"/group_label");
    strcpy(file_path, argv[1]);
    int label=-99;

    std::ifstream infile;
	std::string str;
	std::stringstream strs;

	int i, j;
    int line_count;

    infile.open(file_path);
	line_count = 0;
    for(int i=0; i<855;i++)
    {
        str.clear();
        strs.clear();
        getline(infile, str);
        
        strs << str;

        strs >> hdf5_path;
        // std::cout<<hdf5_path<<" "<<label<<std::endl;

        read_h5(hdf5_path, set_name, &label);
        std::cout<<hdf5_path<<" "<<label<<std::endl;
        // std::cout << str << std::endl;
        for(j=0;j<8;j++)
        {strs >> hdf5_path;}
        // strs.clear();
    }
    infile.close();

    
    
    return 0;
}