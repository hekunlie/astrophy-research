#include<functions_expo_wise.h>

int main(int argc, char **argv)
{
    char file_path[400];
    strcpy(file_path, argv[1]);
    line_count(file_path);
    return 0;
}