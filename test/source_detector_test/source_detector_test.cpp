#include <FQlib.h>
//#include <mpi.h>
#include<hk_iolib.h>


#define MY_FLOAT float
int main(int argc, char **argv)
{
    // int rank, numprocs, namelen;
	// char processor_name[MPI_MAX_PROCESSOR_NAME];

	// MPI_Init(&argc, &argv);
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	// MPI_Get_processor_name(processor_name, &namelen);


    char data_path[200],mask_path[200], set_name[40];
    std::string detect_info;
    MY_FLOAT *image;
    int *mask;
    MY_FLOAT *stamp;

    int i,j,k, nx, ny, pixel_num;
    int tag_s, tag_e;
    int image_label;
    int detect_label;

    fq_paras_float all_paras;
    
    nx = 2048;
    ny = 1489;
    pixel_num = nx*ny;

    int elem_unit = 8;
    int source_num;
    int *source_x = new int[pixel_num] {};
	int *source_y = new int[pixel_num] {};

	all_paras.gal_noise_sig = 0.02507446;
	all_paras.psf_noise_sig = 0;
	all_paras.max_source = 5000;
	all_paras.area_thresh = 5;
	all_paras.detect_thresh = 0.02507446*1.5;
	all_paras.img_x = 2048;
	all_paras.img_y = 1489;
	all_paras.max_distance = 6; 

    image = new MY_FLOAT[nx*ny];
    mask = new int[nx*ny];
    initialize_arr(mask, pixel_num, 0);

    MY_FLOAT *source_data;
    MY_FLOAT *source_para = new MY_FLOAT[elem_unit*all_paras.max_source]{}; // determined by 'max_sources' in paras.
    
    sprintf(data_path,"frame-0.fits");
    read_fits(data_path, image);
	
	source_detector(image, source_x, source_y, source_para, &all_paras, true, source_num, detect_info);

    for (i = 0; i < source_num; i++)
	{
		// start point of source_y(x) of i'th source
		if (i > 0)
		{
			tag_s += source_para[(i - 1)*elem_unit];
		}
		else
		{
			tag_s = 0;
		}
		for (j = tag_s; j < tag_s + source_para[i * elem_unit]; j++)
		{
			// detection mask
			mask[source_y[j] * nx + source_x[j]] = 1;
		}
	
	}
 
    source_data = new MY_FLOAT[source_num*6];
    for (i=0;i<source_num;i++)
    {
        source_data[i*6] = source_para[i*elem_unit + 1];
        source_data[i*6 + 1] = source_para[i*elem_unit + 2];
        source_data[i*6 + 2] = source_para[i*elem_unit + 4];
        source_data[i*6 + 3] = source_para[i*elem_unit];
        source_data[i*6 + 4] = source_para[i*elem_unit + 6];
        source_data[i*6 + 5] = source_para[i*elem_unit + 5];
    }
    sprintf(data_path,"data.hdf5");
    sprintf(set_name, "/data");
    write_h5(data_path, set_name, source_data, source_num, 6, true);
    sprintf(data_path,"!mask.fits");
    write_fits(data_path, mask, ny, nx);
    
    return 0;
}