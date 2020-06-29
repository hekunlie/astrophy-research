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


    char data_path[200],mask_path[200], set_name[30];
    std::string detect_info;
    MY_FLOAT *image;
    int *big_mask;
    MY_FLOAT *stamp;
    int *mask;
    int i,j,k, nx, ny;
    int stamp_tag;
    int stamp_size;
    int image_label;
    int detect_label;

    fq_paras_float all_paras;
    

    nx = 100;
    ny = 100;
    stamp_size = 44;

	all_paras.gal_noise_sig = 0.02507446;
	all_paras.psf_noise_sig = 0;
	all_paras.stamp_size = stamp_size;
	all_paras.max_source = 5000;
	all_paras.area_thresh = 5;
	all_paras.detect_thresh = 0.02507446*1.5;
	all_paras.img_x = 2048;
	all_paras.img_y = 1489;
	all_paras.max_distance = 6; 

    image = new MY_FLOAT[1489*2048];
    mask = new int[1489*2048];

    sprintf(data_path,"frame-0.fits");
    read_fits(data_path, image);
    galaxy_finder(image, mask, &all_paras, false, detect_label, detect_info);
    sprintf(data_path,"!mask.fits");
    write_fits(data_path, mask, stamp_size, stamp_size);
    std::cout<<all_paras.gal_size<<" "<<all_paras.gal_osnr<<std::endl;
/*
    image = new MY_FLOAT[nx*ny*stamp_size*stamp_size];
    stamp = new MY_FLOAT[stamp_size*stamp_size];
    mask = new int[stamp_size*stamp_size];
    big_mask = new int[nx*ny*stamp_size*stamp_size];


    sprintf(data_path,"/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/galsim_dimmer/0/gal_chip_0499.fits");
    read_fits(data_path, image);

    initialize_arr(big_mask, nx*ny*stamp_size*stamp_size, 0);

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            initialize_arr(mask, stamp_size*stamp_size, -1);
            initialize_arr(stamp, stamp_size*stamp_size, 0);


            stamp_tag = i*nx + j;
            segment(image, stamp, stamp_tag, stamp_size, nx, ny);
            if(stamp_tag > -1)
            {
                galaxy_finder(stamp, mask, &all_paras, false, detect_label,detect_info);

            }
            mask[0] = detect_label;
            mask[1] = i;
            mask[2] = j;

            stack(big_mask,mask,stamp_tag, stamp_size, ny, nx);

        }
    }

    sprintf(data_path,"!/mnt/ddnfs/data_users/hkli/selection_bias/gal_0499_mask.fits");
    write_fits(data_path, big_mask, ny*stamp_size, nx*stamp_size);

    sprintf(data_path,"!/mnt/ddnfs/data_users/hkli/selection_bias/gal_0499.fits");
    write_fits(data_path, image, ny*stamp_size, nx*stamp_size);

    delete[] mask;
    delete[] big_mask;
    delete[] stamp;
    delete[] image;
*/
   // MPI_Finalize();

    return 0;
}