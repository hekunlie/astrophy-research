//#include <stdafx.h>
//#include<algorithm> // sort(), std::max()
//#include<functional> // std::less, std::greater..
//#include<ciso646> // for "and, not, or, ..."
//#include<iostream>
#include<FQlib.h>

void find_block_(const pts_info *infos, const double radius, const double *bound_y, const double *bound_x, int *block_mask)
{
	int i, j, k, lb, lby, lb_d, lb_d_, seq = 0;
	int idy = infos->idy;// block id of the point
	int idx = infos->idx;
	double y = infos->y;// the coordinates of the point
	double x = infos->x;
	double scale = infos->scale;// the length of the side of the square blocks
	int ny = infos->ny; // the number of blocks along each axis 
	int nx = infos->nx;
	int num = infos->blocks_num;// the numbers of total blocks

	int nx_left, nx_right, ny_up, ny_down; // the max blocks number in each direction
												// away from the point within radius_e
	int nx_s, nx_e, ny_s, ny_e;

	double re_sq = radius * radius;
	double dx, dy;
	double distance[4];
	// "distance" stores the distance of the four vertexes of each block,
	// the sequence of the vertexes£º
	// | (y1,x1), (y1,x2) |
	// | (y2,x1), (y2,x2) |

	// find the minimum square that contains the target blocks
	nx_left = (int)((radius - x + bound_x[0]) / scale) + idx + 1;
	nx_right = (int)((radius + x - bound_x[0]) / scale) - idx;
	ny_up = (int)((radius + y - bound_y[0]) / scale) - idy;
	ny_down = (int)((radius - y + bound_y[0]) / scale) + idy + 1;

	nx_s = std::max(idx - nx_left, 0);
	nx_e = std::min(idx + nx_right + 1, nx);
	ny_s = std::max(idy - ny_down, 0);
	ny_e = std::min(idy + ny_up + 1, ny);

	// initialiize the mask
	for (i = 0; i < num; i++)
	{
		block_mask[i] = -1;
	}

	for (i = ny_s; i < ny_e; i++)
	{
		lby = i * nx; // for speed
		for (j = nx_s; j < nx_e; j++)
		{
			lb = lby + j; // label of the block
			lb_d = lb * 4;
			for (k = 0; k < 4; k++)
			{
				lb_d_ = lb_d + k;
				dy = bound_y[lb_d_] - y;
				dx = bound_x[lb_d_] - x;
				distance[k] = dy * dy + dx * dx;
			}
			sort_arr(distance, 4, 1); // ascending order
			if (distance[0] > re_sq and i not_eq idy and j not_eq idx)
			{
				//	"distance[0] > re_sq and i not_eq idy and j not_eq idx":  
				//  if the minimum distance between the vertexes of a block ,
				//	which is not on the cross centers on  (idy, idx) , is larger 
				//	than radius_e, it is not the target one. 
				//	The blocks on the cross centers on  (idy, idx) are the 
				//	targets.
				continue;
			}
			else
			{
				block_mask[seq] = lb;
				seq++;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	int ny, nx;
	int i, j, k;
	int tag, bin_label;
	double scale = 0.15;
	double x, y;
	double rad_s, rad_e;
	int shape[2];
	int area_id = 1, num;
	int ra_bin_num, dec_bin_num;

	char path[50], data_path[200];
	char buf[100], set_name[50], attr_name[50];


	// sky area, the label of the foreground galaxy and the radius for searching
	area_id = std::atoi(argv[1]);
	tag = std::atoi(argv[2]);
	rad_e = std::atof(argv[3]);


	// for the mask of needed grid
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_grid.hdf5");
	sprintf(set_name, "/background/w_%d", area_id);
	sprintf(attr_name, "grid_shape");
	read_h5_attrs(data_path, set_name, attr_name, shape, "g");
	ny = shape[0];
	nx = shape[1];
	int *mask = new int[ny*nx]{};
	std::cout << "Grid: " << ny << " " << nx << " " << scale << std::endl;
	// the boundary of blocks
	sprintf(set_name, "/background/w_%d/block_boundx", area_id);
	sprintf(attr_name, "shape");	
	read_h5_attrs(data_path, set_name, attr_name, shape, "d");
	double *boundx = new double[shape[0]*shape[1]]{};
	read_h5(data_path, set_name, boundx);

	sprintf(set_name, "/background/w_%d/block_boundy", area_id);
	sprintf(attr_name, "shape");
	read_h5_attrs(data_path, set_name, attr_name, shape, "d");
	double *boundy = new double[shape[0] * shape[1]]{};
	read_h5(data_path, set_name, boundy);

	sprintf(set_name, "/background/w_%d/RA_bin", area_id);
	sprintf(attr_name, "shape");
	read_h5_attrs(data_path, set_name, attr_name, shape, "d");
	double *ra_bin = new double[shape[0] * shape[1]]{};
	read_h5(data_path, set_name, ra_bin);
	ra_bin_num = shape[0] - 1;

	sprintf(set_name, "/background/w_%d/DEC_bin", area_id);
	sprintf(attr_name, "shape");
	read_h5_attrs(data_path, set_name, attr_name, shape, "d");
	double *dec_bin = new double[shape[0] * shape[1]]{};
	read_h5(data_path, set_name, dec_bin);
	dec_bin_num = shape[0] - 1;

	
	// the foreground 
	sprintf(set_name, "/foreground/w_%d/RA", area_id);
	sprintf(attr_name, "shape");
	read_h5_attrs(data_path, set_name, attr_name, shape, "d");
	double *ra = new double[shape[0]]{};
	read_h5(data_path, set_name, ra);

	sprintf(set_name, "/foreground/w_%d/DEC", area_id);
	sprintf(attr_name, "shape");
	read_h5_attrs(data_path, set_name, attr_name, shape, "d");
	double *dec = new double[shape[0]]{};
	read_h5(data_path, set_name, dec);

	sprintf(set_name, "/foreground/w_%d/COS_DEC", area_id);
	sprintf(attr_name, "shape");
	read_h5_attrs(data_path, set_name, attr_name, shape, "d");
	double *cos_dec = new double[shape[0]]{};
	read_h5(data_path, set_name, cos_dec);

	num = shape[0];

	if (tag >= num)
	{
		tag = num - 1;
	}

	pts_info pts;

	// the information of the foreground galaxy
	histogram2d_s(dec[tag], ra[tag], dec_bin, ra_bin, dec_bin_num, ra_bin_num, bin_label);
	std::cout << dec_bin_num << " " << ra_bin_num <<" "<<bin_label<< std::endl;
	show_arr(dec_bin, 1, dec_bin_num + 1);
	std::cout << std::endl;
	show_arr(ra_bin, 1, ra_bin_num + 1);
	std::cout << std::endl;

	pts.idx = bin_label % nx;
	pts.idy = bin_label / nx;
	pts.nx = nx;
	pts.ny = ny;
	pts.scale = scale;
	pts.blocks_num = ny * nx;
	pts.x = ra[tag];
	pts.y = dec[tag];
	std::cout << "Position: "<<dec[tag]<<" "<<ra[tag] <<" Grid: "<<pts.idy<< " " << pts.idx<< " (" << pts.ny<<" "<<pts.nx <<")"<< std::endl;

	find_block(&pts, rad_e, boundy, boundx, mask);

	std::cout << std::endl;
	sprintf(buf, "PTS: %d, %d", pts.idy, pts.idx);
	std::cout << buf << std::endl;
	std::cout << scale << std::endl;
	std::cout << pts.y << ", " << pts.x << std::endl;


	sprintf(path, "!mask.fits");
	write_fits(path, mask, ny, nx);

	delete[] boundx;
	delete[] boundy;
	delete[] mask;
	delete[] ra;
	delete[] dec;
	delete[] cos_dec;
	gsl_free();
	return 0;
}