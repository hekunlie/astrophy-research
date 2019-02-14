//#include <stdafx.h>
//#include<algorithm> // sort(), std::max()
//#include<functional> // std::less, std::greater..
//#include<ciso646> // for "and, not, or, ..."
//#include<iostream>
#include<FQlib.h>
void sort_arr_(double* arr, int size, int order = 1)
{
	if (order == 1)
	{
		std::sort(arr, arr + size, std::less<double>());
	}
	else
	{
		std::sort(arr, arr + size, std::greater<double>());
	}
}

struct pts_info_
{
	int idy, idx; // block id of the point
	double y, x; // the coordinates of the point

	double scale; // the length of the side of the square blocks
	int ny, nx; // the number of blocks along each axis 
	int blocks_num; // the numbers of total blocks
};

void block_bound_(const double scale, const int ny, const int nx, double *bound_y, double *bound_x)
{
	int i, j, k;
	for (i = 0; i < ny; i++)
	{
		for (j = 0; j < nx; j++)
		{
			// the sequence of the vertexes of the blocks£º
			// | (y1,x1), (y1,x2) |
			// | (y2,x1), (y2,x2) |
			bound_y[ 4*(i * nx + j)] = i*scale;
			bound_y[4 * (i * nx + j)+1] = i * scale;
			bound_y[4 * (i * nx + j)+2] = (i+1) * scale;
			bound_y[4 * (i * nx + j)+3] = (i+1) * scale;

			bound_x[4 * (i * nx + j)] = j * scale;
			bound_x[4 * (i * nx + j) + 1] = (j + 1) * scale;
			bound_x[4 * (i * nx + j) + 2] = j* scale;
			bound_x[4 * (i * nx + j) + 3] = (j + 1) * scale;
		}
	}
}

void find_block_(const pts_info *infos, const double radius_s, const double radius_e, const double *bound_y, const double *bound_x, int *block_mask)
{
	int m;
	int i, j, k, lb, lby, lb_d, lb_d_, seq=0;
	int idy = infos->idy;// block id of the point
	int idx = infos->idx;
	double y = infos->y;// the coordinates of the point
	double x = infos->x;
	double scale = infos->scale;// the length of the side of the square blocks
	int ny = infos->ny; // the number of blocks along each axis 
	int nx = infos->nx; 
	int num = infos->blocks_num;// the numbers of total blocks

	int nx_left, nx_right, ny_up; // the max blocks number in each direction
												// away from the point within radius_e
	int nx_s, nx_e, ny_e;

	double rs_sq = radius_s * radius_s; 
	double re_sq = radius_e * radius_e;
	double dx, dy;
	double distance[4];
	// "distance" stores the distance of the four vertexes of each block,
	// the sequence of the vertexes£º
	// | (y1,x1), (y1,x2) |
	// | (y2,x1), (y2,x2) |
	char buf[100];
	// find the minimum square that contains the target blocks
	nx_left = int(((radius_e - x) / scale) + idx + 1);
	nx_right = int(((radius_e + x) / scale) - idx );
	ny_up = int(((radius_e + y) / scale) - idy);
	std::cout << std::endl;
	sprintf(buf, "%d, %d, %d", nx_left, nx_right, ny_up);
	std::cout << buf << std::endl;
	nx_s = std::max(idx - nx_left, 0);
	nx_e = std::min(idx + nx_right + 1, nx);
	ny_e = std::min(idy + ny_up + 1, ny);
	
	std::cout << std::endl;
	sprintf(buf, "%d, %d, %d", nx_s, nx_e, ny_e);
	std::cout << buf << std::endl;
	// initialiize the mask
	for (i = 0; i < num; i++)
	{
		block_mask[i] = -1;
	}
	
	for (i = idy; i < ny_e; i++)
	{	
		lby = i * nx; // for speed
		for (j = nx_s; j < nx_e; j++)
		{
			lb = lby + j; // label of the block
			lb_d = lb * 4;
			for (k = 0; k < 4; k++)
			{	
				lb_d_ = lb_d +k;
				dy = bound_y[lb_d_] - y;
				dx = bound_x[lb_d_] - x;
				distance[k] = dy*dy + dx*dx;				
			}

			sort_arr(distance, 4, 1); // ascending order
			//for (m = 0; m < 4; m++)
			//{
			//	std::cout << distance[m] << ", ";
			//}
			//std::cout << std::endl;

			if (distance[3] < rs_sq or (distance[0] > re_sq and i not_eq idy and j not_eq idx))
			{
				//"distance[3] < rs_sq": 
				// if the max distance between the vertex and the point is smaller than radius_s,
				//	it is not the target block.
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
				if (i > idy or (i == idy and j >= idx))
				{
					block_mask[seq] = lb;
					seq++;
				}
			}
		}
	}
}

int main(int argc, char* argv[])
{

	int ny = 30, nx = 30;
	int i, j, k;
	int tag;
	double scale = 3;
	double x, y;
	double rad_s, rad_e;
	double *boundx = new double[4 * nx*ny]{};
	double *boundy = new double[4 * nx*ny]{};
	int *mask = new int[ny*nx]{};
	char path[50];
	char buf[100];
	tag = std::atoi(argv[1]);
	rad_s = std::atof(argv[2]);
	rad_e = std::atof(argv[3]);
	int seed = 100;
	gsl_rng_initialize(seed+tag);
	block_bound(scale, ny, nx, boundy, boundx);
	
	std::cout << "X: " << std::endl;
	for (i = 0; i < 4; i++)
	{
		std::cout << boundx[4 * tag + i] << ",  ";
	}
	std::cout << std::endl;
	std::cout << "Y: " << std::endl;
	for (i = 0; i < 4; i++)
	{
		std::cout << boundy[4 * tag + i] << ",  ";
	}
	
	pts_info pts;
	pts.idx = tag % nx;
	pts.idy = tag / nx;
	pts.nx = nx;
	pts.ny = ny;
	pts.scale = scale;
	pts.blocks_num = ny * nx;
	pts.x = rand_uniform(boundx[4 * tag], boundx[4 * tag + 1]);
	pts.y = rand_uniform(boundy[4 * tag], boundy[4 * tag + 2]);
	find_block(&pts, rad_s, rad_e, boundy, boundx, mask);
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
	gsl_rng_free();
	return 0;
}