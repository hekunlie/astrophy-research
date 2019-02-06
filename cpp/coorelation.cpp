#include<FQlib.h>

struct pts_info
{
	int idy, idx; // block id of the point
	double y, x; // the coordinates of the point

	double scale; // the length of the side of the square blocks
	int ny, nx; // the number of blocks along each axis 
	int blocks_num; // the numbers of total blocks
};

void find_block(pts_info *infos, const double radius_s, const double radius_e, const double *bound_y, const double *bound_x, int *block_mask)
{
	int i, j, k, lb, lby, lb_d;
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
	double distance[4];
	// "distance" stores the distance of the four vertexes of each block,
	// the sequence of the vertexes£º
	// | (y1,x1), (y1,x2) |
	// | (y2,x1), (y2,x2) |

	// find the minimum square that contains the target blocks
	nx_left = (int)((radius_e - x) / scale) + idx + 1;
	nx_right = (int)((radius_e + x) / scale) - idx + 1;
	ny_up = (int)((radius_e + y) / scale) - idy;

	nx_s = std::max(idx - nx_left, 0);
	nx_e = std::min(idx + nx_right + 1, nx);
	ny_e = std::min(idy + ny_up + 1, ny);

	// initialiize the mask
	for (i = 0; i < num; i++)
	{
		block_mask[i] = 0;
	}
	
	for (i = idy; i < ny_e; i++)
	{	
		lby = i * ny; // for speed
		for (j = nx_s; j < nx_e; j++)
		{
			lb = lby + j;
			lb_d = lb * 4;
			for (k = 0; k < 4; k++)
			{
				distance[k] = (bound_y[lb_d + k] - y)*(bound_y[lb_d + k] - y) + (bound_x[lb_d + k] - x)*(bound_x[lb_d + k] - x);
			}
			sort_arr(distance, 4, 1);// ascending order
			if (distance[3] < rs_sq or (distance[0] > re_sq and i not_eq idy and j not_eq idx))
			{
				/* "distance[3] < rs_sq": if the max distance between the vertex and the point is smaller than radius_s,
														it is not the target block.
					"distance[0] > re_sq and i not_eq idy and j not_eq idx":  if the 

				*/
				continue;
			}
			else
			{
				if (i > idy or (i == idy and j >= idx))
				{
					block_mask[lb] = 1;
				}
			}
		}
	}

	


}
main(int argc, char* argv[])
{

	return 0;
}