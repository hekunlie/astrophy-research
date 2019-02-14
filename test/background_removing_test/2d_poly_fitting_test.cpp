#include"FQlib.h"
#include "Python.h"

void cov_martix_2d_(const double *x, const double *y, const double *fxy, const int data_num, const int order, double *cov_matrix, double *f_vector)
{
	// order >= 1;
	if (order < 1)
	{
		std::cout << "Order < 1 !!!" << std::endl;
		exit(0);
	}

	int terms = (order + 1)*(order + 2) / 2;
	// rows * data number, "row+1" represent the power of x(y)
	int size = 2 * order * data_num;
	// the length of the square matrix of powers of x and y
	int mask_size = (2 * order + 1)*(2 * order + 1);
	// the number of possible "x^n*y^m" term
	int xys_size = (2 * order - 1)*order * data_num, xy_seq;
	int i, j, k, m, n;
	int y_odr, x_odr, xy_count = 0;
	double dn_sum = 0;

	// the powers of x and y of each turn in the polynomial
	int *pow_y = new int[terms] {};
	int *pow_x = new int[terms] {};
	double *ys = new double[size] {};
	double *xs = new double[size] {};

	// xy_pow_mask stores the location of "x^n*y^m" in the array "xys"
	// xy_pow_mask is used as 2 dimessional array , (y order, x order)
	int *xy_pow_mask = new int[mask_size] {};
	double *xys = new double[xys_size] {};

	for (i = 0; i < mask_size; i++)
	{
		xy_pow_mask[i] = -1;
	}

	k = 0;
	for (i = 0; i < order + 1; i++)
	{
		for (j = 0; j < i + 1; j++)
		{
			pow_x[k] = i - j;
			pow_y[k] = j;
			k++;
		}
	}

	// calculate each order of x and y for the covariance matrix
	// to x^(order*order), y^(order*order)
	for (i = 0; i < 2 * order; i++)
	{

		if (0 == i) // the first order
		{
			for (k = 0; k < data_num; k++)
			{
				ys[k] = y[k];
				xs[k] = x[k];
			}
		}
		else // the higher order  
		{
			m = i * data_num;
			n = (i - 1)*data_num;
			for (k = 0; k < data_num; k++)
			{
				ys[m + k] = ys[n + k] * y[k];
				xs[m + k] = xs[n + k] * x[k];
			}
		}
	}

	// calculate the covariance matrix		
	//row
	for (i = 0; i < terms; i++)
	{
		//column
		for (j = 0; j < terms; j++)
		{
			y_odr = pow_y[j] + pow_y[i];
			x_odr = pow_x[j] + pow_x[i];

			std::cout << x_odr << " " << y_odr << ", ";
			if (j == terms - 1)
			{
				std::cout << std::endl;
			}

			// the terms without x^n
			if (0 == x_odr)
			{
				//the cov[0,0] = number of data points
				if (0 == y_odr)
				{
					dn_sum = data_num;
				}
				// the y^n terms
				else
				{
					sum_arr(ys, size, (y_odr - 1)*data_num, y_odr*data_num, dn_sum);
				}
			}
			// the terms with x^n
			else
			{
				// the x^n terms
				if (0 == y_odr)
				{
					sum_arr(xs, size, (x_odr - 1)*data_num, x_odr*data_num, dn_sum);
				}
				// the x^n*y^m terms
				else
				{
					// if this term has been gotten
					if (xy_pow_mask[y_odr*(2*order + 1) + x_odr] > -1)
					{
						xy_seq = xy_pow_mask[y_odr*(2 * order + 1) + x_odr];
						sum_arr(xys, xys_size, xy_seq*data_num, (xy_seq + 1)*data_num, dn_sum);
						if (x_odr == 5 && y_odr == 1)
						{
							std::cout << " GOT|"<<i<<"|"<<j<<"| " << x_odr << "|" << y_odr << "| " << xy_seq << "| ";
						}
					}
					// if not
					else
					{
						xy_pow_mask[y_odr*(2 * order + 1) + x_odr] = xy_count;
						if (x_odr == 5 && y_odr == 1)
						{
							std::cout <<" |"<< i << "|" << j << "| " << " |" << x_odr << "|" << y_odr << "| " << xy_count << "| ";
						}
						for (k = 0; k < data_num; k++)
						{
							xys[xy_count*data_num + k] = xs[(x_odr - 1)*data_num + k] * ys[(y_odr - 1)*data_num + k];
						}
						sum_arr(xys, xys_size, xy_count*data_num, (xy_count + 1)*data_num, dn_sum);
						xy_count++;
					}
				}
			}
			cov_matrix[i*terms + j] = dn_sum;
		}
	}

	// the vector of the right side
	for (i = 0; i < terms; i++)
	{
		dn_sum = 0;
		y_odr = pow_y[i];
		x_odr = pow_x[i];

		//std::cout << x_odr << " " << y_odr << ", ";
		// the terms without x^n
		if (0 == x_odr)
		{
			if (0 == y_odr)
			{
				sum_arr(fxy, data_num, 0, data_num, dn_sum);
			}
			// the f(x,y)*y^n terms
			else
			{
				for (j = 0; j < data_num; j++)
				{
					dn_sum += fxy[j] * ys[(y_odr - 1)*data_num + j];
				}
			}
		}
		// the terms with x^n
		else
		{
			// the f(x,y)*x^n terms
			if (0 == y_odr)
			{
				for (j = 0; j < data_num; j++)
				{
					dn_sum += fxy[j] * xs[(x_odr - 1)*data_num + j];
				}
			}
			// the f(x,y)* x^n*y^m terms
			else
			{
				xy_seq = xy_pow_mask[y_odr*(2*order + 1) + x_odr];
				for (j = 0; j < data_num; j++)
				{
					dn_sum += fxy[j] * xys[xy_seq*data_num + j];
				}
			}
		}
		f_vector[i] = dn_sum;
	}

	//std::cout <<"beign cov:   "<<terms<< std::endl;
	//for (i = 0; i < terms; i++)
	//{
	//	for (j = 0; j < terms; j++)
	//	{
	//		std::cout << cov_matrix[i*terms + j]<<", ";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << "beign vector:   " << terms << std::endl;
	//for (j = 0; j < terms; j++)
	//{
	//	std::cout << f_vector[j] << ", ";
	//}
	//std::cout << std::endl;

	delete[] xys;
	delete[] pow_x;
	delete[] pow_y;
	delete[] ys;
	delete[] xs;
	delete[] xy_pow_mask;
}

void poly_fit_2d_(const double *x, const double *y, const double *fxy, const int data_num, const int order, double *coeffs)//bug!!!
{
	// order >= 1
	if (order < 1)
	{
		std::cout << "Order < 1 !!!" << std::endl;
		exit(0);
	}

	int terms = (order + 1)*(order + 2) / 2;
	int i, s, j, k;
	double temp_scale = 1;

	double *cov_matrix = new double[terms*terms]{};
	double *f_vector = new double[terms] {};
	double *check = new double[terms*terms]{};

	cov_martix_2d_(x, y, fxy, data_num, order, cov_matrix, f_vector);


	gsl_matrix *cov_mat = gsl_matrix_alloc(terms, terms);
	gsl_vector * vect_b = gsl_vector_alloc(terms);
	gsl_matrix *mat_inv = gsl_matrix_alloc(terms, terms);
	gsl_permutation *permu = gsl_permutation_alloc(terms);
	gsl_vector *pamer = gsl_vector_alloc(terms);


	std::cout << std::endl;
	std::cout << "beign vector:   " << terms << std::endl;
	for (j = 0; j < terms; j++)
	{
		gsl_vector_set(vect_b, j, f_vector[j] / temp_scale);
		std::cout << gsl_vector_get(vect_b, j) << ", ";
		//std::cout << gsl_vector_get(&vect_b.vector, j)<<"("<< f_vector[j]<<")" << ", ";
	}
	std::cout << std::endl;

	std::cout << "beign cov:   " << terms << std::endl;
	for (i = 0; i < terms; i++)
	{
		for (j = 0; j < terms; j++)
		{
			gsl_matrix_set(cov_mat, i, j, cov_matrix[i*terms + j] / temp_scale);
			std::cout << gsl_matrix_get(cov_mat, i, j) << ", ";
			//std::cout << gsl_matrix_get(&mat.matrix, i, j) << "(" << cov_matrix[i*terms+j] << ")" << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	gsl_linalg_LU_decomp(cov_mat, permu, &s);
	gsl_linalg_LU_invert(cov_mat, permu, mat_inv);
	gsl_linalg_LU_solve(cov_mat, permu, vect_b, pamer);

	std::cout << "Inv: " << std::endl;
	for (i = 0; i < terms; i++)
	{
		for (j = 0; j < terms; j++)
		{
			std::cout << gsl_matrix_get(mat_inv, i, j) << ", ";
			//std::cout << gsl_matrix_get(&mat.matrix, i, j) << "(" << cov_matrix[i*terms+j] << ")" << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for (i = 0; i < terms; i++)
	{
		coeffs[i] = gsl_vector_get(pamer, i);
	}

	gsl_matrix_free(cov_mat);
	gsl_vector_free(vect_b);
	gsl_matrix_free(mat_inv);
	gsl_vector_free(pamer);
	gsl_permutation_free(permu);

	delete[] cov_matrix;
	delete[] f_vector;
	delete[] check;
}


main(int argc, char* argv[])
{
	int i, j, k,m,n;
	int order = std::atoi(argv[1]);
	int img_num = std::atoi(argv[2]);
	std::cout << order << ", " << img_num << std::endl;
	int num = std::atoi(argv[3]);
	int terms = (order + 1)*(order + 2) / 2;
	double *x = new double[num] {};
	double *y = new double[num] {};
	double *fxy = new double[num] {};
	double *coeff = new double[terms] {};
	double *cov_cal = new double[terms*terms]{};
	double *f_vt_cal = new double[terms] {};
	char path[150], set1[20], set2[20], set3[20];
	char set_cov1[20], set_cov2[20];

	double temp = 0;
	double t1, t2;
	int s,s_b;
	double dx=0, scales=1, yscale = 1;
	double *covs = new double[terms*terms]{ };
	double *vbs = new double[terms]{ };


	/*for (k = 0; k < 1; k++)
	{
		gsl_matrix *mat = gsl_matrix_alloc(terms, terms);
		gsl_matrix* invs = gsl_matrix_alloc(terms, terms);
		gsl_vector *vb = gsl_vector_alloc(terms);
		gsl_vector *vp = gsl_vector_alloc(terms);
		gsl_permutation *permu = gsl_permutation_alloc(terms);

		for (i = 0; i < terms; i++)
		{
			gsl_vector_set(vb, i, vbs[i]);
			for (j = 0; j < terms; j++)
			{
				gsl_matrix_set(mat, i, j, covs[i*terms + j]);
			}
		}
		for (i = 0; i < terms; i++)
		{
			for (j = 0; j < terms; j++)
			{
				std::cout << gsl_matrix_get(mat, i, j) << ", ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		gsl_linalg_LU_decomp(mat, permu, &s);
		gsl_linalg_LU_invert(mat, permu, invs);
		gsl_linalg_LU_solve(mat, permu, vb, vp);
		for (i = 0; i < terms; i++)
		{
			for (j = 0; j < terms; j++)
			{
				std::cout << gsl_matrix_get(invs, i, j) << ", ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
		for (j = 0; j < terms; j++)
		{
			std::cout << gsl_vector_get(vp, j) << ", ";
		}
		std::cout << std::endl;
		std::cout << std::endl;

		gsl_matrix_free(mat);
		gsl_matrix_free(invs);
		gsl_vector_free(vb);
		gsl_vector_free(vp);
		gsl_permutation_free(permu);
	}*/

	for (i = 0; i < img_num; i++)
	{
		initialize_arr(x, num);
		initialize_arr(y, num);
		initialize_arr(fxy, num);
		initialize_arr(coeff, terms);		
	

		sprintf(path, "/home/hkli/work/cpp/test/data/%d.hdf5", i);
		sprintf(set1, "/x");
		sprintf(set2, "/y");
		sprintf(set3, "/fxy");
		sprintf(set_cov1, "/cov");
		sprintf(set_cov2, "/vect");
		read_h5(path, set1, x);
		read_h5(path, set2, y);
		read_h5(path, set3, fxy);
		read_h5(path, set_cov1, covs);
		read_h5(path, set_cov2, vbs);
		temp = 0;

		arr_rescale(x, dx, scales, num);
		arr_rescale(y, dx, scales, num);
		arr_rescale(fxy, 0, yscale, num);

		t1 = clock();		
		poly_fit_2d(x, y, fxy, num, order, coeff);
		t2 = clock();
		temp += t2 - t1;

		//cov_martix_2d_(x, y, fxy, num, order, cov_cal, f_vt_cal);
		// print the C-result
		for (j = 0; j < terms; j++)
		{
			for (k = 0; k < terms; k++)
			{
				std::cout << cov_cal[j*terms + k]<< ", ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		// print the difference between the python- and C- result
		for (j = 0; j< terms; j++)
		{
			for (k = 0; k < terms; k++)
			{
				std::cout << cov_cal[j*terms + k] - covs[j*terms + k] << ", ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		for (j = 0; j < terms; j++)
		{
			std::cout << coeff[j] << ",  ";
		}
		std::cout << std::endl;
		std::cout << "The estimated: " << std::endl;
		std::cout << fval_at_xy(0, 0, order, coeff) << ", " << fval_at_xy(55, 20, order, coeff) << std::endl;
		std::cout << fval_at_xy(dx*scales, dx*scales, order, coeff)/ yscale << ", " << fval_at_xy((dx+55)*scales, (dx+20)*scales, order, coeff)/ yscale << std::endl;

	}
	std::cout << temp / CLOCKS_PER_SEC;

	delete[] x;
	delete[] y;
	delete[] fxy;
	delete[] coeff;
	delete[] covs;
	delete[] vbs;
	delete[] cov_cal;
	delete[] f_vt_cal;
 	return 0;
}

/*terms = 3;
	gsl_matrix *mat = gsl_matrix_alloc(terms, terms);
	gsl_vector *v = gsl_vector_alloc(terms);
	gsl_permutation *permu = gsl_permutation_alloc(terms);
	gsl_vector * pamer = gsl_vector_alloc(terms);

	gsl_matrix *mat_balance = gsl_matrix_alloc(terms, terms);
	gsl_vector *v_balance = gsl_vector_alloc(terms);
	gsl_vector *vd = gsl_vector_alloc(terms);
	gsl_permutation *permu_balance = gsl_permutation_alloc(terms);
	gsl_vector * pamer_balance = gsl_vector_alloc(terms);
	
	std::cout << "Finish" << std::endl;

	for (i = 0; i < terms; i++)
	{	
		gsl_vector_set(v, i, a2[i]);
		gsl_vector_set(v_balance, i, a2[i]);

		gsl_vector_set(vd, i, 0);
		gsl_vector_set(pamer, i, 0);
		gsl_vector_set(pamer_balance, i, 0);
		
		for (j = 0; j < terms; j++)
		{
			gsl_matrix_set(mat, i, j, a1[i*terms + j]);
			gsl_matrix_set(mat_balance, i, j, a1[i *terms + j]);
		}
	}
	
	gsl_linalg_balance_matrix(mat_balance, vd);
	for (i = 0; i < terms; i++)
	{
		std::cout << gsl_vector_get(vd, i) << ", ";
		gsl_vector_set(v_balance, i, a2[i]/ gsl_vector_get(vd, i));
	}
	std::cout << std::endl;
	std::cout << std::endl;
	for (i = 0; i < terms; i++)
	{
		for (j = 0; j < terms; j++)
		{
			std::cout << gsl_matrix_get(mat_balance, i, j) << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	gsl_linalg_LU_decomp(mat_balance, permu_balance, &s_b);
	gsl_linalg_LU_solve(mat_balance, permu_balance, v_balance, pamer_balance);

	gsl_linalg_LU_decomp(mat, permu, &s);
	gsl_linalg_LU_solve(mat, permu, v, pamer);

	for (i = 0; i < terms; i++)
	{
		std::cout << gsl_vector_get(pamer, i) << ", " << gsl_vector_get(pamer_balance, i) << std::endl;
	}

	gsl_matrix_free(mat_balance);
	gsl_vector_free(v_balance);
	gsl_vector_free(vd);
	gsl_matrix_free(mat);
	gsl_vector_free(v);
	gsl_permutation_free(permu);
	gsl_vector_free(pamer);
	gsl_permutation_free(permu_balance);
	gsl_vector_free(pamer_balance);*/