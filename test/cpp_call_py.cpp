#include<Python.h>
#include<numpy/arrayobject.h>
#include<iostream>

main(int argc, char *argv[])
{
	// this memory block, allocated by Py_DecodeLocale(), must be freed by PyMem_RawFree();
	wchar_t *program = Py_DecodeLocale(argv[0], NULL);
	if (program == NULL) {
		fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
		exit(1);
	}
	char home_path[50];
	sprintf(home_path, "/home/hkli/anaconda3");
	wchar_t *home = Py_DecodeLocale(home_path, NULL);
	Py_SetProgramName(program);
	Py_SetPythonHome(home);
	Py_Initialize();
	import_array();
	double *array_1 = new double[6]{ 2,5,6,5,6,5 };
	npy_intp dims_1[] = { 2,3 };
	PyObject *mat_1 = PyArray_SimpleNewFromData(2, dims_1, NPY_DOUBLE, array_1);
	
	double array_2[3][4] = { { 1,3,0,4 },{ 2,2,5,3 },{ 1,2,1,4 } };
	npy_intp dims_2[] = { 3,4 };
	PyObject *mat_2 = PyArray_SimpleNewFromData(2, dims_2, NPY_DOUBLE, array_2);

	PyObject *prod = PyArray_MatrixProduct(mat_1, mat_2);

	PyArrayObject *mat_3;
	PyArray_OutputConverter(prod, &mat_3);
	npy_intp *shape = PyArray_SHAPE(mat_3);
	double *array_3 = (double*)PyArray_DATA(mat_3);

	std::cout << "numpy result:\n";
	for (int i = 0; i < shape[0]; i++)
	{
		for (int j = 0; j < shape[1]; j++)
		{
			std::cout << array_3[i*shape[1] + j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << "\nC result:\n";
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			double t = 0;
			for (int k = 0; k < 3; k++)
				t += array_1[i*3+k] * array_2[k][j];
			std::cout << t << "\t";
		}
		std::cout << std::endl;
	}

	if (Py_FinalizeEx() < 0)
	{
		exit(120);
	}
	PyMem_RawFree(home);
	PyMem_RawFree(program);
	delete[] array_1;
	getchar();
	return 0;
}
