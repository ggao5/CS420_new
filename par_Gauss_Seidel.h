#include "matrix.h"
void par_omp_gauss_seidel();
double  transfer_boundary(matrix_t* matrix, int d);
void calculate(matrix_t* matrix, int ifirst, int ilast);
double combine_matrix(matrix_t* matrix);
void par_mpi_gauss_seidel();
