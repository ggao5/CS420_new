#pragma once 
#define min(a,b) ((a)<(b) ? a:b)
#define max(a,b) ((a)>(b) ? a:b)

typedef struct{
  int* Ifirst;//array of the ifirst, which have 2n-5 element
  int* Ilast;//array of the ilast, which have 2n-5 element
  int n;//total size of the array
  int process;
} matrix_t;

matrix_t* generate_matrix(int n);

void destroy_matrix(matrix_t* matrix);
