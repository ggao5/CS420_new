#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include "par_Gauss_Seidel.h"
#include "matrix.h"

#define CACHE_SIZE 12582912

double array[N][N];

// Used for testing.
double correct_results[N][N];

void basic_gauss_seidel();
void par_omp_gauss_seidel();
// Initialize array to zeros.
void reset_array()
{
  memset(array, 0, sizeof(array));
}

void unit_test()
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      if (fabs(correct_results[i][j] - array[i][j]) > ERROR_THRESHOLD)
      {
        printf("Incorrect result.\n");
        return;
      }
    }
  }
  printf("GOOD -- Test passed.\n");
}

void clear_cache() {
  size_t array_size = CACHE_SIZE * 2.5;
  double temp = 0.0;
  double* dummy = (double*)malloc(sizeof(double) * array_size);

  for (size_t i = 0; i < array_size; ++i) {
    dummy[i] = i / array_size;
  }

  for (size_t j = 0; j < 10; ++j) {
    for (size_t i = 0; i < array_size; ++i) {
      temp += dummy[i] / j;
    }
  }

  fprintf(stderr, "Cleared cache (ignore this value: %f)\n", temp);
}	
int main (int argc, char** argv)
{

  MPI_Init(&argc, &argv);

  int size, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  reset_array();

  printf("Running basic sequential version to save the results.\n");


  double exe_time = MPI_Wtime();
  {
    basic_gauss_seidel();
  }
  exe_time = MPI_Wtime()-exe_time;
  printf("execution time for basic version:%f\n", exe_time);

  memcpy(correct_results, array, sizeof(array));

  reset_array();
  #pragma omp parallel
  {
    omp_set_num_threads(NUM_THREADS);
    if (omp_get_thread_num() == 0)
      printf("Running openmp only with %d threads.\n", omp_get_num_threads());
  }
  clear_cache();
  exe_time = MPI_Wtime();
  par_omp_gauss_seidel();
  exe_time = MPI_Wtime()-exe_time;
  printf("execution time for openmp only version:%f\n", exe_time);
  unit_test();


  reset_array();
  clear_cache();
  MPI_Barrier(MPI_COMM_WORLD);//as it should be sequantial in the d index for both MPI and openmp
  printf("Running openmp with mpi\n");
  par_mpi_gauss_seidel();

  unit_test();
  MPI_Finalize();
}
