extern double array[N][N];
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"


//just for two process

double send_recv2p(matrix_t* matrix,int* LABELS, int d){
  double exe_time = MPI_Wtime();
  int i, j;
  int ifirst, ilast;
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(d%2==0){
  if(rank==0){
   i = matrix->Ilast[d];
   j = d+2-i;
   MPI_Send(&array[i][j],1,MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
  }
  else{
   i = matrix->Ifirst[d]+1;
   j = d+2-i;
   MPI_Recv(&array[i][j],1,MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  }
  else if(d%2){
  if(rank==1){
   i = matrix->Ifirst[d];
   j = d+2-i;
   MPI_Send(&array[i][j],1,MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  else{
   i = matrix->Ilast[d]-1;
   j = d+2-i;
   MPI_Recv(&array[i][j],1,MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  }
  exe_time = MPI_Wtime()-exe_time;
  return exe_time;
}

void compute_labels(int* labels, matrix_t* matrix){
  int i, j, d;
  int ifirst, ilast;
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for(i=0;i<N;i++)
	  for(j=0;j<N;j++){
	  labels[i*N+j]=0;
	  }
  for(d=0;d<N+N-5;d++){
      ifirst=matrix->Ifirst[d];
      ilast=matrix->Ilast[d];
      if(ifirst+1!=0){
      for( i=ifirst; i>=ilast;i-- ){
      j=d+2-i;
      //labels[i*N+j]=labels[i*N+j]+rank;
      labels[i*N+j]=rank;
      //printf("OK i:%d, j%d,%d\n",i, j, labels[i*N+j]);
      }
      }
  }
}
void calculate(matrix_t* matrix, int ifirst, int ilast, int d){
//the function to calculate element in the matrix
  int i;
  int IFIRST, ILAST;
  IFIRST=(d<N-2)?d+1:N-2;
  ILAST=(d<N-2)?1:d-N+4;
  if(ifirst>=ILAST) {
  if (ifirst-ilast>LENGTH_OMP){//decide whether or not use openmpi
  #pragma omp parallel for
  for (i = ifirst; i >= ilast; i--)
  {
      int j=d+2-i;
      array[i][j]=0.25*(array[i-1][j]+array[i+1][j]+array[i][j-1]+array[i][j+1]);
  }
  }
  else
  for (i =ifirst; i >= ilast; i--)
  {
      int j=d+2-i;
      array[i][j]=0.25*(array[i-1][j]+array[i+1][j]+array[i][j-1]+array[i][j+1]);
  }
}
}

double transfer_d_row(matrix_t* matrix, int d){
//not use MPI VECTOR, just transfer the elements one by one to the first
//process
  MPI_Barrier(MPI_COMM_WORLD);//as it should be sequantial in the d index for both MPI and openmp
  double exe_time = MPI_Wtime();
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int i, j;
  int k;
  int ifirst, ilast;
  //get the datatype
  MPI_Datatype New_Type;
  int IFIRST=(d<N-2)?d+1:N-2;//THE i index for the entire matrix
  int ILAST=(d<N-2)?1:d-N+4;
  MPI_Type_vector(IFIRST-ILAST+1, 1, N-1, MPI_DOUBLE, &New_Type);
  MPI_Type_commit(&New_Type);
  
   if(rank ==0 ){
     for( k=1; k<size;k++){
       MPI_Recv(&(ifirst),1,MPI_INT, k,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       MPI_Recv(&(ilast),1,MPI_INT, k,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       if(ifirst+1!=0)
       for(i=ifirst;i>=ilast;i--) {
       j=d+2-i;
       MPI_Recv(&(array[i][j]),1,MPI_DOUBLE, k , 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       }
     }
   }
   else{
     ifirst= matrix->Ifirst[d];
     ilast= matrix->Ilast[d];
     MPI_Send(&(ifirst),1,MPI_INT, 0,0,MPI_COMM_WORLD);
     MPI_Send(&(ilast),1,MPI_INT, 0,0,MPI_COMM_WORLD);
     if(ifirst+1!=0)
     for(i=ifirst;i>=ilast;i--) {
     j=d+2-i;
     MPI_Send(&(array[i][j]),1,MPI_DOUBLE, 0 ,0, MPI_COMM_WORLD);
     }
   }
  //MPI_Bcast(&array, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&array[ILAST][d+2-ILAST],1 , New_Type, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  exe_time = MPI_Wtime()-exe_time;
  return exe_time;
}

void par_mpi_gauss_seidel(int argc, char** argv)
// the function to use MPI for Gauss-Seidel iteration
{
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int i, j, d;
  int iter;
  int ifirst, ilast;
  int IFIRST, ILAST;
  double total_time= MPI_Wtime();
  double inner_time=0;
  int* labels=(int*)malloc(sizeof(int)*N*N);
  int* LABELS=(int*)malloc(sizeof(int)*N*N);
  matrix_t* matrix = generate_matrix(N);

//  if(rank==0){
//  for(d=0;d<2*N-5;d++)
//  printf("d:%d, Ifirst:%d, Ilast:%d\n",d,matrix->Ifirst[d],matrix->Ilast[d]);
//  printf("\n");
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//  if(rank==1){
//  for(d=0;d<2*N-5;d++)
//  printf("d:%d, Ifirst:%d, Ilast:%d\n",d,matrix->Ifirst[d],matrix->Ilast[d]);
//  printf("\n");
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//  if(rank==2){
//  for(d=0;d<2*N-5;d++)
//  printf("d:%d, Ifirst:%d, Ilast:%d\n",d,matrix->Ifirst[d],matrix->Ilast[d]);
//  printf("\n");
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//  if(rank==3){
//  for(d=0;d<2*N-5;d++)
//  printf("d:%d, Ifirst:%d, Ilast:%d\n",d,matrix->Ifirst[d],matrix->Ilast[d]);
//  printf("\n");
//  }
  //initialize the boundary value
  for(j=0;j<N;j++)
    array[0][j]=BOTTOM_BOUNDARY_VALUE;
  for(j=0;j<N;j++)
    array[N-1][j]=TOP_BOUNDARY_VALUE;
  for(i=1;i<N-1;i++)
    array[i][0]=LEFT_BOUNDARY_VALUE;
  for(i=1;i<N-1;i++)
    array[i][N-1]=RIGHT_BOUNDARY_VALUE;
  
  for (iter = 0; iter < ITERS ; iter++)
  {
    for(d=0;d<N+N-5;d++){
      ifirst = matrix->Ifirst[d];
      ilast = matrix->Ilast[d];
      calculate(matrix, ifirst, ilast,d);
      if(size!=0)
      inner_time+=send_recv2p(matrix,LABELS, d);//update the whole matrix once the element has changed
      }
  }
  

  if(size!=0)
  for(d=0;d<N+N-5;d++){
  inner_time+=transfer_d_row(matrix, d);//update the whole matrix once the element has changed
  }

  MPI_Barrier(MPI_COMM_WORLD);
  total_time = MPI_Wtime()-total_time;

  destroy_matrix(matrix);
  printf("execution time for hybrid:%f\ntime for transfer in the loop: %f\n", total_time,inner_time);
}
