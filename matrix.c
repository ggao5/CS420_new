#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "matrix.h"

#include <mpi.h>

////generate matrix for the size n//////
matrix_t* generate_matrix(int n){

  matrix_t* matrix = (matrix_t*)malloc(sizeof(matrix_t));
  int size, rank;
  int i, j, d;
  int ifirst, ilast;//THE ORIGINAL ifirst and ilast
  int i_interval;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //matrix->array=(double*)malloc(n*n*sizeof(double));
  matrix->Ifirst=(int*)malloc((2*n-5)*sizeof(int));
  matrix->Ilast=(int*)malloc((2*n-5)*sizeof(int));
  matrix->process = rank;
  //initialize the element of array in the boundary
  //get the Ifirst and Ilast
  for(d=0;d<n+n-5;d++){
      ifirst=(d<n-2)?d+1:n-2;
      ilast=(d<n-2)?1:d-n+4;
      i_interval=(ifirst-ilast+1)%size==0?max((int)((ifirst-ilast+1)/size),1):max((int)((ifirst-ilast+1)/size),1)+1;
      if (ifirst-ilast+1>=size*LENGTH_MPI)//decide whether to transfer
        {
        if((int)(ifirst-(rank+1)*i_interval+1)>=ilast){
          matrix->Ifirst[d]=(int)(ifirst-rank*i_interval);
          matrix->Ilast[d]=(int)(ifirst-(rank+1)*i_interval+1);
        }
        else if((int)(ifirst-(rank)*i_interval)>=ilast){
          matrix->Ifirst[d]=(int)(ifirst-rank*i_interval);
          matrix->Ilast[d]=ilast;
        }
        else
        {
           matrix->Ifirst[d]=ilast-1;
           matrix->Ilast[d]=ilast-1;
	   //printf("rank%d, d:%d, ifirst %d, ilast%d", rank, d, ilast-1,ilast-1);
           //matrix->Ilast[d]=-1;
           //matrix->Ifirst[d]=-1;
        }
        }
        else
        {
        if(rank==0){
         matrix->Ifirst[d]=ifirst;
         matrix->Ilast[d]=ilast;
        }
        else{
           matrix->Ifirst[d]=ilast-1;
           matrix->Ifirst[d]=ilast-1;
           //matrix->Ilast[d]=-1;
           //matrix->Ifirst[d]=-1;
        }
        }
	}
return matrix;
}

//destroy matrix free space
void destroy_matrix(matrix_t* matrix){
  free(matrix->Ifirst);
  free(matrix->Ilast);
  free(matrix);
}
