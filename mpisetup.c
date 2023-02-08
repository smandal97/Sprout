#include "defs.h"

int mpiSetup( struct domain * theDomain , int argc , char * argv[] ){

   MPI_Comm_size(MPI_COMM_WORLD,&(theDomain->size));
   int size = theDomain->size;
   int error;
   int * dim_rank   = theDomain->dim_rank;
   int * dim_size   = theDomain->dim_size;
   int * left_rank  = theDomain->left_rank;
   int * right_rank = theDomain->right_rank;

   int NUM_X = theDomain->theParList.Num_x;
   int NUM_Y = theDomain->theParList.Num_y;
   int NUM_Z = theDomain->theParList.Num_z;

   //Dimensions and other arguments for making new MPI_COMM
   int wraparound[3] = {1,1,1};
   dim_size[0] = 0;
   dim_size[1] = 0;
   dim_size[2] = 0;
   dim_rank[0] = 0;
   dim_rank[1] = 0;
   dim_rank[2] = 0;
   int reorder = 1, ndims = 3;
 
   //Conditional for 1D and 2D domains
   if( NUM_X == 1 ){  ndims = 2; dim_size[0] = 1; };
   if( NUM_Y == 1 ){  ndims = 2; dim_size[1] = 1; };
   if( NUM_Z == 1 ){  ndims = 2; dim_size[2] = 1; };
  
   error = MPI_Dims_create(size,ndims,dim_size);
   if( error != MPI_SUCCESS ) return(1);
   
   MPI_Cart_create(MPI_COMM_WORLD,ndims,dim_size,wraparound,reorder,&(theDomain->theComm));
   MPI_Comm_rank(theDomain->theComm,&(theDomain->rank));
   int rank = theDomain->rank;
   MPI_Cart_coords(theDomain->theComm,rank,ndims,dim_rank);

   if(rank==0){printf("dim1 = %d, dim2 = %d, dim3 = %d\n",dim_size[0],dim_size[1],dim_size[2]);}
   
   //Convert topological ranks to my left and right in all 3 dimensions to process ranks
   int next_rank[3];
   int i;
   for( i=0 ; i<3 ; ++i ){
      next_rank[i] = dim_rank[i];
   }
   for( i=0 ; i<3 ; ++i ){
      next_rank[i] += 1;
      MPI_Cart_rank(theDomain->theComm,next_rank,&(right_rank[i]));
      next_rank[i] -= 2;
      MPI_Cart_rank(theDomain->theComm,next_rank,&(left_rank[i]));
      next_rank[i] += 1;
   }   
   return(0);
}
