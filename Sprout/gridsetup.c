
#include "defs.h"

int getN0( int drank , int dsize , int dnum ){
   int N0 = (dnum*drank)/dsize;
   return(N0);
}

void setupGrid( struct domain * theDomain ){

   int Ng = NUM_G;
   theDomain->Ng = Ng;
   int rank = theDomain->rank;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   
   int Nx = theDomain->theParList.Num_x;
   int Ny = theDomain->theParList.Num_y;
   int Nz = theDomain->theParList.Num_z;

   double Lx = theDomain->theParList.Lx;
   double Ly = theDomain->theParList.Ly;
   double Lz = theDomain->theParList.Lz;

   double dx = 0.0;
   double dy = 0.0;
   double dz = 0.0;

   if( rank==0 ){
      dx = Lx/Nx;
      dy = Ly/Ny;
      dz = Lz/Nz;
   }

   MPI_Allreduce( MPI_IN_PLACE , &dx , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , &dy , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , &dz , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );
   MPI_Barrier( theDomain->theComm );

   theDomain->dx = dx;
   theDomain->dy = dy;
   theDomain->dz = dz;

   int N0x = getN0( dim_rank[0]   , dim_size[0] , Nx );
   int N0y = getN0( dim_rank[1]   , dim_size[1] , Ny );
   int N0z = getN0( dim_rank[2]   , dim_size[2] , Nz );

   int N1x = getN0( dim_rank[0]+1 , dim_size[0] , Nx );
   int N1y = getN0( dim_rank[1]+1 , dim_size[1] , Ny );
   int N1z = getN0( dim_rank[2]+1 , dim_size[2] , Nz );

   Nx = N1x - N0x;
   Ny = N1y - N0y;
   Nz = N1z - N0z;

   theDomain->Nx = Nx;
   theDomain->Ny = Ny;
   theDomain->Nz = Nz;
 
   int domain_size = 1;
   if( theDomain->theParList.Num_x != 1 ) domain_size *= Nx+2*Ng;
   if( theDomain->theParList.Num_y != 1 ) domain_size *= Ny+2*Ng;
   if( theDomain->theParList.Num_z != 1 ) domain_size *= Nz+2*Ng;

   theDomain->theCells = (struct cell *) malloc( domain_size*sizeof(struct cell) );
   printf("dim_rank = [%i %i %i], Nx = %i, Ny = %i, Nz = %i\n", dim_rank[0],dim_rank[1],dim_rank[2],Nx,Ny,Nz);
   int i,j,k,ijk;
   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Ny ; ++j ){
         for( i=0 ; i<Nx ; ++i ){
            ijk  = i+Ng;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*(j+Ng);
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*(k+Ng);
            struct cell * c = theDomain->theCells+ijk;
            //set indices
            c->xi[0] = ((double)(i+N0x)+0.5)*dx;
            c->xi[1] = ((double)(j+N0y)+0.5)*dy;
            c->xi[2] = ((double)(k+N0z)+0.5)*dz;
         }
      }
   }
   


}
