
#include "defs.h"

void report( struct domain * theDomain , double t , double dt ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   double W  = theDomain->W;
   MPI_Comm grid_comm = theDomain->theComm;

   int i,j,k,ijk;
   double ESum = 0.0;
   double MSum = 0.0;

   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Ny ; ++j ){
         for( i=0 ; i<Nx ; ++i ){
            ijk  = i+Ng;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*(j+Ng);
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*(k+Ng);
            struct cell * c = theCells+ijk;
            double E = c->cons[TAU];
            double M = c->cons[DEN];
            ESum += E;
            MSum += M;
         }
      }
   }

   MPI_Allreduce( MPI_IN_PLACE , &ESum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &MSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );

   if( rank==0 ){
      FILE * rFile = fopen("report.dat","a");
      fprintf(rFile,"%e %e %e %e %e %e %e %e\n", t, dt, W , dx, dy, dz, ESum, MSum);
      fclose(rFile);
   }

}
