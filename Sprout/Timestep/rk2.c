
#include "../defs.h"

void onestep( struct domain * , double , double , int , int );
void exchangeData( struct domain * );
void boundary( struct domain * );


void timestep( struct domain * theDomain , double dt ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   
   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;
   
   int i,j,k,ijk;
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            ijk  = i;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
            struct cell * c = theCells+ijk;
            memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
            memcpy( c->RKxi   , c->xi   ,     3*sizeof(double) );
         }
      }
   }

   onestep( theDomain , 0.0 ,     dt , 1 , 0 );
   theDomain->t += dt;
   exchangeData( theDomain );
   boundary( theDomain );

   onestep( theDomain , 0.5 , 0.5*dt , 0 , 1 );
   exchangeData( theDomain );
   boundary( theDomain );

   theDomain->count_steps++;



}
