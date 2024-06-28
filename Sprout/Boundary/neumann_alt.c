
#include "../defs.h"

void initial( double * , double * , double );

void mimic_cells( struct cell * cR , struct cell * cS , double dl , int theDIM ){

   memcpy( cR->prim   , cS->prim   , NUM_Q*sizeof(double) ); 
   memcpy( cR->cons   , cS->cons   , NUM_Q*sizeof(double) ); 
   memcpy( cR->RKcons , cS->RKcons , NUM_Q*sizeof(double) );
   memcpy( cR->xi     , cS->xi     ,     3*sizeof(double) );
   cR->xi[theDIM] += dl;

}

void outerboundary_x( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int i,j,k,i_end,j_end,k_end,ijk_send,ijk_recv;

   if( dim_rank[0] == 0 ){
      i_end = Ng;      if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ny+2*Ng; if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Nz+2*Ng; if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=0 ; k<k_end ; ++k ){
         for( j=0 ; j<j_end ; ++j ){
            for( i=0 ; i<i_end ; ++i ){
               ijk_send = Ng;
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * csend = theCells+ijk_send;
               struct cell * crecv = theCells+ijk_recv;
               mimic_cells( crecv , csend , (double)(i-Ng)*dx , 0 );
            }
         }
      }
   }
   
   if( dim_rank[0] == dim_size[0]-1 ){
      i_end = Nx+2*Ng; if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ny+2*Ng; if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Nz+2*Ng; if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=0 ; k<k_end ; ++k ){
         for( j=0 ; j<j_end ; ++j ){
            for( i=Nx+Ng ; i<i_end ; ++i ){
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_send = Nx+Ng-1;
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * crecv = theCells+ijk_recv;
               struct cell * csend = theCells+ijk_send;
               mimic_cells( crecv , csend , (double)(i-Nx-Ng+1)*dx , 0 );
            }
         }
      }
   }

}

void outerboundary_y( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dy = theDomain->dy;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int i,j,k,i_end,j_end,k_end,ijk_send,ijk_recv;

   if( dim_rank[1] == 0 ){
      i_end = Nx+2*Ng; if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ng;      if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Nz+2*Ng; if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=0 ; k<k_end ; ++k ){
         for( j=0 ; j<j_end ; ++j ){
            for( i=0 ; i<i_end ; ++i ){
               ijk_send = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*Ng;
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * csend = theCells+ijk_send;
               struct cell * crecv = theCells+ijk_recv;
               mimic_cells( crecv , csend , (double)(j-Ng)*dy , 1 );
            }
         }
      }
   }

   if( dim_rank[1] == dim_size[1]-1 ){
      i_end = Nx+2*Ng; if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ny+2*Ng; if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Nz+2*Ng; if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=0 ; k<k_end ; ++k ){
         for( j=Ny+Ng ; j<j_end ; ++j ){
            for( i=0 ; i<i_end ; ++i ){
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_send = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*(Ny+Ng-1);
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * crecv = theCells+ijk_recv;
               struct cell * csend = theCells+ijk_send;
               mimic_cells( crecv , csend , (double)(j-Ny-Ng+1)*dy , 1 );
            }
         }
      }
   }

}

void outerboundary_z( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dz = theDomain->dz;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int i,j,k,i_end,j_end,k_end,ijk_send,ijk_recv;

   if( dim_rank[1] == 0 ){
      i_end = Nx+2*Ng; if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ny+2*Ng; if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Ng;      if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=0 ; k<k_end ; ++k ){
         for( j=0 ; j<j_end ; ++j ){
            for( i=0 ; i<i_end ; ++i ){
               ijk_send = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*Ng;
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * csend = theCells+ijk_send;
               struct cell * crecv = theCells+ijk_recv;
               mimic_cells( crecv , csend , (double)(k-Ng)*dz , 2 );
            }
         }
      }
   }

   if( dim_rank[1] == dim_size[1]-1 ){
      i_end = Nx+2*Ng; if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ny+2*Ng; if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Nz+2*Ng; if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=Nz+Ng ; k<k_end ; ++k ){
         for( j=0 ; j<j_end ; ++j ){
            for( i=0 ; i<i_end ; ++i ){
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_send = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*(Nz+Ng-1);
               struct cell * crecv = theCells+ijk_recv;
               struct cell * csend = theCells+ijk_send;
               mimic_cells( crecv , csend , (double)(k-Nz-Ng+1)*dz , 2 );
            }
         }
      }
   }

}

void boundary( struct domain * theDomain ){

   int Nxtot = theDomain->theParList.Num_x;
   int Nytot = theDomain->theParList.Num_y;
   int Nztot = theDomain->theParList.Num_z;

   if( Nxtot > 1 )outerboundary_x( theDomain );
   if( Nytot > 1 )outerboundary_y( theDomain );
   if( Nztot > 1 )outerboundary_z( theDomain );

}
