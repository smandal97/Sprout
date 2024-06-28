
#include "../defs.h"

void prim2cons( double * , double * , double * , double );

void set_ghostcells( struct cell * cR , struct cell * cS0 , struct cell * cS1 , double dl , double dV , int theDIM ){

   memcpy( cR->xi , cS0->xi , 3*sizeof(double) );
   cR->xi[theDIM] += dl;
   int var;
   for( var=0 ; var<NUM_Q ; ++var ){
      cR->prim[var] = 2.*cS0->prim[var] - cS1->prim[var];
   }
   prim2cons( cR->prim , cR->cons , cR->xi , dV );

}

void boundary1D_neumann_local( struct domain * theDomain , int theDIM ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   double dl;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int i,j,k,i_end,j_end,k_end,ijk_sen0,ijk_sen1,ijk_recv;

   int n[3]   = {0};
   n[theDIM]  = 1;
   int cn[3]; cn[0] = 1; cn[1] = 1; cn[2] = 1;
   cn[theDIM] = 0;

   if( dim_rank[theDIM] == 0 ){
      i_end = Ng-1+cn[0]*(Nx+Ng); if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ng-1+cn[1]*(Ny+Ng); if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Ng-1+cn[2]*(Nz+Ng); if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=k_end ; k>=0 ; --k ){
         for( j=j_end ; j>=0 ; --j ){
            for( i=i_end ; i>=0 ; --i ){
               ijk_sen0 = i + n[0]*1;
               ijk_sen1 = i + n[0]*2;
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ){
                  ijk_sen0 += (Nx+2*Ng)*( j + n[1]*1 );
                  ijk_sen1 += (Nx+2*Ng)*( j + n[1]*2 );
                  ijk_recv += (Nx+2*Ng)*j;
               }
               if( theDomain->theParList.Num_z != 1 ){
                  ijk_sen0 += (Nx+2*Ng)*(Ny+2*Ng)*( k + n[2]*1 );
                  ijk_sen1 += (Nx+2*Ng)*(Ny+2*Ng)*( k + n[2]*2 );
                  ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               }
               struct cell * csen0 = theCells+ijk_sen0;
               struct cell * csen1 = theCells+ijk_sen1;
               struct cell * crecv = theCells+ijk_recv;
               dl = -(double)(n[0])*dx - (double)(n[1])*dy - (double)(n[2])*dz;
               set_ghostcells( crecv , csen0 , csen1 , dl , dx*dy*dz , theDIM );
            }
         }
      }
   }
   
   if( dim_rank[theDIM] == dim_size[theDIM]-1 ){
      i_end = Nx+2*Ng; if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ny+2*Ng; if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Nz+2*Ng; if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=n[2]*(Nz+Ng) ; k<k_end ; ++k ){
         for( j=n[1]*(Ny+Ng) ; j<j_end ; ++j ){
            for( i=n[0]*(Nx+Ng) ; i<i_end ; ++i ){
               ijk_sen0 = i - n[0]*1;
               ijk_sen1 = i - n[0]*2;
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ){
                  ijk_sen0 += (Nx+2*Ng)*( j - n[1]*1 );
                  ijk_sen1 += (Nx+2*Ng)*( j - n[1]*2 );
                  ijk_recv += (Nx+2*Ng)*j;
               }
               if( theDomain->theParList.Num_z != 1 ){
                  ijk_sen0 += (Nx+2*Ng)*(Ny+2*Ng)*( k - n[2]*1 );
                  ijk_sen1 += (Nx+2*Ng)*(Ny+2*Ng)*( k - n[2]*2 );
                  ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               }
               struct cell * csen0 = theCells+ijk_sen0;
               struct cell * csen1 = theCells+ijk_sen1;
               struct cell * crecv = theCells+ijk_recv;
               dl = (double)(n[0])*dx + (double)(n[1])*dy + (double)(n[2])*dz;
               set_ghostcells( crecv , csen0 , csen1 , dl , dx*dy*dz , theDIM );
            }
         }
      }
   }

}

void boundary( struct domain * theDomain ){

   int Nxtot = theDomain->theParList.Num_x;
   int Nytot = theDomain->theParList.Num_y;
   int Nztot = theDomain->theParList.Num_z;

   if( Nxtot > 1 )boundary1D_neumann_local( theDomain , 0 );
   if( Nytot > 1 )boundary1D_neumann_local( theDomain , 1 );
   if( Nztot > 1 )boundary1D_neumann_local( theDomain , 2 );

}
