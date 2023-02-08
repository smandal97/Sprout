
#include "../defs.h"

/* These are useful for situations involving mirror symmetry about a corner, e.g.,
   an octant of a 3D explosion. The left hand boundaries incorporate a Neumann
   BC with zero derivative with the additional condition that the velocity
   component of the concerned dimension be mirrored in the ghost zones. The right
   hand side boundaries use a Dirichlet BC. */

void initial( double * , double * , double );
void prim2cons( double * , double * , double * , double );

void set_cells_left( struct cell * cR , struct cell * cS , double dl , double t , int theDIM ){

   memcpy( cR->prim , cS->prim , NUM_Q*sizeof(double) );
   cR->prim[UU1+theDIM] = -1.*cS->prim[UU1+theDIM];
   memcpy( cR->xi , cS->xi , 3*sizeof(double) );
   cR->xi[theDIM] += dl;
   //make pressure dirichlet
   /*double buff[NUM_Q];
   initial(buff,cR->xi,t);
   cR->prim[PPP] = buff[PPP];*/

}

void set_cells_right( struct cell * cR , struct cell * cS , double dl , double t , int theDIM ){

   memcpy( cR->xi , cS->xi , 3*sizeof(double) );
   cR->xi[theDIM] += dl;
   initial( cR->prim , cR->xi , t );
   //make density neumann
   /*double buff[NUM_Q];
   memcpy( buff , cS->prim , NUM_Q*sizeof(double) );
   cR->prim[RHO] = buff[RHO]; cR->prim[PPP] = buff[PPP];*/
}

void boundary1D_leftZN_rightD( struct domain * theDomain , int theDIM ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   double t  = theDomain->t;
   double dl;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int i,j,k,i_end,j_end,k_end,ijk_send,ijk_recv;

   int n[3]   = {0};
   n[theDIM]  = 1;
   int cn[3]; cn[0] = 1; cn[1] = 1; cn[2] = 1;
   cn[theDIM] = 0;

   if( dim_rank[theDIM] == 0 ){
      i_end = Ng+cn[0]*(Nx+Ng); if( theDomain->theParList.Num_x == 1 ) i_end = 1;
      j_end = Ng+cn[1]*(Ny+Ng); if( theDomain->theParList.Num_y == 1 ) j_end = 1;
      k_end = Ng+cn[2]*(Nz+Ng); if( theDomain->theParList.Num_z == 1 ) k_end = 1;

      for( k=0 ; k<k_end ; ++k ){
         for( j=0 ; j<j_end ; ++j ){
            for( i=0 ; i<i_end ; ++i ){
               ijk_send = cn[0]*i + n[0]*(Ng+1-i);
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ){
                  ijk_send += (Nx+2*Ng)*( cn[1]*j + n[1]*(Ng+1-j) );
                  ijk_recv += (Nx+2*Ng)*j;
               }
               if( theDomain->theParList.Num_z != 1 ){
                  ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*( cn[2]*k + n[2]*(Ng+1-k) );
                  ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               }
               struct cell * csend = theCells+ijk_send;
               struct cell * crecv = theCells+ijk_recv;
               dl = (double)(n[0]*(2*i-Ng-1))*dx + (double)(n[1]*(2*j-Ng-1))*dy + (double)(n[2]*(2*k-Ng-1))*dz;
               set_cells_left( crecv , csend , dl , t , theDIM );
               prim2cons( crecv->prim , crecv->cons , crecv->xi , dx*dy*dz );
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
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_send = cn[0]*i + n[0]*(2*Nx+2*Ng-1-i);
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*(cn[1]*j + n[1]*(2*Ny+2*Ng-1-j));
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*(cn[2]*k + n[2]*(2*Nz+2*Ng-1-k));
               struct cell * crecv = theCells+ijk_recv;
               struct cell * csend = theCells+ijk_send;
               dl = (double)(n[0]*(2*(i-Nx-Ng)+1))*dx + (double)(n[1]*(2*(j-Ny-Ng)+1))*dy + (double)(n[2]*(2*(k-Nz-Ng)+1))*dz;
               set_cells_right( crecv , csend , dl , t , theDIM );
               prim2cons( crecv->prim , crecv->cons , crecv->xi , dx*dy*dz );
            }
         }
      }
   }

}

void boundary( struct domain * theDomain ){

   int Nxtot = theDomain->theParList.Num_x;
   int Nytot = theDomain->theParList.Num_y;
   int Nztot = theDomain->theParList.Num_z;

   if( Nxtot > 1 )boundary1D_leftZN_rightD( theDomain , 0 );
   if( Nytot > 1 )boundary1D_leftZN_rightD( theDomain , 1 );
   if( Nztot > 1 )boundary1D_leftZN_rightD( theDomain , 2 );

}