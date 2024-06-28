
#include "../defs.h"

void initial( double * , double * , double );
void prim2cons( double * , double * , double * , double );

void set_coords( struct cell * cR , struct cell * cS , double dl , int theDIM ){
   memcpy( cR->xi , cS->xi , 3*sizeof(double) );
   cR->xi[theDIM] += dl;
}


void boundary1D_periodic( struct domain * theDomain , int theDIM ){

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
               ijk_send = cn[0]*i + n[0]*Ng;
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*(cn[1]*j + n[1]*Ng);
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*(cn[2]*k + n[2]*Ng);
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * csend = theCells+ijk_send;
               struct cell * crecv = theCells+ijk_recv;
               dl = (double)(n[0]*(i-Ng))*dx + (double)(n[1]*(j-Ng))*dy + (double)(n[2]*(k-Ng))*dz;
               set_coords( crecv , csend , dl , theDIM );
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
               ijk_send = cn[0]*i + n[0]*(Nx+Ng-1);
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*(cn[1]*j + n[1]*(Ny+Ng-1));
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*(cn[2]*k + n[2]*(Nz+Ng-1));
               struct cell * crecv = theCells+ijk_recv;
               struct cell * csend = theCells+ijk_send;
               dl = (double)(n[0]*(i-Nx-Ng+1))*dx + (double)(n[1]*(j-Ny-Ng+1))*dy + (double)(n[2]*(k-Nz-Ng+1))*dz;
               set_coords( crecv , csend , dl , theDIM );
               prim2cons( crecv->prim , crecv->cons , crecv->xi , dx*dy*dz );
            }
         }
      }
   }

}


void boundary1D_dirichlet( struct domain * theDomain , int theDIM ){

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
               ijk_send = cn[0]*i + n[0]*Ng;
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*(cn[1]*j + n[1]*Ng);
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*(cn[2]*k + n[2]*Ng);
               ijk_recv = i;
               if( theDomain->theParList.Num_y != 1 ) ijk_recv += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk_recv += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * csend = theCells+ijk_send;
               struct cell * crecv = theCells+ijk_recv;
               dl = (double)(n[0]*(i-Ng))*dx + (double)(n[1]*(j-Ng))*dy + (double)(n[2]*(k-Ng))*dz;
               set_coords( crecv , csend , dl , theDIM );
               initial( crecv->prim , crecv->xi , t );
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
               ijk_send = cn[0]*i + n[0]*(Nx+Ng-1);
               if( theDomain->theParList.Num_y != 1 ) ijk_send += (Nx+2*Ng)*(cn[1]*j + n[1]*(Ny+Ng-1));
               if( theDomain->theParList.Num_z != 1 ) ijk_send += (Nx+2*Ng)*(Ny+2*Ng)*(cn[2]*k + n[2]*(Nz+Ng-1));
               struct cell * crecv = theCells+ijk_recv;
               struct cell * csend = theCells+ijk_send;
               dl = (double)(n[0]*(i-Nx-Ng+1))*dx + (double)(n[1]*(j-Ny-Ng+1))*dy + (double)(n[2]*(k-Nz-Ng+1))*dz;
               set_coords( crecv , csend , dl , theDIM );
               initial( crecv->prim , crecv->xi , t );
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

   if( Nxtot > 1 )boundary1D_periodic( theDomain , 0 );
   if( Nytot > 1 )boundary1D_dirichlet( theDomain , 1 );
   if( Nztot > 1 )boundary1D_periodic( theDomain , 2 );

}
