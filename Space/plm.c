
#include "../defs.h"

double minmod( double x , double y , double z ){
   double a = x;
   if( x*y < 0.0 ) a = 0.0;
   if( fabs(y) < fabs(a) ) a = y;
   if( y*z < 0.0 ) a = 0.0;
   if( fabs(z) < fabs(a) ) a = z;
   return(a);
}

void space_recon1D( struct domain * theDomain , int theDIM ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   double PLM = theDomain->theParList.PLM;

   int n[3]  = {0};
   n[theDIM] = 1;
   double pL,pC,pR,sL,sR,sC;
   double dl = (double)n[0]*dx + (double)n[1]*dy + (double)n[2]*dz;

   int i_end = Nx+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_end = 1;
   int j_end = Ny+2*Ng;
   if( theDomain->theParList.Num_y == 1 ) j_end = 1;
   int k_end = Nz+2*Ng;
   if( theDomain->theParList.Num_z == 1 ) k_end = 1;
   int i,j,k,q,ijk,ijk_l,ijk_r;

   for( k=0 ; k<k_end ; ++k ){
      for( j=0 ; j<j_end ; ++j ){
         for( i=0 ; i<i_end ; ++i ){
            ijk   = i;
            if( i==0 ) ijk_l = i;
            else ijk_l = i-n[0];
            if( i==i_end-1 ) ijk_r = i;
            else ijk_r = i+n[0];
            if( theDomain->theParList.Num_y != 1 ){
               ijk   += (Nx+2*Ng)*j;
               if( j==0 ) ijk_l += (Nx+2*Ng)*j;
               else ijk_l += (Nx+2*Ng)*(j-n[1]);
               if( j==j_end-1 ) ijk_r += (Nx+2*Ng)*j;
               else ijk_r += (Nx+2*Ng)*(j+n[1]);                     
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk   += (Nx+2*Ng)*(Ny+2*Ng)*k;
               if( k==0 ) ijk_l += (Nx+2*Ng)*(Ny+2*Ng)*k;
               else ijk_l += (Nx+2*Ng)*(Ny+2*Ng)*(k-n[2]);
               if( k==k_end-1 ) ijk_r += (Nx+2*Ng)*(Ny+2*Ng)*k;
               else ijk_r += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]);                             
            }
            struct cell * c  = theCells+ijk;
            struct cell * cL = theCells+ijk_l;
            struct cell * cR = theCells+ijk_r;
            for( q=0 ; q<NUM_Q ; ++q ){
               pL = cL->prim[q];
               pC =  c->prim[q];
               pR = cR->prim[q];
               sL = (pC - pL)/dl;
               sR = (pR - pC)/dl;
               sC = (pR - pL)/dl;
               if( theDIM==0 ) c->gradx[q] = minmod( PLM*sL , sC , PLM*sR );
               if( theDIM==1 ) c->grady[q] = minmod( PLM*sL , sC , PLM*sR );
               if( theDIM==2 ) c->gradz[q] = minmod( PLM*sL , sC , PLM*sR );
           }
         }
      }
   }

}
