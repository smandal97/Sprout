enum{MH,PH};

#include "../defs.h"

#define eps 1e-40

int set_accuracy(void){
   return 3;
}


double set_face_values( double um1 , double u0 , double up1 , double dl ){
   double dl2,b0,b1,tau,a0,a1,w0,w1,p0,p1;

   dl2 = dl*dl;

   b0 = up1-u0;
   b1 = um1-u0;
   b0 = b0*b0;
   b1 = b1*b1;

   tau = up1 - 2.*u0 + um1;
   tau = tau*tau;

   a0 = 2.*( 1.+tau/(b0+dl2) );
   a1 = 1.*( 1.+tau/(b1+dl2) );

   w0 = a0/(a0+a1);
   w1 = a1/(a0+a1);

   p0 = (u0+up1)/2.;
   p1 = (3.*u0-um1)/2.;

   return(w0*p0+w1*p1);

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


   int n[3]  = {0};
   n[theDIM] = 1;
   double pm2,pm1,p00,pp1,pmh,pph;
   double dl = (double)n[0]*dx + (double)n[1]*dy + (double)n[2]*dz;

   int i_end = Nx+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_end = 1;
   int j_end = Ny+2*Ng;
   if( theDomain->theParList.Num_y == 1 ) j_end = 1;
   int k_end = Nz+2*Ng;
   if( theDomain->theParList.Num_z == 1 ) k_end = 1;
   int i,j,k,q,ijkm2,ijkm1,ijk0,ijkp1;

   for( k=0 ; k<k_end ; ++k ){
      for( j=0 ; j<j_end ; ++j ){
         for( i=0 ; i<i_end ; ++i ){
            ijkm2 = i-2*n[0];
            ijkm1 = i-1*n[0];
            ijk0  = i;
            ijkp1 = i+1*n[0];
            if( i==0 ){
               ijkm2 = i;
               ijkm1 = i;
            }else if( i==1 ){
               ijkm2 = i;
            }else if( i==i_end-1 ){
               ijkp1 = i_end-1;
            }
            if( theDomain->theParList.Num_y != 1 ){
               ijk0 += (Nx+2*Ng)*j;
               if( j==0 || j==1 ){
                  ijkp1 += (Nx+2*Ng)*(j+1*n[1]);
               }else if( j==j_end-1 ){
                  ijkm2 += (Nx+2*Ng)*(j-2*n[1]);
                  ijkm1 += (Nx+2*Ng)*(j-1*n[1]);
                  ijkp1 += (Nx+2*Ng)*(j_end-1);
               }else{
                  ijkm2 += (Nx+2*Ng)*(j-2*n[1]);
                  ijkm1 += (Nx+2*Ng)*(j-1*n[1]);
                  ijkp1 += (Nx+2*Ng)*(j+1*n[1]);
               }
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk0 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               if( k==0 || k==1 ){
                  ijkp1 += (Nx+2*Ng)*(Ny+2*Ng)*(k+1*n[2]);
               }else if( k==k_end-1 ){
                  ijkm2 += (Nx+2*Ng)*(Ny+2*Ng)*(k-2*n[2]);
                  ijkm1 += (Nx+2*Ng)*(Ny+2*Ng)*(k-1*n[2]);
                  ijkp1 += (Nx+2*Ng)*(Ny+2*Ng)*(k_end-1);
               }else{
                  ijkm2 += (Nx+2*Ng)*(Ny+2*Ng)*(k-2*n[2]);
                  ijkm1 += (Nx+2*Ng)*(Ny+2*Ng)*(k-1*n[2]);
                  ijkp1 += (Nx+2*Ng)*(Ny+2*Ng)*(k+1*n[2]);
               }
            }

            struct cell * cm2 = theCells+ijkm2;
            struct cell * cm1 = theCells+ijkm1;
            struct cell * c00 = theCells+ijk0;
            struct cell * cp1 = theCells+ijkp1;

            for( q=0 ; q<NUM_Q ; ++q ){
               pm2 = cm2->prim[q];
               pm1 = cm1->prim[q];
               p00 = c00->prim[q];
               pp1 = cp1->prim[q];

               pmh = set_face_values(pm2,pm1,p00,dl);
               pph = set_face_values(pm1,p00,pp1,dl);

               if( theDIM==0 ){
                  c00->gradx[q] = (pph-pmh)/dl;
                  c00->pblax[q] = (pph+pmh)/2.-p00;
               }
               if( theDIM==1 ){
                  c00->grady[q] = (pph-pmh)/dl;
                  c00->pblay[q] = (pph+pmh)/2.-p00;
               }
               if( theDIM==2 ){
                  c00->gradz[q] = (pph-pmh)/dl;
                  c00->pblaz[q] = (pph+pmh)/2.-p00;
               }


            }
           
         }
      }
   }

}
