enum{MH,PH};

#include "../defs.h"

#define eps 1e-40

int set_accuracy(void){
   return 3;
}

void set_betas( double um1 , double u0 , double up1 , double * beta_1 , double * beta_2 , int mode ){
   if(mode==MH){
      *beta_1 = 13.*um1*um1/12. - u0*um1/6. + u0*u0/12.;
      *beta_2 = u0*u0/12. - u0*up1/6. + 13.*up1*up1/12.;
   }else{
      *beta_1 = um1*um1/12. - u0*um1/6. + 13.*u0*u0/12.;
      *beta_2 = 13.*u0*u0/12. - u0*up1/6. + up1*up1/12.;
   }
}

void set_face_values( double um1 , double u0 , double up1 , double beta_1 , double beta_2 , double * uh , int mode ){
   double w1,w2,u1_h,u2_h;
   if(mode==MH){
      w1   = 1./( 1. + (eps+beta_1)*(eps+beta_1)/2./(eps+beta_2)/(eps+beta_2)  );
      w2   = 1./( 1. + (eps+beta_2)*(eps+beta_2)*2./(eps+beta_1)/(eps+beta_1)  );
      u1_h = (u0+um1)/2.;
      u2_h = (3.*u0-up1)/2.;     
   }else{
      w1   = 1./( 1. + (eps+beta_1)*(eps+beta_1)*2./(eps+beta_2)/(eps+beta_2)  );
      w2   = 1./( 1. + (eps+beta_2)*(eps+beta_2)/2./(eps+beta_1)/(eps+beta_1)  );
      u1_h = (3.*u0-um1)/2.;
      u2_h = (u0+up1)/2.;
   }
   *uh  = w1*u1_h + w2*u2_h; 
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
   double pL,pC,pR;
   double beta_1,beta_2,pimh,piph;
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
               set_betas(pL,pC,pR,&beta_1,&beta_2,MH);
               set_face_values(pL,pC,pR,beta_1,beta_2,&pimh,MH);
               set_betas(pL,pC,pR,&beta_1,&beta_2,PH);
               set_face_values(pL,pC,pR,beta_1,beta_2,&piph,PH);
               if( theDIM==0 ){
                  c->gradx[q] = (piph-pimh)/dl;
                  c->pblax[q] = (piph+pimh)/2.-pC;
               }
               if( theDIM==1 ){
                  c->grady[q] = (piph-pimh)/dl;
                  c->pblay[q] = (piph+pimh)/2.-pC;
               }
               if( theDIM==2 ){
                  c->gradz[q] = (piph-pimh)/dl;
                  c->pblaz[q] = (piph+pimh)/2.-pC;
               }
           }
         }
      }
   }

}
