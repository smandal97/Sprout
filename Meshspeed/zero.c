
#include "../defs.h"

static double x_cen = 0.0;
static double y_cen = 0.0;
static double z_cen = 0.0;

void setMeshMotionParams( struct domain * theDomain ){

   x_cen = theDomain->theParList.MM_x0;
   y_cen = theDomain->theParList.MM_y0;
   z_cen = theDomain->theParList.MM_z0;
    
}

void set_W( struct domain * theDomain , int reset ){

   theDomain->W = 0.0;

}

double get_wn( double * xl , double dx , double dy , double dz , double W , int theDIM ){

   double wn = W;
   if( theDIM==0 ) wn *= ( xl[0] + dx/2. - x_cen );
   if( theDIM==1 ) wn *= ( xl[1] + dy/2. - y_cen );
   if( theDIM==2 ) wn *= ( xl[2] + dz/2. - z_cen );
   return(wn);

}

void calc_dxs( struct domain * theDomain , double dt ){

   double W  = theDomain->W;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   dx *= ( 1. + W*dt );
   dy *= ( 1. + W*dt );
   dz *= ( 1. + W*dt );

   theDomain->dx = dx;
   theDomain->dy = dy;
   theDomain->dz = dz; 

}

void regrid( struct domain * theDomain , double dt ){

   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double W  = theDomain->W;

   int i,j,k,ijk;
   for( k=0 ; k<Nz+2*Ng ; ++k ){
      for( j=0 ; j<Ny+2*Ng ; ++j ){
         for( i=0 ; i<Nx+2*Ng ; ++i ){
            ijk  = i;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
            struct cell * c = theDomain->theCells+ijk;
            c->xi[0] = c->xi[0] * ( 1. + W*dt ) - x_cen*W*dt;
            c->xi[1] = c->xi[1] * ( 1. + W*dt ) - y_cen*W*dt;
            c->xi[2] = c->xi[2] * ( 1. + W*dt ) - z_cen*W*dt;
         }
      }
   }

}
