
#include "../defs.h"

static double x_speed = 0.0;
static double y_speed = 0.0;
static double z_speed = 0.0;

void setMeshMotionParams( struct domain * theDomain ){

   x_speed = theDomain->theParList.Lx*0.;
   y_speed = theDomain->theParList.Ly*0.;
   z_speed = theDomain->theParList.Lz*0.;
    
}

void set_W( struct domain * theDomain , int reset ){

}

double get_wn( double * xl , double dx , double dy , double dz , double W , int theDIM ){

   double wn;
   if( theDIM==0 ) wn = x_speed;
   if( theDIM==1 ) wn = y_speed;
   if( theDIM==2 ) wn = z_speed;
   return(wn);

}

void calc_dxs( struct domain * theDomain , double dt ){

}

void regrid( struct domain * theDomain , double dt ){

   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;

   int i,j,k,ijk;
   for( k=0 ; k<Nz+2*Ng ; ++k ){
      for( j=0 ; j<Ny+2*Ng ; ++j ){
         for( i=0 ; i<Nx+2*Ng ; ++i ){
            ijk  = i;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
            struct cell * c = theDomain->theCells+ijk;
            c->xi[0] += x_speed * dt;
            c->xi[1] += y_speed * dt;
            c->xi[2] += z_speed * dt;
         }
      }
   }

}
