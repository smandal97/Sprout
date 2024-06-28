
#include "../defs.h"

static double Lx    = 0.0;
static double Ly    = 0.0;
static double d     = 0.0;
static double x_cen = 0.0;
static double y_cen = 0.0;
static double z_cen = 0.0;


void setICParams( struct domain * theDomain ){
   Lx = theDomain->theParList.Lx;
   Ly = theDomain->theParList.Ly;
   x_cen = theDomain->theParList.MM_x0 * theDomain->theParList.Lx;
   y_cen = theDomain->theParList.MM_y0 * theDomain->theParList.Ly;
   z_cen = theDomain->theParList.MM_z0 * theDomain->theParList.Lz;

   if( theDomain->theParList.Num_x!=1 ) ++d;
   if( theDomain->theParList.Num_y!=1 ) ++d;
   if( theDomain->theParList.Num_z!=1 ) ++d;
}

void initial( double * prim , double * xi , double t ){

   double x   = xi[0] - x_cen;
   double y   = 0.;
   if(d>1.) y = xi[1] - y_cen;
   double z   = 0.;
   if(d>2.) z = xi[2] - z_cen;
   double r   = sqrt( x*x + y*y +z*z );
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;
   prim[RHO] = 1e0;
   prim[PPP] = 1e-5;
   prim[XXX] = 1e-6;
 
}
