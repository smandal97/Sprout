
#include "../defs.h"

static double g_acc = 0.0;
static double x_len = 0.0;
static double y_len = 0.0;
static double z_len = 0.0;



void setICParams( struct domain * theDomain ){
   y_interface = theDomain->theParList.Ly/2.;
   x_len = theDomain->theParList.Lx;
   y_len = theDomain->theParList.Ly;
   z_len = theDomain->theParList.Lz;
   g_acc  = -1.0 * theDomain->theParList.Central_Mass * theDomain->theParList.Grav_G;
}

double f( double x , double y ){
   return( ( 1.-cos(2.*M_PI*x/x_len) )*( 1. - cos(2.*M_PI*y/y_len) ) );
}

void initial( double * prim , double * xi , double t ){
   
   double x = xi[0];
   double y = xi[1];
   double z = xi[2];
   prim[RHO] = 1.0;
   prim[PPP] = 2.5 + prim[RHO] * y;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;

}
