
#include "../defs.h"

static double g_acc = 0.0;
static double y_interface = 0.0;
static double x_len = 0.0;
static double y_len = 0.0;
static double z_len = 0.0;

//t = 20.0, dirichlet_periodic BC, Lx = 0.5, Ly = 1.5, g = 0.1

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
   double w0 = 0.0025;
   if( y<=y_interface ){
      prim[RHO] = 1.0;
   }else{
      prim[RHO] = 2.0;
   }
   prim[PPP] = 5.0 + g_acc * prim[RHO] * (y-y_interface);
   prim[UU1] = 0.0;
   prim[UU2] = w0 * f(x,y);
   prim[UU3] = 0.0;

}
