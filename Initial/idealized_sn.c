
#include "../defs.h"

static double R = 0.0;
static double vmax = 0.0;

void setICParams( struct domain * theDomain ){
   R = 0.05;
   vmax = 2.; //sqrt(10./3.);
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0];
   double y = xi[1];
   double z = xi[2]*0.;
   double r = sqrt( x*x + y*y + z*z );
   if( r<=R ){
      prim[RHO] = 1./(M_PI*R*R); //1./( 4.*M_PI*R*R*R/3 );
      prim[UU1] = x/R * vmax;
      prim[UU2] = y/R * vmax;
      prim[UU3] = z/R * vmax;
      prim[XXX] = 1.;
   }else{
      prim[RHO] = 1.;
      prim[UU1] = 0.;
      prim[UU2] = 0.;
      prim[UU3] = 0.;
      prim[XXX] = 0.;
   }
   prim[PPP] = 1e-5*vmax*vmax*prim[RHO];
}
