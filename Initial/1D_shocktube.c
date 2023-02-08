
#include "../defs.h"

static double x0 = 0.0;

// x0 = 0.5, gamma = 1.66666667, t = 0.1

void setICParams( struct domain * theDomain ){
   x0 = theDomain->theParList.Lx/2.;
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0];
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;
   if( x<x0 ){
      prim[RHO] = 1.0;
      prim[PPP] = 1.0;
   }else if( x>=x0 ){
      prim[RHO] = 0.125;
      prim[PPP] = 0.1;
   }

}
