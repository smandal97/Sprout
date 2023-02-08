
#include "../defs.h"

static double x1 = 0.0;
static double x2 = 0.0;

// x1 = 0.1, x2 = 0.9, gamma = 1.4, t = 0.038

void setICParams( struct domain * theDomain ){
   x1 = theDomain->theParList.Lx * 0.1;
   x2 = theDomain->theParList.Lx * 0.9;
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0];
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;
   if( x<x1 ){      
      prim[PPP] = 1.0;
      prim[RHO] = 1.0;
   }else if( x>=x1 && x<x2 ){
      prim[PPP] = 0.01;
      prim[RHO] = 0.01;
   }else{
      prim[PPP] = 1.0;
      prim[RHO] = 1.0;
   }

}
