
#include "../defs.h"

static double r0 = 0.0;

void setICParams( struct domain * theDomain ){
   r0 = theDomain->theParList.Lx/2.;
}

void initial( double * prim , double * xi ){

   double x = xi[0];
   double y = xi[1];
   double r = sqrt(x*x+y*y);
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;
   if( r<r0 ){
      prim[RHO] = 1.0;
      prim[PPP] = 1.0;
   }else{
      prim[RHO] = 0.125;
      prim[PPP] = 0.1;
   }

}
