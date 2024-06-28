
#include "../defs.h"

static double r0 = 0.0;
static double p0 = 0.0;
static double xo = 0.0;
static double yo = 0.0;
static double zo = 0.0;

//domain size = 2.5e-2 for w=1

void setICParams( struct domain * theDomain ){
   r0 = 0.005;
   xo = theDomain->theParList.Lx/2.;
   yo = theDomain->theParList.Ly/2.;
   p0 = 3.5e5;
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0]-xo;
   double y = xi[1]-yo;
   double r = sqrt(x*x+y*y);
   if( r<=r0 )
      prim[PPP] = p0;
   else{
      prim[PPP] = 1e-6;
   }
   prim[RHO] = 1.0;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;

}
