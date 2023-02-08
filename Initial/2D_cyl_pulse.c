
#include "../defs.h"

static double xo  = 0.0;
static double yo  = 0.0;
static double gam = 0.0;

// x0 = 1, y0 = 1, r = 0.5, gamma = 5/3, K = 1.0, t = 3.0

void setICParams( struct domain * theDomain ){
   xo  = theDomain->theParList.Lx/2.*0.;
   yo  = theDomain->theParList.Ly/2.*0.;
   gam = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0];
   double y = xi[1];
   double K = 1.;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0; 
   prim[RHO] = 1.0 + 3.0*exp( -80.0*((x-xo)*(x-xo)+(y-yo)*(y-yo)) );
   prim[PPP] = K * pow( prim[RHO] , gam );

}
