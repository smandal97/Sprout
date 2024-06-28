
#include "../defs.h"

static double xo  = 0.0;
static double yo  = 0.0;
static double zo  = 0.0;
static double gam = 0.0;

// r = 0.5, gamma = 5/3, K = 1.0, t = 0.1
// use H0 or W0 = 4.0,t0 = 1.0, and r = 0.3 for MM
// or maybe follow shock

void setICParams( struct domain * theDomain ){
   xo  = theDomain->theParList.Lx/2.*0.;
   yo  = theDomain->theParList.Ly/2.*0.;
   zo  = theDomain->theParList.Lz/2.*0.;
   gam = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0];
   double y = xi[1];
   double z = xi[2];
   double K = 1.;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0; 
   prim[RHO] = 1.0 + 3.0*exp( -80.0*((x-xo)*(x-xo)+(y-yo)*(y-yo)+(z-zo)*(z-zo)) );
   prim[PPP] = K * pow( prim[RHO] , gam );

}