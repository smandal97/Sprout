
#include "../defs.h"

static double x0  = 0.0;
static double gam = 0.0;

// L = 0.3, gamma = 5/3, t = 0.8

void setICParams( struct domain * theDomain ){
   x0  = theDomain->theParList.Lx * 0.;
   gam = theDomain->theParList.Adiabatic_Index;
}


double func( double x ){
   double val = 0.;
   //val = pow( x*x/L/L-1. , 2. );
   val = 3.0 * exp( -80.*(x-x0)*(x-x0) );
   return val;
}

void initial( double * prim , double * xi , double t ){

   double x    = xi[0];
   double K    = 1.;
   double rho  = 1. + func(x);
   double P    = K * pow( rho , gam );

   prim[RHO] = rho;
   prim[PPP] = P;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;
   

}
