
#include "../defs.h"

static double gam = 0.0;

// Lx = 1, Ly = 1, gamma = 5/3, K = 1.0, t = 1.0

void setICParams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0];
   double y = xi[1];
   double K = 1.0;
   double rho_ref = 1.0;
   prim[RHO] = rho_ref * ( 1.0 + sin(M_PI*(x+y))*sin(M_PI*(x+y)) );  
   prim[PPP] = K * pow( prim[RHO] , gam );
   double csref = sqrt( gam*K*pow( rho_ref,gam-1 ) );
   double cs = sqrt( gam*prim[PPP]/prim[RHO] );
   prim[UU1] = 1.0/sqrt(2.0)/(gam-1.0) * ( cs-csref );
   prim[UU2] = 1.0/sqrt(2.0)/(gam-1.0) * ( cs-csref );;
   prim[UU3] = 0.0;

}
