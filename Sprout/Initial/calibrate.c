
#include "../defs.h"

static double x_0 = 0.0;
static double y_0 = 0.0;
static double z_0 = 0.0;


void setICParams( struct domain * theDomain ){
   x_0 = theDomain->theParList.Lx/2.*0.;
   y_0 = theDomain->theParList.Ly/2.*0.;
   z_0 = theDomain->theParList.Lz/2.*0.;
}

void initial( double * prim , double * xi , double t ){

   double x  = xi[0]-x_0;
   double y  = xi[1]-y_0;
   double z  = xi[2]-z_0;
   double r  = sqrt( x*x + y*y + z*z );
   double r0 = 1e9;
   double theta = acos(z/r); double cth = z/r;
   double phi = atan(y/x);
   //prim[RHO] = 0.;
   //prim[PPP] = 0.;
   //prim[UU1] = 0.;
   //prim[UU2] = 0.;
   //prim[UU3] = 0.;
   if( r<r0 || 1 ){
      prim[RHO] = exp(-r);//0.4886025119*z/r + 
5.*0.37317633259*(5.*z*z*z-3.*z*r*r)/(r*r*r) + 0.06828427691*(429.*pow(cth,7.)-693.*pow(cth,5.)+315.*pow(cth,3.)-35.*cth);
      prim[PPP] = 1.0;
      prim[UU1] = 2.0;
      prim[UU2] = 3.0;
      prim[UU3] = 4.0;
   }


}
