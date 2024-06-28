
#include "../defs.h"

static double r0 = 0.0;
static double p0 = 0.0;
static double xo = 0.0;
static double yo = 0.0;
static double zo = 0.0;

//domain size = 2.5e-2 for w=1

void setICParams( struct domain * theDomain ){
   r0 = 0.02;
   xo = theDomain->theParList.Lx/2.*0.;
   yo = theDomain->theParList.Ly/2.*0.;
   zo = theDomain->theParList.Lz/2.*0.;
   double E = 1.0;
   double rhoe = E/( 4.*M_PI*r0*r0*r0/3. );
   p0 = rhoe*(theDomain->theParList.Adiabatic_Index - 1.);
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0]-xo;
   double y = xi[1]-yo;
   double z = xi[2]-zo;
   double r = sqrt(x*x+y*y+z*z);
   if( r<=r0 )
      prim[PPP] = p0;
   else{
      prim[PPP] = 1e-4;
   }
   prim[RHO] = 1.0 * pow(r,-1.);
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;

}
