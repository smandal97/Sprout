
#include "../defs.h"

static double L, h0, nx, gam;


void setNozzleParams( struct domain * theDomain ){
   double nh = 4.;
   h0 = nh*theDomain->dx;
   nx = theDomain->theParList.Num_x;

   L = theDomain->theParList.Nozzle_pow;
   gam = theDomain->theParList.Adiabatic_Index;
      

}


void nozz_src( double * prim , double * cons , double * x , double dx , double dy , double dz , double dt , double t ){

   double dV = dx*dy*dz;
   if(x[1]<=h0) cons[TAU] += dt * dV * L/nh/nx;

}
