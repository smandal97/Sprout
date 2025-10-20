
#include "../defs.h"

static double xo, yo, A0, B0, R;

// Lx = 2, Ly = 1, R = 0.3

void setICParams( struct domain * theDomain ){
   xo  = theDomain->theParList.Lx/4.;
   yo  = theDomain->theParList.Ly/2.;
   R   = 0.25;
   A0  = 1e-3;
   B0  = 1e-4;
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0];
   double y = xi[1];
   double r = sqrt((x-xo)*(x-xo)+(y-yo)*(y-yo));
   double Bphi = B0*pow(sin(M_PI*r/R) , 2.)*sqrt(2.*r/R);
   
   prim[RHO] = 1.;
   prim[PPP] = 1.;
   prim[UU1] = 1.;
   prim[UU2] = 0.;
   prim[UU3] = 0.;
   if(r<R){
      prim[BB1] = -Bphi * (y-yo)/r;
      prim[BB2] = Bphi * (x-xo)/r;
   }else{
      prim[BB1] = 0.;          
      prim[BB2] = 0.;
   }
   prim[BB3] = 0.;

}


/* Gardiner and Stone 2005; Athena; A0 = 1e-3
   prim[RHO] = 1.;
   prim[PPP] = 1.;
   prim[UU1] = 2.;
   prim[UU2] = 1.;
   prim[UU3] = 0.0;
   if(r<R && r!=0.){
      prim[BB1] = A0*(R-y/r);
      prim[BB2] = A0*(x/r-R);
   }else{
      prim[BB1] = 0.;          
      prim[BB2] = 0.;
   }
   prim[BB3] = 0.;
*/
