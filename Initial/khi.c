
#include "../defs.h"

static double y_one  = 0.0;
static double y_two  = 0.0;

// t = 0.0-2.0, rsolution=64, dirichlet_periodic BC

void setICParams( struct domain * theDomain ){
   y_one  = theDomain->theParList.Ly*.25;
   y_two  = theDomain->theParList.Ly*.75;
}

double f( double y , double ya, double yb ){
   double twice_sigma_squared = 0.05 * 0.05;
   return( exp( -(y-ya)*(y-ya)/twice_sigma_squared ) + exp( -(y-yb)*(y-yb)/twice_sigma_squared ) );
}

void initial( double * prim , double * xi , double t ){

   double x  = xi[0];
   double y  = xi[1];
   double w0 = 0.1;
   prim[PPP] = 2.5;
   prim[UU3] = 0.0;
   prim[UU2] = w0 * sin(4.*M_PI*x) * f(y,y_one,y_two);
   if( y>=y_one && y<=y_two ){
      prim[UU1] = 0.5;
      prim[RHO] = 2.0;
      prim[XXX] = 1.0;
   }else{
      prim[UU1] = -0.5;
      prim[RHO] = 1.0;
      prim[XXX] = 0.0;
   }

}
