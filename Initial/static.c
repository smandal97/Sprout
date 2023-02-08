
#include "../defs.h"

static double Lx = 0.0;
static double Ly = 0.0;

double f_2dgrid( double x , double y ){
   if( x>Lx/3. && x<Lx*2./3. && y>Ly/3. && y<Ly*2./3. ) return 1.;
   else return 1e-5;
}

void setICParams( struct domain * theDomain ){
   Lx = theDomain->theParList.Lx;
   Ly = theDomain->theParList.Ly;
}

void initial( double * prim , double * xi , double t ){

   double x = xi[0];
   //double y = xi[1];
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;
   prim[RHO] = 1.0; //f_2dgrid(x,y);
   prim[PPP] = 1e-4;
   //prim[XXX] = f_2dgrid(x,y);
 
}
