
#include "../defs.h"

static double t_min  = 0.0;
static double x_zero = 0.0;
static double y_zero = 0.0;
static double z_zero = 0.0;
static double gam    = 0.0;
static double A      = 0.0;
static double D      = 0.0;
static double vt     = 0.0;

void setICParams( struct domain * theDomain ){
   t_min  = theDomain->theParList.t_min;
   x_zero = theDomain->theParList.Lx/2.;
   y_zero = theDomain->theParList.Ly/2.;
   z_zero = theDomain->theParList.Lz/2.;
   if( theDomain->theParList.Num_x!=1 ) D += 1.;
   if( theDomain->theParList.Num_y!=1 ) D += 1.;
   if( theDomain->theParList.Num_z!=1 ) D += 1.;
   gam = theDomain->theParList.Adiabatic_Index;
   vt = 1.0;
   A = 1e0;
}

double v_dependence( double v , double vt ){
   return A;
}



void initial( double * prim , double * xi , double t ){
   
   double x = xi[0]-x_zero;
   double y = 0.;
   double z = 0.;
   if(D>1.) y = xi[1]-y_zero;
   if(D>2.) z = xi[2]-z_zero;

   double vx = x/t;
   double vy = y/t;
   double vz = z/t;
   double v  = sqrt(vx*vx+vy*vy+vz*vz);

   double phi = acos(vx/v);
   double A   = 4e-2;
   double k   = 43.7;
   double ptb = A*sin(k*phi); 
   

   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz;
   prim[RHO] = pow( t/t_min , -D ) * v_dependence(v,vt) * (1.+ptb);
   prim[PPP] = 1e-6 * pow( t/t_min , -D*gam ) * v_dependence(v,vt);
   prim[XXX] = 1.;

}
