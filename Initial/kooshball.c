
#include "../defs.h"

static double D, t0, x_zero, y_zero, z_zero;

void setICParams( struct domain * theDomain ){
   D = 0.;
   if( theDomain->theParList.Num_x!=1 ) D += 1.;
   if( theDomain->theParList.Num_y!=1 ) D += 1.;
   if( theDomain->theParList.Num_z!=1 ) D += 1.;
   t0 = theDomain->theParList.t_min;
   x_zero = theDomain->theParList.MM_x0 * theDomain->theParList.Lx;
   y_zero = theDomain->theParList.MM_y0 * theDomain->theParList.Ly;
   z_zero = theDomain->theParList.MM_z0 * theDomain->theParList.Lz;
}

void initial(  double * prim , double * xi , double t ){
   
   double x = xi[0]-x_zero;
   double y = 0.;
   double z = 0.;
   if(D>1.) y = xi[1]-y_zero;
   if(D>2.) z = xi[2]-z_zero;
   double r  = sqrt(x*x+y*y+z*z);

   double E = 1.0;
   double M = 1.0;

   double q = 8e-1;

   double Mshell = M*q;
   double Rshell = 2.;
   double sigma  = 0.02;

   double v0 = sqrt(E/6./M);
   double R0 = v0*t0;
   double rho0 = M/(8.*M_PI*R0*R0*R0);

   double rho_vacuum = 1e-8; 
   double rho_csm = rho_vacuum + Mshell/4./M_PI/Rshell/Rshell/sqrt(2.*M_PI)/sigma*exp(-pow(r-Rshell,2.)/2./sigma/sigma);
   double rho_ej = rho0*exp(-r/R0);

   double rho = rho_csm + rho_ej;
   double XX  = rho_ej/rho;
   double vx  = XX*x/t;
   double vy  = XX*y/t;
   double vz  = XX*z/t;
   double Pmin = 1e-6*rho;

   prim[RHO] = rho;
   prim[PPP] = Pmin;
   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz;
   prim[XXX] = XX;


}