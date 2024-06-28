
#include "../defs.h"

static double x_zero, y_zero, z_zero;
static double t0, gam, D, rho_0, rho_csm, v_0, r_0;

void setICParams( struct domain * theDomain ){
   t0     = theDomain->t_init;
   x_zero = theDomain->theParList.Lx/2.*1.;
   y_zero = theDomain->theParList.Ly/2.*1.;
   z_zero = theDomain->theParList.Lz/2.*1.;
   D = 0.;
   if( theDomain->theParList.Num_x!=1 ) D += 1.;
   if( theDomain->theParList.Num_y!=1 ) D += 1.;
   if( theDomain->theParList.Num_z!=1 ) D += 1.;
   gam = theDomain->theParList.Adiabatic_Index;

   double M = 1.;
   double E = 1.;
   v_0 = sqrt(E/6./M);
   r_0  = v_0*t0;
   rho_0  = M/(8.*M_PI*r_0*r_0*r_0);
   rho_csm = 1e0;//1e-2*rho_0;
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
   double r  = sqrt(x*x+y*y+z*z);
   double v  = r/t;

   double phi = atan2(vy,vx);
   double th  = acos(vz/v);
   double amp = 1e-2;
   double k   = 77.;
   double ptb = amp*sin(k*phi)*sin(k*th);
   
   double rho_ej  = rho_0*exp(-r/r_0); //pow(r/r_0,-9.);
   double rho_tot = rho_ej;
   double XX = 1.0;
   if(rho_ej<=rho_csm){
      XX = 0.0;
      rho_tot = rho_csm*(1.+ptb);
   }

   prim[RHO] = rho_tot;
   prim[PPP] = rho_tot*1e-6;
   prim[UU1] = vx*XX;
   prim[UU2] = vy*XX;
   prim[UU3] = vz*XX;
   prim[XXX] = XX;

}
