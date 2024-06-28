
#include "../defs.h"

//Joggerst et al. 2014 Journal of Computational Physics

static double x_zero, y_zero, z_zero, t0, gam, D;
static double R1, R2, R3, Rp, n, A, rho_0, u_b;

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

   rho_0 = 0.05;
   R1 = 10.; 
   R2 = 12.; 
   R3 = 15.;
   Rp = R1/10.;
   n  = 0.;
   A  = rho_0*(3.-n)/(3.*pow(Rp/R1,n) - n*pow(Rp/R1,3.));
}




void initial( double * prim , double * xi , double t ){
   
   double x = xi[0]-x_zero;
   double y = 0.;
   double z = 0.;
   if(D>1.) y = xi[1]-y_zero;
   if(D>2.) z = xi[2]-z_zero;
   double r  = sqrt(x*x+y*y+z*z);

   double rho, P, vx, vy, vz, ps;
   vx = 0.; vy = 0.; vz = 0.;

   double phi = atan2(y,x);
   double th  = acos(z/r);
   double amp = 1e-2;
   double k   = 77.;
   double ptb = amp*sin(k*phi)*sin(k*th);
   
   if(r<=R1*(1.+ptb)){
      rho = rho_0;
      P   = 0.1;
      ps  = 0.;
   }else if(r>R1*(1.+ptb) && r<=R2){
      rho = rho_0*20.; 
      P   = 0.1;
      ps  = 10.;
   }else if(r>R2){
      rho = rho_0*2.;
      //P   = 10.;
      ps  = 5.;
      if(t<0.5) P = 10.;
      else if(t>=0.5 && t<=3.) P = 11.98 - 3.96*t;
      //vx = -x/(R2/0.24-t);
      //vy = -y/(R2/0.24-t);
      //vz = -z/(R2/0.24-t);
   }

   prim[RHO] = rho;
   prim[PPP] = P;
   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz;
   prim[XXX] = ps;

}
