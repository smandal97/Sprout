
#include "../defs.h"
//Lx = 8e-5 for m=0, eta_on=1.3 (same for m=2 for now)

static double t_min, x_zero, y_zero, z_zero, E_ej, M_ej;
static double gam, D, m, n, rho_t, vt, zeta_m, zeta_e;

void setICParams( struct domain * theDomain ){
   t_min  = theDomain->theParList.t_min;
   x_zero = theDomain->theParList.Lx/2.;
   y_zero = theDomain->theParList.Ly/2.;
   z_zero = theDomain->theParList.Lz/2.;
   D = 0.;
   if( theDomain->theParList.Num_x!=1 ) D += 1.;
   if( theDomain->theParList.Num_y!=1 ) D += 1.;
   if( theDomain->theParList.Num_z!=1 ) D += 1.;
   gam = theDomain->theParList.Adiabatic_Index;


   m = 0.0;
   n = 9.0;
   E_ej = 1e0; //2.76e13;  //these values make rho_t=vt=1 for m=0
   M_ej = 1e0; //3.85e8;
   double vt2  = E_ej/M_ej * 2.*(5.-m)*(n-5.)/((3.-m)*(n-3.));
   double zeta = (3.-m)*(n-3.)/(n-m)/4./M_PI;
   vt    = sqrt(vt2);
   rho_t = zeta*M_ej/pow(vt,3.);
   //zeta_m = (3.-m)*(n-3.)/(n-m)/4./M_PI;
   //zeta_e = (5.-m)*(n-5.)/(n-m)/2./M_PI;
   //vt = pow(zeta_e/zeta_m, 0.5) * pow(E_ej/M_ej, 0.5);
   //rho_t = pow(zeta_m*M_ej, 2.5) * pow(zeta_e*E_ej, -1.5);
   printf("vt = %e, rho_t_0 = %e, A = %e\n",vt, rho_t*pow(t_min,-D), rho_t*pow(vt,m));
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

   double phi = atan2(vy,vx);
   double th  = acos(vz/v);
   double amp = 1e-2*0.;
   double k   = 77.;
   double ptb = amp*sin(k*phi)*sin(k*th);
   

   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz;
   if(v<=vt || 1){
      prim[RHO] = rho_t * pow(v/vt , -m) * pow(t,-D) * (1.+ptb);
      prim[PPP] = pow(v/vt , -m) * pow(t,-D*gam) * 1e-4 * pow(rho_t,gam);
      prim[XXX] = 0.0;
   }
   else{
      prim[RHO] = rho_t * pow(v/vt , -n) * pow(t,-D) * (1.+ptb);
      prim[PPP] = pow(v/vt , -n) * pow(t,-D*gam) * 1e-4 * pow(rho_t,gam);
      prim[XXX] = 0.0;
   }


}
