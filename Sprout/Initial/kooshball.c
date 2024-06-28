
#include "../defs.h"

static double d = 0.0;
static double x_cen, y_cen, z_cen, t0;
//star and environment parameters
static double M,E,eps,q,s,M_sh,v_e,r_e,v_w,t_w,r_w,sig,Mdot,rho_0,rho_w0,rho_vac;

void setICParams( struct domain * theDomain ){
   t0 = theDomain->t_init;
   if( theDomain->theParList.Num_x!=1 ) ++d;
   if( theDomain->theParList.Num_y!=1 ) ++d;
   if( theDomain->theParList.Num_z!=1 ) ++d;
   x_cen = theDomain->theParList.MM_x0 * theDomain->theParList.Lx;
   y_cen = theDomain->theParList.MM_y0 * theDomain->theParList.Ly;
   z_cen = theDomain->theParList.MM_z0 * theDomain->theParList.Lz;
   //set star and env parameters
   M       = 1.;
   E       = 1.;
   s       = 2.;
   q       = 1e-1;
   t_w     = 1e-1;
   v_w     = 1e-1;
   M_sh    = q*M;
   v_e     = sqrt(E/6./M);
   r_e     = v_e*t0;
   r_w     = v_w*t_w;
   sig     = 2.*r_w;
   Mdot    = M_sh/t_w;
   rho_0   = M/(8.*M_PI*r_e*r_e*r_e);
   rho_w0  = Mdot/4./M_PI/v_w;//M_sh*(3.-s)/4./M_PI;
   rho_vac = 1e-8;
}



void initial( double * prim , double * xi , double t ){

   double x   = xi[0] - x_cen;
   double y   = 0.;
   double z   = 0.;
   if(d>1.) y = xi[1] - y_cen;
   if(d>2.) z = xi[2] - z_cen;
   double r   = sqrt( x*x + y*y +z*z );

   double phi = acos(x/sqrt(x*x+y*y));
   double th  = acos(z/r);
   double A   = 1e-2;
   double k   = 100.;
   double ptb = A * sin( k*phi ) * sin( k*log(r) );
   if(d>2.) ptb *= sin( k*th );

   double rho_ej  = rho_0*exp(-r/r_e) * (1.+ptb);
   double rho_csm = rho_w0*pow(r,-s); 
   //if(rho_csm>rho_ej) rho_csm = 0.;
   rho_csm *= exp(-(r-r_w)*(r-r_w)/2./sig/sig);
   double rho_tot = rho_ej + rho_csm + rho_vac;

   double XX  = rho_ej/rho_tot;
   //if(r>r_st) XX = 0.0;
   double vx  = XX*x/t;
   double vy  = XX*y/t;
   double vz  = XX*z/t;
   double Pmin = 1e-6*rho_tot;


   prim[RHO] = rho_tot;
   prim[PPP] = Pmin;
   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz;
   prim[XXX] = XX;





}