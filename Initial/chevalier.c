
#include "../defs.h"

static double t0      = 0.0;
static double R0      = 0.0;
static double Rcutoff = 0.0;
static double rho_0   = 0.0;
static double d       = 0.0;
static double x_cen   = 0.0;
static double y_cen   = 0.0;
static double z_cen   = 0.0;

void setICParams( struct domain * theDomain ){
   t0 = theDomain->t_init;
   R0 = theDomain->theParList.Lx/2.;
   rho_0 = 1.0;
   Rcutoff = theDomain->dx*4.;
   if( theDomain->theParList.Num_x!=1 ) ++d;
   if( theDomain->theParList.Num_y!=1 ) ++d;
   if( theDomain->theParList.Num_z!=1 ) ++d;
   x_cen = theDomain->theParList.MM_x0;
   y_cen = theDomain->theParList.MM_y0;
   z_cen = theDomain->theParList.MM_z0;
}


void initial( double * prim , double * xi , double t ){

   double n = 10.0;
   double s = 2.0;
   double A = 1e-2;
   double k = 10.5;

   double x   = xi[0] - x_cen;
   double y   = 0.;
   if(d>1.) y = xi[1] - y_cen;
   double z   = 0.;
   if(d>2.) z = xi[2] - z_cen;
   double r   = sqrt( x*x + y*y +z*z );
   double phi = acos(x/sqrt(x*x+y*y));
   double th  = acos(z/r);
   double Rc  = R0 * pow( t/t0 , (n-d)/(n-s) );
   double ptb = A * sin( k*phi ) * sin( k*log(r) );
   if(d>2.) ptb *= sin( k*th );

   double rho,vx,vy,vz,ps;

   if( r<Rc ){
      rho = rho_0 / ( pow(R0/Rcutoff,-n) + pow(r/R0,n) ) * pow( t/t0 , -1.*d );
      vx  = x/t;
      vy  = y/t;
      vz  = z/t;
      ps  = 1.;
   }else if( r>=Rc ){
      rho = rho_0 * pow( R0/r , s ) * ( 1.+ptb );
      vx  = 0.;
      vy  = 0.;
      vz  = 0.;
      ps  = 1e-12;
   }

   prim[RHO] = rho;
   prim[PPP] = 1e-6*rho;
   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz; 
   prim[XXX] = ps;


}
