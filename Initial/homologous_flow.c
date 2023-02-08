
#include "../defs.h"

static double t_min  = 0.0;
static double x_zero = 0.0;
static double y_zero = 0.0;
static double z_zero = 0.0;
static double gam    = 0.0;
static double D      = 0.0;
static double r0     = 0.0;

void setICParams( struct domain * theDomain ){
   t_min  = theDomain->theParList.t_min;
   x_zero = theDomain->theParList.Lx/2.*0.;
   y_zero = theDomain->theParList.Ly/2.*0.;
   z_zero = theDomain->theParList.Lz/2.*0.;
   if( theDomain->theParList.Num_x!=1 ) D += 1.;
   if( theDomain->theParList.Num_y!=1 ) D += 1.;
   if( theDomain->theParList.Num_z!=1 ) D += 1.;
   gam = theDomain->theParList.Adiabatic_Index;
   r0 = 0.1;//theDomain->dx*4.;
}

double f(double r){
   double val = 1.0;
   //if(r>r0) val = 1e0 * pow(r/r0, -2.);
   return val;
}

/*
double blob(double x , double y , double z, double r){
   double a = 0.;
   double rl = sqrt(x*x+y*y+z*z);
   if(rl<r) a = 1.;
   return a;
}

double psvals(double x , double y , double z){
   double psval = 0.0;
   double b = 0.18;
   double r = 0.07;
   double xc,yc,zc;
   if(D>0.) psval += blob(x-b,y,z,r) + blob(x+b,y,z,r);
   if(D>1.) psval += blob(x,y-b,z,r) + blob(x,y+b,z,r);
   if(D>2.) psval += blob(x,y,z-b,r) + blob(x,y,z+b,r);
   return psval;

}
*/

double psvals(double x , double y , double z){
   double psval, r_e, r_ul, r_ll, r2_1, r2_2, r2_3, r2_4;
   double xe1, ye1, xe2, ye2, xu, yu, xl, yl;
   psval = 0.0;

   r_e   = x_zero*0.15;
   xe1   = -x_zero*0.25;
   ye1   = y_zero*0.25;
   xe2   = x_zero*0.25;
   ye2   = y_zero*0.25;

   r_ul  = x_zero*1.1;
   r_ll  = x_zero*0.5;
   xu    = 0.;
   yu    = y_zero*1.0;
   xl    = 0.;
   yl    = y_zero*0.1;

   r2_1  = (x-xe1)*(x-xe1) + (y-ye1)*(y-ye1);
   r2_2  = (x-xe2)*(x-xe2) + (y-ye2)*(y-ye2);
   r2_3  = (x-xu)*(x-xu) + (y-yu)*(y-yu);
   r2_4  = (x-xl)*(x-xl) + (y-yl)*(y-yl);

   //printf("%e\n",x_zero);
   if( r2_1<r_e*r_e || r2_2<r_e*r_e || (r2_3>r_ul*r_ul && r2_4<r_ll*r_ll) ){
   //if( r2_4<r_ll*r_ll ){
      psval = 1.0;
   }

   return psval*0.+1.;

}


void initial( double * prim , double * xi , double t ){
   
   double x = xi[0]-x_zero;
   double y = 0.;
   double z = 0.;
   if(D>1.) y = xi[1]-y_zero;
   if(D>2.) z = xi[2]-z_zero;
   double r = sqrt(x*x+y*y+z*z);
   double m = 1.;//0.7777777777777777;

   prim[RHO] = f(r) * pow( t/t_min , -D*m );
   prim[PPP] = f(r) * 1e-4 * pow( t/t_min , -D*gam*m );
   prim[UU1] = x/t*m;
   prim[UU2] = y/t*m;
   prim[UU3] = z/t*m;
   prim[XXX] = psvals(x,y,z);

}