
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
   vt = 6e1;
   A = pow(t_min,D);
   printf("vt = %e, A = %e\n",vt, A*pow(vt,2.)/3.);
}

double v_dependence( double v , double vt ){
   double m  = 2.;
   double Am = A * (3.-m)/3. * pow(vt,m);  //constant inner ejecta mass
   return Am*pow(v,-m);
   //return Am*pow(vt,-m)/( pow(v/vt,m) + pow(v/vt,9.) );
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
   double amp = 1e-2;
   double k   = 77.;
   double ptb = amp*sin(k*phi); 
   

   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz;
   if(v<=vt || 1){
      prim[RHO] = pow( t , -D ) * v_dependence(v,vt) * (1.+ptb);
      prim[PPP] = 1e-8 * pow( t , -D*gam ) * v_dependence(v,vt);
      prim[XXX] = 1e2;
   }
   else{
      prim[RHO] = A * pow(v/vt , -9.) * pow( t/t_min , -D ) * (1.+ptb);
      prim[PPP] = 1e-4 * A * pow(v/vt , -9.) * pow( t/t_min , -D );
      prim[XXX] = 1e-6;
   }
   

}
