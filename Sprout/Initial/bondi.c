
#include "../defs.h"

static double M_cen = 0.0;
static double G     = 0.0;
static double gam   = 0.0;
static double rmax  = 0.0;
static double Mdot  = 0.0;
static double eta   = 0.0;


void setICParams( struct domain * theDomain ){
   M_cen = theDomain->theParList.Central_Mass;
   G     = theDomain->theParList.Grav_G;
   gam   = theDomain->theParList.Adiabatic_Index;
   rmax  = theDomain->theParList.rmax;
   Mdot  = theDomain->theParList.Mdot;
   eta   = theDomain->theParList.eta_P;
}

void initial( double * prim , double * xi , double t  ){
 
   double K   = eta * pow( (4.*M_PI/Mdot) , (gam-1.) );
   K *= pow( (2.*G*M_cen) , ((gam+1.)/2.) );
   K *= pow( rmax , ((3.*gam-5.)/2.) );
   double v   = -sqrt( 2.*G*M_cen/r );
   double rho = Mdot/4./M_PI/r/r/fabs(v);
   prim[VRR]  = v;
   prim[RHO]  = rho;
   prim[PPP]  = K*pow( rho,gam );

}
