
#include "../defs.h"

static double M      = 0.0;
static double G      = 0.0;
static int G_mode    = 0;

void setGravityParams( struct domain * theDomain ){
   M      = theDomain->theParList.Central_Mass;
   G      = theDomain->theParList.Grav_G;
   G_mode = theDomain->theParList.Gravity_Switch;
}


void grav_src( double * prim , double * cons , double * x , double dVdt , double t ){

   if( G_mode==1 ){
      double g_acc[3];
      g_acc[0] = -G * M * 1.;
      g_acc[1] = -G * M * 0.;
      g_acc[2] =  G * M * 0.;
      double rho = prim[RHO]; 
      cons[TAU] += dVdt * rho * ( g_acc[0]*prim[UU1] + g_acc[1]*prim[UU2] + g_acc[2]*prim[UU3] );
      cons[SS1] += dVdt * rho * g_acc[0];
      cons[SS2] += dVdt * rho * g_acc[1];
      cons[SS3] += dVdt * rho * g_acc[2];
   }

}
