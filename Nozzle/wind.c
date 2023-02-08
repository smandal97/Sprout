
#include "../defs.h"

static double L   = 0.0;
static double vw  = 0.0;
static double D   = 0.0;
static double nr  = 0.0;
static double r0  = 0.0;
static double nx0 = 0.0;
static double ny0 = 0.0;
static double nz0 = 0.0;


void setNozzleParams( struct domain * theDomain ){
   nr = 10.;
   r0 = nr*theDomain->dx;
   double t=theDomain->t;

   L = 1e0;
   vw = 1e2;
   //if(vw<r0/t) vw = r0/t;

   nx0 = theDomain->theParList.Nozzle_x0;
   ny0 = theDomain->theParList.Nozzle_y0;
   nz0 = theDomain->theParList.Nozzle_z0;

   if( theDomain->theParList.Num_x!=1 ) D += 1.;
   if( theDomain->theParList.Num_y!=1 ) D += 1.;
   if( theDomain->theParList.Num_z!=1 ) D += 1.;
      

}


void nozz_src( double * prim , double * cons , double * x , double dx , double dy , double dz , double dt , double t ){

   if(t<1e9){
      r0 = dx*nr;
      double xl,yl,zl,rl;
      double dV, Vol, f, SE, SM, SS;
      //stuff to set location and dV of wind
      if( D==1. ){
         xl = x[0]-nx0;
         yl = 0.;
         zl = 0.;
         dV = dx;
      }else if( D==2. ){
         xl = x[0]-nx0;
         yl = x[1]-ny0;
         zl = 0.;
         dV = dx*dy;        
      }else{
         xl = x[0]-nx0;
         yl = x[1]-ny0;
         zl = x[2]-nz0;
         dV = dx*dy*dz;
      }

      Vol = pow( r0*sqrt(M_PI) , D );
      rl  = sqrt( xl*xl + yl*yl + zl*zl );
      f   = exp( -(rl*rl)/2./r0/r0 )/Vol;
      SE  = L*f;
      SS  = SE/.5/vw*.9;
      SM  = SS/vw;
      cons[XXX] += dt * dV * cons[DEN] * 1.;
      cons[DEN] += dt * dV * SM;
      cons[TAU] += dt * dV * SE;
      cons[SS1] += dt * dV * SS * xl/rl;
      cons[SS2] += dt * dV * SS * yl/rl;
      cons[SS3] += dt * dV * SS * zl/rl;
         
   }

}
