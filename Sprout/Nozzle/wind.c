
#include "../defs.h"

static double L, vw, D, nr, r0;
static double nx0, ny0, nz0, K, gam;


void setNozzleParams( struct domain * theDomain ){
   nr = 4.;
   r0 = nr*theDomain->theParList.Lx/(double)(theDomain->theParList.Num_x);  //theDomain->dx;
   double t=theDomain->t;

   L = theDomain->theParList.Nozzle_pow;
   vw = theDomain->theParList.Nozzle_v; //8 for m=0
   //printf("vw=%e, r0/t=%e\n",vw, (r0/t));
   //if(vw<4e1*r0/t) vw = 4e1*r0/t; 

   double Ms = 10.0;
   gam = theDomain->theParList.Adiabatic_Index;
   K = 0.5 + 1./(gam*(gam-1)*Ms*Ms);

   nx0 = theDomain->theParList.Nozzle_x0 * theDomain->theParList.Lx;
   ny0 = theDomain->theParList.Nozzle_y0 * theDomain->theParList.Ly;
   nz0 = theDomain->theParList.Nozzle_z0 * theDomain->theParList.Lz;

   D = 0.;
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

      Vol = pow( 2.*M_PI*r0*r0 , D/2. );
      rl  = sqrt( xl*xl + yl*yl + zl*zl );
      f   = rl*rl/r0/r0*exp( -(rl*rl)/2./r0/r0 )/Vol;
      //f = 1.;
      //if(r0<rl) f = 0.;
      //if(rl<r0)
         //printf("rho*vx=%e\n",cons[SS1]);
      SE  = L*f;
      SS  = SE/vw/K;
      SM  = SS/vw;
      cons[XXX] += dt * dV * SM * 1e2;
      cons[DEN] += dt * dV * SM;
      cons[TAU] += dt * dV * SE;
      if( rl!=0. ){
         cons[SS1] += dt * dV * SS * xl/rl;
         cons[SS2] += dt * dV * SS * yl/rl;
         cons[SS3] += dt * dV * SS * zl/rl;
      }
         
   }

}
