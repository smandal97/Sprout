
#include "../defs.h"

static double GAMMA_LAW = 0.0;
static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;

void setHydroParams( struct domain * theDomain ){

   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;

}

double get_entropy( double * prim ){
   return( log( prim[PPP] / pow( prim[RHO] , GAMMA_LAW ) ) );
}


void prim2cons( double * prim , double * cons , double * x , double dV ){
   
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double vz  = prim[UU3];
   double v2  = vx*vx + vy*vy + vz*vz;
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   double phi = 0.;

   cons[DEN] = rho*dV;
   cons[SS1] = rho*vx*dV;
   cons[SS2] = rho*vy*dV;
   cons[SS3] = rho*vz*dV;
   cons[TAU] = (.5*rho*v2 + rhoe + rho*phi)*dV;
   
   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q )
      cons[q] = cons[DEN]*prim[q];

}

void cons2prim( double * cons , double * prim , double * x , double dV ){

   double rho = cons[DEN]/dV;
   double vx  = cons[SS1]/rho/dV;
   double vy  = cons[SS2]/rho/dV;
   double vz  = cons[SS3]/rho/dV;
   double E   = cons[TAU]/dV;

   double phi = 0.;

   double v2   = vx*vx + vy*vy + vz*vz;
   double rhoe = E - .5*rho*v2 - rho*phi;
   double gam  = GAMMA_LAW;
   double Pp   = (gam-1.)*rhoe;

   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   if( Pp < PRE_FLOOR*rho ) {
      Pp = PRE_FLOOR*rho; 
      //cons[XXX] += 1e10*rho*dV; 
      //printf("Pfloor at r=%e\n",sqrt(x[0]*x[0]+x[1]*x[1]));
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q )
      prim[q] = cons[q]/cons[DEN];

}

void getUstar( double * prim , double * UStar , double * x , double Sk , double Ss , double * n ){
   
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double vz  = prim[UU3];
   
   double v2  = vx*vx + vy*vy + vz*vz;
   double vn  = vx*n[0] + vy*n[1] + vz*n[2];
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   double phi = 0.;

   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   UStar[DEN] = rhostar;
   UStar[SS1] = rhostar*( vx + (Ss-vn)*n[0] );
   UStar[SS2] = rhostar*( vy + (Ss-vn)*n[1] );
   UStar[SS3] = rhostar*( vz + (Ss-vn)*n[2] );
   UStar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vn) + rhostar*phi + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q )
      UStar[q] = prim[q]*UStar[DEN];
   
}

void flux( double * prim , double * flux , double * x , double * n ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double vz  = prim[UU3];
   double v2  = vx*vx + vy*vy + vz*vz;
   double vn  = vx*n[0] + vy*n[1] + vz*n[2];
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   double phi = 0.; 
   
   flux[DEN] = rho*vn;
   flux[SS1] = rho*vx*vn + Pp*n[0];
   flux[SS2] = rho*vy*vn + Pp*n[1];
   flux[SS3] = rho*vz*vn + Pp*n[2];
   flux[TAU] = (.5*rho*v2 + rhoe + rho*phi + Pp)*vn;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q )
      flux[q] = flux[DEN]*prim[q];

}

void source( double * prim , double * cons , double * x , double dVdt ){
   //Silence is golden.
}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss , double * n ){
   
   double gam = GAMMA_LAW;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vx1  = prim1[UU1];
   double vy1  = prim1[UU2];
   double vz1  = prim1[UU3];
   double vn1  = vx1*n[0]+vy1*n[1]+vz1*n[2];

   double cs1 = sqrt(fabs(gam*P1/rho1));

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vx2  = prim2[UU1];
   double vy2  = prim2[UU2];
   double vz2  = prim2[UU3];
   double vn2  = vx2*n[0]+vy2*n[1]+vz2*n[2];

   double cs2 = sqrt(fabs(gam*P2/rho2));

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;
   
}

double mindt( double * prim , double * w , int * dim_ind , double dx , double dy , double dz ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double vz  = prim[UU3];
   double gam = GAMMA_LAW;

   double cs = sqrt( gam*Pp/rho );

   double maxvx = 0.;
   double maxvy = 0.;
   double maxvz = 0.;
   double dtx = 0.;
   double dty = 0.;
   double dtz = 0.;
   double dt = 0.;

   if( dim_ind[0]!=0 ){
      maxvx = cs + fabs( vx - w[0] );
      dtx = dx/maxvx;
      dt = dtx; //dt = 1./dtx;
   }
   if( dim_ind[1]!=0 ){
      maxvy = cs + fabs( vy - w[1] );
      dty = dy/maxvy;
      if( dt > dty ) dt = dty; //dt += 1./dty;
   }
   if( dim_ind[2]!=0 ){
      maxvz = cs + fabs( vz - w[2] );
      dtz = dz/maxvz;
      if( dt > dtz ) dt = dtz; //dt += 1./dtz;
   }
   return( dt ); //return(1./dt);

}

