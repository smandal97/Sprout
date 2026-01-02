
#include "../defs.h"

static double GAMMA_LAW = 0.0;
static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;

#define FL 1e-10

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

   double Bx  = prim[BB1];
   double By  = prim[BB2];
   double Bz  = prim[BB3];
   double B2  = Bx*Bx+By*By+Bz*Bz;

   cons[DEN] = rho*dV;
   cons[SS1] = rho*vx*dV;
   cons[SS2] = rho*vy*dV;
   cons[SS3] = rho*vz*dV;
   cons[TAU] = (.5*rho*v2 + rhoe + rho*phi + .5*B2)*dV;
   
   cons[BB1] = Bx*dV;
   cons[BB2] = By*dV;
   cons[BB3] = Bz*dV;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q )
      cons[q] = cons[DEN]*prim[q];

}

void cons2prim( double * cons , double * prim , double * x , double dV ){

   double rho = cons[DEN]/dV;
   double vx  = cons[SS1]/rho/dV;
   double vy  = cons[SS2]/rho/dV;
   double vz  = cons[SS3]/rho/dV;
   double E   = cons[TAU]/dV;

   double Bx  = cons[BB1]/dV;
   double By  = cons[BB2]/dV;
   double Bz  = cons[BB3]/dV;
   double B2  = Bx*Bx+By*By+Bz*Bz;

   double phi = 0.;

   double v2   = vx*vx + vy*vy + vz*vz;
   double rhoe = E - .5*rho*v2 - rho*phi - .5*B2;
   double gam  = GAMMA_LAW;
   double Pp   = (gam-1.)*rhoe;

   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   if( Pp < PRE_FLOOR*rho ) {
      Pp = PRE_FLOOR*rho; 
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = vx;
   prim[UU2] = vy;
   prim[UU3] = vz;

   prim[BB1] = Bx;
   prim[BB2] = By;
   prim[BB3] = Bz;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q )
      prim[q] = cons[q]/cons[DEN];

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

   double Bx  = prim[BB1];
   double By  = prim[BB2];
   double Bz  = prim[BB3];
   double B2  = Bx*Bx+By*By+Bz*Bz;
   double Bn  = Bx*n[0] + By*n[1] + Bz*n[2];
   double vB  = vx*Bx + vy*By + vz*Bz;

   double phi = 0.; 
   
   flux[DEN] = rho*vn;
   flux[SS1] = rho*vx*vn + (Pp+.5*B2)*n[0] - Bx*Bn;
   flux[SS2] = rho*vy*vn + (Pp+.5*B2)*n[1] - By*Bn;
   flux[SS3] = rho*vz*vn + (Pp+.5*B2)*n[2] - Bz*Bn;
   flux[TAU] = (.5*rho*v2 + rhoe + rho*phi + Pp + B2)*vn - vB*Bn;

   flux[BB1] = vn*Bx - Bn*vx;
   flux[BB2] = vn*By - Bn*vy;
   flux[BB3] = vn*Bz - Bn*vz;


   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q )
      flux[q] = flux[DEN]*prim[q];

}

void source( double * prim , double * cons , double * x , double dVdt ){
   //Silence is golden.
}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss , double * 
n ){
   
 //only for the HLL(C) solvers
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

void vel_and_Pts( double * primL , double * primR , double * vel , double * Pt_star , double * n ){
   
   double gam = GAMMA_LAW;

   double P_L   = primL[PPP];
   double rho_L = primL[RHO];
   double vn_L  = primL[UU1]*n[0]       + primL[UU2]*n[1]       + primL[UU3]*n[2];
   double Bn_L  = primL[BB1]*n[0]       + primL[BB2]*n[1]       + primL[BB3]*n[2];
   double B2_L  = primL[BB1]*primL[BB1] + primL[BB2]*primL[BB2] + primL[BB3]*primL[BB3];
   double Pt_L  = P_L + .5*B2_L;

   double P_R   = primR[PPP];
   double rho_R = primR[RHO];
   double vn_R  = primR[UU1]*n[0]       + primR[UU2]*n[1]       + primR[UU3]*n[2];
   double Bn_R  = primR[BB1]*n[0]       + primR[BB2]*n[1]       + primR[BB3]*n[2];
   double B2_R  = primR[BB1]*primR[BB1] + primR[BB2]*primR[BB2] + primR[BB3]*primR[BB3];
   double Pt_R  = P_R + .5*B2_R;

   double Bn = .5*(Bn_L+Bn_R);
 
   double cf2_L = .5*(gam*P_L + B2_L + sqrt( fabs( pow( gam*P_L + B2_L , 2. ) - 4.*gam*P_L*Bn_L*Bn_L ) ) )/rho_L;
   double cf2_R = .5*(gam*P_R + B2_R + sqrt( fabs( pow( gam*P_R + B2_R , 2. ) - 4.*gam*P_R*Bn_R*Bn_R ) ) )/rho_R;

   double cf2;
   if( cf2_L > cf2_R ) cf2 = cf2_L; else cf2 = cf2_R;
   double umin,umax;
   if( vn_L < vn_R ){ umin = vn_L; umax = vn_R; }else{ umin = vn_R; umax = vn_L; }

   double S_R = umax + sqrt(fabs(cf2));
   double S_L = umin - sqrt(fabs(cf2));

   double S_M = ( ( S_R - vn_R )*rho_R*vn_R - ( S_L - vn_L )*rho_L*vn_L - Pt_R + Pt_L )/( ( S_R - vn_R )*rho_R - (S_L-vn_L)*rho_L );

   double rhostar_L = rho_L*( S_L - vn_L )/( S_L - S_M );
   double rhostar_R = rho_R*( S_R - vn_R )/( S_R - S_M );

   double Ss_L = S_M - fabs(Bn)/sqrt(rhostar_L);
   double Ss_R = S_M + fabs(Bn)/sqrt(rhostar_R);

   vel[0] = S_L;
   vel[1] = Ss_L;
   vel[2] = S_M;
   vel[3] = Ss_R;
   vel[4] = S_R;
   
   double num   = ( S_R - vn_R )*rho_R*Pt_L - ( S_L - vn_L )*rho_L*Pt_R + rho_L*rho_R*( S_R - vn_R )*( S_L - vn_L )*( vn_R - vn_L );
   double denom = ( S_R - vn_R )*rho_R - ( S_L - vn_L )*rho_L;

   *Pt_star = num/denom;
   
}


void get_single_star( double * prim , double * Ustar , double S_K , double S_M , double Pt_star , double * n ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double vz  = prim[UU3];
   double Bx  = prim[BB1];
   double By  = prim[BB2];
   double Bz  = prim[BB3];

   double u_K   = vx*n[0] + vy*n[1] + vz*n[2];
   double Bn    = Bx*n[0] + By*n[1] + Bz*n[2];
   double vdotB = vx*Bx + vy*By + vz*Bz;
   double v2    = vx*vx + vy*vy + vz*vz;
   double B2    = Bx*Bx + By*By + Bz*Bz;

   double Pt = Pp + .5*B2;
 
   double gam = GAMMA_LAW;
   double e = .5*rho*v2 + Pp/(gam-1.) + .5*B2;

   double rhostar = rho*( S_K - u_K )/( S_K - S_M );

   double denom = rho*(S_K-u_K)*(S_K-S_M) - Bn*Bn ;
   if( fabs(denom) < FL ) denom = FL;

   double vx_star = vx - Bn*Bx*( S_M - u_K )/denom;
   double vy_star = vy - Bn*By*( S_M - u_K )/denom;
   double vz_star = vz - Bn*Bz*( S_M - u_K )/denom;
 
   double Bx_star = Bx*( rho*pow(S_K-u_K,2.) - Bn*Bn )/denom;
   double By_star = By*( rho*pow(S_K-u_K,2.) - Bn*Bn )/denom;
   double Bz_star = Bz*( rho*pow(S_K-u_K,2.) - Bn*Bn )/denom;

   double vn_star = vx_star*n[0] + vy_star*n[1] + vz_star*n[2];
   double Bn_star = Bx_star*n[0] + By_star*n[1] + Bz_star*n[2];

   vx_star += ( S_M - vn_star )*n[0];
   vy_star += ( S_M - vn_star )*n[1];
   vz_star += ( S_M - vn_star )*n[2];

   Bx_star += ( Bn  - Bn_star )*n[0];
   By_star += ( Bn  - Bn_star )*n[1];
   Bz_star += ( Bn  - Bn_star )*n[2];

   double vBstar = vx_star*Bx_star + vy_star*By_star + vz_star*Bz_star;
   double e_star = ( (S_K-u_K)*e - Pt*u_K + Pt_star*S_M + Bn*( vdotB - vBstar ) )/( S_K - S_M );

   Ustar[DEN] = rhostar;
   Ustar[SS1] = rhostar*vx_star;
   Ustar[SS2] = rhostar*vy_star;
   Ustar[SS3] = rhostar*vz_star;
   Ustar[TAU] = e_star;
   Ustar[BB1] = Bx_star;
   Ustar[BB2] = By_star;
   Ustar[BB3] = Bz_star;

}


void get_double_star( double * UsL , double * UsR , double * Uss , double * vel , double * n , int LR ){

   double rho_L = UsL[RHO];
   double rho_R = UsR[RHO];
   double rhostar_K;
   if( LR==0 ) rhostar_K = rho_L; else rhostar_K = rho_R;
   double rrh_L = sqrt( rho_L );
   double rrh_R = sqrt( rho_R );
   double rrh_K = sqrt( rhostar_K );
  
   double Sx_L = UsL[SS1];
   double Sy_L = UsL[SS2];
   double Sz_L = UsL[SS3];
   double Bx_L = UsL[BB1];
   double By_L = UsL[BB2];
   double Bz_L = UsL[BB3];

   double vx_L = Sx_L/rho_L;
   double vy_L = Sy_L/rho_L;
   double vz_L = Sz_L/rho_L;
 
   double Sx_R = UsR[SS1];
   double Sy_R = UsR[SS2];
   double Sz_R = UsR[SS3];
   double Bx_R = UsR[BB1];
   double By_R = UsR[BB2];
   double Bz_R = UsR[BB3];

   double vx_R = Sx_R/rho_R;
   double vy_R = Sy_R/rho_R;
   double vz_R = Sz_R/rho_R;
 
   double Bn_L = Bx_L*n[0] + By_L*n[1] + Bz_L*n[2];
   double Bn_R = Bx_R*n[0] + By_R*n[1] + Bz_R*n[2];
   double Bn = .5*(Bn_L+Bn_R);
   double signBn;
   if( Bn>0.0 ) signBn = 1.0; else signBn = -1.0;
 
   double denom = rrh_L+rrh_R;

   double vx_ss = ( rrh_L*vx_L + rrh_R*vx_R + ( Bx_R - Bx_L )*signBn )/denom;
   double vy_ss = ( rrh_L*vy_L + rrh_R*vy_R + ( By_R - By_L )*signBn )/denom;
   double vz_ss = ( rrh_L*vz_L + rrh_R*vz_R + ( Bz_R - Bz_L )*signBn )/denom;

   double vn_ss = vx_ss*n[0]+vy_ss*n[1]+vz_ss*n[2];
   vx_ss += (vel[2]-vn_ss)*n[0];
   vy_ss += (vel[2]-vn_ss)*n[1];
   vz_ss += (vel[2]-vn_ss)*n[2];

   double Bx_ss = ( rrh_L*Bx_R + rrh_R*Bx_L + rrh_L*rrh_R*( vx_R - vx_L )*signBn )/denom;
   double By_ss = ( rrh_L*By_R + rrh_R*By_L + rrh_L*rrh_R*( vy_R - vy_L )*signBn )/denom;
   double Bz_ss = ( rrh_L*Bz_R + rrh_R*Bz_L + rrh_L*rrh_R*( vz_R - vz_L )*signBn )/denom;

   double Bn_ss = Bx_ss*n[0]+By_ss*n[1]+Bz_ss*n[2];
   Bx_ss += (Bn-Bn_ss)*n[0];
   By_ss += (Bn-Bn_ss)*n[1];
   Bz_ss += (Bn-Bn_ss)*n[2];

   double vB_ss;
   vB_ss = vx_ss*Bx_ss + vy_ss*By_ss + vz_ss*Bz_ss;

   double vBstar;
   double e_star;
   if( LR==0 ){
      vBstar = vx_L*Bx_L + vy_L*By_L + vz_L*Bz_L;
      e_star = UsL[TAU];
   }else{
      vBstar = vx_R*Bx_R + vy_R*By_R + vz_R*Bz_R;
      e_star = UsR[TAU];
   }
   double plusminus = -1.0;
   if( LR==1 ) plusminus = 1.0;

   double e_ss = e_star + plusminus*rrh_K*( vBstar - vB_ss )*signBn;

   Uss[DEN] = rhostar_K;
   Uss[SS1] = rhostar_K*vx_ss;
   Uss[SS2] = rhostar_K*vy_ss;
   Uss[SS3] = rhostar_K*vz_ss;
   Uss[TAU] = e_ss;
   Uss[BB1] = Bx_ss;
   Uss[BB2] = By_ss;
   Uss[BB3] = Bz_ss;

}

/*
void get_Ustar_HLLD( double w , double * primL , double * primR , double * F , double * U , double * n . double * x , double dV){

   double vel[5];
   double Pt_star;
   vel_and_Pts( primL , primR , vel , &Pt_star , n );

   int q;

   if( w < vel[2] ){
      //Left Side
      double UL[NUM_Q];
      flux( primL , F , x , n );
      prim2cons( primL , UL , x , dV );
      if( w < vel[0] ){
         //Upwind Left
         prim2cons( primL , U , x , dV );
      }else if( w < vel[1] ){
         //Single Star Left
         get_single_star( primL , U , vel[0] , vel[2] , Pt_star , n );
         for( q=0 ; q<NUM_Q ; ++q ){
            F[q] += vel[0]*( U[q] - UL[q] );
         }
      }else{
         //Double Star Left
         double UsL[NUM_Q];
         double UsR[NUM_Q];
         get_single_star( primL , UsL , vel[0] , vel[2] , Pt_star , n );
         get_single_star( primR , UsR , vel[4] , vel[2] , Pt_star , n );
         get_double_star( UsL , UsR , U , vel , n , 0 );
         for( q=0 ; q<NUM_Q ; ++q ){
            F[q] += vel[1]*( U[q] - UsL[q] ) + vel[0]*( UsL[q] - UL[q] );
         }
      }
   }else{
      //Right Side
      double UR[NUM_Q];
      flux( primR , F , x , n );
      prim2cons( primR , UR , x , dV );
      if( w > vel[4] ){
         //Upwind Right
         prim2cons( primR , U , x , dV );
      }else if( w > vel[3] ){
         //Single Star Right
         get_single_star( primR , U , vel[4] , vel[2] , Pt_star , n );
         for( q=0 ; q<NUM_Q ; ++q ){
            F[q] += vel[4]*( U[q] - UR[q] );
         }
      }else{
         //Double Star Right
         double UsL[NUM_Q];
         double UsR[NUM_Q];
         get_single_star( primL , UsL , vel[0] , vel[2] , Pt_star , n );
         get_single_star( primR , UsR , vel[4] , vel[2] , Pt_star , n );
         get_double_star( UsL , UsR , U , vel , n , 1 );
         for( q=0 ; q<NUM_Q ; ++q ){
            F[q] += vel[3]*( U[q] - UsR[q] ) + vel[4]*( UsR[q] - UR[q] );
         }
      }
   }


}
*/


double mindt( double * prim , double * w , int * dim_ind , double dx , double dy , double dz ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double vz  = prim[UU3];
   double Bx  = prim[BB1];
   double By  = prim[BB2];
   double Bz  = prim[BB3];
   double gam = GAMMA_LAW;

   double B2  = Bx*Bx + By*By + Bz*Bz;

   double cf = sqrt( (gam*Pp+B2)/rho );

   double maxvx = 0.;
   double maxvy = 0.;
   double maxvz = 0.;
   double dtx = 0.;
   double dty = 0.;
   double dtz = 0.;
   double dt = 0.;

   if( dim_ind[0]!=0 ){
      maxvx = cf + fabs( vx - w[0] );
      dtx = dx/maxvx;
      dt = dtx; //dt = 1./dtx;
   }
   if( dim_ind[1]!=0 ){
      maxvy = cf + fabs( vy - w[1] );
      dty = dy/maxvy;
      if( dt > dty ) dt = dty; //dt += 1./dty;
   }
   if( dim_ind[2]!=0 ){
      maxvz = cf + fabs( vz - w[2] );
      dtz = dz/maxvz;
      if( dt > dtz ) dt = dtz; //dt += 1./dtz;
   }
   return( dt ); //return(1./dt);

}

