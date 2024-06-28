
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
   double ux  = prim[UU1];
   double uy  = prim[UU2];
   double uz  = prim[UU3];
   double u0  = sqrt( 1. + ux*ux + uy*uy + uz*uz );
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);
   double rhoh = rho + rhoe + Pp;

   cons[DEN] = rho*u0*dV;
   cons[SS1] = rhoh*u0*ux*dV;
   cons[SS2] = rhoh*u0*uy*dV;
   cons[SS3] = rhoh*u0*uz*dV;
   cons[TAU] = (rhoh*u0*u0 - Pp - rho*u0)*dV;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = cons[DEN]*prim[q];
   }

}

void newt_f( double p , double * f , double * df , double D , double S2 , double E ){
  
   double G = GAMMA_LAW; 
   double v2  = S2/pow( E + p , 2. );
   double rhoe = E*(1.-v2) - D*sqrt(fabs(1.-v2)) - p*v2;
   double Pnew = (G-1.)*rhoe;

   double oe = 1./(sqrt(fabs(1.-v2))*E/D - 1. - p/D*v2/sqrt(fabs(1.-v2)) );
   double c2 = (G-1.)/(1.+oe/G);
   
   *f  = Pnew - p;
   *df = v2*c2 - 1.;

}

void cons2prim( double * cons , double * prim , double * x , double dV ){

   double D   = cons[DEN]/dV;
   double Sx  = cons[SS1]/dV;
   double Sy  = cons[SS2]/dV;
   double Sz  = cons[SS3]/dV;
   double tau = cons[TAU]/dV;
   double E = tau+D;
   double S2 = Sx*Sx + Sy*Sy + Sz*Sz;

   double Pguess = prim[PPP];
   double f,dfdp;
   newt_f( Pguess , &f , &dfdp , D , S2 , E );
   int stop=0;
   while( fabs(f/Pguess) > 1e-8 && stop<10 ){
      Pguess -= f/dfdp;
      newt_f( Pguess , &f , &dfdp , D , S2 , E );
      ++stop;
   }
   
   double Pp  = Pguess;
   double v2  = S2/pow( E + Pp , 2. );
   double rho = D*sqrt(fabs(1.-v2));
   double vx  = Sx/(E+Pp);
   double vy  = Sy/(E+Pp);
   double vz  = Sz/(E+Pp);
   double ux  = vx/sqrt(fabs(1.-v2));
   double uy  = vy/sqrt(fabs(1.-v2));
   double uz  = vz/sqrt(fabs(1.-v2));

   if( rho< RHO_FLOOR ) rho=RHO_FLOOR;
   if( Pp < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = ux;
   prim[UU2] = uy;
   prim[UU3] = uz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DEN];
   }

}

void getUstar( double * prim , double * Ustar , double * x , double Sk , double Ss , double * n ){
   
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double ux  = prim[UU1];
   double uy  = prim[UU2];
   double uz  = prim[UU3];
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);
   double rhoh = rho+rhoe+Pp;

   double un = ux*n[0] + uy*n[1] + uz*n[2];
   double u0  = sqrt( 1. + ux*ux + uy*uy + uz*uz );
   double vn = un/u0;

   double kappa = (Sk-vn)/(Sk-Ss);
   double mn = rhoh*u0*un;
   double E  = rhoh*u0*u0 - Pp;

   double Ps = ( Pp - ( Sk + Ss - vn )*mn + Sk*Ss*E )/( 1.0 - Sk*Ss );

   double alpha1 = ( Ps - Pp )/( Sk - Ss );
   double alpha2 = ( Ss*Ps - vn*Pp )/( Sk - Ss );

   Ustar[DEN] = kappa*rho*u0;
   Ustar[SS1] = kappa*rhoh*u0*ux + alpha1*n[0];
   Ustar[SS2] = kappa*rhoh*u0*uy + alpha1*n[1];
   Ustar[SS3] = kappa*rhoh*u0*uz + alpha1*n[2];
   Ustar[TAU] = kappa*E + alpha2 - kappa*rho*u0;
   
   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DEN];
   }
   
}

void flux( double * prim , double * flux , double * x , double * n ){
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double ux  = prim[UU1];
   double uy  = prim[UU2];
   double uz  = prim[UU3];
   double u0  = sqrt( 1. + ux*ux + uy*uy + uz*uz );
   double un  = ux*n[0] + uy*n[1] + uz*n[2];
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);
   double rhoh = rho+rhoe+Pp;
   
   flux[DEN] = rho*un;
   flux[SS1] = rhoh*ux*un + Pp*n[0];
   flux[SS2] = rhoh*uy*un + Pp*n[1];
   flux[SS3] = rhoh*uz*un + Pp*n[2];
   flux[TAU] = (rhoh*u0-rho)*un;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
      flux[q] = flux[DEN]*prim[q];
   }

}

void source( double * prim , double * cons , double * x , double dVdt ){
   //Silence is golden.
}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss , double * n ){
   
   double gam = GAMMA_LAW;

   double ux1  = prim1[UU1];
   double uy1  = prim1[UU2];
   double uz1  = prim1[UU3];
   double un1  = ux1*n[0]+uy1*n[1]+uz1*n[2];
   double W1   = sqrt( 1. + ux1*ux1 + uy1*uy1 + uz1*uz1 );
   double vn1  = un1/W1;

   double rho  = prim1[RHO];
   double Pp1  = prim1[PPP];
   double rhoe = Pp1/(gam-1.);
   double rhoh1 = rho + rhoe + Pp1;
   double cs = sqrt(fabs( gam*Pp1/rhoh1 ) );
 
   double sigmas = cs*cs/W1/W1/(1.0-cs*cs);
   *Sl = (vn1 - sqrt( sigmas*(1.0 - vn1*vn1 + sigmas) ) )/( 1.0 + sigmas );
   *Sr = (vn1 + sqrt( sigmas*(1.0 - vn1*vn1 + sigmas) ) )/( 1.0 + sigmas );

   double ux2  = prim2[UU1];
   double uy2  = prim2[UU2];
   double uz2  = prim2[UU3];
   double un2  = ux2*n[0]+uy2*n[1]+uz2*n[2];
   double W2   = sqrt( 1. + ux2*ux2 + uy2*uy2 + uz2*uz2 );
   double vn2  = un2/W2;

   rho  = prim2[RHO];
   double Pp2  = prim2[PPP];
   rhoe = Pp2/(gam-1.);
   double rhoh2 = rho + rhoe + Pp2;
   cs = sqrt(fabs( gam*Pp2/rhoh2 ) );
 
   sigmas = cs*cs/W2/W2/(1.0-cs*cs);
   double sl = (vn2 - sqrt( sigmas*(1.0 - vn2*vn2 + sigmas) ) )/( 1.0 + sigmas );
   double sr = (vn2 + sqrt( sigmas*(1.0 - vn2*vn2 + sigmas) ) )/( 1.0 + sigmas );

   if( *Sl > sl ) *Sl = sl;
   if( *Sr < sr ) *Sr = sr;

   double El = rhoh1*W1*W1 - Pp1;
   double Er = rhoh2*W2*W2 - Pp2;
   double Ml = rhoh1*W1*un1;
   double Mr = rhoh2*W2*un2;
   double Fl = rhoh1*un1*un1 + Pp1;
   double Fr = rhoh2*un2*un2 + Pp2;

   double aL = *Sr;
   double aR = -*Sl;

   double FE = ( aL*Ml + aR*Mr + aL*aR*( El - Er ) )/(aL + aR); 
   double UE = ( aR*El + aL*Er + Ml - Mr )/(aL + aR); 
   double FM = ( aL*Fl + aR*Fr + aL*aR*( Ml - Mr ) )/(aL + aR); 
   double UM = ( aR*Ml + aL*Mr + Fl - Fr )/(aL + aR); 

   if( fabs(FE*UM/pow(UE+FM,2.0)) < 1e-10 ){
      *Ss = UM/(UE + FM); 
   }else{
      *Ss = ( UE + FM - sqrt(fabs((UE + FM)*(UE + FM) - 4.0*FE*UM )) )/( 2.0*FE );
   }   
   
}

double get_maxv( double * prim , double w , int dim ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double ux  = prim[UU1];
   double uy  = prim[UU2];
   double uz  = prim[UU3];
   double W2  = 1. + ux*ux + uy*uy + uz*uz;
   double vx  = ux/sqrt(W2);
   double vy  = uy/sqrt(W2);
   double vz  = uz/sqrt(W2);
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);
   double rhoh = rho + rhoe + Pp;

   double cs2 = gam*Pp/rhoh;
   double sigmas = cs2/W2/(1.0-cs2);

   double vn = vx;
   if( dim==1 ) vn = vy;
   if( dim==2 ) vn = vz;

   double sl = (vn - sqrt( sigmas*(1.0 - vn*vn + sigmas) ) )/( 1.0 + sigmas ) - w;
   double sr = (vn - sqrt( sigmas*(1.0 - vn*vn + sigmas) ) )/( 1.0 + sigmas ) - w;

   double maxv = fabs(sl);
   if( maxv < fabs(sr) ) maxv = fabs(sr);

   return(maxv);
}

double mindt( double * prim , double * w , int * dim_ind , double dx , double dy , double dz ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double ux  = prim[UU1];
   double uy  = prim[UU2];
   double uz  = prim[UU3];
   double W2  = 1. + ux*ux + uy*uy + uz*uz;
   double vx  = ux/sqrt(W2);
   double vy  = uy/sqrt(W2);
   double vz  = uz/sqrt(W2);
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);
   double rhoh = rho + rhoe + Pp;

   double cs = sqrt( gam*Pp/rhoh );

   double maxvx = get_maxv(prim,w[0],0);
   double maxvy = get_maxv(prim,w[1],1);
   double maxvz = get_maxv(prim,w[2],2);
   double dtx = 0.;
   double dty = 0.;
   double dtz = 0.;
   double dt = 0.;

   if( dim_ind[0]!=0 ){
      maxvx = cs + fabs( vx - w[0] );
      dtx = dx/maxvx;
      dt = 1./dtx; //dt = dtx;
   }
   if( dim_ind[1]!=0 ){
      maxvy = cs + fabs( vy - w[1] );
      dty = dy/maxvy;
      dt += 1./dty; //if( dt > dty ) dt = dty;
   }
   if( dim_ind[2]!=0 ){
      maxvz = cs + fabs( vz - w[2] );
      dtz = dz/maxvz;
      dt += 1./dtz; //if( dt > dtz ) dt = dtz;
   }
   return(1./dt);//return( dt );

}
