#include "../defs.h"

void prim2cons( double * , double * , double * , double );
void flux( double * , double * , double * , double * );
void vel( double * , double * , double * , double * , double * , double * );
void getUstar( double * , double * , double * , double , double , double * );

double get_wn( double * , double , double , double , double , int );

void get_flux_coefficients( int no_of_dims , int first_step , int last_step , double W , double dt , double * C_F , double * C_U ){

   if( first_step==1 ){
      if( no_of_dims==1 ){
         *C_F = 1.;
         *C_U = 1.;
      }
      if( no_of_dims==2 ){
         *C_F = 1. + W*dt/2.;
         *C_U = 1. + W*dt/2.;
      }
      if( no_of_dims==3 ){
         *C_F = 1. + W*dt + W*W*dt*dt/3.;
         *C_U = 1. + W*dt + W*W*dt*dt/3.;
      }   
   }

   if( first_step==0 ){
      if( no_of_dims==1 ){
         *C_F = 1.;
         *C_U = 1./(1. + W*dt*2.);
      }
      if( no_of_dims==2 ){ 
         *C_F = (1. + W*dt)/(1. + W*dt*2.);
         *C_U = (1. + W*dt)/(1. + W*dt*2.)/(1. + W*dt*2.);
      }
      if( no_of_dims==3 ){ 
         *C_F = (1. + W*dt*2. + W*W*dt*dt*4./3.)/(1. + W*dt*2.)/(1. + W*dt*2.);
         *C_U = (1. + W*dt*2. + W*W*dt*dt*4./3.)/(1. + W*dt*2.)/(1. + W*dt*2.)/(1. + W*dt*2.);
      } 
   }
   
}

void riemann1D( struct cell * cL , struct cell * cR , double dx , double dy , double dz , double dt , double W , int no_of_dims , int theDIM , int first_step , int last_step ){
   
   double primL[NUM_Q];
   double primR[NUM_Q];

   double n[3] = {0.0};
   n[theDIM] = 1.0;

   double * xl , * xr;
   xl = cL->xi;
   xr = cR->xi;

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + .5 * (cL->gradx[q]*dx*n[0] + cL->grady[q]*dy*n[1] + cL->gradz[q]*dz*n[2]);
      primR[q] = cR->prim[q] - .5 * (cR->gradx[q]*dx*n[0] + cR->grady[q]*dy*n[1] + cR->gradz[q]*dz*n[2]);
   }

   double Sl,Sr,Ss;

   vel( primL , primR , &Sl , &Sr , &Ss , n );

   double Fl[NUM_Q];
   double Ul[NUM_Q];
   double Fr[NUM_Q];
   double Ur[NUM_Q];

   double Flux[NUM_Q];
   
   double wn = get_wn( xl , dx , dy , dz , W , theDIM );

   double C_F, C_U;
   get_flux_coefficients( no_of_dims , first_step , last_step , W , dt , &C_F , &C_U );

   double Mach_L, Mach_R, Mach_local, Mach_limit;
   Mach_L = fabs( (cL->prim[UU1+theDIM]-wn)/sqrt(5./3.*cL->prim[PPP]/cL->prim[RHO]) );
   Mach_R = fabs( (cR->prim[UU1+theDIM]-wn)/sqrt(5./3.*cR->prim[PPP]/cR->prim[RHO]) );
   Mach_local = Mach_L;
   if(Mach_R>Mach_L) Mach_local = Mach_R;
   Mach_limit = 0.5;

   double phi = 1.0;
   if( Mach_local<Mach_limit ) phi = sin(Mach_local/Mach_limit*M_PI/2.);

   double Sl_LM = Sl*phi;
   double Sr_LM = Sr*phi;


   if( wn < Sl_LM ){
      flux( primL , Fl , xl , n );
      prim2cons( primL , Ul , xl , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = C_F*Fl[q] - C_U*wn*Ul[q];
      }
   }else if( wn > Sr_LM ){
      flux( primR , Fr , xr , n );
      prim2cons( primR , Ur , xr , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = C_F*Fr[q] - C_U*wn*Ur[q];
      }
   }else{
      double UstarL[NUM_Q];
      double UstarR[NUM_Q];
      prim2cons( primL , Ul , xl , 1.0 );
      prim2cons( primR , Ur , xr , 1.0 );
      getUstar( primL , UstarL , xl , Sl , Ss , n );
      getUstar( primR , UstarR , xr , Sr , Ss , n );
      flux( primL , Fl , xl , n );
      flux( primR , Fr , xr , n );

      if( wn < Ss ){
         for( q=0 ; q<NUM_Q ; ++q )
            Flux[q] = C_F*0.5*( Fl[q] + Fr[q] + Sl_LM*(UstarL[q]-Ul[q]) + Sr_LM*(UstarR[q]-Ur[q]) + Ss*(UstarL[q]-UstarR[q]) ) - C_U*wn*UstarL[q];
            //Flux[q] = C_F*( Fl[q] + Sl*( UstarL[q] - Ul[q] ) ) - C_U*wn*UstarL[q];
      }else{
         for( q=0 ; q<NUM_Q ; ++q )
            Flux[q] = C_F*0.5*( Fl[q] + Fr[q] + Sl_LM*(UstarL[q]-Ul[q]) + Sr_LM*(UstarR[q]-Ur[q]) - Ss*(UstarL[q]-UstarR[q]) ) - C_U*wn*UstarR[q];
            //Flux[q] = C_F*( Fr[q] + Sr*( UstarR[q] - Ur[q] ) ) - C_U*wn*UstarR[q];
      }   

   }

   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];
   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dt*dA;
      cR->cons[q] += Flux[q]*dt*dA;
   }

}
