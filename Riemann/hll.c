
#include "../defs.h"

void prim2cons( double * , double * , double * , double );
void flux( double * , double * , double * , double * );
void vel( double * , double * , double * , double * , double * , double * );

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
         *C_U = 1.;
      }
      if( no_of_dims==2 ){ 
         *C_F = (1. + W*dt) /(1. + W*dt*2.);
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
   double Fr[NUM_Q];
   double Ul[NUM_Q];
   double Ur[NUM_Q];

   double Flux[NUM_Q];

   double wn = get_wn( xl , dx , dy , dz , W , theDIM );

   double C_F, C_U;
   get_flux_coefficients( no_of_dims , first_step , last_step , W , dt , &C_F , &C_U );

   if( wn < Sl ){
      flux( primL , Fl , xl , n );
      prim2cons( primL , Ul , xl , 1.0 );
      
      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = C_F*Fl[q] - C_U*wn*Ul[q];
      }
   }else if( wn > Sr ){
      flux( primR , Fr , xr , n );
      prim2cons( primR , Ur , xr , 1.0 );
      
      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = C_F*Fr[q] - C_U*wn*Ur[q];
      }
   }else{
      
      double Fstar;
      double Ustar;

      double aL =  Sr;
      double aR = -Sl;
      
      prim2cons( primL , Ul , xl , 1.0 );
      prim2cons( primR , Ur , xr , 1.0 );
      flux( primL , Fl , xl , n );
      flux( primR , Fr , xr , n );
      
      for( q=0 ; q<NUM_Q ; ++q ){

         Fstar = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
         Ustar = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );

         Flux[q] = C_F*Fstar - C_U*wn*Ustar;
      }

   }


   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];
   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dt*dA;
      cR->cons[q] += Flux[q]*dt*dA;
   }
 
   //printf("%e %e\n", xl[0]+dx/2-.5,Flux[SS1]); //for file

   //printf("x = %e, Delta_F = %e\n", xl[0]+dx/2-.5, Flux[DEN]); //mass

   //printf("x = %e, Delta_F = %e\n", xl[0]+dx/2-.5, -Flux[SS1] + (xl[0]+dx/2-.5)*(xl[0]+dx/2-.5) + 30. - (xl[0]+dx/2) );
   //printf("x = %e, Delta_F = %e\n", xl[0]+dx/2-.5,-Flux[SS1]+primL[PPP]+(xl[0]+dx/2-.5)*(xl[0]+dx/2-.5)); //momentum

   //printf("x = %e, Delta_F = %e\n", xl[0]+dx/2-.5,-Flux[TAU]+0.5*pow((xl[0]+dx/2-.5),3.) + 2.5*(30.-(xl[0]+dx/2))*(xl[0]+dx/2-.5));  //energy
}
