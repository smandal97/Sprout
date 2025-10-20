#include "../defs.h"

void prim2cons( double * , double * , double * , double );
void flux( double * , double * , double * , double * );
void vel_and_Pts( double * , double * , double * , double * , double * );
void get_single_star( double * , double * , double , double , double , double * );
void get_double_star( double * , double * , double * , double * , double * , int );


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

void riemann1D( struct cell * cL , struct cell * cR , double dx , double dy , double dz , double dt , double W , int no_of_dims , int theDIM , int first_step , int last_step , int CT ){
   
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

   double wn = get_wn( xl , dx , dy , dz , W , theDIM );
   double C_F, C_U;
   get_flux_coefficients( no_of_dims , first_step , last_step , W , dt , &C_F , &C_U );

   double vel[5];
   double Pt_star;
   vel_and_Pts( primL , primR , vel , &Pt_star , n );

   double F[NUM_Q];
   double U[NUM_Q];
   double Flux[NUM_Q];



   if( wn < vel[2] ){
      //Left Side
      double UL[NUM_Q];
      flux( primL , F , xl , n );
      prim2cons( primL , UL , xl , 1.0 );
      if( wn < vel[0] ){
         //Upwind Left
         prim2cons( primL , U , xl , 1.0 );
      }else if( wn < vel[1] ){
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
      flux( primR , F , xr , n );
      prim2cons( primR , UR , xr , 1.0 );
      if( wn > vel[4] ){
         //Upwind Right
         prim2cons( primR , U , xr , 1.0 );
      }else if( wn > vel[3] ){
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

   
   for( q=0 ; q<NUM_Q ; ++q )
      Flux[q] = C_F*F[q] - C_U*wn*U[q];
   


   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];
   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dt*dA;
      cR->cons[q] += Flux[q]*dt*dA;
      //if(q==7 && Flux[q]>1e-10) printf("B Riemann flux earlier = %e\n", Flux[q]);
   }

   //if(Flux[7]>1e-10) printf("Riemann flux test HLLD: %e\n", Flux[7]);


   if( NUM_M!=0 && CT==1 ){
      int DIM_p1 = (theDIM+1)%3;
      int DIM_p2 = (theDIM+2)%3;
      cR->Phi_R[theDIM*(NUM_M-1)+0] =  Flux[DIM_p2+NUM_C+NUM_N];
      cR->Phi_R[theDIM*(NUM_M-1)+1] = -Flux[DIM_p1+NUM_C+NUM_N];
      //if(theDIM==1 && fabs(Flux[BB1])>1e-10) printf("N1, N2 = %i, %i, x = %e, B Riemann flux later= %e\n", BB1, DIM_p2+NUM_C+NUM_N, cL->xi[0], Flux[BB1]);

      cR->CMode[theDIM] = Flux[RHO];
   }

}
