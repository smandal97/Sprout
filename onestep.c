
#include "defs.h"

void exchangeData( struct domain * );
void boundary( struct domain * );

void adjust_rk_cons( struct domain * theDomain , double RK ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;

   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;
   
   int i,j,k,ijk,q;
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            ijk  = i;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
            struct cell * c = theCells+ijk;
            for( q=0 ; q<NUM_Q ; ++q )
               c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
         }
      }   
   }

}

void space_recon1D( struct domain * , int );
void riemann1D( struct cell * , struct cell * , double , double , double , double , double , int , int , int , int );

void flux1D( struct domain * theDomain , int theDIM , double dt , int first_step , int last_step ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   double W = theDomain->W;
   
   int n[3]  = {0};
   n[theDIM] = 1;

   int no_of_dims = 3;
   int i_index = Nx+2*Ng-1;
   int j_index = Ny+2*Ng-1;
   int k_index = Nz+2*Ng-1;
   if( theDomain->theParList.Num_x == 1 ){
      i_index = 1;
      no_of_dims -= 1;
   }
   if( theDomain->theParList.Num_y == 1 ){
      j_index = 1;
      no_of_dims -= 1;
   }
   if( theDomain->theParList.Num_z == 1 ){
      k_index = 1;
      no_of_dims -= 1;
   }


   space_recon1D( theDomain , theDIM );

   int i,j,k,ijk_L,ijk_R;
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            ijk_L = i;
            ijk_R = i+n[0];
            if( theDomain->theParList.Num_y != 1 ){
               ijk_L += (Nx+2*Ng)*j;
               ijk_R += (Nx+2*Ng)*(j+n[1]);
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk_L += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_R += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]);
            }
            struct cell * cL = theCells+ijk_L;
            struct cell * cR = theCells+ijk_R;
            riemann1D( cL , cR , dx , dy , dz , dt , W , no_of_dims , theDIM , first_step , last_step );
         }
      }   
   }

}

void add_flux( struct domain * theDomain , double dt , int first_step , int last_step ){

   if( theDomain->theParList.Num_x != 1 ) flux1D( theDomain , 0 , dt , first_step , last_step );
   if( theDomain->theParList.Num_y != 1 ) flux1D( theDomain , 1 , dt , first_step , last_step );
   if( theDomain->theParList.Num_z != 1 ) flux1D( theDomain , 2 , dt , first_step , last_step );

}

void source( double * , double * , double * , double );
void grav_src( double * , double * , double * , double , double );
void nozz_src( double * , double * , double * , double , double , double , double , double );

double get_source_coefficient( int no_of_dims , int first_step , int last_step , double W , double dt ){

   double C_S;
   if( first_step==1 ){
      if( no_of_dims==1 ) C_S = 1. + W*dt/2.;
      if( no_of_dims==2 ) C_S = 1. + W*dt + W*W*dt*dt/3.;
      if( no_of_dims==3 ) C_S = 1. + 3./2.*W*dt + W*W*dt*dt + W*W*W*dt*dt*dt/4.; 
   }

   if( first_step==0 ){
      if( no_of_dims==1 ) C_S = (1. + W*dt)/(1. + W*dt*2.);
      if( no_of_dims==2 ) C_S = (1. + W*dt*2. + W*W*dt*dt*4./3.)/(1. + W*dt*2.)/(1. + W*dt*2.);
      if( no_of_dims==3 ) C_S = (1. + 3.*W*dt + 4.*W*W*dt*dt + W*W*W*dt*dt*dt*2.)/(1. + W*dt*2.)/(1. + W*dt*2.)/(1. + W*dt*2.); 
   }
   
   return C_S;

}

void add_source( struct domain * theDomain , double dt , int first_step , int last_step ){

   struct cell * theCells = theDomain->theCells;
   int grav_switch = theDomain->theParList.Gravity_Switch;
   int nozz_switch = theDomain->theParList.Nozzle_Switch;
   
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   double W = theDomain->W;
   double t = theDomain->t;
   
   int no_of_dims = 3;
   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ){
      no_of_dims -= 1;
      i_index = 1;
   }
   if( theDomain->theParList.Num_y == 1 ){ 
      no_of_dims -= 1;
      j_index = 1;
   }
   if( theDomain->theParList.Num_z == 1 ){
      no_of_dims -= 1;
      k_index = 1;
   }

   double C_S;
   C_S = get_source_coefficient( no_of_dims , first_step , last_step , W , dt );
   
   int i,j,k,ijk;
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            ijk  = i;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
            struct cell * c = theCells+ijk;
            source( c->prim , c->cons , c->xi , dx*dy*dz*dt*C_S );
            if( grav_switch ) grav_src( c->prim , c->cons , c->xi , dx*dy*dz*dt*C_S , t );
            if( nozz_switch ) nozz_src( c->prim , c->cons , c->xi , dx, dy, dz, dt*C_S , t );
            
         }
      }
   }

}

void set_W( struct domain * , int );
void calc_dxs( struct domain * , double );
void regrid( struct domain * , double );

void cons2prim( double * , double * , double * , double );

void calc_prim( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   
   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;
   
   int i,j,k,ijk;
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            ijk  = i;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
            struct cell * c = theCells+ijk;
            cons2prim( c->cons , c->prim , c->xi , dx*dy*dz );
            //printf("x = %e, check = %e\n", c->xi[0], -c->cons[TAU]/0.01 + 0.5/1.000660897*pow( (c->xi[0]-.5)/1.000660897,2. ) + 1.5 *( c->prim[PPP] ) );
         }
      }
   }

}

void onestep( struct domain * theDomain , double RK , double dt , int first_step , int last_step ){
   
   if( first_step ) set_W( theDomain , 1 );
   
   adjust_rk_cons( theDomain , RK );
   
   add_flux( theDomain , dt , first_step , last_step );
   add_source( theDomain , dt , first_step , last_step );

   if( first_step ) calc_dxs( theDomain , dt );
   calc_prim( theDomain );

   if( first_step ) regrid( theDomain , dt );
}
