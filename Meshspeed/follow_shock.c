
#include "../defs.h"

static int D = 0;
static double x_cen = 0.0;
static double y_cen = 0.0;
static double z_cen = 0.0;
static double t_min  = 0.0;
static double eta_on = 0.0;
static double eps = 0.0;
static double shock_pos = 0.0;



void setMeshMotionParams( struct domain * theDomain ){

   x_cen = theDomain->theParList.MM_x0 * theDomain->theParList.Lx;
   y_cen = theDomain->theParList.MM_y0 * theDomain->theParList.Ly;
   z_cen = theDomain->theParList.MM_z0 * theDomain->theParList.Lz;
   if( theDomain->theParList.Num_x!=1 ) D += 1;
   if( theDomain->theParList.Num_y!=1 ) D += 1;
   if( theDomain->theParList.Num_z!=1 ) D += 1;
   t_min  = theDomain->t_init;
   eta_on = theDomain->theParList.eta_on; 
   shock_pos = 0.85;
   eps = 0.25;
    
}

double marker( struct cell * c , int dim , double x0 ){
   //return ( c->prim[UU1+dim]*c->prim[RHO] );
   //return ( c->prim[UU1+dim]*c->prim[PPP]*(c->xi[dim]-x0) );
   return ( c->prim[RHO] );
   //return ( (1.-c->prim[XXX])*c->prim[RHO] );
}


double clip_v( double L , double L0 , double x ){
   double L_SP = shock_pos*L + (1.-shock_pos)*L0;
   return (1.+log(x/L_SP))/log(L/L_SP);
}

void clip_r( double L , double L0 , double * x ){
   if( *x>L0+shock_pos*L ) *x = L0+shock_pos*L;
   //if( *x>L0+shock_pos*L/2. ) *x = L0+shock_pos*L/2.;
   //if( *x<L0-shock_pos*L/2. ) *x = L0-shock_pos*L/2.;
}


void set_W( struct domain * theDomain , int reset ){

   double t = theDomain->t;

   if( t<t_min*eta_on ){
      theDomain->W = 0.;
   }else{
      int i,j,k,ijk,dim,fastdim;
      double r_fast,v_fast,pv;
      double pv_fast  = 0.;
      double r_from_c = 0.;

      int Nx = theDomain->Nx;
      int Ny = theDomain->Ny;
      int Nz = theDomain->Nz;
      int Ng = theDomain->Ng;
      double dx = theDomain->dx;
      double dy = theDomain->dy;
      double dz = theDomain->dz;
      int Num_x = theDomain->theParList.Num_x;
      int Num_y = theDomain->theParList.Num_y;
      int Num_z = theDomain->theParList.Num_z;
      

      double Ls[3],L0s[3],Cs[3];
      Ls[0]  = (double)Num_x * dx;
      Ls[1]  = (double)Num_y * dy;
      Ls[2]  = (double)Num_z * dz;
      L0s[0] = theDomain->theParList.Lx;
      L0s[1] = theDomain->theParList.Ly;
      L0s[2] = theDomain->theParList.Lz;
      Cs[0]  = x_cen;
      Cs[1]  = y_cen;
      Cs[2]  = z_cen;
 
      int i0=0 ; int j0=0 ; int k0=0; 
      int i1=1 ; int j1=1 ; int k1=1;

      if( theDomain->theParList.Num_x != 1 ){
         i0 = Ng; 
         i1 = Nx+Ng;
      }
      if( theDomain->theParList.Num_y != 1 ){
         j0 = Ng; 
         j1 = Ny+Ng;
      }
      if( theDomain->theParList.Num_z != 1 ){
         k0 = Ng; 
         k1 = Nz+Ng;
      }

   
      for( k=k0 ; k<k1 ; ++k ){
         for( j=j0 ; j<j1 ; ++j ){
            for( i=i0 ; i<i1 ; ++i ){
               ijk  = i;
               if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * c = theDomain->theCells+ijk;
               for( dim=0 ; dim<1 ; ++dim ){
                  pv = marker( c, dim, Cs[dim] );
                  if( pv_fast<pv ){
                     pv_fast  = pv;
                     //v_fast   = c->prim[UU1+dim];
                     //r_fast   = c->xi[dim];
                     fastdim  = dim;
                  }
               }
               
            }
         }
      }

      MPI_Allreduce( MPI_IN_PLACE , &pv_fast , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );

      //if(fabs(t/1.6e-4-1.)<.2) printf("CHECK: t=%.4e, marker=%.4e, r=%.4e\n",t,pv_fast,r_fast);    
 
      
      for( k=k0 ; k<k1 ; ++k ){
         for( j=j0 ; j<j1 ; ++j ){
            for( i=i0 ; i<i1 ; ++i ){
               ijk  = i;
               if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * c = theDomain->theCells+ijk;
               pv = marker( c, fastdim, Cs[fastdim] );
               if( pv>eps*pv_fast && fabs(c->xi[fastdim]-Cs[fastdim])>r_from_c ){ //
                  r_from_c = fabs(c->xi[fastdim]-Cs[fastdim]);
                  r_fast   = c->xi[fastdim];
                  v_fast   = c->prim[UU1+fastdim];
               }
               /*for( dim=0 ; dim<D ; ++dim ){
                  pv = fabs(c->prim[UU1+dim]*c->prim[PPP]);
                  if( pv>0.01*pv_fast && fabs(c->xi[dim]-Cs[dim])>r_from_c ){
                     r_from_c = fabs(c->xi[dim]-Cs[dim]);
                     r_fast   = c->xi[fastdim];
                  }
               }*/
            }
         }
      }

            
      MPI_Allreduce( MPI_IN_PLACE , &v_fast , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );
      MPI_Allreduce( MPI_IN_PLACE , &r_fast , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );

      //if(fabs(t/1.6e-4-1.)<.2) printf("CHECK: t=%.4e, marker=%.4e, r=%.4e\n",t,v_fast,r_fast);



      double W_local = 0.;
      //v_fast *= clip_v( Ls[fastdim] , Cs[fastdim] , r_fast );
      //clip_r( Ls[fastdim] , Cs[fastdim] , &r_fast );
      W_local = fabs( v_fast/(r_fast - Cs[fastdim]) );
      //MPI_Allreduce( MPI_IN_PLACE , &W_local , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );
      theDomain->W = W_local;
   }


}

double get_wn( double * xl , double dx , double dy , double dz , double W , int theDIM ){

   double wn = W;
   if( theDIM==0 ) wn *= ( xl[0] + dx/2. - x_cen );
   if( theDIM==1 ) wn *= ( xl[1] + dy/2. - y_cen );
   if( theDIM==2 ) wn *= ( xl[2] + dz/2. - z_cen );
   return(wn);

}

void calc_dxs( struct domain * theDomain , double dt ){

   double W  = theDomain->W;

   if( theDomain->theParList.Num_x != 1 ){
      theDomain->dx *= ( 1. + W*dt );
      theDomain->theParList.Lx *= ( 1. + W*dt );
   }
   if( theDomain->theParList.Num_y != 1 ){
      theDomain->dy *= ( 1. + W*dt );
      theDomain->theParList.Ly *= ( 1. + W*dt );
   }
   if( theDomain->theParList.Num_z != 1 ){
      theDomain->dz *= ( 1. + W*dt );
      theDomain->theParList.Lz *= ( 1. + W*dt );
   }

}

void regrid( struct domain * theDomain , double dt ){

   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double W  = theDomain->W;
 
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
            struct cell * c = theDomain->theCells+ijk;
            if( theDomain->theParList.Num_x != 1 )
               c->xi[0] = c->xi[0] * ( 1. + W*dt ) - x_cen*W*dt;
            if( theDomain->theParList.Num_y != 1 )
               c->xi[1] = c->xi[1] * ( 1. + W*dt ) - y_cen*W*dt;
            if( theDomain->theParList.Num_z != 1 )
               c->xi[2] = c->xi[2] * ( 1. + W*dt ) - z_cen*W*dt;
         }
      }
   }

}
