
#include "../defs.h"

static double x_cen = 0.0;
static double y_cen = 0.0;
static double z_cen = 0.0;
static double shock_pos = 0.9;

void setMeshMotionParams( struct domain * theDomain ){

   x_cen = theDomain->theParList.MM_x0;
   y_cen = theDomain->theParList.MM_y0;
   z_cen = theDomain->theParList.MM_z0;
   shock_pos = 0.60;
    
}

double time_to_boundary( double L , double x , double v ){
   if( v>0. ) return ( fabs((L-x)/v) );
   if( v<0. ) return ( fabs(-1.*x/v) );
   if( v==0. ) return ( 1e10 );
}

double get_r_fast_zero( int dim ){
   if( dim==0 ) return x_cen;
   if( dim==1 ) return y_cen;
   if( dim==2 ) return z_cen;
}

double get_L( double dx , double dy , double dz , int Nx , int Ny , int Nz , int dim ){
   if( dim==0 ) return( (double)Nx * dx );
   if( dim==1 ) return( (double)Ny * dy );
   if( dim==2 ) return( (double)Nz * dz );
}

void set_W( struct domain * theDomain , int reset ){

   int i,j,k,ijk,dim;
   double v_fast , r_fast , r_fast_0; 
   double L,t_local,t_min = 1e30; 

   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   int Num_x = theDomain->theParList.Num_x;
   int Num_y = theDomain->theParList.Num_y;
   int Num_z = theDomain->theParList.Num_z;
   int xsize = theDomain->dim_size[0];
   int ysize = theDomain->dim_size[1];
   int zsize = theDomain->dim_size[2];
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
 
   int i0=0 ; int j0=0 ; int k0=0; 
   int i1=1 ; int j1=1 ; int k1=1;
   int maxdim  = 0;
   if( theDomain->theParList.Num_x != 1 ){
      i0 = Ng; i1 = Nx+Ng;
      maxdim += 1;
   }
   if( theDomain->theParList.Num_y == 1 ){
      j0 = Ng; j1 = Ny+Ng;
      maxdim += 1;
   }
   if( theDomain->theParList.Num_z == 1 ){
      k0 = Ng; k1 = Nz+Ng;
      maxdim += 1;
   }

   
   for( k=k0 ; k<k1 ; ++k ){
      for( j=j0 ; j<j1 ; ++j ){
         for( i=i0 ; i<i1 ; ++i ){
            ijk  = i;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
            struct cell * c = theDomain->theCells+ijk;
            for( dim=0 ; dim<maxdim ; ++dim  ){
               if( dim==0 ) L = (double)Nx * (double)xsize * dx;
               if( dim==1 ) L = (double)Ny * (double)ysize * dy;
               if( dim==2 ) L = (double)Nz * (double)zsize * dz;
               t_local = time_to_boundary( L , c->xi[dim] , c->prim[UU1+dim] );
               if( t_local<t_min ){
                  t_min    = t_local;
                  v_fast   = c->prim[UU1+dim];
                  r_fast   = c->xi[dim];
                  r_fast_0 = get_r_fast_zero(dim);
                  L        = get_L( dx,dy,dz,Num_x,Num_y,Num_z,dim );
               }   
            }
         }
      }
   }

   double W_local = 0.;
   //printf("Lshock = %e, rfast = %e, vfast = %e\n",L*shock_pos, r_fast, v_fast);
   if( r_fast>= shock_pos*L )
      W_local = fabs( v_fast/(shock_pos*L - r_fast_0) );
   MPI_Allreduce( MPI_IN_PLACE , &W_local , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );
   theDomain->W = W_local;

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
