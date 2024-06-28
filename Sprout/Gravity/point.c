
#include "../defs.h"

static double M         = 0.0;
static double G         = 0.0;
static int G_mode       = 0;
static double eps       = 0.0;
static double v_star    = 0.0;
static double t0        = 0.0;
static double dims      = 0.0;
static double r_ref[3]  = {0.0};

void setGravityParams( struct domain * theDomain ){
   M      = theDomain->theParList.Central_Mass;
   G      = theDomain->theParList.Grav_G;
   G_mode = theDomain->theParList.Gravity_Switch;
   eps    = theDomain->dx*2.;
   v_star = theDomain->theParList.Ly*.5/(theDomain->t_fin - theDomain->t_init);
   t0     = theDomain->t_init;
   r_ref[0]  = theDomain->theParList.Lx/2.;
   r_ref[1]  = theDomain->theParList.Ly/4.;
   r_ref[2]  = theDomain->theParList.Lz/2.;
   if( theDomain->theParList.Num_x != 1 ) dims += 1.;
   if( theDomain->theParList.Num_y != 1 ) dims += 1.;
   if( theDomain->theParList.Num_z != 1 ) dims += 1.;
}

double get_pot( double x_cell , double y_cell , double z_cell , double * r_0 ){
   double x   = x_cell - r_0[0];
   double y   = y_cell - r_0[1];
   double z   = z_cell - r_0[2];
   double r   = sqrt( x*x + y*y + z*z + eps*eps );
   double phi = -G*M/r;
   return( phi );
}


/*//////////////////////////////////////////////////////
 **This is how the phi indices are arranged (C=center)**


        2***********************3          y
      * *                    *  *          ^
    *   *                  *    *          *
 6***********************7      *          *
 *      *          C     *      *          *
 *      *                *      *          * * * * * >x
 *      0***********************1        *
 *     *                 *     *       *
 *   *                   *   *       *
 * *                     * *       z
 4***********************5

/////////////////////////////////////////////////////*/


void get_g( double * g , double * r , double * r_0 , double dx , double dy , double dz ){
   
   g[0] = 0.; g[1] = 0.; g[2] = 0.;
   double phi[8] = {0.0};
   if( dims==1. ){
      phi[0] = get_pot( r[0]-dx/2. , r[1] , r[2] , r_0 );
      phi[1] = get_pot( r[0]+dx/2. , r[1] , r[2] , r_0 );
      g[0] = ( phi[0] - phi[1] )/dx;
   } else if( dims==2. ){
      phi[0] = get_pot( r[0]-dx/2. , r[1]-dy/2. , r[2] , r_0 );
      phi[1] = get_pot( r[0]+dx/2. , r[1]-dy/2. , r[2] , r_0 );
      phi[2] = get_pot( r[0]-dx/2. , r[1]+dy/2. , r[2] , r_0 );
      phi[3] = get_pot( r[0]+dx/2. , r[1]+dy/2. , r[2] , r_0 );
      g[0] = ( phi[0] + phi[2] - phi[1] - phi[3] )/dx/2.;
      g[1] = ( phi[0] + phi[1] - phi[2] - phi[3] )/dy/2.;
   } else if( dims==3. ){
      phi[0] = get_pot( r[0]-dx/2. , r[1]-dy/2. , r[2]-dz/2. , r_0 );
      phi[1] = get_pot( r[0]+dx/2. , r[1]-dy/2. , r[2]-dz/2. , r_0 );
      phi[2] = get_pot( r[0]-dx/2. , r[1]+dy/2. , r[2]-dz/2. , r_0 );
      phi[3] = get_pot( r[0]+dx/2. , r[1]+dy/2. , r[2]-dz/2. , r_0 );
      phi[4] = get_pot( r[0]-dx/2. , r[1]-dy/2. , r[2]+dz/2. , r_0 );
      phi[5] = get_pot( r[0]+dx/2. , r[1]-dy/2. , r[2]+dz/2. , r_0 );
      phi[6] = get_pot( r[0]-dx/2. , r[1]+dy/2. , r[2]+dz/2. , r_0 );
      phi[7] = get_pot( r[0]+dx/2. , r[1]+dy/2. , r[2]+dz/2. , r_0 );
      g[0] = ( phi[0] + phi[2] + phi[4] + phi[6] - phi[1] - phi[3] - phi[5] - phi[7] )/dx/4.;
      g[1] = ( phi[0] + phi[1] + phi[4] + phi[5] - phi[2] - phi[3] - phi[6] - phi[7] )/dy/4.;
      g[2] = ( phi[0] + phi[1] + phi[2] + phi[3] - phi[4] - phi[5] - phi[6] - phi[7] )/dz/4.;
   }
   
}

void grav_src( double * prim , double * cons , double * r , double dx , double dy , double dz , double dt , double t ){

   if( G_mode ){
      double g_acc[3] , r_0[3];
      r_0[0] = r_ref[0];
      r_0[1] = r_ref[1] + v_star*(t-t0);
      r_0[2] = r_ref[2];

      get_g( g_acc , r , r_0 , dx , dy , dz );

      double rho  = prim[RHO];
      double dVdt = dx*dy*dz*dt;

      cons[TAU] += dVdt * rho * ( g_acc[0]*prim[UU1] + g_acc[1]*prim[UU2] + g_acc[2]*prim[UU3] );
      cons[SS1] += dVdt * rho * g_acc[0];
      cons[SS2] += dVdt * rho * g_acc[1];
      cons[SS3] += dVdt * rho * g_acc[2];
   }

}
