
#include "defs.h"

void setupDomain( struct domain * theDomain ){
     
   theDomain->t      = theDomain->theParList.t_min;
   theDomain->t_init = theDomain->theParList.t_min;
   theDomain->t_fin  = theDomain->theParList.t_max;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;

   theDomain->final_step = 0;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;

   theDomain->count_steps = 0;
   
}

void restart( struct domain * );
void exchange1D( struct domain * , int );
void initial( double * , double * , double );
void prim2cons( double * , double * , double * , double );

int getN0( int , int , int );

void setupCells( struct domain * theDomain ){

   int restart_flag = theDomain->theParList.restart_flag;
   if( restart_flag ) restart( theDomain );
      
   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   double t  = theDomain->t;


   int i,j,k,q,ijk,ijk_0,ijk_m1;
   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Ny ; ++j ){
         for( i=0 ; i<Nx ; ++i ){
            ijk  = i+Ng;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*(j+Ng);
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*(k+Ng);
            struct cell * c = theCells+ijk;
            //set prims
            if(!restart_flag) initial( c->prim , c->xi , t );
            prim2cons( c->prim , c->cons , c->xi , dx*dy*dz );
            //set gradients
            for( q=0 ; q<NUM_Q ; ++q ){
               c->gradx[q] = 0.0;
               c->grady[q] = 0.0;
               c->gradz[q] = 0.0;
            }
         }
      }
   }

   //MHD; set B-fluxes at face centers
   
   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;
   int i_n;       //index and max_index in the direction of traversal
   int n[3] = {0};
   double dA;


   for( q=0 ; q<NUM_M ; ++q ){

      n[q] = 1;
      double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];
       
      for( k=0 ; k<k_index ; ++k ){
         for( j=0 ; j<j_index ; ++j ){
            for( i=0 ; i<i_index ; ++i ){
               i_n = n[0]*i+n[1]*j+n[2]*k;

               ijk_0 = i;
               if( i_n==0 ) ijk_m1 = i;
               else ijk_m1 = i-n[0];
               if( theDomain->theParList.Num_y != 1 ){
                  ijk_0 += (Nx+2*Ng)*j;
                  if( i_n==0 ) ijk_m1 += (Nx+2*Ng)*j;
                  else ijk_m1 += (Nx+2*Ng)*(j-n[1]);
               }
               if( theDomain->theParList.Num_z != 1 ){
                  ijk_0 += (Nx+2*Ng)*(Ny+2*Ng)*k;
                  if( i_n==0 ) ijk_m1 += (Nx+2*Ng)*(Ny+2*Ng)*k;
                  else ijk_m1 += (Nx+2*Ng)*(Ny+2*Ng)*(k-n[2]);
               }
               struct cell * c0  = theCells+ijk_0;
               struct cell * cm1 = theCells+ijk_m1;
               c0->Phi_B[q] = (c0->prim[BB1+q]+cm1->prim[BB1+q])/2.*dA;
            }
         }   
      }

      n[q] = 0;

   }


}


void freeDomain( struct domain * theDomain ){
   free( theDomain->theCells );
}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final = 0;
   int rank = theDomain->rank;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }
   
   if( rank==0 ){
      FILE * abort = NULL;
      abort = fopen("abort","r");
      if( abort ){ final = 1; fclose(abort); }
   }

   MPI_Allreduce( MPI_IN_PLACE , &final , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   if( final ) theDomain->final_step = 1;   

}

double get_wn( double * , double , double , double , double , int );
double mindt( double * , double * , int * , double , double , double );

double getmindt( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double W  = theDomain->W;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   double dt = 1e100;
   double wcell[3];
   int dim_ind[3] = {1};
   int i,j,k,ijk,dim;

   if( theDomain->theParList.Num_x == 1 ) dim_ind[0] = 0;
   if( theDomain->theParList.Num_y == 1 ) dim_ind[1] = 0;
   if( theDomain->theParList.Num_z == 1 ) dim_ind[2] = 0;

   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Ny ; ++j ){
         for( i=0 ; i<Nx ; ++i ){
            ijk  = i+Ng;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*(j+Ng);
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*(k+Ng);
            struct cell * c = theCells+ijk;
            for( dim = 0 ; dim<3 ; ++dim )
               wcell[dim] = get_wn( c->xi , dx , dy , dz , W , dim );
            double dt_temp = mindt( c->prim , wcell , dim_ind , dx , dy , dz );
            if( dt > dt_temp ) dt = dt_temp;
         }
      }
   }
   
   dt *= theDomain->theParList.CFL; 
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , theDomain->theComm );

   return( dt );

}

void report( struct domain * , double , double);
void snapshot( struct domain * , char * );
void output( struct domain * , char * );

void possiblyOutput( struct domain * theDomain , double dt , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   double Nchk = theDomain->N_chk;
   double Nsnp = theDomain->N_snp;
   int LogOut = theDomain->theParList.Out_LogTime;
   int n0;

   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      report( theDomain , t , dt );
      if( theDomain->rank==0 ) printf("t = %.3e\n",t);
   }

   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         if(theDomain->rank==0) printf("Creating Checkpoint #%04d...\n",n0);
         sprintf(filename,"checkpoint_%04d",n0);
         output( theDomain , filename );
      }else{
         if(theDomain->rank==0) printf("Creating Final Checkpoint...\n");
         output( theDomain , "output" );
      }
   }

   n0 = (int)( t*Nsnp/t_fin );
   if( LogOut ) n0 = (int)( Nsnp*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nsnp < n0 && Nsnp>0) || override ){
      theDomain->nsnp = n0;
      char filename[256];
      if(!override) sprintf( filename , "map_%04d.dat" , n0 );
      else sprintf( filename , "map.dat" );
      snapshot( theDomain , filename );
   }

}
