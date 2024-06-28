
#include "defs.h"

void slice( struct domain * theDomain , char * filestart );
void sphere_surface( struct domain * theDomain , char * filestart );

void snapshot( struct domain * theDomain , char * filestart ){
   slice(theDomain,filestart);
   //sphere_surface(theDomain,filestart);
}


void slice( struct domain * theDomain , char * filestart ){
//get planar slice
//is not foolproof in 2D

   int * dim_rank = theDomain->dim_rank;
   int Num_x = theDomain->theParList.Num_x;
   int Num_y = theDomain->theParList.Num_y;
   double dx = theDomain->dx;
   double t  = theDomain->t;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   int i,j,ijk,ij_sl;
   int i0=0; 
   int j0=0;
   int i1=1; 
   int j1=1;
   int k0 = 0; //choose slice for snapshot


   if( Num_x != 1 ){
      i0 = Ng; 
      i1 = Nx+Ng;
   }
   if( Num_y != 1 ){
      j0 = Ng; 
      j1 = Ny+Ng;
   }


   double * slice = (double *)calloc( Num_x*Num_y,sizeof(double) );


   for( j=j0 ; j<j1 ; ++j ){
      for( i=i0 ; i<i1 ; ++i ){
         ijk   = i + (Nx+2*Ng)*j + (Nx+2*Ng)*(Ny+2*Ng)*k0;
         ij_sl = (i-Ng+dim_rank[0]*Nx) + (j-Ng+dim_rank[1]*Ny)*Num_x; 
         struct cell * c = theDomain->theCells+ijk;
         slice[ij_sl] = c->prim[RHO];  //choose variable for snapshot
      }
   }
   
   MPI_Allreduce( MPI_IN_PLACE , slice , Num_x*Num_y , MPI_DOUBLE , MPI_SUM , theDomain->theComm );

   if( theDomain->rank==0 ){
      char fname1[256];
      FILE * FP;
      strcpy( fname1 , "slice" );
      strcat( fname1 , filestart );
      FP = fopen( fname1 , "w" );
      fprintf(FP, "# %e %e\n", dx, t);
      for( j=0 ; j<Num_y ; ++j ){
       for( i=0 ; i<Num_x ; ++i ){
          ij_sl = i+j*Num_x;
          fprintf(FP, "%e ", slice[ij_sl] );   
       }
       fprintf(FP, "\n");
      }

      fclose(FP);

   }


   free(slice);

}





void sphere_surface( struct domain * theDomain , char * filestart ){
//calculate radially integrated spherical surface maps
//ONLY WORKS AS INTENDED WITH 3D

   //set snapshot parameters; DO AUTO LATER
   int Nt,Np;
   double rmin,rmax,tmin,pmin,tmax,pmax;
   double xo,yo,zo,eta_min,eta_max;
   Nt = 200;
   Np = 200;
   tmin = 0.0;
   pmin = 0.0;
   tmax = M_PI/2.;
   pmax = M_PI/2.;
   xo  = theDomain->theParList.MM_x0 * theDomain->theParList.Lx;
   yo  = theDomain->theParList.MM_y0 * theDomain->theParList.Ly;
   zo  = theDomain->theParList.MM_z0 * theDomain->theParList.Lz;
   eta_min = 0.35;
   eta_max = 0.99;
   rmin = theDomain->dx*(double)(theDomain->theParList.Num_x) * eta_min;
   rmax = theDomain->dx*(double)(theDomain->theParList.Num_x) * eta_max;


   int p,q,pq;
   double vavg=0., davg=0.;
   double *vmap, *dmap, *map_wt;
   vmap   = (double *)calloc( Nt*Np,sizeof(double) );
   dmap   = (double *)calloc( Nt*Np,sizeof(double) );
   map_wt = (double *)calloc( Nt*Np,sizeof(double) );


   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   int i,j,k,ijk;
   int i0=0 ; int j0=0 ; int k0=0; 
   int i1=1 ; int j1=1 ; int k1=1;
   double x,y,z,r,th,ph,rho,ppp,vr;

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
            x   = c->xi[0] - xo;
            y   = c->xi[1] - yo;
            z   = c->xi[2] - zo;
            r   = sqrt(x*x+y*y+z*z); //printf("r = %e,%e,%e\n",r,rmin,rmax);
            th  = acos(z/r);
            ph  = atan2(y,x) + M_PI;
            rho = c->prim[RHO];
            ppp = c->prim[PPP];
            vr  = (c->prim[UU1]*x + c->prim[UU2]*y + c->prim[UU3]*z)/r;
            p   = floor( ph/((pmax-pmin)/(float)Np) ); 
            q   = floor( th/((tmax-tmin)/(float)Nt) ); 
            pq  = p + Np*q;
            if( r>rmin && r<rmax ){
            	vmap[pq]   += ppp*vr; //printf("p,q = %i,%i\n",p,q);
            	dmap[pq]   += ppp*rho;
            	map_wt[pq] += ppp;
            	vavg       += vr;
            	davg       += rho;
            }

         }
      }
   }

   vavg /= (double)(Nt*Np);
   davg /= (double)(Nt*Np);

   MPI_Allreduce( MPI_IN_PLACE , vmap   , Nt*Np , MPI_DOUBLE , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , dmap   , Nt*Np , MPI_DOUBLE , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , map_wt , Nt*Np , MPI_DOUBLE , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , &vavg  , 1     , MPI_DOUBLE , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , &davg  , 1     , MPI_DOUBLE , MPI_SUM , theDomain->theComm );

   if( theDomain->rank==0 ){
      char fname1[256],fname2[256];
      FILE * FP1, * FP2;
      strcpy( fname1 , "v" );
      strcat( fname1 , filestart );
      FP1 = fopen( fname1 , "w" );
      strcpy( fname2 , "d" );
      strcat( fname2 , filestart );
      FP2 = fopen( fname2 , "w" );
      for( q=0 ; q<Nt ; ++q ){
       for( p=0 ; p<Np ; ++p ){
          pq  = p + Np*q;
          fprintf(FP1, "%e\n", vmap[pq]/map_wt[pq]/vavg );
          fprintf(FP2, "%e\n", dmap[pq]/map_wt[pq]/davg );   
       }
      }

      fclose(FP2);
      fclose(FP1);

   }



   free(vmap);
   free(dmap);
   free(map_wt);


}
