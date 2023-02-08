enum{SND,RCV};

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

struct cell_lite{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double xi[3];
   double RKxi[3];
};

void generate_mpi_cell( MPI_Datatype * cell_mpi ){

   struct cell_lite test;
   int count = 5;
   int blocksize[]      = {NUM_Q,NUM_Q,NUM_Q,3,3};
   MPI_Datatype types[] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
   MPI_Aint offsets[5];

   offsets[0] = (char *)&(test.prim)   - (char *)(&test);
   offsets[1] = (char *)&(test.cons)   - (char *)(&test);
   offsets[2] = (char *)&(test.RKcons) - (char *)(&test);
   offsets[3] = (char *)&(test.xi)     - (char *)(&test);
   offsets[4] = (char *)&(test.RKxi)   - (char *)(&test);

   MPI_Type_create_struct( count , blocksize , offsets , types , cell_mpi );
   MPI_Type_commit( cell_mpi );

}

void copy_cell_to_lite( struct cell * c , struct cell_lite * cl ){
  
   memcpy( cl->prim   , c->prim   , NUM_Q*sizeof(double) ); 
   memcpy( cl->cons   , c->cons   , NUM_Q*sizeof(double) ); 
   memcpy( cl->RKcons , c->RKcons , NUM_Q*sizeof(double) );
   memcpy( cl->xi     , c->xi     ,     3*sizeof(double) ); 
   memcpy( cl->RKxi   , c->RKxi   ,     3*sizeof(double) );

}

void copy_lite_to_cell( struct cell_lite * cl , struct cell * c ){

   memcpy( c->prim   , cl->prim   , NUM_Q*sizeof(double) ); 
   memcpy( c->cons   , cl->cons   , NUM_Q*sizeof(double) ); 
   memcpy( c->RKcons , cl->RKcons , NUM_Q*sizeof(double) );
   memcpy( c->xi     , cl->xi     ,     3*sizeof(double) ); 
   memcpy( c->RKxi   , cl->RKxi   ,     3*sizeof(double) );

}

void choosecells( struct domain * theDomain , struct cell_lite * pl , struct cell_lite * pr , int theDIM , int mode ){
   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   int i,j,k,i1,j1,k1,ijk,ct;
   
   int n[3]   = {0};
   n[theDIM]  = 1;
   int cn[3]; cn[0] = 1; cn[1] = 1; cn[2] = 1;
   cn[theDIM] = 0;
   
   if( mode==SND ){
      //LEFT SLAB
      ct = 0;
      i1 = Nx*cn[0]+2*Ng; if( theDomain->theParList.Num_x == 1 ) i1 = 1;
      j1 = Ny*cn[1]+2*Ng; if( theDomain->theParList.Num_y == 1 ) j1 = 1;
      k1 = Nz*cn[2]+2*Ng; if( theDomain->theParList.Num_z == 1 ) k1 = 1;

      for( k=Ng*n[2] ; k<k1 ; ++k ){
         for( j=Ng*n[1] ; j<j1 ; ++j ){
            for( i=Ng*n[0] ; i<i1 ; ++i ){
               ijk  = i;
               if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * cL = theCells+ijk;
               copy_cell_to_lite( cL , pl+ct );
               ++ct;
            }
         }
      }
      
      //RIGHT SLAB
      ct = 0;
      i1 = Nx+(1+cn[0])*Ng; if( theDomain->theParList.Num_x == 1 ) i1 = 1;
      j1 = Ny+(1+cn[1])*Ng; if( theDomain->theParList.Num_y == 1 ) j1 = 1;
      k1 = Nz+(1+cn[2])*Ng; if( theDomain->theParList.Num_z == 1 ) k1 = 1;

      for( k=Nz*n[2] ; k<k1 ; ++k ){
         for( j=Ny*n[1] ; j<j1 ; ++j ){
            for( i=Nx*n[0] ; i<i1 ; ++i ){
               ijk  = i;
               if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * cR = theCells+ijk;
               copy_cell_to_lite( cR , pr+ct );
               ++ct;
            }
         }
      }
      
   }else if( mode==RCV ){
      //LEFT SLAB
      ct=0;
      i1 = (Nx+Ng)*cn[0]+Ng; if( theDomain->theParList.Num_x == 1 ) i1 = 1;
      j1 = (Ny+Ng)*cn[1]+Ng; if( theDomain->theParList.Num_y == 1 ) j1 = 1;
      k1 = (Nz+Ng)*cn[2]+Ng; if( theDomain->theParList.Num_z == 1 ) k1 = 1;

      for( k=0 ; k<k1 ; ++k ){
         for( j=0 ; j<j1 ; ++j ){
            for( i=0 ; i<i1 ; ++i ){
               ijk  = i;
               if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * cL = theCells+ijk;
               copy_lite_to_cell( pl+ct , cL );
               ++ct;
            }
         }
      }
      
      //RIGHT SLAB
      ct = 0;
      i1 = Nx+2*Ng; if( theDomain->theParList.Num_x == 1 ) i1 = 1;
      j1 = Ny+2*Ng; if( theDomain->theParList.Num_y == 1 ) j1 = 1;
      k1 = Nz+2*Ng; if( theDomain->theParList.Num_z == 1 ) k1 = 1;

      for( k=(Nz+Ng)*n[2] ; k<k1 ; ++k ){
         for( j=(Ny+Ng)*n[1] ; j<j1 ; ++j ){
            for( i=(Nx+Ng)*n[0] ; i<i1 ; ++i ){
               ijk  = i;
               if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
               if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
               struct cell * cR = theCells+ijk;
               copy_lite_to_cell( pr+ct , cR );
               ++ct;
            }
         }
      }
      
   }
}

void exchange1D( struct domain * theDomain , int theDIM ){
   
   MPI_Datatype cell_mpi = {0}; 
   generate_mpi_cell( &cell_mpi );

   MPI_Comm grid_comm = theDomain->theComm;
   int * left_rank  = theDomain->left_rank;
   int * right_rank = theDomain->right_rank;
   int Ng = theDomain->Ng;
   
   int tag = 0;
   MPI_Status status;
   
   //Set packet dimensions for sending
   int i,send_num[3]; 
   send_num[0] = theDomain->Nx + Ng; 
   send_num[1] = theDomain->Ny + Ng;
   send_num[2] = theDomain->Nz + Ng;
   send_num[theDIM] = 0; 
   for( i=0 ; i<3 ; ++i ) send_num[i] += Ng;
   
   int send_size = 1;
   if( theDomain->theParList.Num_x != 1 ) send_size *= send_num[0];
   if( theDomain->theParList.Num_y != 1 ) send_size *= send_num[1];
   if( theDomain->theParList.Num_z != 1 ) send_size *= send_num[2];
   
   ////////////////////////
   //Send and receive cells
   ////////////////////////

   struct cell_lite * pl_send , * pr_send , * pl_recv , * pr_recv;
   pl_send = (struct cell_lite *) malloc( send_size*sizeof(struct cell_lite) );
   pr_send = (struct cell_lite *) malloc( send_size*sizeof(struct cell_lite) );
   pl_recv = (struct cell_lite *) malloc( send_size*sizeof(struct cell_lite) );
   pr_recv = (struct cell_lite *) malloc( send_size*sizeof(struct cell_lite) );
   
   //Build up list of cells to send...
   choosecells( theDomain , pl_send , pr_send , theDIM , SND ); 
   //Send!
   MPI_Sendrecv( pl_send , send_size , cell_mpi ,  left_rank[theDIM] , tag+2,
                 pr_recv , send_size , cell_mpi , right_rank[theDIM] , tag+2, grid_comm , &status);
   MPI_Sendrecv( pr_send , send_size , cell_mpi , right_rank[theDIM] , tag+3,
                 pl_recv , send_size , cell_mpi ,  left_rank[theDIM] , tag+3, grid_comm , &status);
   //Now take the list of cells and put them into the appropriate locations...
   choosecells( theDomain , pl_recv , pr_recv , theDIM , RCV );
   
   free(pl_send);
   free(pr_send);
   free(pl_recv);
   free(pr_recv);
   MPI_Type_free( &cell_mpi );
     
}

void exchangeData( struct domain * theDomain ){

   if( theDomain->theParList.Num_x != 1 )exchange1D( theDomain , 0 );
   if( theDomain->theParList.Num_y != 1 )exchange1D( theDomain , 1 );
   if( theDomain->theParList.Num_z != 1 )exchange1D( theDomain , 2 );
 
}
