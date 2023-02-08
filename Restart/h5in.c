
#include "../defs.h"
#include <hdf5.h>
#include <string.h>

void getH5dims( char * file , char * group , char * dset , hsize_t * dims ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );
   hid_t h5spc = H5Dget_space( h5dst );

   H5Sget_simple_extent_dims( h5spc , dims , NULL);

   H5Sclose( h5spc );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readPatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];

   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = loc_size[d];
      fdims[d] = glo_size[d];

      fstart[d]  = start[d];
      fstride[d] = 1;
      fcount[d]  = loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}


void Doub2Cell( double * Q , struct cell * c , int v , int count ){
   if( v<NUM_Q )
      c->prim[v] = Q[count];
   else
      c->xi[v-NUM_Q] = Q[count];
}

int getN0( int , int , int );
void freeDomain( struct domain * );

void restart( struct domain * theDomain ){

   freeDomain( theDomain );

   int Ng = NUM_G;
   int Ndoub = NUM_Q;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   char filename[256];
   strcpy(filename,"input.h5");
   char group1[256];
   strcpy(group1,"Grid");
   char group2[256];
   strcpy(group2,"Data");

   if( rank==0 ) printf("Restarting from file...\n");

   int Ntots[3];
   double dxs[3];
   double Oxs[3]={0};
   double tstart;

   int start_grid[1]    = {0};
   int loc_size_grid[1] = {3};
   int glo_size_grid[1] = {3};
   
   if(rank==0){
      readSimple( filename , group1 , "T" , &tstart , H5T_NATIVE_DOUBLE );
      readPatch( filename , group1 , "Ntots"  , Ntots , H5T_NATIVE_INT    , 1 , start_grid , loc_size_grid , glo_size_grid );
      readPatch( filename , group1 , "dxs"    , dxs   , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );
      readPatch( filename , group1 , "Origin" , Oxs   , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );
   }
   MPI_Bcast( Ntots   , 3 , MPI_INT    , 0 , theDomain->theComm );
   MPI_Bcast( dxs     , 3 , MPI_DOUBLE , 0 , theDomain->theComm );
   MPI_Bcast( Oxs     , 3 , MPI_DOUBLE , 0 , theDomain->theComm );
   MPI_Bcast( &tstart , 1 , MPI_DOUBLE , 0 , theDomain->theComm );
   theDomain->theParList.Num_x = Ntots[0];
   theDomain->theParList.Num_y = Ntots[1];
   theDomain->theParList.Num_z = Ntots[2];
   theDomain->dx = dxs[0];
   theDomain->dy = dxs[1];
   theDomain->dz = dxs[2];
   theDomain->t = tstart;
 
   int Nxtot = Ntots[0];
   int Nytot = Ntots[1];
   int Nztot = Ntots[2];
   //The following is very similar to equivalent code in gridsetup.c
   //Now you're just doing the process over because you're restarting
   //from file. 
   int N0x = getN0( dim_rank[0]   , dim_size[0] , Nxtot );
   int N0y = getN0( dim_rank[1]   , dim_size[1] , Nytot );
   int N0z = getN0( dim_rank[2]   , dim_size[2] , Nztot );

   int N1x = getN0( dim_rank[0]+1 , dim_size[0] , Nxtot );
   int N1y = getN0( dim_rank[1]+1 , dim_size[1] , Nytot );
   int N1z = getN0( dim_rank[2]+1 , dim_size[2] , Nztot );

   int Nx = N1x - N0x;
   int Ny = N1y - N0y;
   int Nz = N1z - N0z;

   theDomain->Nx = Nx;
   theDomain->Ny = Ny;
   theDomain->Nz = Nz;
 
   int domain_size = 1;
   if( theDomain->theParList.Num_x != 1 ) domain_size *= Nx+2*Ng;
   if( theDomain->theParList.Num_y != 1 ) domain_size *= Ny+2*Ng;
   if( theDomain->theParList.Num_z != 1 ) domain_size *= Nz+2*Ng;

   theDomain->theCells = (struct cell *) malloc( domain_size*sizeof(struct cell) );
   int v,i,j,k,ijk,rnk,count=0;
   
   double * Qread = (double *)malloc( Nx*Ny*Nz*Ndoub*sizeof(double) );

   for( rnk=0 ; rnk<size ; ++rnk ){
      if( rnk == rank ){

         int start[4]    = {0,N0z,N0y,N0x};
         int loc_size[4] = {Ndoub,Nz,Ny,Nx};
         int glo_size[4] = {Ndoub,Nztot,Nytot,Nxtot};

         readPatch( filename , group2 , "Cells" , Qread , H5T_NATIVE_DOUBLE , 4 , start , loc_size , glo_size );

         for( v=0 ; v<Ndoub ; ++v ){
            for( k=0 ; k<Nz ; ++k ){
               for( j=0 ; j<Ny ; ++j ){
                  for( i=0 ; i<Nx ; ++i ){
                     ijk  = i+Ng;
                     if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*(j+Ng);
                     if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*(k+Ng);
                     struct cell * c = theDomain->theCells + ijk;
                     Doub2Cell( Qread , c , v , count );
                     c->xi[0] = Oxs[0]+dxs[0]*(double)(i+N0x);
                     c->xi[1] = Oxs[1]+dxs[1]*(double)(j+N0y);
                     c->xi[2] = Oxs[2]+dxs[2]*(double)(k+N0z);
                     ++count;
                  }
               }
            }
         }
         

      }

      MPI_Barrier(theDomain->theComm);

   }
 
   free(Qread);
   
}

