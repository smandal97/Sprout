
#include "../defs.h"
#include <hdf5.h>

void createFile( char * fname ){
   hid_t h5file = H5Fcreate( fname , H5F_ACC_TRUNC , H5P_DEFAULT , H5P_DEFAULT );
   H5Fclose( h5file );
}

void createGroup( char * fname , char * gname ){
   hid_t h5file = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5group = H5Gcreate2( h5file , gname , H5P_DEFAULT , H5P_DEFAULT , H5P_DEFAULT );
   H5Gclose( h5group );
   H5Fclose( h5file );
}

void createDataset( char * fname , char * gname , char * dname , int dim , hsize_t * fdims , hid_t type ){
   hid_t h5file  = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5group = H5Gopen2( h5file , gname , H5P_DEFAULT );
   hid_t fspace  = H5Screate_simple(dim,fdims,NULL);
   hid_t h5dset  = H5Dcreate2( h5group , dname , type , fspace , H5P_DEFAULT , H5P_DEFAULT , H5P_DEFAULT);
   H5Sclose( fspace );
   H5Dclose( h5dset );
   H5Gclose( h5group );
   H5Fclose( h5file );
}

void writeSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dwrite( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void writePatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
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

   H5Dwrite( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void Cell2Doub( struct cell * c , double * Q , int v , int count ){
   if( v<NUM_Q )
      Q[count] = c->prim[v];
   else
      Q[count] = c->xi[v-NUM_Q];
}

int getN0( int , int , int );

void output( struct domain * theDomain , char * filestart ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;   
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   char filename[256];
   sprintf(filename,"%s.h5",filestart);

   int Nxtot = theDomain->theParList.Num_x + dim_size[0]*2*Ng;
   int Nytot = theDomain->theParList.Num_y + dim_size[1]*2*Ng;
   int Nztot = theDomain->theParList.Num_z + dim_size[2]*2*Ng;
   int Ndoub = NUM_Q+3;  

   hsize_t fdims[1];
   hsize_t datdims[4];
   if( rank==0 ){
      printf("Writing checkpoint...\n");
      
      createFile(filename);
      createGroup(filename,"Grid");

      fdims[0] = 1;
      createDataset(filename,"Grid","T"    ,1,fdims,H5T_NATIVE_DOUBLE);
      fdims[0] = 3;
      createDataset(filename,"Grid","Ntots",1,fdims,   H5T_NATIVE_INT);
      fdims[0] = 3;
      createDataset(filename,"Grid","dxs"  ,1,fdims,H5T_NATIVE_DOUBLE);
      fdims[0] = 3;
      createDataset(filename,"Grid","Origin",1,fdims,H5T_NATIVE_DOUBLE);
      
      createGroup(filename,"Data");
      
      datdims[3] = Nxtot;
      datdims[2] = Nytot;
      datdims[1] = Nztot;
      datdims[0] = Ndoub;
      
      createDataset(filename,"Data","Cells",4,datdims,H5T_NATIVE_DOUBLE);
   }
   MPI_Barrier( theDomain->theComm ); 
   if( rank==0 ){
      //write time
      writeSimple(filename,"Grid","T",&(theDomain->t),H5T_NATIVE_DOUBLE);

      int Ntots[3]; 
      Ntots[0] = Nxtot; 
      Ntots[1] = Nytot; 
      Ntots[2] = Nztot;

      double dxs[3]; 
      dxs[0] = dx; 
      dxs[1] = dy; 
      dxs[2] = dz;

      double Oxs[3];
      int ijk_origin = Ng;
      if( theDomain->theParList.Num_y != 1 ) ijk_origin += (Nx+2*Ng)*Ng;
      if( theDomain->theParList.Num_z != 1 ) ijk_origin += (Nx+2*Ng)*(Ny+2*Ng)*Ng;
      struct cell * c_origin = theCells + ijk_origin;
      Oxs[0] = c_origin->xi[0];
      Oxs[1] = c_origin->xi[1];
      Oxs[2] = c_origin->xi[2];

      int start_grid[1]    = {0};
      int loc_size_grid[1] = {3};
      int glo_size_grid[1] = {3};
      
      //write grid data
      writePatch( filename , "Grid" , "Ntots" , Ntots , H5T_NATIVE_INT , 1 , start_grid , loc_size_grid , glo_size_grid );
      writePatch( filename , "Grid" , "dxs"  , dxs , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );
      writePatch( filename , "Grid" , "Origin" , Oxs   , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );
   }
   
   //write cell data
   int v,i,j,k,ijk,rnk,count=0;
   double * Qwrite = (double *)malloc( (Nx+2*Ng)*(Ny+2*Ng)*(Nz+2*Ng)*Ndoub*sizeof(double) );

   for( rnk=0 ; rnk<size ; ++rnk ){
      if( rnk == rank ){

         int start[4]    = {0,0,0,0};
         int loc_size[4] = {Ndoub,Nz+2*Ng,Ny+2*Ng,Nx+2*Ng};
         int glo_size[4] = {Ndoub,Nztot,  Nytot,  Nxtot};

         for( v=0 ; v<Ndoub ; ++v ){
            for( k=0 ; k<Nz+2*Ng ; ++k ){
               for( j=0 ; j<Ny+2*Ng ; ++j ){
                  for( i=0 ; i<Nx+2*Ng ; ++i ){
                     ijk  = i;
                     if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
                     if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
                     struct cell * c = theCells + ijk;
                     Cell2Doub( c , Qwrite , v , count );
                     ++count;
                  }
               }
            }
         }

         start[1] = getN0( dim_rank[2] , dim_size[2] , Nztot-dim_size[2]*2*Ng ) + dim_rank[2]*2*Ng; //z_offset
         start[2] = getN0( dim_rank[1] , dim_size[1] , Nytot-dim_size[1]*2*Ng ) + dim_rank[1]*2*Ng; //y_offset
         start[3] = getN0( dim_rank[0] , dim_size[0] , Nxtot-dim_size[0]*2*Ng ) + dim_rank[0]*2*Ng; //x_offset
         
         writePatch( filename , "Data" , "Cells" , Qwrite , H5T_NATIVE_DOUBLE , 4 , start , loc_size , glo_size );

      }
    
      MPI_Barrier(theDomain->theComm);    

   }

   free(Qwrite);

}





