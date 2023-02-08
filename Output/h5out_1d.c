
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

int Cell2Doub( struct cell * c , double * Q , int mode ){
   if( mode==0 ) return(NUM_Q+3); else{
      int q;
      for( q=0 ; q<NUM_Q ; ++q ) Q[q] = c->prim[q];
      for( q=0 ; q<3 ; ++q ) Q[NUM_Q+q] = c->xi[q];
      return(0);
   }
}

void output( struct domain * theDomain , char * filestart ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;   
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   char filename[256];
   sprintf(filename,"%s.h5",filestart);

   int Nxtot = theDomain->theParList.Num_x;
   int Nytot = theDomain->theParList.Num_y;
   int Nztot = theDomain->theParList.Num_z;
   int Ntot  = Nxtot * Nytot * Nztot;
   int Ndoub = Cell2Doub(NULL,NULL,0);  

   hsize_t fdims[1];
   hsize_t datdims[2];
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
      
      createGroup(filename,"Data");
      
      datdims[0] = Ntot;
      datdims[1] = Ndoub;
     
      createDataset(filename,"Data","Cells",2,datdims,H5T_NATIVE_DOUBLE);
   }
   MPI_Barrier( theDomain->theComm ); 
   if( rank==0 ){
      //write time
      writeSimple(filename,"Grid","T",&(theDomain->t),H5T_NATIVE_DOUBLE);

      int Ntots[3]; Ntots[0] = Nxtot; Ntots[1] = Nytot; Ntots[2] = Nztot;
      double dxs[3]; dxs[0] = dx; dxs[1] = dy; dxs[2] = dz;

      int start_grid[1]    = {0};
      int loc_size_grid[1] = {3};
      int glo_size_grid[1] = {3};
      
      //write grid data
      writePatch( filename , "Grid" , "Ntots" , Ntots , H5T_NATIVE_INT , 1 , start_grid , loc_size_grid , glo_size_grid );
      writePatch( filename , "Grid" , "dxs"  , dxs , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );
       
   }
   
   //write cell data
   int i,j,k,ijk,rnk;
   double * Qwrite = (double *)malloc( Nx*Ndoub*sizeof(double) );

   for( rnk=0 ; rnk<size ; ++rnk ){
      if( rnk == rank ){

         int start[2]    = {0,0};
         int loc_size[2] = {Nx  ,Ndoub};
         int glo_size[2] = {Ntot,Ndoub};
         int x_offset,y_offset,z_offset;

         for( k=0 ; k<Nz ; ++k ){
            for( j=0 ; j<Ny ; ++j ){
               for( i=0 ; i<Nx ; ++i ){
                  ijk  = i+Ng;
                  if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*(j+Ng);
                  if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*(k+Ng);
                  struct cell * c = theCells + ijk;
                  Cell2Doub( c , Qwrite+Ndoub*i , 1 );
               }
               z_offset = (dim_rank[2] + k) * Nxtot * Nytot;
               y_offset = (dim_rank[1] * Ny + j) * Nxtot;
               x_offset =  dim_rank[0] * Nx;
               start[0] = x_offset + y_offset + z_offset;
               writePatch( filename , "Data" , "Cells" , Qwrite , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
            }
      
         }

      }
    
      MPI_Barrier(theDomain->theComm);    

   }

   free(Qwrite);
   

}
