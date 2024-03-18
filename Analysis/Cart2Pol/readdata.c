
#include "defs.h"
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

void readdims( char * filename , struct domain * theDomain ){

   char group1[256];
   strcpy(group1,"Grid");

   int Ntots[3];
   double tstart;

   int start_grid[1]    = {0};
   int loc_size_grid[1] = {3};
   int glo_size_grid[1] = {3};
   
   readSimple( filename , group1 , "T" , &tstart , H5T_NATIVE_DOUBLE );
   readPatch( filename , group1 , "Ntots" , Ntots , H5T_NATIVE_INT , 1 , start_grid , loc_size_grid , glo_size_grid );
    
   theDomain->Nx = Ntots[0];
   theDomain->Ny = Ntots[1];
   theDomain->Nz = Ntots[2];
  
}

void readdxs( char * filename , struct domain * theDomain ){

   char group1[256];
   strcpy(group1,"Grid");

   double dxs[3];
   double tstart;

   int start_grid[1]    = {0};
   int loc_size_grid[1] = {3};
   int glo_size_grid[1] = {3};
   
   readSimple( filename , group1 , "T" , &tstart , H5T_NATIVE_DOUBLE );
   readPatch( filename , group1 , "dxs" , dxs , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );
    
   theDomain->dx = dxs[0];
   theDomain->dy = dxs[1];
   theDomain->dz = dxs[2];
}

void readmincoords( char * filename , struct domain * theDomain ){

   char group1[256];
   strcpy(group1,"Grid");

   double mins[3];
   double tstart;

   int start_grid[1]    = {0};
   int loc_size_grid[1] = {3};
   int glo_size_grid[1] = {3};
   
   readSimple( filename , group1 , "T" , &tstart , H5T_NATIVE_DOUBLE );
   readPatch( filename , group1 , "Origin" , mins , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );
    
   theDomain->xmin = mins[0];
   theDomain->ymin = mins[1];
   theDomain->zmin = mins[2];
}

void readdata( char * filename , struct domain * theDomain ){

   readdxs(  filename , theDomain );
   readdims( filename , theDomain );
   readmincoords( filename , theDomain );
   int Ndoub = NUM_Q;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;

   theDomain->cartdata = (double *)malloc( Ndoub*Nx*Ny*Nz*sizeof(double) );

   char group2[256];
   strcpy(group2,"Data");

   int start[4]    = {0,0,0,0};
   int loc_size[4] = {Ndoub,Nz,Ny,Nx};
   int glo_size[4] = {Ndoub,Nz,Ny,Nx};

   readPatch( filename , group2 , "Cells" , theDomain->cartdata , H5T_NATIVE_DOUBLE , 4 , start , loc_size , glo_size );

}
