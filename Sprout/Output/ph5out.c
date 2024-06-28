
#include "../defs.h"
#include <hdf5.h>

char * change_ext( char * mystr );
void makexdmf( struct domain * theDomain , char * filename , double * Oxs );



void createDataset( hid_t group_id , char * dname , int dim , hsize_t * fdims , hid_t type ){
   hid_t fspace  = H5Screate_simple(dim,fdims,NULL);
   hid_t h5dset  = H5Dcreate2( group_id , dname , type , fspace , H5P_DEFAULT , H5P_DEFAULT , H5P_DEFAULT);
   H5Sclose( fspace );
   H5Dclose( h5dset );
}

void writeSimple( hid_t h5grp , char * dset , void * data , hid_t type ){
   hid_t h5dst   = H5Dopen1( h5grp , dset );
   //set up the collective transfer properties list
   //hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
   //H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
   hid_t dxpl_id = H5P_DEFAULT; //use for independent writes 
   H5Dwrite( h5dst , type , H5S_ALL , H5S_ALL , dxpl_id , data );
   H5Pclose( dxpl_id );
   H5Dclose( h5dst );
}

void writePatch( hid_t h5grp , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
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

   //set up the collective transfer properties list
   //hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER); //use for collective writes
   //H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
   hid_t dxpl_id = H5P_DEFAULT; //use for independent writes

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );
   H5Dwrite( h5dst , type , mspace , fspace , dxpl_id , data );

   H5Pclose( dxpl_id );
   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
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



   //prepare grid data
   int Nxtot = theDomain->theParList.Num_x;
   int Nytot = theDomain->theParList.Num_y;
   int Nztot = theDomain->theParList.Num_z;
   int Ndoub = NUM_Q;
   int Ntots[3]; 
   Ntots[0] = Nxtot; 
   Ntots[1] = Nytot; 
   Ntots[2] = Nztot;
   double dxs[3]; 
   dxs[0] = theDomain->dx;
   dxs[1] = theDomain->dy; 
   dxs[2] = theDomain->dz;
   double Oxs[3];
   int ijk_origin = Ng;
   if( theDomain->theParList.Num_y != 1 ) ijk_origin += (Nx+2*Ng)*Ng;
   if( theDomain->theParList.Num_z != 1 ) ijk_origin += (Nx+2*Ng)*(Ny+2*Ng)*Ng;
   struct cell * c_origin = theCells + ijk_origin;
   Oxs[0] = c_origin->xi[0];
   Oxs[1] = c_origin->xi[1];
   Oxs[2] = c_origin->xi[2];



   //define file variables
   char filename[256];
   sprintf(filename,"%s.h5",filestart);

   int ierr;
   MPI_Info FILE_INFO;
   hid_t FILE_ID;
   hid_t FILE_ACCESS;

   //START WRITING!
   time_t t_b, t_e;
   t_b = time(NULL);
   if( rank==0 ) printf("Writing checkpoint...\n");
   //if(Nxtot>1 && Nytot>1 && Nztot>1 && rank==0) makexdmf( theDomain,filename,Oxs );

   //set file access template for parallel I/O access
   FILE_ACCESS = H5Pcreate(H5P_FILE_ACCESS);
   //ensure metadata is accessed collectively for all namespace objects
   //H5Pset_all_coll_metadata_ops(FILE_ACCESS, 1);
   //H5Pset_coll_metadata_write(FILE_ACCESS, 1);


   //-------------------------------------------
   //PLATFORM-DEPENDENCE START
   //-------------------------------------------
   //tune access template for filesystem 

   //create MPI_INFO object
   ierr = MPI_Info_create(&FILE_INFO);

   ierr = H5Pset_sieve_buf_size(FILE_ACCESS, 262144); 
   ierr = H5Pset_alignment(FILE_ACCESS, 524288, 262144);
  
   ierr = MPI_Info_set(FILE_INFO, "access_style", "write_once");
   ierr = MPI_Info_set(FILE_INFO, "collective_buffering", "true");
   ierr = MPI_Info_set(FILE_INFO, "cb_block_size", "1048576");
   ierr = MPI_Info_set(FILE_INFO, "cb_buffer_size", "4194304");
   //ierr = MPI_Info_set(FILE_INFO, "cb_nodes, 3072");


   //tell the HDF5 library that we want to use MPI-IO to do the writing
   ierr = H5Pset_fapl_mpio(FILE_ACCESS, theDomain->theComm, FILE_INFO);

   //-------------------------------------------
   //PLATFORM-DEPENDENCE END
   //-------------------------------------------
   //create file and set file identifier
   FILE_ID = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, FILE_ACCESS);

   //release the file access template
   ierr = H5Pclose(FILE_ACCESS);
   ierr = MPI_Info_free(&FILE_INFO);


   //-------------------------------------------
   //create GRID group and all relevant datasets
   //-------------------------------------------
   hid_t GROUP_G = H5Gcreate2( FILE_ID , "Grid" , H5P_DEFAULT , H5P_DEFAULT , H5P_DEFAULT );

   hsize_t fdims[1];
   //create datasets
   fdims[0] = 1;
   createDataset(GROUP_G,"T"     ,1,fdims,H5T_NATIVE_DOUBLE);
   fdims[0] = 3;
   createDataset(GROUP_G,"Ntots" ,1,fdims,   H5T_NATIVE_INT);
   fdims[0] = 3;
   createDataset(GROUP_G,"dxs"   ,1,fdims,H5T_NATIVE_DOUBLE);
   fdims[0] = 3;
   createDataset(GROUP_G,"Origin",1,fdims,H5T_NATIVE_DOUBLE);

   //write datasets
   if(rank==0){
      //write time
      writeSimple(GROUP_G,"T",&(theDomain->t),H5T_NATIVE_DOUBLE);

      int start_grid[1]    = {0};
      int loc_size_grid[1] = {3};
      int glo_size_grid[1] = {3};
      
      //write grid data
      writePatch( GROUP_G , "Ntots"  , Ntots , H5T_NATIVE_INT    , 1 , start_grid , loc_size_grid , glo_size_grid );
      writePatch( GROUP_G , "dxs"    , dxs   , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );
      writePatch( GROUP_G , "Origin" , Oxs   , H5T_NATIVE_DOUBLE , 1 , start_grid , loc_size_grid , glo_size_grid );

   }
   
   //close group
   H5Gclose(GROUP_G);


   //-------------------------------------------
   //create DATA group and all relevant datasets
   //-------------------------------------------
   hid_t GROUP_D = H5Gcreate2( FILE_ID , "Data" , H5P_DEFAULT , H5P_DEFAULT , H5P_DEFAULT );
   
   hsize_t datdims[4];
   datdims[3] = Nxtot;
   datdims[2] = Nytot;
   datdims[1] = Nztot;
   datdims[0] = Ndoub;
   //create dataset
   createDataset(GROUP_D,"Cells",4,datdims,H5T_NATIVE_DOUBLE);
   /*
   //prep data for writing
   int v,i,j,k,ijk,rnk,count=0;
   double * Qwrite = (double *)malloc( Nx*Ny*Nz*Ndoub*sizeof(double) );

   int start[4]    = {0,0,0,0};
   int loc_size[4] = {Ndoub,Nz,Ny,Nx};
   int glo_size[4] = {Ndoub,Nztot,Nytot,Nxtot};

   for( v=0 ; v<Ndoub ; ++v ){
      for( k=0 ; k<Nz ; ++k ){
         for( j=0 ; j<Ny ; ++j ){
            for( i=0 ; i<Nx ; ++i ){
               ijk  = i+Ng;
               if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*(j+Ng);
               if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*(k+Ng);
               struct cell * c = theCells + ijk;
               Cell2Doub( c , Qwrite , v , count );
               ++count;
            }
         }
      }
   }

   start[1] = getN0( dim_rank[2] , dim_size[2] , Nztot );
   start[2] = getN0( dim_rank[1] , dim_size[1] , Nytot );
   start[3] = getN0( dim_rank[0] , dim_size[0] , Nxtot );

   //write data and free buffer array Qwrite
   writePatch( GROUP_D , "Cells" , Qwrite , H5T_NATIVE_DOUBLE , 4 , start , loc_size , glo_size );
   free(Qwrite);
   */
   //close group
   H5Gclose(GROUP_D);
   //close file
   H5Fclose(FILE_ID);


   t_e = time(NULL);
   if(rank==0) printf("Time taken = %ld seconds\n", t_e-t_b);

}










//function to make xdmf files for VisIt; change_ext() needs work

char * change_ext( char * mystr ) {

   char * tmp = (char *)malloc( (strlen(mystr)-2)*sizeof(char) );
   unsigned int x;

   for (x = 0; x < (strlen(mystr) - 2); x++)
      tmp[x] = mystr[x];
   
   char ext[] = "xmf";
   strcat(tmp,ext);

   return tmp;
}


void makexdmf( struct domain * theDomain , char * filename , double * Oxs ){

   int Nx = theDomain->theParList.Num_x;
   int Ny = theDomain->theParList.Num_y;
   int Nz = theDomain->theParList.Num_z;
   int Nv = NUM_Q;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   FILE *xmf = 0;
 
    /*
     * Open the file and write the XML description of the mesh..
     */
    char * xmf_fname = change_ext(filename);

    xmf = fopen(xmf_fname, "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%d %d %d\"/>\n", Nz+1, Ny+1, Nx+1);
    fprintf(xmf, "     <Geometry GeometryType=\"Origin_DxDyDz\">\n");
    fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
    fprintf(xmf, "        %e %e %e\n", Oxs[0], Oxs[1], Oxs[2] );
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
    fprintf(xmf, "        %e %e %e\n", dx, dy, dz );
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"rho\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n", Nz, Ny, Nx);
    fprintf(xmf, "         <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
    fprintf(xmf, "           0 0 0 0\n");
    fprintf(xmf, "           1 1 1 1\n");
    fprintf(xmf, "           1 %d %d %d\n", Nz, Ny, Nx);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d %d %d %d\" Format=\"HDF\">\n", Nv, Nz, Ny, Nx);
    fprintf(xmf, "           %.*s:/Data/Cells\n", 1+(int)sizeof filename, filename);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    
    fprintf(xmf, "     <Attribute Name=\"p\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n", Nz, Ny, Nx);
    fprintf(xmf, "         <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
    fprintf(xmf, "           1 0 0 0\n");
    fprintf(xmf, "           1 1 1 1\n");
    fprintf(xmf, "           2 %d %d %d\n", Nz, Ny, Nx);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d %d %d %d\" Format=\"HDF\">\n", Nv, Nz, Ny, Nx);
    fprintf(xmf, "           %.*s:/Data/Cells\n", 1+(int)sizeof filename, filename);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    fprintf(xmf, "     <Attribute Name=\"vx\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n", Nz, Ny, Nx);
    fprintf(xmf, "         <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
    fprintf(xmf, "           2 0 0 0\n");
    fprintf(xmf, "           1 1 1 1\n");
    fprintf(xmf, "           3 %d %d %d\n", Nz, Ny, Nx);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d %d %d %d\" Format=\"HDF\">\n", Nv, Nz, Ny, Nx);
    fprintf(xmf, "           %.*s:/Data/Cells\n", 1+(int)sizeof filename, filename);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    fprintf(xmf, "     <Attribute Name=\"vy\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n", Nz, Ny, Nx);
    fprintf(xmf, "         <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
    fprintf(xmf, "           3 0 0 0\n");
    fprintf(xmf, "           1 1 1 1\n");
    fprintf(xmf, "           4 %d %d %d\n", Nz, Ny, Nx);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d %d %d %d\" Format=\"HDF\">\n", Nv, Nz, Ny, Nx);
    fprintf(xmf, "           %.*s:/Data/Cells\n", 1+(int)sizeof filename, filename);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    fprintf(xmf, "     <Attribute Name=\"vz\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n", Nz, Ny, Nx);
    fprintf(xmf, "         <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
    fprintf(xmf, "           4 0 0 0\n");
    fprintf(xmf, "           1 1 1 1\n");
    fprintf(xmf, "           5 %d %d %d\n", Nz, Ny, Nx);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d %d %d %d\" Format=\"HDF\">\n", Nv, Nz, Ny, Nx);
    fprintf(xmf, "           %.*s:/Data/Cells\n", 1+(int)sizeof filename, filename);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    fprintf(xmf, "     <Attribute Name=\"xx\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n", Nz, Ny, Nx);
    fprintf(xmf, "         <DataItem Dimensions=\"3 4\" Format=\"XML\">\n");
    fprintf(xmf, "           5 0 0 0\n");
    fprintf(xmf, "           1 1 1 1\n");
    fprintf(xmf, "           6 %d %d %d\n", Nz, Ny, Nx);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d %d %d %d\" Format=\"HDF\">\n", Nv, Nz, Ny, Nx);
    fprintf(xmf, "           %.*s:/Data/Cells\n", 1+(int)sizeof filename, filename);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);

}



