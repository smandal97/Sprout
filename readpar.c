enum{VAR_INT,VAR_DOUB,VAR_STR};

#include "defs.h"
#include <string.h>

int readvar( char * filename , char * varname , int vartype , void * ptr ){

   FILE * inFile = fopen( filename , "r" );
   char s[512];
   char nm[512];
   char s1[512];
   int found = 0;
   
   while( (fgets(s,512,inFile) != NULL) && found==0 ){
      sscanf(s,"%s ",nm);
      if( strcmp(nm,varname)==0 ){
         strcpy(s1,s);
         found=1;
      }
   }
   
   fclose( inFile );
   if( found==0 ){ printf("%s not found\n",varname); return(1); }

   char * s2 = s1+strlen(nm)+strspn(s1+strlen(nm),"\t :=>_");

   double temp;
   char stringval[256];

   sscanf(s2,"%lf",&temp);
   sscanf(s2,"%256s",stringval);

   if( vartype == VAR_INT ){
      *((int *)   ptr) = (int)temp;
   }else if( vartype == VAR_DOUB ){
      *((double *)ptr) = (double)temp;
   }else{
      strcpy( ptr , stringval );
   }

   return(0);
}

int read_par_file( struct domain * theDomain ){

   MPI_Comm_rank(MPI_COMM_WORLD,&(theDomain->rank));
   MPI_Comm_size(MPI_COMM_WORLD,&(theDomain->size));

   int rank = theDomain->rank;
   int size = theDomain->size;

   struct param_list * theList = &( theDomain->theParList );

   char pfile[] = "in.par";

   int err=0;  

   int nrank;
   for( nrank=0 ; nrank<size ; ++nrank ){
      if( rank==nrank ){
         err += readvar( pfile , "Num_X"           , VAR_INT  , &(theList->Num_x)           );
         err += readvar( pfile , "Num_Y"           , VAR_INT  , &(theList->Num_y)           );
         err += readvar( pfile , "Num_Z"           , VAR_INT  , &(theList->Num_z)           );
         err += readvar( pfile , "Num_Reports"     , VAR_INT  , &(theList->NumRepts)        );
         err += readvar( pfile , "Num_Snapshots"   , VAR_INT  , &(theList->NumSnaps)        );
         err += readvar( pfile , "Num_Checkpoints" , VAR_INT  , &(theList->NumChecks)       );
         err += readvar( pfile , "Use_Logtime"     , VAR_INT  , &(theList->Out_LogTime)     );
         err += readvar( pfile , "T_Start"         , VAR_DOUB , &(theList->t_min)           );
         err += readvar( pfile , "T_End"           , VAR_DOUB , &(theList->t_max)           );
         err += readvar( pfile , "Lx"              , VAR_DOUB , &(theList->Lx)              );
         err += readvar( pfile , "Ly"              , VAR_DOUB , &(theList->Ly)              );
         err += readvar( pfile , "Lz"              , VAR_DOUB , &(theList->Lz)              );
         err += readvar( pfile , "CFL"             , VAR_DOUB , &(theList->CFL)             );
         err += readvar( pfile , "PLM"             , VAR_DOUB , &(theList->PLM)             );
         err += readvar( pfile , "W0"              , VAR_DOUB , &(theList->W0)              );
         err += readvar( pfile , "W_frac"          , VAR_DOUB , &(theList->W_frac)          );
         err += readvar( pfile , "Mesh_motion_x0"  , VAR_DOUB , &(theList->MM_x0)           );
         err += readvar( pfile , "Mesh_motion_y0"  , VAR_DOUB , &(theList->MM_y0)           );
         err += readvar( pfile , "Mesh_motion_z0"  , VAR_DOUB , &(theList->MM_z0)           );
         err += readvar( pfile , "Adiabatic_Index" , VAR_DOUB , &(theList->Adiabatic_Index) );
         err += readvar( pfile , "Density_Floor"   , VAR_DOUB , &(theList->Density_Floor)   );
         err += readvar( pfile , "Pressure_Floor"  , VAR_DOUB , &(theList->Pressure_Floor)  );
         err += readvar( pfile , "Gravity_Switch"  , VAR_INT  , &(theList->Gravity_Switch)  );
         err += readvar( pfile , "Central_Mass"    , VAR_DOUB , &(theList->Central_Mass)    );
         err += readvar( pfile , "Grav_G"          , VAR_DOUB , &(theList->Grav_G)          );
         err += readvar( pfile , "Bondi_Mdot"      , VAR_DOUB , &(theList->Mdot)            );
         err += readvar( pfile , "Bondi_P_coeff"   , VAR_DOUB , &(theList->eta_P)           );
         err += readvar( pfile , "Restart"         , VAR_INT  , &(theList->restart_flag)    );
         err += readvar( pfile , "Nozzle_Switch"   , VAR_INT  , &(theList->Nozzle_Switch)   );
         err += readvar( pfile , "Nozzle_x0"       , VAR_DOUB , &(theList->Nozzle_x0)       );
         err += readvar( pfile , "Nozzle_y0"       , VAR_DOUB , &(theList->Nozzle_y0)       );
         err += readvar( pfile , "Nozzle_z0"       , VAR_DOUB , &(theList->Nozzle_z0)       );
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   int errtot;
   MPI_Allreduce( &err , &errtot , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );

   if( errtot > 0 ){
      printf("Read Failed, err = %d\n",err);
      return(1);
   }

   return(0);

}
