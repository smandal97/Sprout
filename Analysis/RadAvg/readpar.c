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

int read_par_file( struct domain * Data ){

   struct parlist * theList = &( Data->theParList );

   char pfile[] = "in.par";

   int err=0;  

   err += readvar( pfile , "Num_r"     , VAR_INT  , &(theList->Nr)      );
   err += readvar( pfile , "Num_theta" , VAR_INT  , &(theList->Nt)      );
   err += readvar( pfile , "Num_phi"   , VAR_INT  , &(theList->Np)      );
   err += readvar( pfile , "Theta_min" , VAR_DOUB , &(theList->th_min)  );
   err += readvar( pfile , "Theta_max" , VAR_DOUB , &(theList->th_max)  );
   err += readvar( pfile , "Phi_min"   , VAR_DOUB , &(theList->ph_min)  );
   err += readvar( pfile , "Phi_max"   , VAR_DOUB , &(theList->ph_max)  );
   err += readvar( pfile , "x0"        , VAR_DOUB , &(theList->x0)      );
   err += readvar( pfile , "y0"        , VAR_DOUB , &(theList->y0)      );
   err += readvar( pfile , "z0"        , VAR_DOUB , &(theList->z0)      );

   if( err > 0 ){
      printf("Read Failed, err = %d\n",err);
      return(1);
   }

   return(0);

}
