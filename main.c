
#include "defs.h"


int read_par_file( struct domain * );

int mpiSetup( struct domain * , int , char *[] );
void setupGrid( struct domain * );
void setupDomain( struct domain * );
void setupCells( struct domain * );
void setICParams( struct domain * );
void setHydroParams( struct domain * );
void setGravityParams( struct domain * );
void setNozzleParams( struct domain * );
void setMeshMotionParams( struct domain * );

void exchangeData( struct domain * );
void boundary( struct domain * );

void check_dt( struct domain * , double * );
double getmindt( struct domain * );
void timestep( struct domain * , double );
void printCells( struct domain * );
void output( struct domain * , char * );
void possiblyOutput( struct domain * , double , int );
void freeDomain( struct domain * );

void start_clock( struct domain * );
void generate_log( struct domain * );

int main( int argc , char * argv[] ){

   MPI_Init( &argc , &argv );
   struct domain theDomain = {0};
   start_clock( &theDomain );
   read_par_file( &theDomain );

   int error = mpiSetup(&theDomain,argc,argv);
   if( error ) return(0);

   if(theDomain.rank==0) remove("abort");

   setupGrid( &theDomain ); 
   setupDomain( &theDomain ); 
   setICParams( &theDomain );
   setHydroParams( &theDomain );
   setGravityParams( &theDomain );
   setNozzleParams( &theDomain );
   setMeshMotionParams( &theDomain );
   setupCells( &theDomain );
   
   exchangeData( &theDomain );
   boundary( &theDomain );
   
   double dtf;
     
   while( !(theDomain.final_step) ){

      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      timestep( &theDomain , dt );
      possiblyOutput( &theDomain , dt , 0 );
      dtf = dt;
   }
   /*
   int g;
   for( g=0 ; g<0 ; ++g ){
      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      timestep( &theDomain , dt );
      possiblyOutput( &theDomain , dt , 0 );
      dtf = dt;
   }*/
   
   possiblyOutput( &theDomain , dtf , 1 );
   generate_log( &theDomain );
   MPI_Barrier( theDomain.theComm );
   freeDomain( &theDomain );
   MPI_Finalize();

   return(0);

}
