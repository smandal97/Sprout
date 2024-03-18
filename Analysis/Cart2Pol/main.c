
#include "defs.h"

int read_par_file( struct domain * );
void readdata( char * , struct domain * );
void makepoldata( char * , struct domain * );


int main( int argc, char **argv ){

   if( argc < 2 ){
      printf("Please specify the input file.\n");
      exit(1);
   }
   char infile[256];
   if( argv[1] ){
      strcpy( infile , argv[1] );
   }

   struct domain theDomain = {0};

   read_par_file( &theDomain );
   readdata( infile , &theDomain );
   makepoldata( infile , &theDomain );
   return (0);

}
