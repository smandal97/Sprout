
#include "../defs.h"

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
   int size = theDomain->size;
   if( size!=1 ){ printf("Failure! Exiting..."); exit(1); }
   
   double rho,p,v1,x;
   int i,j,k,ijk;
   char filename[256];
   sprintf(filename,"%s.dat",filestart);
   FILE * op;
   op = fopen(filename,"w");
   fprintf(op,"#t = %.5e\n", theDomain->t);
   fprintf(op,"#x rho p v1\n");
   fclose(op);
   
   printf("Writing checkpoint...\n");
      
   FILE * pFile = fopen( filename , "a" );

   for( k=0 ; k<Nz ; ++k ){
     for( j=0 ; j<Ny ; ++j ){
        for( i=0 ; i<Nx ; ++i ){
           ijk  = i+Ng;
           if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*(j+Ng);
           if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*(k+Ng);
           struct cell * c = theCells + ijk;
           rho = c->prim[RHO]; 
           p = c->prim[PPP]; 
           v1 = c->prim[UU1];
           x = c->xi[0]; 
           fprintf(pFile,"%e %e %e %e\n",x,rho,p,v1);
        }
     }
   }

   fclose(pFile);
         
   

}
