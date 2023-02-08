
#include "defs.h"

void start_clock( struct domain * theDomain ){
	theDomain->Wallt_init = time(NULL);
}


void generate_log( struct domain * theDomain ){
	time_t endtime = time(NULL);
	int seconds = (int) (endtime-theDomain->Wallt_init);

	/*int Ng = theDomain->Ng;
	int Nx = theDomain->Nx;
	int Ny = theDomain->Ny;
	int Nz = theDomain->Nz;
	int Nc = 1;
	if(theDomain->theParList.Num_x!=1) Nc *= Nx; //+2*Ng;
	if(theDomain->theParList.Num_y!=1) Nc *= Ny; //+2*Ng;
	if(theDomain->theParList.Num_z!=1) Nc *= Nz; //+2*Ng;*/

	int Nc = theDomain->theParList.Num_x*theDomain->theParList.Num_y*theDomain->theParList.Num_z;
	int Nt = theDomain->count_steps;

	double avgdt = (double)seconds/2./(double)Nc/(double)Nt;

	int size = theDomain->size;

    if( theDomain->rank==0 ){
      FILE * logfile = fopen("times.log","w");
      fprintf(logfile,"Run using %d MPI process",size);
      if( theDomain->size > 1 ) fprintf(logfile,"es");
      fprintf(logfile,".\n");
      fprintf(logfile,"Total time = %d sec\n",seconds);
      fprintf(logfile,"Number of cells = %d\n",Nc);
      fprintf(logfile,"Number of timesteps = %d (x%d)\n",Nt,2);
      fprintf(logfile,"Megazones per second = %.2e\n",1./(avgdt*1e6));
      fprintf(logfile,"Megazones per CPU second = %.2e\n",1./(avgdt*1e6*size));
      fprintf(logfile,"Time/zone/step = %.2e microseconds\n",(avgdt*1e6));
      fclose(logfile);
    }
}