enum{RHO,PPP,UU1,UU2,UU3,XXX,XPS,YPS,ZPS};

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define NUM_C 5
#define NUM_N 1
#define NUM_Q (NUM_C+NUM_N)
#define NUM_G 2

struct parlist{
   int Nr, Nt, Np;
   double th_min, th_max;
   double ph_min, ph_max;
   double x0, y0, z0;
   double eta_min, eta_max;
};


struct domain{
   double * polgrid;
   double * cartdata;
   double * polmap1;
   double * polmap2;
   
   double dx, dy, dz;
   int Nx, Ny, Nz;
   double xmin, ymin, zmin;

   struct parlist theParList;
};
