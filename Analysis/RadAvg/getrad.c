
#include "defs.h"



void get_rad_avg( char * filename , struct domain * theDomain ){
   
   int Nt,Np,Nr,Nx,Ny,Nz;
   double tmin,tmax,pmin,pmax;
   double r, ph, th, rho, P, vr, KE_turb;

   
   Nx = theDomain->Nx;
   Ny = theDomain->Ny;
   Nz = theDomain->Nz;
   Nr = theDomain->theParList.Nr;
   Nt = theDomain->theParList.Nt;
   Np = theDomain->theParList.Np;
   tmin = theDomain->theParList.th_min;
   tmax = theDomain->theParList.th_max;
   pmin = theDomain->theParList.ph_min;
   pmax = theDomain->theParList.ph_max;


   double i0, j0, k0;
   i0 = (int)(theDomain->Nx * theDomain->theParList.x0);
   j0 = (int)(theDomain->Ny * theDomain->theParList.y0);
   k0 = (int)(theDomain->Nz * theDomain->theParList.z0);

   FILE * fp;
   fp = fopen("rad_profiles.dat", "w");
   fprintf(fp, "#r rho P vr KE_turb");

   int l,m,n;
   int q,i,j,k,qijk;
   for( l=0 ; l<Nr-1 ; ++l ){
      rho     = 0.;
      P       = 0.;
      vr      = 0.;
      KE_turb = 0.;
      for( m=0 ; m<Np ; ++m ){
         for( n=0 ; n<Nt ; ++n ){
            for( q=0 ; q<NUM_Q ; ++q ){
               r  = (double)l*(double)(Nx-i0)/(double)(Nr-1);
               ph = (double)m * (pmax-pmin) / (double)(Np-1);
               th = (double)n * (tmax-tmin) / (double)(Nt-1);

               //convert to i,j,k and adjust for origin of spherical grid w.r.t. Cartesian grid
               i = (int)(r*sin(th)*cos(ph) + i0);
               j = (int)(r*sin(th)*sin(ph) + j0);
               k = (int)(r*cos(th) + k0);

               //get values of variables from Cartesian grid and populate radial data arrays
               qijk = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*q;
               if(q==RHO)
                  rho += theDomain->cart_data[qijk];
               else if(q==PPP)
                  P += theDomain->cart_data[qijk];
               else if(q==UU1){
                  vr      += theDomain->cart_data[qijk] * sin(th) * cos(ph);
                  KE_turb += theDomain->cart_data[qijk] * theDomain->cart_data[qijk];
               }
               else if(q==UU2){
                  vr      += theDomain->cart_data[qijk] * sin(th) * sin(ph);
                  KE_turb += theDomain->cart_data[qijk] * theDomain->cart_data[qijk];
               }
               else if(q==UU3){
                  vr      += theDomain->cart_data[qijk] * cos(th);
                  KE_turb += theDomain->cart_data[qijk] * theDomain->cart_data[qijk];
               }
               

            }
         }
      }
      r *= theDomain->dx;
      rho /= Nt*Np;
      P   /= Nt*Np;
      vr  /= Nt*Np;
      KE_turb /= Nt*Np;
      KE_turb -= vr*vr;
      //write to file
      fprintf(fp, "%e %e %e %e %e\n", r, rho, P, vr, KE_turb);
   }

   fclose(fp);
   free(theDomain->cart_data);



}

