
#include "defs.h"

double getmin( double a , double b , double c , int Na , int Nb , int Nc ){
   a *= (double)Na;
   b *= (double)Nb;
   c *= (double)Nc;
   if(a>b)a=b;
   if(a>c)a=c;
   return(a);
}

void custom_vals( struct domain * theDomain , int i , int j , int k , double th, double ph , double * op1 , double * op2 , double * op3 ){
   double * cartdata = theDomain->cartdata;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int v1,v2,v3,v4,vijk;

   //this is where you CHOOSE WHAT POLMAP STORES
   v1 = 2;
   v2 = 3;
   v3 = 4;
   v4 = 1;

   double vx,vy,vz,P,rho;
   if(i<0 || j<0 || k<0 || i>(Nx-1) || j>(Ny-1) || k>(Nz-1)){
      vx  = 0.;
      vy  = 0.;
      vz  = 0.;
      P   = 0.;
      rho = 0.;
   }
   else{
      vijk = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*v1;
      vx = cartdata[vijk];
      vijk = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*v2;
      vy = cartdata[vijk];
      vijk = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*v3;
      vz = cartdata[vijk];
      vijk = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*v4;
      P = cartdata[vijk];
      vijk = i + Nx*j + Nx*Ny*k;
      rho = cartdata[vijk];
   }
   
   
   *op1 = P*rho; //( vx*sin(th)*cos(ph) + vy*sin(th)*sin(ph) + vz*cos(th) );
   *op2 = P;
   *op3 = (*op1)/(*op2);
}

/*void change_extension( char * myFilename , char * newFilename , int flag ){
   char * ptrFile = strrchr(myFilename, '/');
   ptrFile = (ptrFile) ? myFilename : ptrFile+1;

   char ext[8];
   if(flag==0)
      ext = "_mtd.dat";
   else if(flag==1)
      ext = "_grd.dat";
   else
      ext = "_map.dat";

   char * ptrExt = strrchr(ptrFile, '.');
   if (ptrExt != NULL)
      strcpy(ptrExt, ext);
   else
      strcat(ptrFile, ext);
}*/


void makepoldata( char * filename , struct domain * theDomain ){
   
   int Nt,Np,Nr,Nx,Ny,Nz;
   double x_0,y_0,z_0,dx,dy,dz,xmin,ymin,zmin;
   double rmin,rmax,tmin,tmax,pmin,pmax,eta_min,eta_max;

   dx = theDomain->dx;
   dy = theDomain->dy;
   dz = theDomain->dz;
   Nx = theDomain->Nx;
   Ny = theDomain->Ny;
   Nz = theDomain->Nz;
   xmin = theDomain->xmin;
   ymin = theDomain->ymin;
   zmin = theDomain->zmin;
   Nr = theDomain->theParList.Nr;
   Nt = theDomain->theParList.Nt;
   Np = theDomain->theParList.Np;
   x_0  = theDomain->theParList.x0;
   y_0  = theDomain->theParList.y0;
   z_0  = theDomain->theParList.z0;
   tmin = theDomain->theParList.th_min;
   tmax = theDomain->theParList.th_max;
   pmin = theDomain->theParList.ph_min;
   pmax = theDomain->theParList.ph_max;
   eta_min = theDomain->theParList.eta_min;
   eta_max = theDomain->theParList.eta_max;
   rmin = getmin(dx,dy,dz,1,1,1);
   rmax = getmin(dx,dy,dz,Nx,Ny,Nz) - x_0*0.;

   printf("pmin = %e, pmax = %e\n",pmin,pmax);
   printf("tmin = %e, tmax = %e\n",tmin,tmax);
   int zflag = 1;	
   if(Nz==1){
      zflag=0;
      Nt=1;
   }

    theDomain->polgrid = (double *)malloc( Nr*Nt*Np*sizeof(double) );
    theDomain->polmap1 = (double *)calloc( Nt*Np,sizeof(double) );
    theDomain->polmap2 = (double *)calloc( Nt*Np,sizeof(double) );

    /*int a1;
    FILE * fsl = fopen("slice.dat", "w");
    for( a1=67108864 ; a1<67371008 ; ++a1 )
       fprintf(fsl, "%e\n", theDomain->cartdata[a1]);
    fclose(fsl);*/

    int l,p,q,pq,lpq,v,i,j,k,vijk;
    double r, th, ph, x, y, z, op1, op2, op3, norm;
    norm = 0.0;
    for( q=0 ; q<Nt ; ++q ){
       for( p=0 ; p<Np ; ++p ){
          for( l=0 ; l<Nr ; ++l ){
             pq  = p + Np*q;
             lpq = l + Nr*p + Nr*Np*q;
             r   = rmin + (double)l*(rmax-rmin)/(double)(Nr-1);
             ph  = pmin + (double)p*(pmax-pmin)/(double)(Np-1);
             th  = tmin + (double)q*(tmax-tmin)/(double)(Nt-1);
             if(!zflag) th = M_PI/2.;

             //convert to x,y,z
             x = r*sin(th)*cos(ph);
             y = r*sin(th)*sin(ph);
             z = r*cos(th);

             //adjust for mesh motion and chosen center
             x -= xmin+x_0;
             y -= ymin+y_0;
             z -= zmin+z_0;

             //change to i,j,k with attention to zflag
             i = (int)(x/dx);
             j = (int)(y/dy);
             k = 0;
             if(zflag) k = (int)(z/dz);
            
             
             //CHOOSEE PRIMITIVE VARIABLE FOR POLGRID
             v = 0;
             vijk = i + Nx*j + Nx*Ny*k + Nx*Ny*Nz*v;
             
             //populate polgrid
             if(i<0 || j<0 || k<0 || i>(Nx-1) || j>(Ny-1) || k>(Nz-1))
                theDomain->polgrid[lpq] = 0.;   //if i,j,k out of bounds
             else
                theDomain->polgrid[lpq] = theDomain->cartdata[vijk];

             //populate polmaps
             custom_vals( theDomain , i , j , k , th , ph , &op1 , &op2 , &op3 );
             //if( r>rmin+eta_min*(rmax-rmin) || r<rmin+eta_max*(rmax-rmin) ){
             if( r>eta_min*rmax || r<eta_max*rmax ){
                theDomain->polmap1[pq] += op1;
                theDomain->polmap2[pq] += op2;
                norm += op3;
             }
          }
       }
    }

    norm /= (double)(Nt*Np);

    //char newfilename[strlen(filename)+8];
    FILE * fp1, * fp2;
    //write metadata
    fp1 = fopen("metadata.dat", "w");
    fprintf(fp1, "%i %e %e\n%i %e %e\n%i %e %e", Nr, rmin, rmax, Np, pmin, pmax, Nt, tmin, tmax);
    fclose(fp1);
    //write polgrid
    fp1 = fopen("polgrid.dat", "w");
    fp2 = fopen("dmap_.dat", "w");
    for( q=0 ; q<Nt ; ++q ){
       for( p=0 ; p<Np ; ++p ){
          pq  = p + Np*q;
          for( l=0 ; l<Nr ; ++l ){
             lpq = l + Nr*p + Nr*Np*q;
             fprintf(fp1, "%e\n", theDomain->polgrid[lpq]);
          }
          fprintf(fp2, "%e\n", theDomain->polmap1[pq]/theDomain->polmap2[pq]/norm );   
       }
    }

    
    fclose(fp1);
    fclose(fp2);

    //free the domain
    free(theDomain->polgrid);
    free(theDomain->cartdata);
    free(theDomain->polmap1);
    free(theDomain->polmap2);

}

