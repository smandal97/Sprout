
#include "../defs.h"

int set_accuracy(void){
   return 3;
}

double get_monotone_slope( double wm1 , double w0 , double wp1 ){
   if( (wp1-w0)*(w0-wm1)<=0. ){
      return 0.;
   }else{
      double del_w0 = 0.5*(wp1+wm1);
      double diff_w_r = 2.*fabs(wp1-w0);
      double diff_w_l = 2.*fabs(w0-wm1);
      double min = fabs(del_w0);
      if(min>diff_w_r) min = diff_w_r;
      if(min>diff_w_l) min = diff_w_l;
      if(del_w0<0.) min *= -1.;
      return min;
   }
}

double second_derivative( double ym1, double yo , double yp1 , double dx ){
   return (yp1+ym1-2.*yo)/6./dx/dx;
}

double steepen_cd( struct cell * cm2 , struct cell * cm1 , struct cell * c0 , struct cell * cp1 , struct cell * cp2 , double gam , double dl ){

   double pm1,pp1;
   double dm2,dm1,dp1,dp2,d0;
   pm1 = cm1->prim[PPP];
   pp1 = cp1->prim[PPP];
   dm2 = cm2->prim[RHO];
   dm1 = cm1->prim[RHO];
   d0  =  c0->prim[RHO];
   dp1 = cp1->prim[RHO];
   dp2 = cp2->prim[RHO];

   double pmin,dmin;
   if(pp1>pm1)
      pmin = pm1;
   else
      pmin = pp1;
   if(pmin==0) pmin = 1e-20;
   if(dp1>dm1)
      dmin = dm1;
   else
      dmin = dp1;
   if(dmin==0) dmin = 1e-20;

   if( fabs(dp1-dm1)/dmin<=0.01 ) return 0.;

   if( fabs(pp1-pm1)/pmin>0.1*gam*fabs(dp1-dm1)/dmin ) return 0.;

   double d2rho_m1 = second_derivative(dm2,dm1,d0,dl);
   double d2rho_p1 = second_derivative(d0,dp1,dp2,dl);
   if( d2rho_p1*d2rho_m1<0. ) return 0.;
   double eta = (d2rho_m1-d2rho_p1)*dl*dl/(dp1-dm1);
   if(eta>0.1) eta=1.;
   else eta=20.*eta-1.;
   if(eta<0.) eta=0.;

   return eta;

}

double flatten_shock( struct cell * cm2 , struct cell * cm1 , struct cell * c0 , struct cell * cp1 , struct cell * cp2 , int theDIM ){

   double pm2,pm1,pp1,pp2,um1,up1;
   double S,F,pmin;
   pm2 = cm2->prim[PPP];
   pm1 = cm1->prim[PPP];
   pp1 = cp1->prim[PPP];
   pp2 = cp2->prim[PPP];
   um1 = cm1->prim[UU1+theDIM];
   up1 = cp1->prim[UU1+theDIM];

   S = 0.;
   if(pp2!=pm2) S = (pp1-pm1)/(pp2-pm2);

   F = 10.*(S-0.75);
   if(F>1.) F=1.;
   if(F<0.) F=0;

   if( pp1>pm1 ) pmin = pm1;
   else pmin = pp1;
   if(pmin==0.) pmin=1e-20;
   
   if( (fabs(pp1-pm1)/pmin<0.333333333333) || (up1>um1) ) F=0.;

   return F;
}


void parabolic_interpolate( double wm2 , double wm1 , double w0 , double wp1 , double wp2 , double * wl , double * wr ){

   double del_wp0 = get_monotone_slope(wm1,w0,wp1);
   double del_wp1 = get_monotone_slope(w0,wp1,wp2);
   double del_wm1 = get_monotone_slope(wm2,wm1,w0);
   *wl = 0.5*(w0+wp1) - (del_wp1-del_wp0)/6.;
   *wr = 0.5*(w0+wm1) - (del_wp0-del_wm1)/6.;

}

void make_monotone( double w0 , double * wl , double * wr ){

   if( (*wr-w0)*(w0-*wl)<=0. ){
      *wl = w0;
      *wr = w0; 
   }else if( (*wr-*wl)*(*wl-(3*w0-*wr*2.))<0. ){
      *wr = w0;
      *wl = 3.*w0-*wr*2.;
   }else if( (*wr-*wl)*((3*w0-*wl*2.)-*wr)<0. ){
      *wl = 3.*w0-*wl*2.;
      *wr = w0;
   }

}




void space_recon1D( struct domain * theDomain , int theDIM ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;
   double gam = theDomain->theParList.Adiabatic_Index;

   int n[3]  = {0};
   n[theDIM] = 1;
   double p1,p2,p3,p4,p5,p6,p7;
   double eta,drho_m1,drho_p1,Fm1,F0,Fp1,wl,wr;
   double dl = (double)n[0]*dx + (double)n[1]*dy + (double)n[2]*dz;

   int i_end = Nx+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_end = 1;
   int j_end = Ny+2*Ng;
   if( theDomain->theParList.Num_y == 1 ) j_end = 1;
   int k_end = Nz+2*Ng;
   if( theDomain->theParList.Num_z == 1 ) k_end = 1;

   //set loop variables
   //ijk4 is the cell we wanna interpolate for
   int i,j,k,q,ijk1,ijk2,ijk3,ijk4,ijk5,ijk6,ijk7;

   for( k=0 ; k<k_end ; ++k ){
      for( j=0 ; j<j_end ; ++j ){
         for( i=0 ; i<i_end ; ++i ){
            ijk1 = i-3*n[0];
            ijk2 = i-2*n[0];
            ijk3 = i-1*n[0];
            ijk4 = i;
            ijk5 = i+1*n[0];
            ijk6 = i+2*n[0];
            ijk7 = i+3*n[0];
            if( i==0 ){
               ijk1 = i;
               ijk2 = i;
               ijk3 = i;
            }else if( i==1 ){
               ijk1 = i;
               ijk2 = i;
            }else if( i==2 ){
               ijk1 = i;
            }else if( i==i_end-3 ){
               ijk7 = i_end-1;
            }else if( i==i_end-2 ){
               ijk6 = i_end-1;
               ijk7 = i_end-1;
            }else if( i==i_end-1 ){
               ijk5 = i_end-1;
               ijk6 = i_end-1;
               ijk7 = i_end-1;
            }
            if( theDomain->theParList.Num_y != 1 ){
               ijk4 += (Nx+2*Ng)*j;
               if( j==0 ){
                  ijk5 += (Nx+2*Ng)*(j+1*n[1]);
                  ijk6 += (Nx+2*Ng)*(j+2*n[1]);
                  ijk7 += (Nx+2*Ng)*(j+3*n[1]);
               }else if( j==1 ){
                  ijk5 += (Nx+2*Ng)*(j+1*n[1]);
                  ijk6 += (Nx+2*Ng)*(j+2*n[1]);
                  ijk7 += (Nx+2*Ng)*(j+3*n[1]);
               }else if( j==2 ){
                  ijk3 += (Nx+2*Ng)*(j-1*n[1]);
                  ijk5 += (Nx+2*Ng)*(j+1*n[1]);
                  ijk6 += (Nx+2*Ng)*(j+2*n[1]);
                  ijk7 += (Nx+2*Ng)*(j+3*n[1]);
               }else if( j==j_end-3 ){
                  ijk1 += (Nx+2*Ng)*(j-3*n[1]);
                  ijk2 += (Nx+2*Ng)*(j-2*n[1]);
                  ijk3 += (Nx+2*Ng)*(j-1*n[1]);
                  ijk5 += (Nx+2*Ng)*(j+1*n[1]);
                  ijk6 += (Nx+2*Ng)*(j+2*n[1]);
                  ijk7 += (Nx+2*Ng)*(j_end-1);
               }else if( j==j_end-2 ){
                  ijk1 += (Nx+2*Ng)*(j-3*n[1]);
                  ijk2 += (Nx+2*Ng)*(j-2*n[1]);
                  ijk3 += (Nx+2*Ng)*(j-1*n[1]);
                  ijk5 += (Nx+2*Ng)*(j+1*n[1]);
                  ijk6 += (Nx+2*Ng)*(j_end-1);
                  ijk7 += (Nx+2*Ng)*(j_end-1);
               }else if( j==j_end-1 ){
                  ijk1 += (Nx+2*Ng)*(j-3*n[1]);
                  ijk2 += (Nx+2*Ng)*(j-2*n[1]);
                  ijk3 += (Nx+2*Ng)*(j-1*n[1]);
                  ijk5 += (Nx+2*Ng)*(j_end-1);
                  ijk6 += (Nx+2*Ng)*(j_end-1);
                  ijk7 += (Nx+2*Ng)*(j_end-1);
               }else{
                  ijk1 += (Nx+2*Ng)*(j-3*n[1]);
                  ijk2 += (Nx+2*Ng)*(j-2*n[1]);
                  ijk3 += (Nx+2*Ng)*(j-1*n[1]);
                  ijk5 += (Nx+2*Ng)*(j+1*n[1]);
                  ijk6 += (Nx+2*Ng)*(j+2*n[1]);
                  ijk7 += (Nx+2*Ng)*(j+3*n[1]);
               }
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk4 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               if( k==0 ){
                  ijk5 += (Nx+2*Ng)*(Ny+2*Ng)*(k+1*n[2]);
                  ijk6 += (Nx+2*Ng)*(Ny+2*Ng)*(k+2*n[2]);
                  ijk7 += (Nx+2*Ng)*(Ny+2*Ng)*(k+3*n[2]);
               }else if( k==1 ){
                  ijk5 += (Nx+2*Ng)*(Ny+2*Ng)*(k+1*n[2]);
                  ijk6 += (Nx+2*Ng)*(Ny+2*Ng)*(k+2*n[2]);
                  ijk7 += (Nx+2*Ng)*(Ny+2*Ng)*(k+3*n[2]);
               }else if( k==2 ){
                  ijk3 += (Nx+2*Ng)*(Ny+2*Ng)*(k-1*n[2]);
                  ijk5 += (Nx+2*Ng)*(Ny+2*Ng)*(k+1*n[2]);
                  ijk6 += (Nx+2*Ng)*(Ny+2*Ng)*(k+2*n[2]);
                  ijk7 += (Nx+2*Ng)*(Ny+2*Ng)*(k+3*n[2]);
               }else if( k==k_end-3 ){
                  ijk1 += (Nx+2*Ng)*(Ny+2*Ng)*(k-3*n[2]);
                  ijk2 += (Nx+2*Ng)*(Ny+2*Ng)*(k-2*n[2]);
                  ijk3 += (Nx+2*Ng)*(Ny+2*Ng)*(k-1*n[2]);
                  ijk5 += (Nx+2*Ng)*(Ny+2*Ng)*(k+1*n[2]);
                  ijk6 += (Nx+2*Ng)*(Ny+2*Ng)*(k+2*n[2]);
                  ijk7 += (Nx+2*Ng)*(Ny+2*Ng)*(k_end-1);
               }else if( k==k_end-2 ){
                  ijk1 += (Nx+2*Ng)*(Ny+2*Ng)*(k-3*n[2]);
                  ijk2 += (Nx+2*Ng)*(Ny+2*Ng)*(k-2*n[2]);
                  ijk3 += (Nx+2*Ng)*(Ny+2*Ng)*(k-1*n[2]);
                  ijk5 += (Nx+2*Ng)*(Ny+2*Ng)*(k+1*n[2]);
                  ijk6 += (Nx+2*Ng)*(Ny+2*Ng)*(k_end-1);
                  ijk7 += (Nx+2*Ng)*(Ny+2*Ng)*(k_end-1);
               }else if( k==k_end-1 ){
                  ijk1 += (Nx+2*Ng)*(Ny+2*Ng)*(k-3*n[2]);
                  ijk2 += (Nx+2*Ng)*(Ny+2*Ng)*(k-2*n[2]);
                  ijk3 += (Nx+2*Ng)*(Ny+2*Ng)*(k-1*n[2]);
                  ijk5 += (Nx+2*Ng)*(Ny+2*Ng)*(k_end-1);
                  ijk6 += (Nx+2*Ng)*(Ny+2*Ng)*(k_end-1);
                  ijk7 += (Nx+2*Ng)*(Ny+2*Ng)*(k_end-1);
               }else{
                  ijk1 += (Nx+2*Ng)*(Ny+2*Ng)*(k-3*n[2]);
                  ijk2 += (Nx+2*Ng)*(Ny+2*Ng)*(k-2*n[2]);
                  ijk3 += (Nx+2*Ng)*(Ny+2*Ng)*(k-1*n[2]);
                  ijk5 += (Nx+2*Ng)*(Ny+2*Ng)*(k+1*n[2]);
                  ijk6 += (Nx+2*Ng)*(Ny+2*Ng)*(k+2*n[2]);
                  ijk7 += (Nx+2*Ng)*(Ny+2*Ng)*(k+3*n[2]);
               }
            }
            struct cell * c1 = theCells+ijk1;
            struct cell * c2 = theCells+ijk2;
            struct cell * c3 = theCells+ijk3;
            struct cell * c4 = theCells+ijk4;
            struct cell * c5 = theCells+ijk5;
            struct cell * c6 = theCells+ijk6;
            struct cell * c7 = theCells+ijk7;
            
            //get shock flattening parameter
            Fm1 = flatten_shock(c1,c2,c3,c4,c5,theDIM);
            F0  = flatten_shock(c2,c3,c4,c5,c6,theDIM);
            Fp1 = flatten_shock(c3,c4,c5,c6,c7,theDIM);
            //printf("Fp1=%e\n",Fp1);
            if( c5->prim[PPP]-c3->prim[PPP]<0 ){
               if(Fp1>F0) F0 = Fp1;
            }else{
               if(Fm1>F0) F0 = Fm1;
            }
            
            //get contact steepening parameter
            drho_m1 = get_monotone_slope(c2->prim[RHO],c3->prim[RHO],c4->prim[RHO]);
            drho_p1 = get_monotone_slope(c4->prim[RHO],c5->prim[RHO],c6->prim[RHO]);
            eta = steepen_cd(c2,c3,c4,c5,c6,gam,dl);

            for( q=0 ; q<NUM_Q ; ++q ){
               p1 = c1->prim[q];
               p2 = c2->prim[q];
               p3 = c3->prim[q];
               p4 = c4->prim[q];
               p5 = c5->prim[q];
               p6 = c6->prim[q];
               p7 = c7->prim[q];

               parabolic_interpolate( p2,p3,p4,p5,p6,&wl,&wr );
               /*//do contact steepening
               if( q==0 ){
                  wl = wl*(1.-eta) + (p3+0.5*drho_m1)*eta;
                  wr = wr*(1.-eta) + (p3-0.5*drho_p1)*eta;
               }*/   
               //shock flattening
               wl = F0*p4 + (1.-F0)*wl;
               wr = F0*p4 + (1.-F0)*wr;
               make_monotone( p4,&wl,&wr );

               if( theDIM==0 ){
                  c4->gradx[q] = (wr-wl)/dl;
                  c4->pblax[q] = (wr+wl-2.*p4)/2.;
               }
               if( theDIM==1 ){
                  c4->grady[q] = (wr-wl)/dl;
                  c4->pblay[q] = (wr+wl-2.*p4)/2.;
               }
               if( theDIM==2 ){
                  c4->gradz[q] = (wr-wl)/dl;
                  c4->pblaz[q] = (wr+wl-2.*p4)/2.;
               }

            }

         }
      }
   }

}
