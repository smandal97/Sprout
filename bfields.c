
#include "defs.h"

//use set_B_flag function to make sure not doing mhd with euler.c

/*in onestep.c, do this:
	if( bflag && NUM_EDGES >= 4 ){
      avg_Efields( theDomain );
      subtract_advective_B_fluxes( theDomain );
      update_B_fluxes( theDomain , dt );
   }
	..............

   if( bflag && theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , 1 );
   }
*/



void get_cell_Efields( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;

   int ndims = 3;
   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ){ 
      i_index = 1;
      ndims -= 1;
   }
   if( theDomain->theParList.Num_y == 1 ){
      j_index = 1;
      ndims -= 1;
   }
   if( theDomain->theParList.Num_z == 1 ){
      k_index = 1;
      ndims -= 1;
   }
   
   int i,j,k,ijk,theDIM;
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            ijk  = i;
            if( theDomain->theParList.Num_y != 1 ) ijk += (Nx+2*Ng)*j;
            if( theDomain->theParList.Num_z != 1 ) ijk += (Nx+2*Ng)*(Ny+2*Ng)*k;
            struct cell * c = theCells+ijk;
            for( theDIM=0 ; theDIM<3 ; ++theDIM ){
               int DIM_p1 = (theDIM+1)%3;
               int DIM_p2 = (theDIM+2)%3;
               //c->E_cntr[theDIM] = -c->prim[UU1+DIM_p1]*c->prim[BB1+DIM_p2] + c->prim[UU1+DIM_p2]*c->prim[BB1+DIM_p1];
               c->E_cntr[DIM_p2] = -c->prim[UU1+theDIM]*c->prim[BB1+DIM_p1] + c->prim[UU1+DIM_p1]*c->prim[BB1+theDIM];
            }
         }
      }   
   }

}



void edge_Efields_1D( struct domain * theDomain , int theDIM ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;

   int i_index = Nx+2*Ng-1;
   int j_index = Ny+2*Ng-1;
   int k_index = Nz+2*Ng-1;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;

   int n[3] = {0}; n[theDIM] = 1;
   int DIM_p1 = (theDIM+1)%3; int n_p1[3] = {0}; n_p1[DIM_p1] = 1;
   int DIM_p2 = (theDIM+2)%3; int n_p2[3] = {0}; n_p2[DIM_p2] = 1;


   int i,j,k,ijk_0,ijk_1,ijk_2,ijk_3,ijk_a,ijk_b,ijk_c,ijk_d,ijk_e,ijk_f;     
   //convention: when traversing i-th direcion, 0:ijk, 1:(i+1)jk, 2:(i+1)(j+1)k, 3:(i+1)j(k+1), and
   //A: i(j+1)k, B: ij(k+1), C: i(j-1)k, D: (i+1)(j-1)k, E: ij(k-1), F: (i+1)j(k-1)
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            ijk_1 = i+n[0];
            ijk_2 = i+n[0]+n_p1[0];
            ijk_3 = i+n[0]+n_p2[0];
            ijk_0 = i;
            ijk_a = i+n_p1[0];
            ijk_b = i+n_p2[0];
            if( n_p1[0]*i+n_p1[1]*j+n_p1[2]*k ){ 
               ijk_c = i-n_p1[0];
               ijk_d = i+n[0]-n_p1[0];
            }else{ 
               ijk_c = i;
               ijk_d = i+n[0];
            }
            if( n_p2[0]*i+n_p2[1]*j+n_p2[2]*k ){ 
               ijk_e = i-n_p2[0];
               ijk_f = i+n[0]-n_p2[0];
            }else{ 
               ijk_e = i;
               ijk_f = i+n[0];
            }
            //if(i==10 && j==5 && k==5) printf("ia = %i, ib = %i, ic = %i, id = %i, ie = %i, if = %i\n", i+n_p1[0], i+n_p2[0], i-n_p1[0], i+n[0]-n_p1[0], i-n_p2[0], i+n[0]-n_p2[0] );
            if( theDomain->theParList.Num_y != 1 ){
               ijk_1 += (Nx+2*Ng)*(j+n[1]);
               ijk_2 += (Nx+2*Ng)*(j+n[1]+n_p1[1]);
               ijk_3 += (Nx+2*Ng)*(j+n[1]+n_p2[1]);
               ijk_0 += (Nx+2*Ng)*j;
               ijk_a += (Nx+2*Ng)*(j+n_p1[1]);
               ijk_b += (Nx+2*Ng)*(j+n_p2[1]);;
               if( n_p1[0]*i+n_p1[1]*j+n_p1[2]*k ){ 
                  ijk_c += (Nx+2*Ng)*(j-n_p1[1]);
                  ijk_d += (Nx+2*Ng)*(j+n[1]-n_p1[1]);
               }else{ 
                  ijk_c += (Nx+2*Ng)*j;
                  ijk_d += (Nx+2*Ng)*(j+n[1]);
               }
               if( n_p2[0]*i+n_p2[1]*j+n_p2[2]*k ){ 
                  ijk_e += (Nx+2*Ng)*(j-n_p2[1]);
                  ijk_f += (Nx+2*Ng)*(j+n[1]-n_p2[1]);
               }else{ 
                  ijk_e += (Nx+2*Ng)*j;
                  ijk_f += (Nx+2*Ng)*(j+n[1]);
               }
               //if(i==10 && j==5 && k==5) printf("ja = %i, jb = %i, jc = %i, jd = %i, je = %i, jf = %i\n", j+n_p1[1], j+n_p2[1], j-n_p1[1], j+n[1]-n_p1[1], j-n_p2[1], j+n[1]-n_p2[1] );
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk_1 += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]);
               ijk_2 += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]+n_p1[2]);
               ijk_3 += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]+n_p2[2]);
               ijk_0 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_a += (Nx+2*Ng)*(Ny+2*Ng)*(k+n_p1[2]);
               ijk_b += (Nx+2*Ng)*(Ny+2*Ng)*(k+n_p2[2]);;
               if( n_p1[0]*i+n_p1[1]*j+n_p1[2]*k ){ 
                  ijk_c += (Nx+2*Ng)*(Ny+2*Ng)*(k-n_p1[2]);
                  ijk_d += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]-n_p1[2]);
               }else{ 
                  ijk_c += (Nx+2*Ng)*(Ny+2*Ng)*k;
                  ijk_d += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]);
               }
               if( n_p2[0]*i+n_p2[1]*j+n_p2[2]*k ){ 
                  ijk_e += (Nx+2*Ng)*(Ny+2*Ng)*(k-n_p2[2]);
                  ijk_f += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]-n_p2[2]);
               }else{ 
                  ijk_e += (Nx+2*Ng)*(Ny+2*Ng)*k;
                  ijk_f += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]);
               }
               //if(i==10 && j==5 && k==5) printf("ka = %i, kb = %i, kc = %i, kd = %i, ke = %i, kf = %i\n", k+n_p1[2], k+n_p2[2], k-n_p1[2], k+n[2]-n_p1[2], k-n_p2[2], k+n[2]-n_p2[2] );
            }
            struct cell * c1 = theCells+ijk_1;
            struct cell * c2 = theCells+ijk_2;
            struct cell * c3 = theCells+ijk_3;
            //get spatial averaged E-fields
            c1->E_edge[DIM_p2] += c1->Phi_R[theDIM*(NUM_M-1)+1]/4.;
            c2->E_edge[DIM_p2] += c1->Phi_R[theDIM*(NUM_M-1)+1]/4.;
            c1->E_edge[DIM_p1] += c1->Phi_R[theDIM*(NUM_M-1)+0]/4.;
            c3->E_edge[DIM_p1] += c1->Phi_R[theDIM*(NUM_M-1)+0]/4.;

            //upwind E-fields            
            struct cell * c0 = theCells+ijk_0;
            struct cell * cA = theCells+ijk_a;
            struct cell * cB = theCells+ijk_b;
            struct cell * cC = theCells+ijk_c;
            struct cell * cD = theCells+ijk_d;
            struct cell * cE = theCells+ijk_e;
            struct cell * cF = theCells+ijk_f;
            //*
            if(c1->CMode[theDIM]>0.){
               c2->E_edge[DIM_p2] += (c0->E_cntr[DIM_p2] - cA->Phi_R[theDIM*(NUM_M-1)+1])/4.;
               c1->E_edge[DIM_p2] += (c0->E_cntr[DIM_p2] - c0->Phi_R[theDIM*(NUM_M-1)+1])/4.;
               c3->E_edge[DIM_p1] += (c0->E_cntr[DIM_p1] - cB->Phi_R[theDIM*(NUM_M-1)+0])/4.;
               c1->E_edge[DIM_p1] += (c0->E_cntr[DIM_p1] - c0->Phi_R[theDIM*(NUM_M-1)+0])/4.;
            }
            else if(c1->CMode[theDIM]<0.){
               c2->E_edge[DIM_p2] += (c1->E_cntr[DIM_p2] - c2->Phi_R[theDIM*(NUM_M-1)+1])/4.;   //same as above
               c1->E_edge[DIM_p2] += (c1->E_cntr[DIM_p2] - c1->Phi_R[theDIM*(NUM_M-1)+1])/4.;   //but for the other contact mode
               c3->E_edge[DIM_p1] += (c1->E_cntr[DIM_p1] - c3->Phi_R[theDIM*(NUM_M-1)+0])/4.;
               c1->E_edge[DIM_p1] += (c1->E_cntr[DIM_p1] - c1->Phi_R[theDIM*(NUM_M-1)+0])/4.;
            }
            else{
               c2->E_edge[DIM_p2] += (c0->E_cntr[DIM_p2] - cA->Phi_R[theDIM*(NUM_M-1)+1] + c1->E_cntr[DIM_p2] - c2->Phi_R[theDIM*(NUM_M-1)+1])/8.;
               c1->E_edge[DIM_p2] += (c0->E_cntr[DIM_p2] - c0->Phi_R[theDIM*(NUM_M-1)+1] + c1->E_cntr[DIM_p2] - c1->Phi_R[theDIM*(NUM_M-1)+1])/8.;
               c3->E_edge[DIM_p1] += (c0->E_cntr[DIM_p1] - cB->Phi_R[theDIM*(NUM_M-1)+0] + c1->E_cntr[DIM_p1] - c3->Phi_R[theDIM*(NUM_M-1)+0])/8.;
               c1->E_edge[DIM_p1] += (c0->E_cntr[DIM_p1] - c0->Phi_R[theDIM*(NUM_M-1)+0] + c1->E_cntr[DIM_p1] - c1->Phi_R[theDIM*(NUM_M-1)+0])/8.;
            }//*/
            
   
         }
      }   
   }

}



void get_edge_Efields( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;

   
   int i1 = Nx+2*Ng;
   int j1 = Ny+2*Ng;
   int k1 = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i1 = 1;
   if( theDomain->theParList.Num_y == 1 ) j1 = 1;
   if( theDomain->theParList.Num_z == 1 ) k1 = 1;

   
   int ijk, dim;
   for( ijk=0 ; ijk<i1*j1*k1 ; ++ijk ){
      struct cell * c = theCells+ijk;
      //set edge E-fields to zero first
      for( dim=0 ; dim<NUM_M+1 ; ++dim ) c->E_edge[dim] = 0.;
   }

   //evaluate edge E-fields now
   
   if( theDomain->theParList.Num_x != 1 ) edge_Efields_1D( theDomain , 0 );
   if( theDomain->theParList.Num_y != 1 ) edge_Efields_1D( theDomain , 1 );
   if( theDomain->theParList.Num_z != 1 ) edge_Efields_1D( theDomain , 2 );  
   /*
   edge_Efields_1D( theDomain , 0 );
   edge_Efields_1D( theDomain , 1 );
   edge_Efields_1D( theDomain , 2 ); 
   */

}





void update_B_fluxes_1D( struct domain * theDomain , int theDIM , double dt , int first_step , int last_step ){
   //get E and dl from faces
   //add E*dl*dt to all Phis
   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   int i_index = Nx+2*Ng-1;
   int j_index = Ny+2*Ng-1;
   int k_index = Nz+2*Ng-1;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;

   int n[3] = {0}; n[theDIM] = 1;
   int DIM_p1 = (theDIM+1)%3; int n_p1[3] = {0}; n_p1[DIM_p1] = 1;
   int DIM_p2 = (theDIM+2)%3; int n_p2[3] = {0}; n_p2[DIM_p2] = 1;
   double dl_p1 = dx*(double)n_p1[0] + dy*(double)n_p1[1] + dz*(double)n_p1[2];
   double dl_p2 = dx*(double)n_p2[0] + dy*(double)n_p2[1] + dz*(double)n_p2[2];


   int i,j,k,ijk_1,ijk_2,ijk_3;     
   //convention: when traversing i-th direcion, 0:ijk, 1:(i+1)jk, 2:(i+1)(j+1)k, 3:(i+1)j(k+1)
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            ijk_1 = i;
            ijk_2 = i+n_p1[0];
            ijk_3 = i+n_p2[0];
            //if(i==10 && j==5 && k==5) printf("i1 = %i, i2 = %i, i3 = %i\n", i, i+n_p1[0], i+n_p2[0] );
            if( theDomain->theParList.Num_y != 1 ){
               ijk_1 += (Nx+2*Ng)*j;
               ijk_2 += (Nx+2*Ng)*(j+n_p1[1]);
               ijk_3 += (Nx+2*Ng)*(j+n_p2[1]);
               //if(i==10 && j==5 && k==5) printf("j1 = %i, j2 = %i, j3 = %i\n", j, j+n_p1[1], j+n_p2[1] );
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk_1 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               ijk_2 += (Nx+2*Ng)*(Ny+2*Ng)*(k+n_p1[2]);
               ijk_3 += (Nx+2*Ng)*(Ny+2*Ng)*(k+n_p2[2]);
               //if(i==10 && j==5 && k==5) printf("k1 = %i, k2 = %i, k3 = %i\n", k, k+n_p1[2], k+n_p2[2] );
            }
            struct cell * c1 = theCells+ijk_1;
            struct cell * c2 = theCells+ijk_2;
            struct cell * c3 = theCells+ijk_3;
            //update c1:flux[DIM]
            c1->Phi_B[theDIM] -= (c1->E_edge[DIM_p1]-c3->E_edge[DIM_p1])*dl_p1*dt + (c2->E_edge[DIM_p2]-c1->E_edge[DIM_p2])*dl_p2*dt;
            //if(theDIM==0 && i==27 && j==5 && k==5) printf("x = %e, j,k = %i %i, E_z = %e, E_z = %e, Delta(Phi_B) = %e\n", c1->xi[0], j,k, c1->E_edge[DIM_p2], c2->E_edge[DIM_p2], c1->Phi_B[theDIM] );
         }
      }   
   }
}


void update__B_fluxes( struct domain * theDomain , double dt , int first_step , int last_step ){
   
   if( theDomain->theParList.Num_x != 1 ) update_B_fluxes_1D( theDomain , 0 , dt , first_step , last_step );
   if( theDomain->theParList.Num_y != 1 ) update_B_fluxes_1D( theDomain , 1 , dt , first_step , last_step );
   if( theDomain->theParList.Num_z != 1 ) update_B_fluxes_1D( theDomain , 2 , dt , first_step , last_step );
   /*
   update_B_fluxes_1D( theDomain , 0 , dt , first_step , last_step );
   update_B_fluxes_1D( theDomain , 1 , dt , first_step , last_step );
   update_B_fluxes_1D( theDomain , 2 , dt , first_step , last_step );
   */

}





void B_faces_to_cells_1D( struct domain * theDomain , int theDIM , double dt , int first_step , int last_step ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;

   int i_n, i_n_max;       //index and max_index in the direction of traversal
   int n[3] = {0}; 
   n[theDIM] = 1;
   i_n_max = n[0]*(i_index-1) + n[1]*(j_index-1) + n[2]*(k_index-1);

   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];

   int i,j,k,ijk_0,ijk_1;     
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            i_n = n[0]*i+n[1]*j+n[2]*k;

            ijk_0 = i;
            if( i_n==i_n_max ) ijk_1 = i;
            else ijk_1 = i+n[0];
            if( theDomain->theParList.Num_y != 1 ){
               ijk_0 += (Nx+2*Ng)*j;
               if( i_n==i_n_max ) ijk_1 += (Nx+2*Ng)*j;
               else ijk_1 += (Nx+2*Ng)*(j+n[1]);
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk_0 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               if( i_n==i_n_max ) ijk_1 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               else ijk_1 += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]);
            }
            struct cell * c0 = theCells+ijk_0;
            struct cell * c1 = theCells+ijk_1;
            c0->prim[NUM_C+NUM_N+theDIM] = (c0->Phi_B[theDIM]+c1->Phi_B[theDIM])/2./dA;       //pretty sure this isn't needed
            //c0->cons[NUM_C+NUM_N+theDIM] = (c0->Phi_B[theDIM]+c1->Phi_B[theDIM])/2./dA * dx*dy*dz;
         }
      }   
   }

}


void B_faces_to_cells( struct domain * theDomain , double dt , int first_step , int last_step ){

   if( theDomain->theParList.Num_x != 1 ) B_faces_to_cells_1D( theDomain , 0 , dt , first_step , last_step );
   if( theDomain->theParList.Num_y != 1 ) B_faces_to_cells_1D( theDomain , 1 , dt , first_step , last_step );
   if( theDomain->theParList.Num_z != 1 ) B_faces_to_cells_1D( theDomain , 2 , dt , first_step , last_step );
   /*
   B_faces_to_cells_1D( theDomain , 0 , dt , first_step , last_step );
   B_faces_to_cells_1D( theDomain , 1 , dt , first_step , last_step );
   B_faces_to_cells_1D( theDomain , 2 , dt , first_step , last_step );
   */

}








void B_cells_to_faces_1D( struct domain * theDomain , int theDIM , double dt , int first_step , int last_step ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;

   int i_n;       //index in the direction of traversal
   int n[3] = {0}; 
   n[theDIM] = 1;

   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];

   int i,j,k,ijk_0,ijk_m1;     
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            i_n = n[0]*i+n[1]*j+n[2]*k;

            ijk_0 = i;
            if( i_n==0 ) ijk_m1 = i;
            else ijk_m1 = i-n[0];
            if( theDomain->theParList.Num_y != 1 ){
               ijk_0 += (Nx+2*Ng)*j;
               if( i_n==0 ) ijk_m1 += (Nx+2*Ng)*j;
               else ijk_m1 += (Nx+2*Ng)*(j-n[1]);
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk_0 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               if( i_n==0 ) ijk_m1 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               else ijk_m1 += (Nx+2*Ng)*(Ny+2*Ng)*(k-n[2]);
            }
            struct cell * c0  = theCells+ijk_0;
            struct cell * cm1 = theCells+ijk_m1;
            c0->Phi_B[theDIM] = (c0->prim[NUM_C+NUM_N+theDIM]+cm1->prim[NUM_C+NUM_N+theDIM])/2.*dA;
         }
      }   
   }

}

/*
void B_cells_to_faces_1D( struct domain * theDomain , int theDIM , double dt , int first_step , int last_step ){

   struct cell * theCells = theDomain->theCells;
   int Nx = theDomain->Nx;
   int Ny = theDomain->Ny;
   int Nz = theDomain->Nz;
   int Ng = theDomain->Ng;
   double dx = theDomain->dx;
   double dy = theDomain->dy;
   double dz = theDomain->dz;

   int i_index = Nx+2*Ng;
   int j_index = Ny+2*Ng;
   int k_index = Nz+2*Ng;
   if( theDomain->theParList.Num_x == 1 ) i_index = 1;
   if( theDomain->theParList.Num_y == 1 ) j_index = 1;
   if( theDomain->theParList.Num_z == 1 ) k_index = 1;

   int i_n, i_n_max;       //index and max_index in the direction of traversal
   int n[3] = {0}; 
   n[theDIM] = 1;
   i_n_max = n[0]*(i_index-1) + n[1]*(j_index-1) + n[2]*(k_index-1);

   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];

   int i,j,k,ijk_0,ijk_1;     
   for( k=0 ; k<k_index ; ++k ){
      for( j=0 ; j<j_index ; ++j ){
         for( i=0 ; i<i_index ; ++i ){
            i_n = n[0]*i+n[1]*j+n[2]*k;

            ijk_0 = i;
            if( i_n==i_n_max ) ijk_1 = i;
            else ijk_1 = i+n[0];
            if( theDomain->theParList.Num_y != 1 ){
               ijk_0 += (Nx+2*Ng)*j;
               if( i_n==i_n_max ) ijk_1 += (Nx+2*Ng)*j;
               else ijk_1 += (Nx+2*Ng)*(j+n[1]);
            }
            if( theDomain->theParList.Num_z != 1 ){
               ijk_0 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               if( i_n==i_n_max ) ijk_1 += (Nx+2*Ng)*(Ny+2*Ng)*k;
               else ijk_1 += (Nx+2*Ng)*(Ny+2*Ng)*(k+n[2]);
            }
            struct cell * c0 = theCells+ijk_0;
            struct cell * c1 = theCells+ijk_1;
            c1->Phi_B[theDIM] = (c0->prim[NUM_C+NUM_N+theDIM]+c1->prim[NUM_C+NUM_N+theDIM])/2.*dA;
         }
      }   
   }

}*/



void B_cells_to_faces( struct domain * theDomain , double dt , int first_step , int last_step ){

   if( theDomain->theParList.Num_x != 1 ) B_cells_to_faces_1D( theDomain , 0 , dt , first_step , last_step );
   if( theDomain->theParList.Num_y != 1 ) B_cells_to_faces_1D( theDomain , 1 , dt , first_step , last_step );
   if( theDomain->theParList.Num_z != 1 ) B_cells_to_faces_1D( theDomain , 2 , dt , first_step , last_step );

}



