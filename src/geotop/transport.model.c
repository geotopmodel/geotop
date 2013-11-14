#include "constants.h"
#include "struct.geotop.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "util_math.h"
#include "transport.model.h"
//#include "water.balance.h"
#include "meteodata.h"
#include <time.h>

// warum short???
short transport_model(double Dt, DOUBLEVECTOR *H, ALLDATA *adt){

	long d, i, l, r, c, ch, Nl=1, sy;
	long n=(Nl+1)*adt->P->total_pixel;

	long il[6] = {0, 1,  0, 0,  0, 0};
	long ir[6] = {0, 0, -1, 1,  0, 0};
	long ic[6] = {0, 0,  0, 0, -1, 1};

	double area, ds, dz;

printf("hallo\n\n");

		//VERTICAL FLUXES
		if( i<=n){//land
			
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			sy=adt->S->type->co[r][c];
			
			ch=adt->C->ch->co[r][c];
			area=ds*ds/cos(adt->T->slope->co[r][c]*Pi/180.);
			if (ch>0) area-=adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			if (l>0) {dz = adt->S->pa->co[sy][jdz][l];}
			else {dz = 0;}



		}else{//channel
			
			l=adt->C->lch->co[i-n][1];
			ch=adt->C->lch->co[i-n][2];
			r=adt->C->r->co[ch];
			c=adt->C->c->co[ch];
			sy=adt->C->soil_type->co[ch];
			area=adt->C->length->co[ch] * adt->P->w_dx * ds;
			
			//vertical hydraulic conductivity
			if (l>0) {dz = adt->S->pa->co[sy][jdz][l];}
			else {dz = 0;}
		}


/*
			//flux from cell below
			if (l<Nl) {
				
				I = i+1;
				
				if(l==0){	//overland flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dzn;
					
					if( H->co[i] < H->co[I] ){	
						//upward flux
						kn = k_from_psi(jKn, H->co[I] - (adt->T->Z->co[l+1][r][c]-adt->P->depr_channel), SC->thi->co[l+1][ch], SC->T->co[l+1][ch], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}else{	
						//downward flow
						kn = k_from_psi(jKn, psisat_from(SC->thi->co[l+1][ch], l+1, adt->S->pa->co[sy]), SC->thi->co[l+1][ch], SC->T->co[l+1][ch], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}
					
				}else{	//subsurface flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dz + 0.5*dzn;
					
					if( H->co[i] < H->co[I] ){	
						//upward flux
						kn = k_from_psi(jKn, H->co[I] - (adt->T->Z->co[l+1][r][c]-adt->P->depr_channel), SC->thi->co[l+1][ch], SC->T->co[l+1][ch], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}else{	
						//downward flow
						kn = k_from_psi(jKn, H->co[i] - (adt->T->Z->co[l][r][c]-adt->P->depr_channel), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}			
					
					
					kmax = k_from_psi(jKn, psisat_from(SC->thi->co[l][ch], l, adt->S->pa->co[sy]), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					kmaxn = k_from_psi(jKn, psisat_from(SC->thi->co[l+1][ch], l+1, adt->S->pa->co[sy]), SC->thi->co[l+1][ch], SC->T->co[l+1][ch], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					kmaxn = Fmin(kmax, kmaxn);
					
					kn = Fmin(kn, kmaxn);
					
				}
				
				cnt++;
				Lx->co[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
			}
		}
		
		//LATERAL FLUXES
		if (i<=n){
			
			//lateral hydraulic conductivity
			if(l>0){
				k = k_from_psi(jKl, H->co[i] - adt->T->Z->co[l][r][c], SL->thi->co[l][j], SL->T->co[l][j], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
				kmax = k_from_psi(jKl, psisat_from(SL->thi->co[l][j], l, adt->S->pa->co[sy]), SL->thi->co[l][l], SL->T->co[l][j], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
				if(adt->T->BC_counter->co[r][c]>0){
					if(adt->T->pixel_type->co[r][c] == 1 || adt->T->pixel_type->co[r][c] == 2 || adt->T->pixel_type->co[r][c] == 11 || adt->T->pixel_type->co[r][c] == 12) Klat->co[adt->T->BC_counter->co[r][c]][l] = k;
				}
			}
			
			//4 neighbouring cells
			R = r-1;
			C = c;
			//1.
			if(R>=1 && R<=Nr && C>=1 && C<=Nc){
				if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type->co[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));//[m]
					
					if(l>0){
						
						//Subsurface Flow
						if (H->co[I] > H->co[i]) {
							kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi->co[l][J], l, adt->S->pa->co[syn]), SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
												
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						//Surface Flow
						if (H->co[I] > H->co[i]) {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[I] - adt->T->Z->co[l][R][C]) / cos(adt->T->slope->co[R][C]*Pi/180.);
						}else {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[i] - adt->T->Z->co[l][r][c]) / cos(adt->T->slope->co[r][c]*Pi/180.);
						}
						
						if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			R = r+1;
			C = c;
			dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));
			//2.
			if(R>=1 && R<=Nr && C>=1 && C<=Nc){
				if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type->co[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));//[m]
					
					if(l>0){
						if (H->co[I] > H->co[i]) {
							kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi->co[l][J], l, adt->S->pa->co[syn]), SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						if (H->co[I] > H->co[i]) {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[I] - adt->T->Z->co[l][R][C]) / cos(adt->T->slope->co[R][C]*Pi/180.);
						}else {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[i] - adt->T->Z->co[l][r][c]) / cos(adt->T->slope->co[r][c]*Pi/180.);
						}
						
						if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
					
				}
			}
			
			R = r;
			C = c-1;
			dn = ds/cos(0.5*atan(adt->T->dzdN->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));
			//3.
			if(R>=1 && R<=Nr && C>=1 && C<=Nc){
				if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type->co[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));//[m]
					
					if(l>0){
						if (H->co[I] > H->co[i]) {
							kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi->co[l][J], l, adt->S->pa->co[syn]), SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						if (H->co[I] > H->co[i]) {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[I] - adt->T->Z->co[l][R][C]) / cos(adt->T->slope->co[R][C]*Pi/180.);
						}else {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[i] - adt->T->Z->co[l][r][c]) / cos(adt->T->slope->co[r][c]*Pi/180.);
						}
						
						if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			R = r;
			C = c+1;
			dn = ds/cos(0.5*atan(adt->T->dzdN->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));
			//4.
			if(R>=1 && R<=Nr && C>=1 && C<=Nc){
				if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type->co[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));//[m]
					
					if(l>0){
						if (H->co[I] > H->co[i]) {
							kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi->co[l][J], l, adt->S->pa->co[syn]), SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						if (H->co[I] > H->co[i]) {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[I] - adt->T->Z->co[l][R][C]) / cos(adt->T->slope->co[R][C]*Pi/180.);
						}else {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[i] - adt->T->Z->co[l][r][c]) / cos(adt->T->slope->co[r][c]*Pi/180.);
						}
						
						if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			//exchange with channels
			if (l>0 && adt->C->ch->co[r][c] > 0) {
				
				ch = adt->C->ch->co[r][c];
				syn = adt->C->soil_type->co[ch];
				I = n + adt->C->ch3[l][ch];
				
				if (H->co[I] > H->co[i]) {
					kn = k_from_psi(jKl, H->co[I] - (adt->T->Z->co[l][r][c]-adt->P->depr_channel), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
				}else {
					kn = k;
				}
				
				kmaxn = k_from_psi(jKl, psisat_from(SC->thi->co[l][ch], l, adt->S->pa->co[syn]), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
				kmaxn = Fmin(kmax, kmaxn);
				kn = Fmin(kn, kmaxn);
				
				//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]	
				//Area[m2] = 2 * channel length * layer thickness = 2 * length[m] * dz[mm] * 1.E-3[m/mm]
				//dD[mm] = 0.25 * ds * (1+w_dx)
				dD = find_3Ddistance(ds * (1.0 + adt->P->w_dx) / 4.0, 1.E-3*adt->P->depr_channel) * 1.E3;//[mm]
				
				cnt++;
				Lx->co[cnt] = -(2.*adt->C->length->co[adt->C->ch->co[r][c]]*1.E-3*dz)*kn/dD;						
				
			}
		}	
	}										
	
*/	
	
	return 0;
	
}



/*
        for (iy = 0; iy < ny; iy++){
            for (ix = 0; ix < nx; ix++){
	            for (iz = 0; iz < nz
	                if ( a.Q_x(iy,ix,iz,it_inner) > 0 )
	                    a.dm_x(iy,ix,iz,it_inner) = a.m_inner(iy,ix,iz,it_inner) / a.V_inner(iy,ix,iz,it_inner) * a.Q_x(iy,ix,iz,it_inner) * a.dt_inner ;
	                else
	                    a.dm_x(iy,ix,iz,it_inner) = a.m_inner(iy,ix+1,iz,it_inner) / a.V_inner(iy,ix+1,iz,it_inner) * a.Q_x(iy,ix,iz,it_inner) * a.dt_inner ;
		    end
                end
            }
        }
        
        for iy = 1 : a.ny - 1
            for ix = 1 : a.nx
	            for iz = 1 : a.nz
	                if ( a.Q_y(iy,ix,iz,it_inner) > 0 )
	                    a.dm_y(iy,ix,iz,it_inner) = a.m_inner(iy,ix,iz,it_inner) / a.V_inner(iy,ix,iz,it_inner) * a.Q_y(iy,ix,iz,it_inner) * a.dt_inner ;
	                else
	                    a.dm_y(iy,ix,iz,it_inner) = a.m_inner(iy+1,ix,iz,it_inner) / a.V_inner(iy+1,ix,iz,it_inner) * a.Q_y(iy,ix,iz,it_inner) * a.dt_inner ;
	                end
                end
            end
        end
        
        for iy = 1 : a.ny
            for ix = 1 : a.nx
	            for iz = 1 : a.nz - 1
	                if ( a.Q_z(iy,ix,iz,it_inner) > 0 )
	                    a.dm_z(iy,ix,iz,it_inner) = a.m_inner(iy,ix,iz,it_inner) / a.V_inner(iy,ix,iz,it_inner) * a.Q_y(iy,ix,iz,it_inner) * a.dt_inner ;
	                else
	                    a.dm_z(iy,ix,iz,it_inner) = a.m_inner(iy,ix,iz+1,it_inner) / a.V_inner(iy,ix,iz+1,it_inner) * a.Q_y(iy,ix,iz,it_inner) * a.dt_inner ;
	                end
                end
            end
        end

        a.dm(2:a.ny-1,2:a.nx-1,it_inner) = a.dm_x(2:a.ny-1,1:a.nx-2,2:a.nz-1,it_inner) - a.dm_x(2:a.ny-1,2:a.nx-1,2:a.nz-1,it_inner) + ...
                                           a.dm_y(1:a.ny-2,2:a.nx-1,2:a.nz-1,it_inner) - a.dm_y(2:a.ny-1,2:a.nx-1,2:a.nz-1,it_inner) + ...
                                           a.dm_z(2:a.ny-1,2:a.nx-1,1:a.nz-2,it_inner) - a.dm_z(2:a.ny-1,2:a.nx-1,2:a.nz-1,it_inner) + ...
                                           a.input_concentration(2:a.ny-1,2:a.nx-1,it) .* a.precipitation(2:a.ny-1,2:a.nx-1,it) * a.dx * a.dy * 					   a.dt_inner ; // nodal mass change is calculated

*/

