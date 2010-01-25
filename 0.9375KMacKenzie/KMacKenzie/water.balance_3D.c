
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/


#include "constant.h"
#include "struct.geotop.09375.h"
#include "water.balance_3D.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "lu.h"

extern T_INIT *UV;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;
extern char *error_file_name;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void water_balance_3D(ALLDATA *adt){

	DOUBLEVECTOR *hc;
	DOUBLEMATRIX *h;
	DOUBLETENSOR *P;
	
	long n, r, c, l, s, ch;
	double Dt0, Dt, tb, te, te0, Loss;
	static double cum_losses;	
	
	FILE *f;

	if(adt->I->time == 0) cum_losses=0.0;

	initialize_doublevector(adt->C->Qsup_spread, 0.0);
	initialize_doublevector(adt->C->Qsub_spread, 0.0);
	
	hc=new_doublevector(adt->C->r->nh);	//channel 
	h=new_doublematrix(Nr, Nc);	//overland flow
	P=new_doubletensor(Nl, Nr, Nc);	//soil pressure
			
	te0=0.0;
	
	do{

		Dt0=adt->P->Dt;
		if(te0+Dt0 > adt->P->Dt) Dt0=adt->P->Dt-te0;
		te0 += Dt0;	
		
		printf("te0:%f Dt0:%f\n",te0,Dt0);
	
		n=1;
		te=0;
	
		do{
	
			tb=te;				
		
			do{
		
				Dt=Dt0/(double)n;
											
				Richards_3D(Dt, P, h, hc, &Loss, adt);

				f=fopen(error_file_name,"a");
				fprintf(f,"n:%ld Dt:%f Dt0:%f te0:%f tb:%f Loss:%e\n\n\n",n,Dt,Dt0,te0,tb,Loss);	
				fclose(f);
				
				printf("te0:%f Dt0:%f\n",te0,Dt0);
						
				if(fabs(Loss) > adt->P->MaxErrWb*(Dt/adt->P->Dt) && Dt>adt->P->DtminWb){
					n*=adt->P->nredDtWb;
					printf("reduced time step: dt:%f\n\n\n",Dt=Dt0/(double)n);
				}
									
			}while( fabs(Loss) > adt->P->MaxErrWb*(Dt/adt->P->Dt) && Dt>adt->P->DtminWb); 
		
			te=tb+Dt;  
			
			for(ch=1;ch<=adt->C->r->nh;ch++){
				for(s=1;s<=adt->C->Qsup_spread->nh;s++){
					adt->C->Qsup_spread->co[s] += (adt->C->Qsup->co[ch])*(adt->C->fraction_spread->co[ch][s])*0.001*UV->U->co[1]*UV->U->co[2];//m3/s
					adt->C->Qsub_spread->co[s] += (adt->C->Qsub->co[ch])*(adt->C->fraction_spread->co[ch][s])*0.001*UV->U->co[1]*UV->U->co[2];//m3/s
				}
			} 
										                                  
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(adt->L->LC->co[r][c]!=NoV){
						adt->W->h_sup->co[r][c] = h->co[r][c];
						//if( adt->C->ch->co[r][c] != 0) adt->C->h_sup->co[ adt->C->ch->co[r][c] ] = hc->co[ adt->C->ch->co[r][c] ];
						for(l=1;l<=Nl;l++){
							adt->S->P->co[l][r][c] = P->co[l][r][c];
						}
					}
				}
			}
		
			adt->W->out2->co[8]+=Loss;
			
			cum_losses+=Loss;
		
		}while(te<Dt0);
		
	}while(te0<adt->P->Dt);

	routing3(adt->C);

	free_doublevector(hc);
	free_doublematrix(h);
	free_doubletensor(P);
	
	f=fopen(error_file_name,"a");
	fprintf(f,"Total Loss:%e\n\n",cum_losses);
	fclose(f);
			
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Richards_3D(double Dt, DOUBLETENSOR *P, DOUBLEMATRIX *h, DOUBLEVECTOR *h_ch, double *loss, ALLDATA *adt){

	DOUBLEVECTOR *H00, *H0, *H1, *dH, *B;
	long i, l, r, c, cont, cont2, iter_tot=0, iter=0;
	short sy;
	double mass0=0.0, mass0h=0.0, mass1, mass1h, massloss0, massloss=1.E99, mass_in=0.0, mass_out, nw, res, psi, h_sup;
	short out;
	
	long n=(Nl+1)*adt->P->total_pixel;
		
	FILE *f;
		
	H00=new_doublevector(n);				
	H0=new_doublevector(n);
	dH=new_doublevector(n);
	H1=new_doublevector(n);
	B=new_doublevector(n);	
		
	for(i=1;i<=n;i++){

		l = adt->T->lrc_cont->co[i][1];
		r = adt->T->lrc_cont->co[i][2];
		c = adt->T->lrc_cont->co[i][3];
		sy = adt->S->type->co[r][c];
		
		if(l==0){ //overland flow
 
			H1->co[i] = adt->W->h_sup->co[r][c] + 1.E3*adt->T->Z0dp->co[r][c];
			mass0h += Fmax(0.0, adt->W->h_sup->co[r][c])/(double)adt->P->total_pixel;
			mass_in += adt->W->Pn->co[r][c]*Dt/(double)adt->P->total_pixel;;
		
		}else{	//subsurface flow
			
			H1->co[i] = adt->S->P->co[l][r][c] + adt->T->Z->co[l][r][c];		
			mass0 += adt->S->pa->co[sy][jdz][l]*theta_from_psi(adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->Esoil)/(double)adt->P->total_pixel;
			mass_in -= adt->S->ET->co[l][r][c]/(double)adt->P->total_pixel;;

		}
		
		H0->co[i] = H1->co[i];
		H00->co[i] = H1->co[i];		
	}
		
	cont=0;
	
	do{
		
		cont++;

		for(i=1;i<=n;i++){	
			
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];	
			
			if(adt->P->UpdateK==1) H00->co[i] = H1->co[i];
			
			H0->co[i] = H1->co[i];
		}
								
		for(i=1;i<=n;i++){
			B->co[i] = Find_b(i, H0, H00, Dt, adt);
		}
		
		iter=tridiag_preconditioned_conjugate_gradient_search( adt->P->TolCG, H1, H0, H00, Dt, B, Solve_Richards_3D_p, adt );
		iter_tot+=iter;
		
		cont2=0.0;
		nw=1.0;
		
		massloss0=massloss;
						
		do{
			
			mass1=0.0;	
			mass1h=0.0;
			mass_out=0.0;
					
			for(i=1;i<=n;i++){							
				l=adt->T->lrc_cont->co[i][1];
				r=adt->T->lrc_cont->co[i][2];
				c=adt->T->lrc_cont->co[i][3];
				sy=adt->S->type->co[r][c];

				dH->co[i] = H1->co[i] - H0->co[i];
				H1->co[i] = H0->co[i] + nw*dH->co[i];
				
				if(H1->co[i] != H1->co[i]) printf("no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);

				if(l==0){	//overland flow
				
					h_sup = Fmax(0.0, H1->co[i]-1.E3*adt->T->Z0dp->co[r][c]);
					mass1h += h_sup/(double)adt->P->total_pixel;
					
					/*if(adt->T->pixel_type->co[r][c] == 11 && h_sup > 0){
						adt->C->Q_sup_s->co[1] += sqrt( g * 1.E-3*h_sup ) * 1.E-3*h_sup * UV->U->co[1];
						mass_out += (sqrt( g * 1.E-3*h_sup ) * h_sup * Dt / UV->U->co[2])/(double)adt->P->total_pixel;
					}*/
												
				}else{
						
					psi = H1->co[i] - adt->T->Z->co[l][r][c];	
					mass1 += adt->S->pa->co[sy][jdz][l]*theta_from_psi(psi, l, r, c, adt->S, adt->P->Esoil)/(double)adt->P->total_pixel;	
						
				}
			}
			
			massloss = (mass0+mass0h-(mass1+mass1h+mass_out-mass_in));
						
			cont2++;
			nw/=adt->P->nredCorrWb;
			
		}while( fabs(massloss)>fabs(massloss0) && cont2<adt->P->MaxiterCorrWb );		
				
		res=0.0;
		for(i=1;i<=n;i++){
			res += fabs(H1->co[i] - H0->co[i])/(double)n;
		}
		
		out=0;
		if(res <= adt->P->TolWb) out=1;
		if(fabs(massloss) > adt->P->MaxErrWb) out=0;	
		if(cont==1) out=0;
		if(cont>=adt->P->MaxiterWb) out=1;
		if(iter==0 || cont2==adt->P->MaxiterCorrWb) out=1;

		printf("iter:%ld/%ld cont:%ld cont2:%ld massloss:%e massloss0:%e \n",iter,iter_tot,cont,cont2,massloss,massloss0);
		printf("mass_in:%e mass_out:%e \n",mass_in,mass_out);
		printf("mass0:%e mass0h:%e\n",mass0,mass0h);
		printf("mass1:%e mass1h:%e dsub:%f sup:%f\n",mass1,mass1h,mass1-mass0,mass1h-mass0h);

		f=fopen(error_file_name,"a");
		fprintf(f,"iter:%ld/%ld cont:%ld cont2:%ld massloss:%e massloss0:%e mass_in:%e mass_out:%e mass0:%e mass0h:%e mass1:%e mass1h:%e ",iter,iter_tot,cont,cont2,massloss,massloss0,mass_in,mass_out,mass0,mass0h,mass1,mass1h);
		fprintf(f,"res:%e out:%ld\n",res,out);
		fclose(f);
		
		printf("res:%e out:%ld cont2:%ld/%ld iter:%ld \n\n\n",res,out,cont2,adt->P->MaxiterCorrWb,iter);

		for(i=1;i<n;i++){
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
				
			if(l==0){
				h->co[r][c] = H1->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
			}else{
				P->co[l][r][c] = H1->co[i] - adt->T->Z->co[l][r][c];
			}
		}
		
		supflow3(Dt, adt->P->Dt, h, adt->T, adt->L, adt->W, adt->C, adt->P);
		
																								
	}while(out==0);
	

	*loss = massloss;
	
	free_doublevector(H00);					
	free_doublevector(H0);
	free_doublevector(dH);
	free_doublevector(H1);
	free_doublevector(B);
	
}
	
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Solve_Richards_3D(long i, DOUBLEVECTOR *H1, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, ALLDATA *adt){

	long l, r, c, landtype, I, R, C, LANDTYPE;
	short sy;
	double a=0.0;
	double dz=0.0, k=0.0, kn, dzn, psi, ds, A=0.0, h, klim;
	
	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];
	landtype=(long)adt->L->LC->co[r][c];	
		
	//hydraulic capacity (diagonal term)
	if(l==0){
		h = H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
		if(h>0) a += (1./Dt)*H1->co[i];
	}else{
		dz = adt->S->pa->co[sy][jdz][l];
		a += H1->co[i]*(dtheta_dpsi_from_psi(H0->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->Esoil)*dz/Dt);
	}
	
	//vertical hydraulic conductivity
	if(l>0) k = k_from_psi(jKv, H00->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
	
	//Vertical fluxes (diagonal and tridiagonal terms)
	if(l<Nl){
		I = adt->T->i_cont[l+1][r][c];	
		
		if(l==0){	//overland flow
			
			if( H0->co[i]< H0->co[I] ){	//upward flux
				kn = k_from_psi(jKv,  H0->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
			}else{	//downward flow
				kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
			}
			dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];
			
		}else{	//subsurface flow
		
			psi = H00->co[I] - adt->T->Z->co[l+1][r][c];		
			kn = k_from_psi(jKv, psi, l+1, r, c, adt->S, adt->P->imp );
			dzn = adt->S->pa->co[sy][jdz][l+1];
			kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
				k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			dzn = 0.5*dz + 0.5*dzn;
			
		}
		
		a += (kn/dzn)*(H1->co[i]-H1->co[I]);
	}
			
	if(l>0){
		I=adt->T->i_cont[l-1][r][c];
		if(l==1){	

			if( H0->co[I] < H0->co[i] ){	//upward flux
				kn = k_from_psi(jKv,  H0->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{	//downward flow
				kn = k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
			}
			dzn = 0.5*adt->S->pa->co[sy][jdz][l];

		}else{
		
			psi = H00->co[I] - adt->T->Z->co[l-1][r][c];		
			kn = k_from_psi(jKv, psi, l-1, r, c, adt->S, adt->P->imp );
			dzn = adt->S->pa->co[sy][jdz][l-1];
			kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
				k_from_psi(jKv,  psisat_from(l-1, r, c, adt->S), l-1, r, c, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			dzn = 0.5*dz + 0.5*dzn;
			
		}
		
		a+=(kn/dzn)*(H1->co[i]-H1->co[I]);
	}
	
	//lateral hydraulic conductivity
	if(l>0) k = k_from_psi(jKh,   H00->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
	
	//lateral fluxes
	//4 neighbouring cells
	R = r+1;
	C = c;
	ds = 1.E3*UV->U->co[2];	//[mm]
	//1.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(adt->L->LC->co[R][C]!=NoV){
			if(l==0){
			
				I=adt->T->i_cont[l][R][C];
				LANDTYPE = (long)adt->L->LC->co[R][C];
				//A = find_k_sup(i, I, adt->T->lrc_cont, adt->T->slope_H, adt->T->Z0dp, H00, adt->L->ty->co[landtype][jcm], adt->L->ty->co[LANDTYPE][jcm], adt->P->gamma_m);
				if(H0->co[i]<1.E3*adt->T->Z0dp->co[r][c] || H0->co[I]<1.E3*adt->T->Z0dp->co[R][C]) A=0.0;
				a += A*(H1->co[i]-H1->co[I])/pow(ds,2.);	//A in [mm2/s], ds in [mm], H in [mm]
			
			}else{
			
				I=adt->T->i_cont[l][R][C];				
				psi = H00->co[I] - adt->T->Z->co[l][R][C];
				kn = k_from_psi(jKh,   psi, l, R, C, adt->S, adt->P->imp );	
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKh,   psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
					k_from_psi(jKh,   psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				a += kn*(H1->co[i]-H1->co[I])/ds;
				
			}
		}
	}

	R = r-1;
	C = c;
	ds = 1.E3*UV->U->co[2];	//[mm]
	//2.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(adt->L->LC->co[R][C]!=NoV){
			if(l==0){
			
				I=adt->T->i_cont[l][R][C];
				LANDTYPE = (long)adt->L->LC->co[R][C];
				//A = find_k_sup(i, I, adt->T->lrc_cont, adt->T->slope_H, adt->T->Z0dp, H00, adt->L->ty->co[landtype][jcm], adt->L->ty->co[LANDTYPE][jcm], adt->P->gamma_m);
				if(H0->co[i]<1.E3*adt->T->Z0dp->co[r][c] || H0->co[I]<1.E3*adt->T->Z0dp->co[R][C]) A=0.0;
				a += A*(H1->co[i]-H1->co[I])/pow(ds,2.);	//A in [mm2/s], ds in [mm], H in [mm]
			
			}else{
			
				I=adt->T->i_cont[l][R][C];				
				psi = H00->co[I] - adt->T->Z->co[l][R][C];
				kn = k_from_psi(jKh,   psi, l, R, C, adt->S, adt->P->imp );	
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKh,   psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
					k_from_psi(jKh,   psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				a += kn*(H1->co[i]-H1->co[I])/ds;
				
			}
		}
	}

	R = r;
	C = c+1;
	ds = 1.E3*UV->U->co[1];	//[mm]
	//3.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(adt->L->LC->co[R][C]!=NoV){
			if(l==0){
			
				I=adt->T->i_cont[l][R][C];
				LANDTYPE = (long)adt->L->LC->co[R][C];
				//A = find_k_sup(i, I, adt->T->lrc_cont, adt->T->slope_H, adt->T->Z0dp, H00, adt->L->ty->co[landtype][jcm], adt->L->ty->co[LANDTYPE][jcm], adt->P->gamma_m);
				if(H0->co[i]<1.E3*adt->T->Z0dp->co[r][c] || H0->co[I]<1.E3*adt->T->Z0dp->co[R][C]) A=0.0;
				a += A*(H1->co[i]-H1->co[I])/pow(ds,2.);	//A in [mm2/s], ds in [mm], H in [mm]
			
			}else{
			
				I=adt->T->i_cont[l][R][C];				
				psi = H00->co[I] - adt->T->Z->co[l][R][C];
				kn = k_from_psi(jKh,   psi, l, R, C, adt->S, adt->P->imp );	
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKh,   psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
					k_from_psi(jKh,   psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				a += kn*(H1->co[i]-H1->co[I])/ds;
				
			}
		}
	}


	R = r;
	C = c-1;
	ds = 1.E3*UV->U->co[1];	//[mm]
	//4.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(adt->L->LC->co[R][C]!=NoV){
			if(l==0){
			
				I=adt->T->i_cont[l][R][C];
				LANDTYPE = (long)adt->L->LC->co[R][C];
				//A = find_k_sup(i, I, adt->T->lrc_cont, adt->T->slope_H, adt->T->Z0dp, H00, adt->L->ty->co[landtype][jcm], adt->L->ty->co[LANDTYPE][jcm], adt->P->gamma_m);
				if(H0->co[i]<1.E3*adt->T->Z0dp->co[r][c] || H0->co[I]<1.E3*adt->T->Z0dp->co[R][C]) A=0.0;
				a += A*(H1->co[i]-H1->co[I])/pow(ds,2.);	//A in [mm2/s], ds in [mm], H in [mm]
			
			}else{
			
				I=adt->T->i_cont[l][R][C];				
				psi = H00->co[I] - adt->T->Z->co[l][R][C];
				kn = k_from_psi(jKh,   psi, l, R, C, adt->S, adt->P->imp );	
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKh,   psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
					k_from_psi(jKh,   psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				a += kn*(H1->co[i]-H1->co[I])/ds;
				
			}
		}
	}
		
	//outlet
	/*if( (adt->T->pixel_type->co[r][c]==1 || adt->T->pixel_type->co[r][c]==11) && l==0){ 
		h = H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c];		
		if( h > 0 ){	
			a += (1.5*pow(g*1.E-3*h, 0.5)/UV->U->co[1])*H1->co[i];
		}
	}*/

	return(a);
	
}

double Solve_Richards_3D_p(long i, DOUBLEVECTOR *H1, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, void *adt) { 
	return Solve_Richards_3D(i, H1, H0, H00, Dt, (ALLDATA *)adt); 
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Find_b(long i, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, ALLDATA *adt){

	long l, r, c;
	short sy;
	double a = 0.0;
	double h, V0, V1, dz;
	
	
	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];

	//hydraulic capacity (diagonal term)
	if(l==0){
	
		h = H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
		if(h>0) a += (1./Dt)*H0->co[i];
	
		V1 = Fmax(0.0, H0->co[i]-1.E3*adt->T->Z0dp->co[r][c]);
		V0 = Fmax(0.0, adt->W->h_sup->co[r][c]);
		a += (V0-V1)/Dt;
		
	}else{
		dz = adt->S->pa->co[sy][jdz][l];
		a += H0->co[i]*(dtheta_dpsi_from_psi(H0->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->Esoil)*dz/Dt);

		V1 = dz * theta_from_psi(H0->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->Esoil);
		V0 = dz * theta_from_psi(adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->Esoil);
		a += (V0-V1)/Dt;
	}
			
	//outlet
	/*if( (adt->T->pixel_type->co[r][c]==1 || adt->T->pixel_type->co[r][c]==11) && l==0){ 
		h = H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c];		
		if( h > 0 ){	
			a += ( (1.5*sqrt(g*1.E-3*h)/UV->U->co[1])*H0->co[i] - (sqrt(g*1.E-3*h)/UV->U->co[1])*h );
		}
	}*/
	
	//evaporation
	if(l>0){
		a -= adt->S->ET->co[l][r][c];
	}else{
		a += adt->W->Pn->co[r][c];
	}
		
	return(a);

}



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_waterbalance3(WATER *wat, SOIL *sl, PAR *par, DOUBLEMATRIX *Z){

	long i,l,r,c;

	for(i=1;i<=par->chkpt->nrh;i++){
		r=par->rc->co[i][1];
		c=par->rc->co[i][2];

		wat->out1->co[5][i]+=wat->Pn->co[r][c]*par->Dt;
		wat->out1->co[6][i]+=(wat->Pn->co[r][c]-sl->Jinf->co[r][c])*par->Dt;
		wat->out1->co[7][i]+=(wat->h_sup->co[r][c] - wat->out1->co[9][i] - (wat->Pn->co[r][c]-sl->Jinf->co[r][c])*par->Dt);    /*positive if entering into the pixel*/
		wat->out1->co[8][i]=0.0;
		wat->out1->co[9][i]=wat->h_sup->co[r][c];

		for(l=1;l<=Nl;l++){
			wat->out1->co[10][i]-=0.0*par->Dt;/*positive if entering into the pixel*/
		}
	}

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(Z->co[r][c]!=NoV){ /*if the pixel is not a novalue*/
				wat->out2->co[3]+=wat->Pn->co[r][c]*par->Dt;								/*[mm]*/
				wat->out2->co[4]+=sl->Jinf->co[r][c]*par->Dt;							/*[mm]*/
			}
		}
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_slope_H(long ***I, DOUBLEVECTOR *H, DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, DOUBLEMATRIX *slope){


	long r, c;

	DOUBLEMATRIX *a;/*contiene, per ogni punto del bacino, il coefficiente angolare della retta
					interpolante il punto a monte, quello a valle e il punto stesso*/
	DOUBLEMATRIX *b;/*contiene, per ogni punto del bacino,il coefficiente angolare della retta
					interpolante il punto a destra, quello a sinistra e il punto stesso*/


	a = new_doublematrix(Nr, Nc);
	initialize_doublematrix(a, NoV);

	b = new_doublematrix(Nr, Nc);
	initialize_doublematrix(b, NoV);

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(LC->co[r][c]!=NoV);

				if(LC->co[r-1][c]!=NoV && LC->co[r+1][c] != NoV){
					a->co[r][c] = atan( ( Fmax(1.E-3*H->co[I[0][r-1][c]],Z->co[r-1][c]) - Fmax(1.E-3*H->co[I[0][r+1][c]],Z->co[r+1][c]) )/(2*UV->U->co[2]) );
				}else if(LC->co[r-1][c]==NoV && LC->co[r+1][c] != NoV){
					a->co[r][c] = atan( ( Fmax(1.E-3*H->co[I[0][r][c]],Z->co[r][c]) - Fmax(1.E-3*H->co[I[0][r+1][c]],Z->co[r+1][c]) )/(UV->U->co[2]) );
				}else if(LC->co[r-1][c]!=NoV && LC->co[r+1][c] == NoV){
					a->co[r][c] = atan( ( Fmax(1.E-3*H->co[I[0][r-1][c]],Z->co[r-1][c]) - Fmax(1.E-3*H->co[I[0][r][c]],Z->co[r][c]) )/(UV->U->co[2]) );
				}

				if(LC->co[r][c-1]!=NoV && LC->co[r][c+1] != NoV){
					b->co[r][c] = atan( ( Fmax(1.E-3*H->co[I[0][r][c-1]],Z->co[r][c-1]) - Fmax(1.E-3*H->co[I[0][r][c+1]],Z->co[r][c+1]) )/(2*UV->U->co[1]) );
				}else if(LC->co[r][c-1]==NoV && LC->co[r][c+1] != NoV){
					b->co[r][c] = atan( ( Fmax(1.E-3*H->co[I[0][r][c]],Z->co[r][c]) - Fmax(1.E-3*H->co[I[0][r][c+1]],Z->co[r][c+1]) )/(UV->U->co[1]) );
				}else if(LC->co[r][c-1]!=NoV && LC->co[r][c+1] == NoV){
					b->co[r][c] = atan( ( Fmax(1.E-3*H->co[I[0][r][c-1]],Z->co[r][c-1]) - Fmax(1.E-3*H->co[I[0][r][c]],Z->co[r][c]) )/(UV->U->co[1]) );
			}
		}
	}

	/*La matrice "slopes" contiene la pendenza di ogni pixel:*/

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(a->co[r][c]!=NoV){
				slope->co[r][c] = tan(acos(cos(fabs(a->co[r][c]))*cos(fabs(b->co[r][c]))));
				slope->co[r][c] = fabs(slope->co[r][c]);
			}else{
				slope->co[r][c] = 0.0;
			}
		}
	}

	free_doublematrix(a);
	free_doublematrix(b);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void supflow3(double Dt, double Dtmax, DOUBLEMATRIX *h, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par)

{
	long r,c,R,C,ch;                                    
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	double dx,dy;                                         // the two dimensions of a pixel
	double Ks;											  // the Strickler's coefficent calculated with a standard deviation
	double b[10];                                         // area perpendicular to the superficial flow divided by h_sup
	double i;											  // hydraulic gradient
	double q,tb,te,dt;
	short lu;
	
if(par->point_sim==0){	//distributed simulations

	initialize_doublevector(cnet->Qsub, 0.0);
	initialize_doublevector(cnet->Qsup, 0.0);
	
	
	dx=UV->U->co[1];                                    
	dy=UV->U->co[2];                                    
 
	b[1]=0.0;  b[2]=dy;             b[3]=dx/2.0+dy/2.0;  b[4]=dx;            b[5]=dx/2.0+dy/2.0;
	b[6]=dy;   b[7]=dx/2.0+dy/2.0;  b[8]=dx;             b[9]=dx/2.0+dy/2.0;

	te=0.0;
	
	do{
	
		tb=te;
		dt=Dt;
						
		//find dt min
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(h->co[r][c]>0 && top->DD->co[r][c]>=1 && top->DD->co[r][c]<=8 && top->pixel_type->co[r][c]==0){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( h->co[r][c] - h->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];				
					//i=top->i_DD->co[r][c];
					if(i<0) i=0;  //correct false slopes
					q=b[top->DD->co[r][c]+1]*Ks*pow(h->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;	//mm/s
					if(q>0){
						if(0.5*h->co[r][c]/q<dt) dt=0.5*h->co[r][c]/q; 
					}
				}
			}
		}
				
		te=tb+dt;
		if(te>Dt){
			te=Dt;
			dt=te-tb;
		}

		
		//printf("tb:%f dt:%f Dt:%f Dtmax:%f te:%f\n",tb,dt,Dt,Dtmax,te);
		
		
		//h(=height of water over the land-surface) is updated adding the new runoff(=precipita-
		//tion not infiltrated of the actual time); then it is calculated q_sub and it is checked
		//that its value is not greater than the avaible water on the land-surface:
		//Remember the units of measure: q_sup=[mm/s],b[m],Ks[m^(1/3)/s],h[mm],dx[m],dy[m]
		initialize_doublematrix(wat->q_sup, 0.0);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(h->co[r][c]>0 && top->DD->co[r][c]>=1 && top->DD->co[r][c]<=8 && top->pixel_type->co[r][c]==0){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];							
					i=0.001*( h->co[r][c] - h->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];				
					//i=top->i_DD->co[r][c];
					if(i<0) i=0.0;
					wat->q_sup->co[r][c]=b[top->DD->co[r][c]+1]*Ks*pow(h->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;

					if(wat->q_sup->co[r][c]!=wat->q_sup->co[r][c]){
						printf("NO VALUE SUP:%ld %ld %ld %ld %f %f\n",r,c,r+r_DD[top->DD->co[r][c]],c+c_DD[top->DD->co[r][c]],h->co[r][c],h->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]]);
						printf("i:%f Dh:%f Ks:%f pow:%f iF:%f\n",i,h->co[r][c] - h->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]],Ks,pow(h->co[r][c]/1000.0,1.0+par->gamma_m),top->i_DD->co[r][c]);
					}
				}
			}
		}

		//After the computation of the surface flow,these flows are moved trough D8 scheme:
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->pixel_type->co[r][c]==0){
				
					R=r+r_DD[top->DD->co[r][c]];
					C=c+c_DD[top->DD->co[r][c]];
					
					h->co[r][c]-=wat->q_sup->co[r][c]*dt;
					h->co[r][c]=Fmax(h->co[r][c], 0.0);
										
					//the superficial flow is added to the land pixels (code 0):    
					if (top->pixel_type->co[R][C]==0){
						h->co[R][C]+=wat->q_sup->co[r][c]*dt;
						
					//the superficial flow is added to flow which flows into the channel pixels (code 10):
					}else if (top->pixel_type->co[R][C]>=10){
						for(ch=1;ch<=cnet->r->nh;ch++){
							if(R==cnet->r->co[ch] && C==cnet->c->co[ch]){
								cnet->Qsup->co[ch] += wat->q_sup->co[r][c]*dt/Dtmax;	
																
								if(cnet->Qsup->co[ch]!=cnet->Qsup->co[ch]){
									printf("qsup no value: r:%ld c:%ld ch:%ld R:%ld C:%ld qsup:%f hsup:%f\n",r,c,ch,R,C,wat->q_sup->co[r][c],h->co[r][c]);
								}
							}
						}
					}
					
				}else if(top->pixel_type->co[r][c]>=10){
					for(ch=1;ch<=cnet->r->nh;ch++){
						if(r==cnet->r->co[ch] && c==cnet->c->co[ch]){
							cnet->Qsub->co[ch] += h->co[r][c]/Dtmax; //[mm/s]
							h->co[r][c]=0.0;
						}
					}
				}
				
			}
		}  
						
	}while(te<Dt);
	
	
}else{	//point simulation  
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				/*q=0.0;
				lu=(short)land->LC->co[r][c];
				Ks=land->ty->co[lu][jcm];
				i=pow(pow(top->dz_dx->co[r][c],2.0)+pow(top->dz_dy->co[r][c],2.0),0.5);			
				if(h->co[r][c]>0) q=Ks*pow(h->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0;	//mm/s
				h->co[r][c]-=q*Dt;
				if(h->co[r][c]<0) h->co[r][c]=0.0;*/
				h->co[r][c]=0.0;
			}
		}
	}
}

							
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void routing3(CHANNEL *cnet){

	long s;
	 
	/*The just calculated Qsup_spread and Qsub_spread are added to the flow already
	present in the channel-network(Q_sup_s and Q_sub_s) and at once it is made a
	traslation of Q_sup_s e Q_sub_s towards the outlet:*/
	for(s=1;s<=cnet->Qsup_spread->nh-1;s++){
		cnet->Q_sub_s->co[s]=cnet->Q_sub_s->co[s+1]+cnet->Qsub_spread->co[s];
		cnet->Q_sup_s->co[s]=cnet->Q_sup_s->co[s+1]+cnet->Qsup_spread->co[s];
	}
	
	cnet->Q_sub_s->co[cnet->Qsub_spread->nh]=cnet->Qsub_spread->co[cnet->Qsub_spread->nh];
	cnet->Q_sup_s->co[cnet->Qsup_spread->nh]=cnet->Qsup_spread->co[cnet->Qsup_spread->nh];
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
