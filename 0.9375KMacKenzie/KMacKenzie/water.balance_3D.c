
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

#define min_slope 0.0001

//[1/s]: Keff sediments / thickness sediments
#define Kb_ch 1.E-2

//channel fraction: channel width/grid with
#define B_dy 0.3

//discharge coefficient of the weir model
#define Cd 0.61

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void water_balance_3D(ALLDATA *adt){

	DOUBLEMATRIX *h;
	DOUBLETENSOR *PSI;

	long n, r, c, l, ch, s;
	double Dt0, Dt, tb, te, te0, Loss;
	static double cum_losses;

	FILE *f;

	cum_losses=0.0;

	h=new_doublematrix(Nr, Nc);
	PSI=new_doubletensor(Nl, Nr, Nc);

	initialize_doublevector(adt->C->Qsup_spread, 0.0);
	initialize_doublevector(adt->C->Qsub_spread, 0.0);

	te0=0.0;

	do{

		Dt0=adt->P->Dt;
		if(te0+Dt0 > adt->P->Dt) Dt0=adt->P->Dt-te0;
		te0 += Dt0;

		f=fopen(error_file_name, "a");
		printf("te0:%f Dt0:%f\n",te0,Dt0);
		fclose(f);

		n=1;
		te=0;

		do{

			tb=te;

			do{

				Dt=Dt0/(double)n;

				Richards_3D(Dt, PSI, h, &Loss, adt);

				f=fopen(error_file_name, "a");
				fprintf(f,"n:%ld Dt:%f Dt0:%f te0:%f tb:%f Loss:%e\n\n\n",n,Dt,Dt0,te0,tb,Loss);
				fclose(f);

				if(fabs(Loss) > adt->P->MaxErrWb && Dt>adt->P->DtminWb) n*=adt->P->nredDtWb;

			}while( fabs(Loss) > adt->P->MaxErrWb && Dt>adt->P->DtminWb);

			te=tb+Dt;

			for(ch=1;ch<=adt->C->r->nh;ch++){
				for(s=1;s<=adt->C->Qsup_spread->nh;s++){
					adt->C->Qsup_spread->co[s] += (adt->C->Qsup->co[ch])*(adt->C->fraction_spread->co[ch][s]);//m3/s
					adt->C->Qsub_spread->co[s] += (adt->C->Qsub->co[ch])*(adt->C->fraction_spread->co[ch][s]);//m3/s
				}
			}

			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(adt->L->LC->co[r][c]!=NoV){
						adt->W->h_sup->co[r][c] = h->co[r][c];
						for(l=1;l<=Nl;l++){
							adt->S->P->co[l][r][c] = PSI->co[l][r][c];
						}
					}
				}
			}

			adt->W->out2->co[8]+=Loss;

			cum_losses+=Loss;

		}while(te<Dt0);

	}while(te0<adt->P->Dt);

	//if(adt->P->channel_network == 1) routing3(adt->C);
	routing3(adt->C);

	free_doublematrix(h);
	free_doubletensor(PSI);

	f=fopen(error_file_name, "a");
	fprintf(f,"Total Loss:%e\n\n",cum_losses);
	fclose(f);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Richards_3D(double Dt, DOUBLETENSOR *P, DOUBLEMATRIX *h, double *loss, ALLDATA *adt){

	DOUBLEVECTOR *H00, *H0, *H1, *dH, *B;
	long i, l, r, c, cont, cont2, landtype, iter_tot=0, iter=0;
	short sy;
	double mass0=0.0, mass1, massloss0, massloss=1.E99, mass_to_channel, nw, res, psi, A, h_sup;
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
		landtype=(long)adt->L->LC->co[r][c];

		if(l==0){	//overland flow
			H1->co[i] = adt->W->h_sup->co[r][c] + 1.E3*adt->T->Z0dp->co[r][c];
			mass0 += Fmax(0.0, adt->W->h_sup->co[r][c]);
		}else{	//subsurface flow
			H1->co[i] = adt->S->P->co[l][r][c] + adt->T->Z->co[l][r][c];
			mass0 += adt->S->pa->co[sy][jdz][l]*theta_from_psi(adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->Esoil);
		}

		H0->co[i] = H1->co[i];
		H00->co[i] = H1->co[i];
	}

	find_slope_H(adt->T->i_cont, H0, adt->T->Z0dp, adt->L->LC, adt->T->slope_H);

	cont=0;

	do{

		cont++;

		for(i=1;i<=n;i++){
			if(adt->P->UpdateK==1) H00->co[i] = H1->co[i];
			H0->co[i] = H1->co[i];
		}

		for(i=1;i<=n;i++){
			B->co[i] = Find_b(i, H0, H00, Dt, adt);
		}

		iter=tridiag_preconditioned_conjugate_gradient_search(adt->P->TolCG, H1, H0, H00, Dt, B, Solve_Richards_3D_p, adt);
		iter_tot+=iter;

		for(i=1;i<=n;i++){
			dH->co[i] = H1->co[i]-H0->co[i];
		}

		cont2=0.0;
		nw=1.0;

		massloss0=massloss;

		do{

			mass1=0.0;
			mass_to_channel=0.0;

			for(i=1;i<=n;i++){
				H1->co[i] = H0->co[i] + nw*dH->co[i];
				if(H1->co[i] != H1->co[i]) printf("no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);

				l=adt->T->lrc_cont->co[i][1];
				r=adt->T->lrc_cont->co[i][2];
				c=adt->T->lrc_cont->co[i][3];
				sy=adt->S->type->co[r][c];

				if(l==0){
					h_sup = Fmax(0.0, H1->co[i]-1.E3*adt->T->Z0dp->co[r][c]);
					mass1 += h_sup;

					if(adt->T->pixel_type->co[r][c] == 10){
						if(h_sup>0){
							A = Cd*(2./3.)*sqrt(2.*g)*2.*UV->U->co[1];	//m^(3/2) s^(-1)
							adt->C->Qsup->co[ adt->C->ch->co[r][c] ] = A * pow( 1.E-3*h_sup , 1.5 ); //[m3/s]
							mass_to_channel += adt->C->Qsup->co[ adt->C->ch->co[r][c] ]*Dt*1.E3/(UV->U->co[1]*UV->U->co[2]);
						}else{
							adt->C->Qsup->co[ adt->C->ch->co[r][c] ] = 0.0;
						}

					}else if(adt->T->pixel_type->co[r][c] == 1){
						if(h_sup>0){
							A = Cd*(2./3.)*sqrt(2.*g)*2.*UV->U->co[1];	//m^(3/2) s^(-1)
							adt->C->Q_sup_s->co[1] = A * pow( 1.E-3*h_sup , 1.5 ); //[m3/s]
							mass_to_channel += adt->C->Q_sup_s->co[1]*Dt*1.E3/(UV->U->co[1]*UV->U->co[2]);
						}else{
							adt->C->Q_sup_s->co[1] = 0.0;
						}
					}

				}else{
					if(H1->co[i] - adt->T->Z->co[l][r][c] < PSImin) H1->co[i] = PSImin + adt->T->Z->co[l][r][c];
					psi = H1->co[i] - adt->T->Z->co[l][r][c];
					mass1 += adt->S->pa->co[sy][jdz][l]*theta_from_psi(psi, l, r, c, adt->S, adt->P->Esoil);

					if(l==1 && adt->T->pixel_type->co[r][c] == 10){
						if( H1->co[i] - 1.E3*adt->T->Z0dp->co[r][c] > 0 ){	//hyporreic flow
							adt->C->Qsub->co[ adt->C->ch->co[r][c] ] = 1.E-3*Kb_ch*B_dy*(H1->co[i] - 1.E3*adt->T->Z0dp->co[r][c])
								*UV->U->co[1]*UV->U->co[2];	//m3/s
							mass_to_channel += adt->C->Qsub->co[ adt->C->ch->co[r][c] ]*Dt*1.E3/(UV->U->co[1]*UV->U->co[2]);
						}
					}
				}
			}

			massloss = (mass0-(mass1+mass_to_channel))/(double)adt->P->total_pixel;

			cont2++;
			nw/=adt->P->nredCorrWb;

			f=fopen(error_file_name, "a");
			fprintf(f,"iter:%ld/%ld cont:%ld cont2:%ld massloss:%e massloss0:%e mass_to_channel:%e mass0:%e mass1:%e\n",iter,iter_tot,cont,cont2,massloss,massloss0,mass_to_channel,mass0,mass1);
			fclose(f);

			//printf("iter:%ld/%ld cont:%ld cont2:%ld massloss:%e massloss0:%e mass_to_channel:%e mass0:%e mass1:%e\n",iter,iter_tot,cont,cont2,massloss,massloss0,mass_to_channel,mass0,mass1);

		}while( fabs(massloss)>fabs(massloss0) && cont2<adt->P->MaxiterCorrWb );

		res=0.0;
		for(i=1;i<=n;i++){
			res+=fabs(H1->co[i] - H0->co[i])/(double)n;
		}

		out=0;
		if(res <= adt->P->TolWb) out=1;
		if(fabs(massloss) > adt->P->MaxErrWb) out=0;
		if(cont == 1) out=0;
		if(cont >= adt->P->MaxiterWb) out=1;
		if(iter==0 || cont >= adt->P->MaxiterWb) out=1;

		f=fopen(error_file_name, "a");
		fprintf(f,"res:%e out:%d\n",res,out);
		fclose(f);

		//printf("res:%e out:%d\n\n\n",res,out);

	}while(out==0);

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
/*!
\author
\date


\param
\param
\param
\param adt - (ALLDATA *)

\brief This routine solves ....

\return the i-th e


*/
	long l, r, c, landtype, I, R, C, LANDTYPE;
	short sy;
	double a = 0.0;
	double HyC, dz, k, kn, dzn, psi, ds, A, h, klim;

	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];
	landtype=(long)adt->L->LC->co[r][c];

	//hydraulic capacity (diagonal term)
	if(l==0){	//overland flow
		psi = H0->co[i]-1.E3*adt->T->Z0dp->co[r][c];
		if(psi>0){
			a += H1->co[i]*(1.0/Dt);
		}else{
			//a += H1->co[i]*(0.0/Dt);
		}
	}else{	//subsurface flow
		psi = H0->co[i] - adt->T->Z->co[l][r][c];
		HyC = dtheta_dpsi_from_psi(psi, l, r, c, adt->S, adt->P->Esoil);
		dz = adt->S->pa->co[sy][jdz][l];
		a += H1->co[i]*(HyC*dz/Dt);
	}

	//vertical hydraulic conductivity
	if(l>0){
		psi = H00->co[i] - adt->T->Z->co[l][r][c];
		k = k_from_psi( jKv, psi, l, r, c, adt->S, adt->P->imp );
	}

	//Vertical fluxes (diagonal and tridiagonal terms)
	if(l<Nl){
		I = adt->T->i_cont[l+1][r][c];
		if(l==0){	//overland flow
			if( (H0->co[i]-1.E3*adt->T->Z0dp->co[r][c]) < (H0->co[I] - adt->T->Z->co[l+1][r][c]) ){	//upward flux
				psi = H0->co[I] - adt->T->Z->co[l+1][r][c];
				kn = k_from_psi( jKv,  psi, l+1, r, c, adt->S, adt->P->imp );
			}else{	//downward flow
				kn = k_from_psi( jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
			}
			dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];
			if(adt->T->pixel_type->co[r][c]==10) kn*=(1.-B_dy);
		}else{	//subsurface flow
			psi = H00->co[I] - adt->T->Z->co[l+1][r][c];
			kn = k_from_psi( jKv, psi, l+1, r, c, adt->S, adt->P->imp );
			dzn = adt->S->pa->co[sy][jdz][l+1];
			kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);

			klim = Harmonic_Mean( dz, dzn, k_from_psi( jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
				k_from_psi( jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;

			dzn = 0.5*dz + 0.5*dzn;
		}
		a += (kn/dzn)*(H1->co[i]-H1->co[I]);
	}

	if(l>0){
		I=adt->T->i_cont[l-1][r][c];
		if(l==1){	//overland flow
			if( (H0->co[I]-1.E3*adt->T->Z0dp->co[r][c]) < (H0->co[i] - adt->T->Z->co[l][r][c]) ){	//upward flux
				psi = H0->co[i] - adt->T->Z->co[l][r][c];
				kn = k_from_psi( jKv,  psi, l, r, c, adt->S, adt->P->imp );
			}else{	//downward flow
				kn = k_from_psi( jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
			}
			dzn = 0.5*adt->S->pa->co[sy][jdz][l];
			if(adt->T->pixel_type->co[r][c]==10) kn*=(1.-B_dy);
		}else{
			psi = H00->co[I] - adt->T->Z->co[l-1][r][c];
			kn = k_from_psi( jKv, psi, l-1, r, c, adt->S, adt->P->imp );
			dzn = adt->S->pa->co[sy][jdz][l-1];
			kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);

			klim = Harmonic_Mean( dz, dzn, k_from_psi( jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
				k_from_psi( jKv,  psisat_from(l-1, r, c, adt->S), l-1, r, c, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;

			dzn = 0.5*dz + 0.5*dzn;
		}
		a+=(kn/dzn)*(H1->co[i]-H1->co[I]);
	}

	//lateral hydraulic conductivity
	if(l>0){
		psi = H00->co[i] - adt->T->Z->co[l][r][c];
		k = k_from_psi( jKh,  psi, l, r, c, adt->S, adt->P->imp );
	}

	//lateral fluxes
	//4 neighboring cells
	R = r+1;
	C = c;
	ds = 1.E3*UV->U->co[2];	//[mm]
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];
			if(l==0){	//overland flow
				LANDTYPE = (long)adt->L->LC->co[R][C];
				A = find_k_sup(i, I, adt->T->lrc_cont, adt->T->slope_H, adt->T->Z0dp, H00, adt->L->ty->co[landtype][jcm], adt->L->ty->co[LANDTYPE][jcm], adt->P->gamma_m);
				a += A*(H1->co[i]-H1->co[I])/pow(ds,2.);	//A in [mm2/s], ds in [mm], H in [mm]
			}else{
				psi = H00->co[I] - adt->T->Z->co[l][R][C];
				kn = k_from_psi( jKh,  psi, l, R, C, adt->S, adt->P->imp );
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);

				klim = Harmonic_Mean( dz, dzn, k_from_psi( jKh,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
					k_from_psi( jKh,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;

				a += kn*(H1->co[i]-H1->co[I])/ds;
			}
		}
	}

	R = r-1;
	C = c;
	ds = 1.E3*UV->U->co[2];	//[mm]
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];
			if(l==0){	//overland flow
				LANDTYPE = (long)adt->L->LC->co[R][C];
				A = find_k_sup(i, I, adt->T->lrc_cont, adt->T->slope_H, adt->T->Z0dp, H00, adt->L->ty->co[landtype][jcm], adt->L->ty->co[LANDTYPE][jcm], adt->P->gamma_m);
				a += A*(H1->co[i]-H1->co[I])/pow(ds,2.);	//A in [mm2/s], ds in [mm], H in [mm]
			}else{
				psi = H00->co[I] - adt->T->Z->co[l][R][C];
				kn = k_from_psi( jKh,  psi, l, R, C, adt->S, adt->P->imp );
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);

				klim = Harmonic_Mean( dz, dzn, k_from_psi( jKh,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
					k_from_psi( jKh,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;

				a += kn*(H1->co[i]-H1->co[I])/ds;
			}
		}
	}

	R = r;
	C = c+1;
	ds = 1.E3*UV->U->co[1];	//[mm]
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];
			if(l==0){	//overland flow
				LANDTYPE = (long)adt->L->LC->co[R][C];
				A = find_k_sup(i, I, adt->T->lrc_cont, adt->T->slope_H, adt->T->Z0dp, H00, adt->L->ty->co[landtype][jcm], adt->L->ty->co[LANDTYPE][jcm], adt->P->gamma_m);
				a += A*(H1->co[i]-H1->co[I])/pow(ds,2.);	//A in [mm2/s], ds in [mm], H in [mm]
			}else{
				psi = H00->co[I] - adt->T->Z->co[l][R][C];
				kn = k_from_psi( jKh,  psi, l, R, C, adt->S, adt->P->imp );
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);

				klim = Harmonic_Mean( dz, dzn, k_from_psi( jKh,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
					k_from_psi( jKh,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;

				a += kn*(H1->co[i]-H1->co[I])/ds;
			}
		}
	}

	R = r;
	C = c-1;
	ds = 1.E3*UV->U->co[1];	//[mm]
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];
			if(l==0){	//overland flow
				LANDTYPE = (long)adt->L->LC->co[R][C];
				A = find_k_sup(i, I, adt->T->lrc_cont, adt->T->slope_H, adt->T->Z0dp, H00, adt->L->ty->co[landtype][jcm], adt->L->ty->co[LANDTYPE][jcm], adt->P->gamma_m);
				a += A*(H1->co[i]-H1->co[I])/pow(ds,2.);	//A in [mm2/s], ds in [mm], H in [mm]
			}else{
				psi = H00->co[I] - adt->T->Z->co[l][R][C];
				kn = k_from_psi( jKh,  psi, l, R, C, adt->S, adt->P->imp );
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);

				klim = Harmonic_Mean( dz, dzn, k_from_psi( jKh,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
					k_from_psi( jKh,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;

				a += kn*(H1->co[i]-H1->co[I])/ds;
			}
		}
	}

	//channel interactions with surface flow
	if( (adt->T->pixel_type->co[r][c]==10 || adt->T->pixel_type->co[r][c]==1) && l==0){ //channel
		if( H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c] > 0 ){

			h = H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c];

			//Q = A*h^(1.5)*dx [m3/s]
			//q = Q/(dx*dy) [m/s]

			A = Cd*(2./3.)*sqrt(2.*g)*2.*UV->U->co[1];	//m^(3/2) s^(-1)
			A /= (UV->U->co[1]*UV->U->co[2]); //m^(-1/2) s^(-1)
			A *= pow(1.E3, -0.5);  //mm^(-1/2) s^(-1)

			a += 1.5 * A * pow( h, 0.5 ) * H1->co[i];	//mm/s

		}
	}

	//channel interaction with subsurface flow
	if(adt->T->pixel_type->co[r][c]==10 && l==1){ //channel
		if( H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c] > 0 ){
			//q = A*(H-zf)
			A = Kb_ch*B_dy; //s^(-1)
			a += A * H1->co[i];
		}
	}

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

	long l, r, c, landtype;
	short sy;
	double a = 0.0;
	double HyC, dz, theta1, theta0, psi, A, h;

	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];
	landtype=(long)adt->L->LC->co[r][c];

	//hydraulic capacity (diagonal term)
	if(l==0){	//overland flow
		psi = H0->co[i]-1.E3*adt->T->Z0dp->co[r][c];
		if(psi>0){
			a += (1.E3*adt->T->Z0dp->co[r][c] + Fmax(adt->W->h_sup->co[r][c], 0.0))/Dt;
		}else{
			a += (Fmax(adt->W->h_sup->co[r][c], 0.0))/Dt;
		}

	}else{	//subsurface flow
		psi = H0->co[i] - adt->T->Z->co[l][r][c];
		HyC = dtheta_dpsi_from_psi(psi, l, r, c, adt->S, adt->P->Esoil);
		theta1 = theta_from_psi(psi, l, r, c, adt->S, adt->P->Esoil);
		theta0 = theta_from_psi(adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->Esoil);
		dz = adt->S->pa->co[sy][jdz][l];
		a += (HyC*H0->co[i] + theta0 - theta1)*dz/Dt;
	}

	//channel interactions with surface flow
	if( (adt->T->pixel_type->co[r][c]==10 || adt->T->pixel_type->co[r][c]==1) && l==0){ //channel
		if( H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c] > 0 ){

			h = H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c];

			//Q = A*h^(1.5)*dx [m3/s]
			//q = Q/(dx*dy) [m/s]

			A = Cd*(2./3.)*sqrt(2.*g)*2.*UV->U->co[1];	//m^(3/2) s^(-1)
			A /= (UV->U->co[1]*UV->U->co[2]); //m^(-1/2) s^(-1)
			A *= pow(1.E3, -0.5);  //mm^(-1/2) s^(-1)

			a += ( 1.5 * A * pow( h , 0.5 ) * H0->co[i] - A * pow( h , 1.5 ) );	//mm/s

		}
	}

	//channel interaction with subsurface flow
	if(adt->T->pixel_type->co[r][c]==10 && l==1){ //channel
		if( H0->co[i] - 1.E3*adt->T->Z0dp->co[r][c] > 0 ){
			//q = A*(H-zf)
			A = Kb_ch*B_dy; //s^(-1)
			a += A * 1.E3*adt->T->Z0dp->co[r][c];
		}
	}

	return(a);

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

double find_k_sup(long i1, long i2, LONGMATRIX *lrc, DOUBLEMATRIX *slope, DOUBLEMATRIX *Z, DOUBLEVECTOR *H, double cm1, double cm2, double gamma){

	double a, h, cm, dHds;
	long r1, c1, r2, c2;

	r1 = lrc->co[i1][2];
	c1 = lrc->co[i1][3];
	if(lrc->co[i1][1]!=0) t_error("Error 1 in find_k_sup");

	r2 = lrc->co[i2][2];
	c2 = lrc->co[i2][3];
	if(lrc->co[i2][1]!=0) t_error("Error 2 in find_k_sup");

	h = Arithmetic_Mean( 1., 1., Fmax(0.0, (1.E-3*H->co[i1] - Z->co[r1][c1])), Fmax(0.0, (1.E-3*H->co[i2] - Z->co[r2][c2])) );	//[m]
	cm = Harmonic_Mean( 1., 1., cm1, cm2 );
	dHds = Fmax ( min_slope , Arithmetic_Mean( 1., 1., slope->co[r1][c1], slope->co[r2][c2] ) );

	a = 1.E6 * pow(h, 1.0+gamma) * cm * pow(dHds, -0.5);	//[mm2/s]

	return(a);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
