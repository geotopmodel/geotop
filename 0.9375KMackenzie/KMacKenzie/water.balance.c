
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


/*Authors: Stefano Endrizzi and Matteo Dall'Amico - October 2008*/
#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "water.balance.h"
#include "pedo.funct.h"`
#include "t_datamanipulation.h"
#include "util_math.h"

void checkErrorSize(char *errfilepath);

extern T_INIT *UV;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;

/*--------------------------------------------*/
void water_balance(TOPO *top, SOIL *sl, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, double time){

	double Dt, DPsiMax, te, tb, dt;
	long i, n, r, c, l;
	FILE *f;
	DOUBLETENSOR *P1;

	P1=new_doubletensor(Nl,Nr,Nc);

	for(i=1;i<=par->nDt_water;i++){

		Dt=par->Dt/(double)par->nDt_water;

		te=0.0; //dt t begin

		initialize_doublevector(cnet->Q_sup_spread,0.0);
		initialize_doublevector(cnet->Q_sub_spread,0.0);

		do{
			tb=te;		//t begin
			te=Dt;		//t end
			dt=te-tb;	//dt integration
			n=1;		//parts in which (te-tb) is divided

			do{

				dt/=(double)n;

				subflow(dt, land, top, sl, par, wat, &DPsiMax, P1);

				if(DPsiMax>par->Dpsi) n++;	//decrease time step

			}while(dt/(double)n>par->dtmin && DPsiMax>par->Dpsi);

			f=fopen(files->co[ferr]+1,"a");
			fprintf(f,"\nWATER BALANCE: tb:%f dt:%f DPsiMax:%f maxadmitted:%f\n",tb,dt,DPsiMax,par->Dpsi);
			fclose(f);

			if(DPsiMax>par->Dpsi){
				printf("NOT ABLE TO REDUCE MAX PSI INCREMENT AFTER SUBSURFACE FLOW BELOW MAX ADMITTED VALUE, with dt:%f s\n",dt);
				printf("MAX PSI increment after lateral subsurface flow: %f, max admitted:%f\n",DPsiMax,par->Dpsi);
				f=fopen(files->co[ferr]+1,"a");
				fprintf(f,"NOT ABLE TO REDUCE MAX PSI INCREMENT AFTER SUBSURFACE FLOW BELOW MAX ADMITTED VALUE, with dt:%f s\n",dt);
				fprintf(f,"MAX PSI increment after lateral subsurface flow: %f, max admitted:%f\n",DPsiMax,par->Dpsi);
				fclose(f);
			}

			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(land->LC->co[r][c]!=NoV){
						for(l=1;l<=Nl;l++){
							sl->P->co[l][r][c]=P1->co[l][r][c];
						}
					}
				}
			}

			subflow_channel(dt, Dt, cnet, sl, par);

			vertical_water_balance(dt, land->LC, sl, wat, time, par);

			supflow(dt, Dt, top, land, wat, cnet, par);

			te=tb+dt;

		}while(te<Dt);

		routing(cnet);

		output_waterbalance(Dt, wat, sl, par, land->LC);

	}

	free_doubletensor(P1);
}


/*--------------------------------------------*/

void vertical_water_balance (double Dt, DOUBLEMATRIX *Z, SOIL *sl, WATER *wat, double time, PAR *par)

{

	long l;					/* the index of depth is l not d because d is used for sl-thickness vector*/
	long r,c;				/* counters of rows and columns*/
	long i;
	DOUBLEVECTOR *psi;
	short sy;
	double masserrorbasin=0.0, masserror;
	double Pnbasin=0.0, Infbasin=0.0;
	double DW, th0, th1, theq, h, psisat, oversat;
	FILE *f;                    /* file which contains the errors of simulation*/

	//long ri=0, ci=0;

	psi=new_doublevector(Nl);

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(Z->co[r][c]!=NoV){

				sy=sl->type->co[r][c];

				/*control of net precipitation:*/
				if (wat->Pn->co[r][c]<-0.01){
					if(par->n_error<par->max_error){
						par->n_error++;
						f=fopen(files->co[ferr]+1,"a");
						fprintf(f,"Pn[%ld,%ld] is negative: %f!\n",r,c,wat->Pn->co[r][c]);
						fclose(f);
						wat->Pn->co[r][c]=0.0;
					}
				}

				/*record psi_old*/
				for(l=1;l<=Nl;l++){
					psi->co[l]=sl->P->co[l][r][c];
				}

				/*solution of Richards' equation*/

				Richards(Dt, r, c, sl, psi, wat->Pn->co[r][c]+wat->h_sup->co[r][c]/Dt, time, par, &masserror);

				masserrorbasin+=masserror/(double)par->total_pixel;

				Pnbasin+=wat->Pn->co[r][c]*Dt/(double)par->total_pixel;

				if(sl->Jinf->co[r][c]!=sl->Jinf->co[r][c]){
					printf("ERROR ON INFILTRATION r:%ld c:%ld at time: %f sec\n",r,c,time);
					stop_execution();
				}

				Infbasin+=sl->Jinf->co[r][c]*Dt/(double)par->total_pixel;

				/*write q_out for specified pixels*/
				for(i=1;i<=par->chkpt->nrh;i++){
					wat->out1->co[28][i]=masserror*3600.0/Dt;
					if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){
						wat->out1->co[15][i]+=sl->J->co[1][r][c];
					}
				}
				
				/*another update of the height of water over the sl-surface:*/
				wat->h_sup->co[r][c]+=(wat->Pn->co[r][c]-sl->Jinf->co[r][c])*Dt; /*in mm*/

				/*update P*/
				for(l=1;l<=Nl;l++){
					if(sl->T->co[l][r][c]<=Tfreezing && (long)sl->pa->co[sy][jsf][l]==1 && par->en_balance==1){

						psisat=psi_saturation(sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
							1-1/sl->pa->co[sy][jns][l]);

						th1=teta_psi(fmin(psi->co[l],psisat),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
							sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);

						theq=teta_psi(Psif2(sl->T->co[l][r][c]),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],
							1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);

						if(th1>theq+1.E-8){

							th0=teta_psi(fmin(sl->P->co[l][r][c],psisat),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);

							oversat=teta_psi(psi->co[l],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil)
								-th1;

							DW=(th1-th0)*(sl->pa->co[sy][jdz][l]);
							h=internal_energy_soil(th0, sl->thice->co[l][r][c], sl->T->co[l][r][c], sl->pa->co[sy][jdz][l], sl->pa->co[sy][jct][l], sl->pa->co[sy][jsat][l]);
							from_internal_soil_energy(r, c, l, h+Lf*DW, &th1, &(sl->thice->co[l][r][c]), &(sl->T->co[l][r][c]), sl->pa->co[sy], par->psimin);

							sl->P->co[l][r][c]=psi_teta(th1+oversat,sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);

						}else{
							sl->P->co[l][r][c]=psi->co[l];
						}

					}else{
						sl->P->co[l][r][c]=psi->co[l];
					}
				}

				//error check
				if (wat->h_sup->co[r][c]<0.0){
					if (wat->h_sup->co[r][c]<-1E-8){
						if(par->n_error<par->max_error){
							par->n_error++;
							f=fopen(files->co[ferr]+1,"a");
							fprintf(f,"h_sup[%ld,%ld], at the end of vertical_water_balance subroutine, is negative: %f!\n",r,c,wat->h_sup->co[r][c]);
							fclose(f);
							printf("h_sup[%ld,%ld], at the end of vertical_water_balance subroutine, is negative: %f!\n",r,c,wat->h_sup->co[r][c]);
							printf("Inf:%f Pnet:%f\n",sl->Jinf->co[r][c],wat->Pn->co[r][c]);
						}
					}
					wat->h_sup->co[r][c]=0.0;
				}

				if (wat->h_sup->co[r][c]!=wat->h_sup->co[r][c]){
					wat->h_sup->co[r][c]=0.0;
					if(par->n_error<par->max_error){
						par->n_error++;
						f=fopen(files->co[ferr]+1,"a");
						fprintf(f,"h_sup[%ld,%ld], at the end of vertical_water_balance subroutine, is not a number!\n",r,c);
						fclose(f);
					}
				}
			}
		}
	}

	free_doublevector(psi);

	wat->out2->co[8]+=masserrorbasin*3600.0/Dt;

	//printf("%f %f\n",Pnbasin,Dt);

	f=fopen(files->co[ferr]+1,"a");
	fprintf(f,"\nERROR MASS BALANCE: %20.18fmm/h - Pnet: %20.18fmm/h - Infiltration: %20.18fmm/h\n",masserrorbasin*3600.0/Dt,Pnbasin*3600.0/Dt,Infbasin*3600.0/Dt);
	fclose(f);

}


/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void Richards(double Dt, long r, long c, SOIL *sl, DOUBLEVECTOR *psi, double Pnet, double t, PAR *par, double *masserrorcum)

{
	double dz,ddw,dup,dzdw,dzup;
	double Inf,theta0,theta1;
	double K1,Kdw,Kup,K1dw,K1up,C1,res,por,adi1,ads1,ad1,b1,b2;
	DOUBLEVECTOR *ad, *adi, *ads, *b, *e0, *e1, *psit;
	long cont,l,n,m;
	double tb,te,dtp,dt;	//tb=time of begin; te=time of end; dtp,dt time steps
	short infcond, sy=sl->type->co[r][c];
	double psisat=par->PsiInf+psi_saturation(sl->thice->co[1][r][c], sl->pa->co[sy][jsat][1], sl->pa->co[sy][jres][1], sl->pa->co[sy][ja][1], sl->pa->co[sy][jns][1], 1-1/sl->pa->co[sy][jns][1]);
	double mass0, mass1, masserror, masserror_tol=par->MaxerrVWB*Dt/3600.0, masserrordisplay;
	long ri=0, ci=0;

	*masserrorcum=0.0;

	ad=new_doublevector(Nl);
	adi=new_doublevector(Nl-1);
	ads=new_doublevector(Nl-1);
	b=new_doublevector(Nl);
	e0=new_doublevector(Nl);
	e1=new_doublevector(Nl);
	psit=new_doublevector(Nl);

	te=0.0; //dt t begin

	//initialization of the fluxes
	for(l=1;l<=Nl;l++){
		sl->J->co[l][r][c]=0.0;
		psit->co[l]=psi->co[l];
		if(r==ri && c==ci) printf("0. l:%ld psit:%f psi:%f %f\n",l,psit->co[l],psi->co[l],psisat);

	}
	sl->Jinf->co[r][c]=0.0;

	//if the first layer is oversaturated due to lateral flow, the excess water is exfiltrated and the pressure of the first layer is set at f
	if(psit->co[1]>psisat){
		//exfiltration occurs
		if(r==ri && c==ci) printf("IN psit:%e psisat:%e Diff:%e ",psit->co[1],psisat,psit->co[1]-psisat);
		theta1=teta_psi(psit->co[1],sl->thice->co[1][r][c],sl->pa->co[sy][jsat][1],sl->pa->co[sy][jres][1],sl->pa->co[sy][ja][1],
					sl->pa->co[sy][jns][1],1-1/sl->pa->co[sy][jns][1],fmin(sl->pa->co[sy][jpsimin][1], Psif(sl->T->co[1][r][c], par->psimin)),par->Esoil);
		por=sl->pa->co[sy][jsat][1]-sl->thice->co[1][r][c];
		if(r==ri && c==ci) printf("theta:%f por:%f ice:%f\n",theta1,por,sl->thice->co[1][r][c]);
		if(theta1>por+1.E-4) sl->Jinf->co[r][c]=-(theta1-por)*sl->pa->co[sy][jdz][1]/Dt;
		//printf("psi:%f psisat:%f theta:%f por:%f inf:%f\n",psit->co[1],psisat,theta1,por,sl->Jinf->co[r][c]);
		psit->co[1]=psisat;
		if(r==ri && c==ci) printf("theta1:%f por:%f\n",theta1,por);
		if(r==ri && c==ci) printf("Jinf0:%e\n ",sl->Jinf->co[r][c]);
	}

	do{
		tb=te;
		te=Dt;	//t end
		dtp=te-tb;	//dt integration
		dt=dtp;
		n=1;		//parts in which dtp is divided

		mass0=0.0;
		for(l=1;l<=Nl;l++){
			mass0+=sl->pa->co[sy][jdz][l]*teta_psi(psit->co[l],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],
				   sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);
		}


		do{

			dt=dt/(double)n;


/*--------------------------------------------------------------------------------------------------------------------*/
/*                                                       INITIALIZATION                                               */
/*--------------------------------------------------------------------------------------------------------------------*/

			for(l=1;l<=Nl;l++){
				e1->co[l]=psit->co[l];
			}

			if(r==ri && c==ci){
				printf("\n\n\n\n\nr:%ld c:%ld tb:%f te:%f dt:%f psimin:%f\n",r,c,tb,te,dt,par->psimin);
				for(l=1;l<=Nl;l++){
					printf("beg. l:%ld psit:%f e1:%f\n",l,psit->co[l],e1->co[l]);
				}
			}

			cont=0;
			Inf=Pnet;	//first attempt - then, if the soil saturates, infiltration is corrected

			do{ //loop of Picard iteration

				for(l=1;l<=Nl;l++){
					e0->co[l]=e1->co[l];
					if(r==ri && c==ci) printf("IN %ld e0:%f\n",l,e0->co[l]);
				}
				
				infcond=0;
				if(e1->co[1]-psisat>-0.001) infcond=1;
				//saturated soil either at the beginning of the time step or as an effect of the iteration (this method prevents numerical instabilities) 
				//the inflitration flux becomes the unknown
				
/*--------------------------------------------------------------------------------------------------------------------*/
/*                                                          BC  BLOCK                                                 */
/*--------------------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------------------------------------------*/
/*                                     MATRIX COEFFIECIENT BUILDING BLOCK                                             */
/*--------------------------------------------------------------------------------------------------------------------*/
				l=1;
				m=l;
				K1=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m], sl->pa->co[sy][ja][m],
					sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
				dz=sl->pa->co[sy][jdz][m];

				m=l+1;
				Kdw=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m], sl->pa->co[sy][ja][m],
					sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m],sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
				ddw=sl->pa->co[sy][jdz][m];

				dzdw=0.5*dz+0.5*ddw;

				K1dw=Mean((short)sl->pa->co[sy][jKav][l], dz, ddw, K1, Kdw);
				C1=dteta_dpsi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], 
					sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil, par->stmin);
				theta0=teta_psi(psit->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);
				theta1=teta_psi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);
										
				if(infcond==1){	
					if(r==ri && c==ci)printf("IRREGULAR Psisat:%f e0:%f D:%f\n",psisat,e0->co[1],psisat-e0->co[1]);					
					ad->co[l]= -1.0;
					ads->co[l]= -K1dw/dzdw;
					b->co[l]= -K1dw - K1dw/dzdw*psisat;
					ad1=C1*dz/dt + K1dw/dzdw;	//keeps memory of the coefficients in case the first soil layer is not saturated
					if(r==ri && c==ci)printf("C1:%e e0:%f ad:%e psimin:%f\n",C1,e0->co[l],ad1, fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)));
					ads1=-K1dw/dzdw;
					b1=(C1*e0->co[l]+theta0-theta1)*dz/dt + (Inf-K1dw);
				}else{							//first layer unsaturated, the suction head remains the unknown
					if(r==ri && c==ci)printf("REGULAR Psisat:%f e0:%f D:%f\n",psisat,e0->co[1],psisat-e0->co[1]);
					ad->co[l]= C1*dz/dt + K1dw/dzdw;
					ads->co[l]=-K1dw/dzdw;
					b->co[l]=(C1*e0->co[l]+theta0-theta1)*dz/dt + (Inf-K1dw);
					if(r==ri && c==ci)printf("C1:%e e0:%f ad:%e psimin:%f\n",C1,e0->co[l],ad1, fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)));
					if(r==ri && c==ci)printf("l:%ld coeff %e %e %e\n",l,ad->co[l],ads->co[l],b->co[l]);
				}

				//middle layers
				for(l=2;l<=Nl-1;l++){

					m=l;
					K1=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
						sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
					dz=sl->pa->co[sy][jdz][m];

					m=l+1;
					Kdw=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
						sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
					ddw=sl->pa->co[sy][jdz][m];

					m=l-1;
					Kup=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
						sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
					dup=sl->pa->co[sy][jdz][m];

					dzdw=0.5*dz+0.5*ddw;
					dzup=0.5*dz+0.5*dup;

					K1dw=Mean((short)sl->pa->co[sy][jKav][l], dz, ddw, K1, Kdw);
					K1up=Mean((short)sl->pa->co[sy][jKav][l], dz, dup, K1, Kup);
					C1=dteta_dpsi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
						sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil, par->stmin);
					theta0=teta_psi(psit->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
						sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);
					theta1=teta_psi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
						sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);

					ad->co[l]=C1*dz/dt + (K1dw/dzdw+K1up/dzup);
					ads->co[l]=-K1dw/dzdw;
					adi->co[l-1]=-K1up/dzup;
					b->co[l]=(C1*e0->co[l]+theta0-theta1)*dz/dt + (K1up-K1dw);
					if(r==ri && c==ci) printf("l:%ld coeff %e %e %e %e\n",l,ad->co[l],ads->co[l],adi->co[l],b->co[l]);		
					
				}
				if(infcond==1){	
				    b2=b->co[2];
					adi1=adi->co[1];
					e0->co[1]=psisat;
					b->co[2]-=adi->co[1]*psisat;
					adi->co[1]=0.0;
				}

				//last layer
				l=Nl;

				m=l;
				K1=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
					sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
				dz=sl->pa->co[sy][jdz][m];

				m=l-1;
				Kup=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
					sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
				dup=sl->pa->co[sy][jdz][m];

				dzup=0.5*dz+0.5*dup;

				K1up=Mean((short)sl->pa->co[sy][jKav][l], dz, dup, K1, Kup);
				C1=dteta_dpsi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil, par->stmin);
				theta0=teta_psi(psit->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);
				theta1=teta_psi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);

				ad->co[l]=C1*dz/dt + (K1up/dzup);
				adi->co[l-1]=-K1up/dzup;
				b->co[l]=(C1*e0->co[l]+theta0-theta1)*dz/dt + (K1up-K1*par->f_bound_Richards);

//---------------------------------------------------
//		solution of the tridiagonal system
//---------------------------------------------------

				tridiag(2,r,c,Nl,adi,ad,ads,b,e1);

				if(r==ri && c==ci){
					for(l=1;l<=Nl;l++){
						if(l<Nl) printf("l: %ld adi:%e ads:%e\n",l,adi->co[l],ads->co[l]);
						printf("l: %ld e1:%e ad:%e b:%e\n",l,e1->co[l],ad->co[l],b->co[l]);
					}
				}

				//saturated soil either at the beginning of the time step or as an effect of the iteration (this method prevents numerical instabilities)
				if(infcond==1){	
					
					theta0=teta_psi(psit->co[1], sl->thice->co[1][r][c], sl->pa->co[sy][jsat][1], sl->pa->co[sy][jres][1], sl->pa->co[sy][ja][1],
						sl->pa->co[sy][jns][1], 1-1/sl->pa->co[sy][jns][1], fmin(sl->pa->co[sy][jpsimin][1], Psif(sl->T->co[1][r][c], par->psimin)), par->Esoil);						
					theta1=sl->pa->co[sy][jsat][1]-sl->thice->co[1][r][c];	//porosity
					if(theta1>theta0){
						Inf=(theta1-theta0)*sl->pa->co[sy][jdz][1]/dt;	//if the soil is not saturated at the beginning of the time step, it is considered saturated now, updating the infilitration
					}else{
						Inf=0.0;
					}
					//printf("Pnet:%e Inf:%e theta0:%f theta1:%f psit:%f psisat:%f\n",Pnet,e1->co[1],theta0,theta1,e0->co[1],psisat);stop_execution();
					if(r==ri && c==ci) printf("CRUCIAL: Pnet:%e Inf:%e e1:%f theta0:%f theta1:%f psit:%f psisat:%f\n",Pnet,Inf,e1->co[1],theta0,theta1,e0->co[1],psisat);
					if(Pnet-Inf>=e1->co[1]){	//there is enough water that can infiltrate (e1 is the infiltration)
						if(r==ri && c==ci) printf("case 1 e1:%f ",e1->co[1]);
						Inf+=e1->co[1];
						if(r==ri && c==ci) printf("Pnet:%e Inf:%e ",Pnet,Inf);
						e1->co[1]=e0->co[1];
						if(r==ri && c==ci) printf("Pnet:%e Inf:%e ",Pnet,Inf);

						if(r==ri && c==ci){
							for(l=1;l<=Nl;l++){
								if(l<Nl) printf("l: %ld adi:%e ads:%e\n",l,adi->co[l],ads->co[l]);
								printf("l: %ld e1:%e ad:%e b:%e\n",l,e1->co[l],ad->co[l],b->co[l]);
							}
						}
						//printf("case1 Pnet:%e Inf:%e\n",Pnet,Inf);
					}else{	//there is not enough water that can infiltrate, the unknown becomes again the suction head
						if(r==ri && c==ci) printf("case2 Pnet:%e Inf:%e e1:%e \n",Pnet,Inf,e1->co[1]);

						l=1;
						ad->co[l]=ad1;
						ads->co[l]=ads1;
						b->co[l]=b1;

						l=2;
						adi->co[1]=ad1;
						b->co[2]=b2;

						Inf=Pnet;

						tridiag(2,r,c,Nl,adi,ad,ads,b,e1);

						if(r==ri && c==ci) printf("case2 Pnet:%e e1:%f \n",Pnet,e1->co[1]);


						if(r==ri && c==ci){
							for(l=1;l<=Nl;l++){
								if(l<Nl) printf("l: %ld adi:%e ads:%e\n",l,adi->co[l],ads->co[l]);
								printf("l: %ld e1:%e ad:%e b:%e\n",l,e1->co[l],ad->co[l],b->co[l]);
							}
						}

					}
				}

				//avoid instability when psi<psimin
				if(r==ri && c==ci)printf("CHECK INST: %f %f\n",e1->co[1],e1->co[2]);
				if(e1->co[1]>1.0E3 || e1->co[2]>1.0E4){					
					theta0=teta_psi(psit->co[1], sl->thice->co[1][r][c], sl->pa->co[sy][jsat][1], sl->pa->co[sy][jres][1], sl->pa->co[sy][ja][1],
						sl->pa->co[sy][jns][1], 1-1/sl->pa->co[sy][jns][1], fmin(sl->pa->co[sy][jpsimin][1], Psif(sl->T->co[1][r][c], par->psimin)), par->Esoil);
							
					if(r==ri && c==ci) printf("INST....theta0:%f theta1:%f e1:%f ",theta0,theta0+Inf*dt/sl->pa->co[sy][jdz][1],e1->co[1]);

					Inf=Pnet;

					e1->co[1]=psi_teta(theta0+Inf*dt/sl->pa->co[sy][jdz][1], sl->thice->co[1][r][c], sl->pa->co[sy][jsat][1], sl->pa->co[sy][jres][1], sl->pa->co[sy][ja][1],
							sl->pa->co[sy][jns][1], 1-1/sl->pa->co[sy][jns][1], fmin(sl->pa->co[sy][jpsimin][1], Psif(sl->T->co[1][r][c], par->psimin)), par->Esoil);

					if(r==ri && c==ci) printf("e1:%f \n",e1->co[1]);

					for(l=2;l<=Nl;l++){
						e1->co[l]=psit->co[l];
					}
				}


				//passage between insaturated and saturated soil (should be avoided)
				infcond=0;
				if(psit->co[1]<psisat && e1->co[1]>=psisat+par->TolPsiInf) infcond=1;
				if(r==ri && c==ci) printf("psit:%f psisat:%f e1:%f Inf:%f Pnet:%f infcond:%ld\n",psit->co[1],psisat,e1->co[1],Inf,Pnet,infcond);
				if(dt/(double)(n+1)<=par->dtminVWB) infcond=0.0;

				cont++;

				res=0.0;
				for(l=1;l<=Nl;l++){
					//if(psi<psimin) psi=psimin
					if(r==ri && c==ci) printf(".%ld...e1:%f psisat:%f",l,e1->co[l],psisat);
					theta1=teta_psi(e1->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
						sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);
					if(r==ri && c==ci) printf(" the1:%f por:%f",theta1,sl->pa->co[sy][jsat][l]-sl->thice->co[l][r][c]);
					e1->co[l]=psi_teta(theta1, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
						sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);
					if(r==ri && c==ci) printf(" e1:%f ",e1->co[l]);

					res+=pow((e1->co[l]-e0->co[l]),2.0);
					if(r==ri && c==ci) printf(" res:%f\n",pow((e1->co[l]-e0->co[l]),2.0));
				}
				res=pow(res,0.5);

				mass1=0.0;
				for(l=1;l<=Nl;l++){
					mass1+=sl->pa->co[sy][jdz][l]*teta_psi(e1->co[l],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],
						sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);
				}

				masserror=fabs(mass1-mass0-Inf*dt);

				if(r==ri && c==ci){
					printf("\nerr:%e res:%e tol:%e dt:%f tb:%f Inf:%e infcond:%ld psisat:%f \n",masserror,res,par->TolVWb,dt,tb,Inf,infcond,psisat);
					for(l=1;l<=Nl;l++){
						printf("end. l:%ld psit:%f e0:%f e1:%f\n",l,psit->co[l],e0->co[l],e1->co[l]);
					}
				}

				/*if(masserror>1){
					printf("masserror:%f\n",masserror);
					stop_execution();
				}*/

			}while(cont<par->MaxiterVWB && (res>par->TolVWb || masserror>masserror_tol*dt/Dt) && infcond==0);

			if(res>par->TolVWb || masserror>masserror_tol*dt/Dt || infcond==1) n++;	//decrease time step

		}while(dt/(double)n>par->dtminVWB && (res>par->TolVWb || masserror>masserror_tol*dt/Dt || infcond==1));

		te=tb+dt;

		*masserrorcum+=masserror;

		masserrordisplay=masserror_tol*dt/Dt;
		if(masserrordisplay<1.0E-6) masserrordisplay=1.0E-6;

		/*if(res>par->TolVWb){
			printf("res %f tol %f mass error %20.18f masserrormax:%20.18f  dt:%20.18f Dt:%f dtp:%20.18f  n:%ld infcond:%ld\n",res,par->TolVWb,masserror,masserror_tol*dt/Dt,dt,Dt,dtp,n,infcond);
			printf("Celia method for Richards equation does not converge (res=%f) even reducing the time step at %f seconds, r:%ld c:%ld! psit:%f e1:%f psisat:%f\n",res,dt,r,c,psit->co[1],e1->co[1],psisat);
			if(r==ri && c==ci) stop_execution();
		}*/

		//Check errors
		/*if(res>par->TolVWb || masserror>masserrordisplay){
			if(res>par->TolVWb){
				printf("res %f tol %f mass error %20.18f masserrormax:%20.18f  dt:%20.18f Dt:%f dtp:%20.18f  n:%ld infcond:%ld\n",res,par->TolVWb,masserror,masserror_tol*dt/Dt,dt,Dt,dtp,n,infcond);
				printf("Celia method for Richards equation does not converge (res=%f) even reducing the time step at %f seconds, r:%ld c:%ld! psit:%f e1:%f psisat:%f\n",res,dt,r,c,psit->co[1],e1->co[1],psisat);
			}
			//stop_execution();
			//if(par->n_error<par->max_error){
				//par->n_error++;
				f=fopen(files->co[ferr]+1,"a");
				fprintf(f,"Celia method for Richards equation does not converge (res=%f) even reducing the time step at %10.18f seconds, Infcond:%ld  %f %f r:%ld c:%ld!\n",res,dt,psit->co[1],e1->co[1],infcond,r,c);
				fprintf(f,"mass error:%f mass0:%f mass1:%f massinf:%f masssink:%f\n",mass1-mass0-Inf*dt,mass0,mass1,Inf*dt,0);
				fprintf(f,"Inf:%f Pnet:%f tb:%f te:%f\n",Inf,Pnet,tb,te);
				for(l=1;l<=Nl;l++){
					fprintf(f,"l:%ld psit:%f e0:%f e1:%f\n",l,psit->co[l],e0->co[l],e1->co[l]);
				}
				fclose(f);
			//}
			if(r==ri && c==ci) stop_execution();
		}*/

		/*Infiltration value update*/
		if(e1->co[1]>psisat){
			//exfiltration occurs
			theta1=teta_psi(e1->co[1],sl->thice->co[1][r][c],sl->pa->co[sy][jsat][1],sl->pa->co[sy][jres][1],sl->pa->co[sy][ja][1],
					sl->pa->co[sy][jns][1],1-1/sl->pa->co[sy][jns][1],fmin(sl->pa->co[sy][jpsimin][1], Psif(sl->T->co[1][r][c], par->psimin)),par->Esoil);

			if(r==ri && c==ci) printf("Inf:%e Infex:%e Infnet:%e\n",Inf,-(theta1-(sl->pa->co[sy][jsat][1]-sl->thice->co[1][r][c]))*sl->pa->co[sy][jdz][1]/dt,Inf-(theta1-(sl->pa->co[sy][jsat][1]-sl->thice->co[1][r][c]))*sl->pa->co[sy][jdz][1]/dt);

			Inf-=(theta1-(sl->pa->co[sy][jsat][1]-sl->thice->co[1][r][c]))*sl->pa->co[sy][jdz][1]/dt;
			e1->co[1]=psisat;
		}

		sl->Jinf->co[r][c]+=Inf*dt/Dt;
		if(r==ri && c==ci) printf("dt:%f Dt:%f Inf:%e %e Inf:%e\n",dt,Dt,Inf, Inf*dt/Dt,sl->Jinf->co[r][c]);

		/*fluxes update*/
		for(l=1;l<=Nl-1;l++){
			m=l;
			K1=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
				sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
			dz=sl->pa->co[sy][jdz][m];

			m=l+1;
			Kdw=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
				sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);
			ddw=sl->pa->co[sy][jdz][m];

			dzdw=0.5*dz+0.5*ddw;
			K1dw=Mean((short)sl->pa->co[sy][jKav][l], dz, ddw, K1, Kdw);

			sl->J->co[l][r][c]+=(K1dw*(e1->co[l]-e1->co[l+1])/dzdw + K1dw)*(dt/Dt);
		}

		l=Nl;
		m=l;
		K1=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
			sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[m][r][c], par->psimin)), sl->T->co[m][r][c]);

		sl->J->co[Nl][r][c]+=par->f_bound_Richards*K1*(dt/Dt);

		if(r==ri && c==ci) printf("END err:%e res:%e dt:%f tb:%f Inf:%e \n",masserror,res,dt,tb,Inf);

		//update psi
		for(l=1;l<=Nl;l++){
			psit->co[l]=e1->co[l];
		}

	}while(te<Dt);

	if(r==ri && c==ci){
		printf("\n\n\nEND TIME STEP Inf:%e Pnet:%e Error:%e\n",sl->Jinf->co[r][c],Pnet,*masserrorcum);
		for(l=1;l<=Nl;l++){
			printf("l:%ld psit:%f e1:%f\n",l,psi->co[l],e1->co[l]);
		}
	}

	/*Assign values*/
	for(l=1;l<=Nl;l++){
		psi->co[l]=e1->co[l];
	}

	/*Deallocate*/
	free_doublevector(ad);
	free_doublevector(adi);
	free_doublevector(ads);
	free_doublevector(b);
	free_doublevector(e0);
	free_doublevector(e1);
	free_doublevector(psit);

}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void supflow(double Dt, double Dtmax, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par)
{
	long r,c,R,C,ch,s;                                    /* counters*/
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; /* differential of number-pixel for rows and*/
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; /* columns, depending on Drainage Directions*/
	double dx,dy;                                         /* the two dimensions of a pixel*/
	double Ks;											  /* the Strickler's coefficent calculated with a standard deviation*/
	double q_sup_max;                                     /* maximum superficial flow (in the first part) or the incoming flow in a channel pixel (in the second part)*/
	double b[10];                                         /* area perpendicular to the superficial flow divided by h_sup*/
	double i;											  /*hydraulic gradient*/
	double q,tb,te,dt;
	short lu;


if(par->point_sim==0){	//distributed simulations

	initialize_doublematrix(wat->q_sup,0.0);
	initialize_doublevector(cnet->Q,0.0);

	dx=UV->U->co[1];                                     /*cell side [m]*/
	dy=UV->U->co[2];                                     /*cell side [m]*/

	b[1]=0.0;  b[2]=dy;             b[3]=dx/2.0+dy/2.0;  b[4]=dx;            b[5]=dx/2.0+dy/2.0;
	b[6]=dy;   b[7]=dx/2.0+dy/2.0;  b[8]=dx;             b[9]=dx/2.0+dy/2.0;

	te=0.0;

	do{

		tb=te;
		dt=Dt;
		//find dt min
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->pixel_type->co[r][c]==0 || top->pixel_type->co[r][c]==10){
					if(top->DD->co[r][c]!=0 && wat->h_sup->co[r][c]>0.1){
						if(top->DD->co[r][c]<=9){
							lu=(short)land->LC->co[r][c];
							Ks=land->ty->co[lu][jcm];
							i=0.001*( wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
							if(i<0) i=0;  //correct false slopes
							q=b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;	//mm/s
							if(q>0){ if(wat->h_sup->co[r][c]/q<dt) dt=wat->h_sup->co[r][c]/q; }
						}
					}
				}
			}
		}
		te=tb+dt;
		if(te>Dt){
			te=Dt;
			dt=te-tb;
		}

		/*wat->h_sup(=height of water over the land-surface) is updated adding the new runoff(=precipita-
		tion not infiltrated of the actual time); then it is calculated q_sub and it is checked
		that its value is not greater than the avaible water on the land-surface:*/
		/*Remember the units of measure: q_sup=[mm/s],b[m],Ks[m^(1/3)/s],wat->h_sup[mm],dx[m],dy[m]*/
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->pixel_type->co[r][c]==0 || top->pixel_type->co[r][c]==10){

					if (top->DD->co[r][c]!=0 && wat->h_sup->co[r][c]>0.1){

						if(top->DD->co[r][c]>9){
							wat->q_sup->co[r][c]=0.99*wat->h_sup->co[r][c]/dt; /*[mm/s]*/
						}else{
							lu=(short)land->LC->co[r][c];
							Ks=land->ty->co[lu][jcm];
							i=0.001*( wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
							if(i<0) i=0;  //correct false slopes
							wat->q_sup->co[r][c]=b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;
							if(wat->q_sup->co[r][c]!=wat->q_sup->co[r][c]){
								printf("NO VALUE SUP:%ld %ld %ld %ld %f %f\n",r,c,r+r_DD[top->DD->co[r][c]],c+c_DD[top->DD->co[r][c]],wat->h_sup->co[r][c],wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]]);
								printf("i:%f Dh:%f Ks:%f pow:%f iF:%f\n",i,wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]],Ks,pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m),top->i_DD->co[r][c]);
							}
						}

					}
				}
			}
		}

		/*After the computation of the surface flow,these flows are moved trough D8 scheme:*/
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->pixel_type->co[r][c]!=9){

					R=r+r_DD[top->DD->co[r][c]];
					C=c+c_DD[top->DD->co[r][c]];

					/*the superficial flow is subtracted to the land pixels (code 0):*/
					if (top->pixel_type->co[r][c]==0){
						wat->h_sup->co[r][c]-=wat->q_sup->co[r][c]*dt;
						//possible correction
						if(wat->h_sup->co[r][c]<0){
							wat->q_sup->co[r][c]+=wat->h_sup->co[r][c]/dt;
							wat->h_sup->co[r][c]=0.0;
						}
					}

					//possible correction for channel pixels (code 10)
					if (top->pixel_type->co[r][c]==10){
						if(wat->h_sup->co[r][c]-wat->q_sup->co[r][c]*dt<0){
							wat->q_sup->co[r][c]+=(wat->h_sup->co[r][c]-wat->q_sup->co[r][c]*dt)/dt;
						}
					}

					/*the superficial flow is added to the land pixels (code 0):*/
					if (top->pixel_type->co[R][C]==0) wat->h_sup->co[R][C]+=wat->q_sup->co[r][c]*dt;

					/*the superficial flow is added to flow which flows into the channel pixels (code 10):*/
					if (top->pixel_type->co[R][C]==10){
						for(ch=1;ch<=cnet->r->nh;ch++){
							if(R==cnet->r->co[ch] && C==cnet->c->co[ch]){
								cnet->Q->co[ch]+=wat->q_sup->co[r][c];
								if(cnet->Q->co[ch]!=cnet->Q->co[ch]){
									printf("qsub no value: r:%ld c:%ld ch:%ld R:%ld C:%ld qsup:%f hsup:%f\n",r,c,ch,R,C,wat->q_sup->co[r][c],wat->h_sup->co[r][c]);
								}
							}
						}
					}
				}
			}
		}

		for(ch=1;ch<=cnet->r->nh;ch++){

			r=cnet->r->co[ch];
			c=cnet->c->co[ch];

			lu=(short)land->LC->co[r][c];
			Ks=land->ty->co[lu][jcm];

			/*computation of value of flow which can go into channel:*/
			q_sup_max=fmin(2.0*b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(top->i_ch->co[r][c])*1000.0/dx/dy,cnet->Q->co[ch]);
			if (top->DD->co[r][c]==10) q_sup_max=cnet->Q->co[ch];	//closure section
			if(q_sup_max!=q_sup_max) printf("q_sup_max no value: qsup:%f r:%ld c:%ld ch:%ld\n",cnet->Q->co[ch],r,c,ch);

			/*wat of channel pixel which becames lateral flow and doesn't go into the channel;
			the variable cnet->q_sup is use in this case improperly:*/
			cnet->Q->co[ch]-=q_sup_max;

			/*the water on surface in channel pixel is updated:*/
			wat->h_sup->co[r][c]-=wat->q_sup->co[r][c]*dt;
			wat->h_sup->co[r][c]+=cnet->Q->co[ch]*dt;

			/*now is assingned to cnet->q_sup is right value:*/
			cnet->Q->co[ch]=q_sup_max;

			/*The new water arrived in the channel-network is spread immediately along a
			virtual coordinate s (=distance from the outlet) so that it can travel to the
			outlet, after simple traslation, compling with the convolution of the De
			Saint Venant equations. This is made either for superficial flow or for the
			subsuperficial flow:*/
			for(s=1;s<=cnet->Q_sup_spread->nh;s++){
				cnet->Q_sup_spread->co[s]+=(cnet->Q->co[ch]*dt/Dtmax)*(cnet->fraction_spread->co[ch][s])*0.001*dx*dy; /*in mc/s*/
			}

		}

	}while(te<Dt);

}else{	//point simulation

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				q=0.0;
				lu=(short)land->LC->co[r][c];
				Ks=land->ty->co[lu][jcm];
				i=pow(pow(top->dz_dx->co[r][c],2.0)+pow(top->dz_dy->co[r][c],2.0),0.5);			
				if(wat->h_sup->co[r][c]>0) q=Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0;	//mm/s
				wat->h_sup->co[r][c]-=q*Dt;
				if(wat->h_sup->co[r][c]<0) wat->h_sup->co[r][c]=0.0;
			}
		}
	}
}


}


/*--------------------------------------------*/
void subflow(double Dt, LAND *land, TOPO *top, SOIL *sl, PAR *par, WATER *wat, double *DeltaPsiMax, DOUBLETENSOR *P1)
{
	long l;  /*the index of depth is l not d because d is used for the vector of the sl-thickness*/
	long r,c,R,C;/*counters of rows and columns*/
	short sy, sy1, flux; //sl type
	double dx,dy; /*the two dimensions of a pixel*/
	double ds,dn;
	double q;    /*auxiliar variable to calculate the subflow for unit of land-surface*/
	double i; //lateral gradient of pressure [-]
	double psisat, psisat1, Kh, Kh1;
	//debug
	long ri=0, ci=0;

	*DeltaPsiMax=0.0;

	dx=UV->U->co[1];	//m
	dy=UV->U->co[2];	//m

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				for(l=1;l<=Nl;l++){
					P1->co[l][r][c]=sl->P->co[l][r][c];
				}
			}
		}
	}


	if(par->point_sim==0){	//2D simulations

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc-1;c++){

				for(l=1;l<=Nl;l++){
					if(sl->P->co[l][r][c]!=sl->P->co[l][r][c]){
						printf("nan %ld %ld\n",r,c);
						stop_execution();
					}
				}

				R=r;
				C=c+1;
				ds=dx;
				dn=dy;

				/*case with two land-pixel:*/
				if((top->pixel_type->co[r][c]==0 || top->pixel_type->co[r][c]==10) && (top->pixel_type->co[R][C]==0 || top->pixel_type->co[R][C]==10)){

					sy=sl->type->co[r][c];
					sy1=sl->type->co[R][C];

					for(l=1;l<=Nl;l++){

						flux=0;
						if( sl->pa->co[sy][jlatfl][l]==1 || sl->pa->co[sy1][jlatfl][l]==1 ) flux=1;
						if( sl->pa->co[sy][jlatfl][l]==2 || sl->pa->co[sy1][jlatfl][l]==2 ) flux=2;

						psisat=-20.0;
						psisat1=-20.0;

						if(r==ri && c==ci) printf("l:%ld Psi:%f Psineigh:%f R:%ld C:%ld h:%f\n",l,sl->P->co[l][r][c],sl->P->co[l][R][C],R,C,wat->h_sup->co[r][c]);

						if( flux==2 || (flux==1 && (sl->P->co[l][r][c]>=psisat || sl->P->co[l][R][C]>=psisat1)) ){
							Kh=K(sl->P->co[l][r][c], sl->pa->co[sy][jKh][l], par->imp, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l],
									sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], sl->pa->co[sy][jv][l],
									fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), sl->T->co[l][r][c]);
							Kh1=K(sl->P->co[l][R][C], sl->pa->co[sy1][jKh][l], par->imp, sl->thice->co[l][R][C], sl->pa->co[sy1][jsat][l],
									sl->pa->co[sy1][jres][l], sl->pa->co[sy1][ja][l], sl->pa->co[sy1][jns][l], 1-1/sl->pa->co[sy1][jns][l], sl->pa->co[sy1][jv][l],
									fmin(sl->pa->co[sy1][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), sl->T->co[l][R][C]);

							i=top->dz_dx->co[r][c]+0.001*(sl->P->co[l][r][c]-sl->P->co[l][R][C])/ds;
							q=Mean(sl->pa->co[sy][jKav][l], ds, ds, Kh, Kh1)*i*sl->pa->co[sy][jdz][l]*dn*1.0E-6;	//[m3/s]
							if(r==ri && c==ci) printf("%ld %ld l:%ld T:%f i:%f if:%f ig:%f P:%f P:%f k:%e q:%f\n",R,C,l,sl->T->co[l][r][c],i,top->dz_dx->co[r][c],0.001*(sl->P->co[l][r][c]-sl->P->co[l][R][C])/ds,sl->P->co[l][r][c],sl->P->co[l][R][C],Mean(sl->pa->co[sy][jKav][l], ds, ds, Kh, Kh1),q);

							set_psi(P1, sl, &q, Dt, l, r, c, R, C, par->psimin, par->Esoil);

							if(r==ri && c==ci) printf("P:%f oldP:%f P:%f oldP:%f q:%e\n",P1->co[l][r][c],sl->P->co[l][r][c],P1->co[l][R][C],sl->P->co[l][R][C],q);
						}
					}
				}
			}
		}

		for(c=1;c<=Nc;c++){
			for(r=1;r<=Nr-1;r++){

				R=r+1;
				C=c;
				ds=dy;
				dn=dx;

				/*case with two land-pixel:*/
				if((top->pixel_type->co[r][c]==0 || top->pixel_type->co[r][c]==10) && (top->pixel_type->co[R][C]==0 || top->pixel_type->co[R][C]==10)){

					sy=sl->type->co[r][c];
					sy1=sl->type->co[R][C];

					for(l=1;l<=Nl;l++){

						flux=0;
						if( sl->pa->co[sy][jlatfl][l]==1 || sl->pa->co[sy1][jlatfl][l]==1 ) flux=1;
						if( sl->pa->co[sy][jlatfl][l]==2 || sl->pa->co[sy1][jlatfl][l]==2 ) flux=2;

						psisat=-20.0;
						psisat1=-20.0;

						if(r==ri && c==ci) printf("l:%ld Psi:%f Psineigh:%f R:%ld C:%ld h:%f\n",l,sl->P->co[l][r][c],sl->P->co[l][R][C],R,C,wat->h_sup->co[r][c]);

						if( flux==2 || (flux==1 && (sl->P->co[l][r][c]>=psisat || sl->P->co[l][R][C]>=psisat1)) ){
							Kh=K(sl->P->co[l][r][c], sl->pa->co[sy][jKh][l], par->imp, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l],
									sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], sl->pa->co[sy][jv][l],
									fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), sl->T->co[l][r][c]);
							Kh1=K(sl->P->co[l][R][C], sl->pa->co[sy1][jKh][l], par->imp, sl->thice->co[l][R][C], sl->pa->co[sy1][jsat][l],
									sl->pa->co[sy1][jres][l], sl->pa->co[sy1][ja][l], sl->pa->co[sy1][jns][l], 1-1/sl->pa->co[sy1][jns][l], sl->pa->co[sy1][jv][l],
									fmin(sl->pa->co[sy1][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), sl->T->co[l][R][C]);

							i=top->dz_dy->co[r][c]+0.001*(sl->P->co[l][r][c]-sl->P->co[l][R][C])/ds;
							q=Mean(sl->pa->co[sy][jKav][l], ds, ds, Kh, Kh1)*i*sl->pa->co[sy][jdz][l]*dn*1.0E-6;	//[m3/s]
							if(r==ri && c==ci) printf("%ld %ld l:%ld T:%f i:%f if:%f ig:%f P:%f P:%f k:%e q:%f\n",R,C,l,sl->T->co[l][r][c],i,top->dz_dx->co[r][c],0.001*(sl->P->co[l][r][c]-sl->P->co[l][R][C])/ds,sl->P->co[l][r][c],sl->P->co[l][R][C],Mean(sl->pa->co[sy][jKav][l], ds, ds, Kh, Kh1),q);

							set_psi(P1, sl, &q, Dt, l, r, c, R, C, par->psimin, par->Esoil);

							if(r==ri && c==ci) printf("P:%f oldP:%f P:%f oldP:%f q:%e\n",P1->co[l][r][c],sl->P->co[l][r][c],P1->co[l][R][C],sl->P->co[l][R][C],q);
						}
					}
				}
			}
		}

		//calculates max increment in Psi, so that it is larger than a prefixed value Dt is reduced
		for(c=1;c<=Nc;c++){
			for(r=1;r<=Nr-1;r++){
				if(land->LC->co[r][c]!=NoV){
					for(l=1;l<=Nl;l++){
						if(*DeltaPsiMax<(P1->co[l][r][c]-sl->P->co[l][r][c])) *DeltaPsiMax=P1->co[l][r][c]-sl->P->co[l][r][c];
					}
				}
			}
		}

	}else{	//1D simulations

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					sy=sl->type->co[r][c];
					for(l=1;l<=Nl;l++){
						Kh=K(sl->P->co[l][r][c], sl->pa->co[sy][jKh][l], par->imp, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l],
								sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], sl->pa->co[sy][jv][l],
								fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), sl->T->co[l][r][c]);
						i=top->dz_dx->co[r][c];	
						q=Kh*i*sl->pa->co[sy][jdz][l]*1.0E-6;	//[m3/s]	
						i=top->dz_dy->co[r][c];	
						q+=Kh*i*sl->pa->co[sy][jdz][l]*1.0E-6;	//[m3/s]
						//printf("l:%ld q:%f i:%f",l,q,i);	
						set_psi_single(P1, sl, &q, Dt, l, r, c, par->psimin, par->Esoil);
						//printf("Pante:%f Ppost:%f\n",sl->P->co[l][r][c],P1->co[l][r][c]);
						//if(fabs(q)>1.0E-6) stop_execution();
						if(*DeltaPsiMax<(P1->co[l][r][c]-sl->P->co[l][r][c])) *DeltaPsiMax=P1->co[l][r][c]-sl->P->co[l][r][c];
					}
				}
			}
		}
	}

}


/*--------------------------------------------*/
void subflow_channel(double Dt, double Dtmax, CHANNEL *cnet, SOIL *sl, PAR *par){

	short sy;
	long ch,l,s,r,c;
	double q, Kh, dP, dpixel=0.5*(UV->U->co[1]+UV->U->co[2]);

	/* Note: cnet->q_sup can be negative, this means that the channel cells are recharging the
	surrounding pixels, but a channel cell cannot recharge more than the amount of water contained in it:
	!no control for this! */
	initialize_doublevector(cnet->Q, 0.0);

	for(ch=1;ch<=cnet->r->nh;ch++){

		for(l=1;l<=Nl;l++){

			r=cnet->r->co[ch];
			c=cnet->c->co[ch];

			sy=sl->type->co[r][c];

			Kh=K(sl->P->co[l][r][c], sl->pa->co[sy][jKh][l], par->imp, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
				sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], sl->pa->co[sy][jv][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),
				sl->T->co[l][r][c]);

			dP=fmax(sl->P->co[l][r][c] - psi_saturation(sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
				sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l]), 0.0);

			//q[m3/s] q=2*Kh*(sin(slope)+4*dP/dn)*D*dpixel
			//slope=slope to the channel (q1)
			//dn=fraction of the pixel covered by land
			//dP=difference f pressure between land and channel
			//parameters q1=slope of the land converging to the channel - q2=fraction of the pixel occupied by the channel
			q=2*Kh*(par->q1+dP/(1.0E3*0.25*dpixel*(1.0-par->q2)))*sl->pa->co[sy][jdz][l]*dpixel*1.0E-6;

			set_psi_single(sl->P, sl, &q, Dt, l, r, c, par->psimin, par->Esoil);

			cnet->Q->co[ch]+=q;
		}
	}

	for(ch=1;ch<=cnet->r->nh;ch++){
		for(s=1;s<=cnet->Q_sub_spread->nh;s++){
			cnet->Q_sub_spread->co[s]+=(cnet->Q->co[ch]*Dt/Dtmax)*(cnet->fraction_spread->co[ch][s]); /*in mc/s*/
		}
	}
}


/*--------------------------------------------*/
void routing(CHANNEL *cnet){

	long s;

	/*The just calculated Q_sup_spread and Q_sub_spread are added to the flow already
	present in the channel-network(Q_sup_s and Q_sub_s) and at once it is made a
	traslation of Q_sup_s e Q_sub_s towards the outlet:*/
	for(s=1;s<=cnet->Q_sub_spread->nh-1;s++){
		cnet->Q_sub_s->co[s]=cnet->Q_sub_s->co[s+1]+cnet->Q_sub_spread->co[s];
		cnet->Q_sup_s->co[s]=cnet->Q_sup_s->co[s+1]+cnet->Q_sup_spread->co[s];
	}
	cnet->Q_sub_s->co[cnet->Q_sub_spread->nh]=cnet->Q_sub_spread->co[cnet->Q_sub_spread->nh];
	cnet->Q_sup_s->co[cnet->Q_sup_spread->nh]=cnet->Q_sup_spread->co[cnet->Q_sup_spread->nh];

}

/*--------------------------------------------*/
void output_waterbalance(double Dt, WATER *wat, SOIL *sl, PAR *par, DOUBLEMATRIX *Z){

	long i,l,r,c;

	for(i=1;i<=par->chkpt->nrh;i++){
		r=par->rc->co[i][1];
		c=par->rc->co[i][2];

		wat->out1->co[5][i]+=wat->Pn->co[r][c]*Dt;
		wat->out1->co[6][i]+=(wat->Pn->co[r][c]-sl->Jinf->co[r][c])*Dt;
		wat->out1->co[7][i]+=(wat->h_sup->co[r][c] - wat->out1->co[9][i] - (wat->Pn->co[r][c]-sl->Jinf->co[r][c])*Dt);    /*positive if entering into the pixel*/
		wat->out1->co[8][i]+=sl->J->co[Nl][r][c]*Dt;
		wat->out1->co[9][i]=wat->h_sup->co[r][c];

		for(l=1;l<=Nl;l++){
			wat->out1->co[10][i]-=0.0*Dt;/*positive if entering into the pixel*/
		}
	}

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(Z->co[r][c]!=NoV){ /*if the pixel is not a novalue*/
				wat->out2->co[3]+=wat->Pn->co[r][c]*Dt;								/*[mm]*/
				wat->out2->co[4]+=sl->Jinf->co[r][c]*Dt;							/*[mm]*/
			}
		}
	}

}

/*--------------------------------------------*/
//Q(from r1,c1 to r2,c2) in [m3/s]
//D in [m]

void set_psi(DOUBLETENSOR *P1, SOIL *sl, double *Q, double dt, long l, long r1, long c1, long r2, long c2, double psimin, double Esoil){

	long rout,cout,rin,cin,r,c;
	double dx,dy,D,V,theta,thetamin,q=*Q;
	short sy,sgnchg=0;

	long ri=0, ci=0;

	dx=UV->U->co[1];	//[m]
	dy=UV->U->co[2];	//[m]
	D=0.001*sl->pa->co[sl->type->co[r1][c1]][jdz][l];	//layer thickness [m]
	V=D*dx*dy;			//[m3]=cell volume

	if(q>0){
		rout=r1;
		cout=c1;
		rin=r2;
		cin=c2;
	}else{
		rout=r2;
		cout=c2;
		rin=r1;
		cin=c1;
		q*=(-1);
		sgnchg=1;
	}

	r=rout;
	c=cout;
	sy=sl->type->co[r][c];

	theta=teta_psi(P1->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), Esoil);

	thetamin=teta_psi(fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), Esoil);

	if(theta-q*dt/V<thetamin){
		q=(theta-thetamin)*V/dt;
	}

	P1->co[l][r][c]=psi_teta(theta-q*dt/V, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), Esoil);

	if(r==ri && c==ci) printf("LIM 1. th:%f thmin:%f q:%f Pin:%f Pf:%f\n",theta,thetamin,q,P1->co[l][r][c],P1->co[l][r][c]);


	r=rin;
	c=cin;
	sy=sl->type->co[r][c];

	theta=teta_psi(P1->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), Esoil);

	P1->co[l][r][c]=psi_teta(theta+q*dt/V, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), Esoil);

	if(r==ri && c==ci) printf("LIM 2. th:%f thmin:%f q:%f Pin:%f Pf:%f\n",theta,thetamin,q,P1->co[l][r][c],P1->co[l][r][c]);

	if(sgnchg==0){
		*Q=q;
	}else{
		*Q=-q;
	}

}

/*--------------------------------------------*/

void set_psi_single(DOUBLETENSOR *psi, SOIL *sl, double *q, double dt, long l, long r, long c, double psimin, double Esoil){

	double dx, dy, D, V, theta, thetamin;
	short sy;

	sy=sl->type->co[r][c];

	dx=UV->U->co[1];	//[m]
	dy=UV->U->co[2];	//[m]
	D=0.001*sl->pa->co[sy][jdz][l];	//layer thickness [m]
	V=D*dx*dy;			//[m3]=cell volume

	theta=teta_psi(psi->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
		1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), Esoil);

	thetamin=teta_psi(fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), Esoil);

	if(theta-(*q)*dt/V<thetamin){
		*q=(theta-thetamin)*V/dt;
		if(*q<0) *q=0.0;
	}
	if(r==25 && c==12) printf("set %ld %f %f %f \n",l,psi->co[l][r][c],theta,*q);


	psi->co[l][r][c]=psi_teta(theta-(*q)*dt/V, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
		1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), Esoil);

	if(r==25 && c==12) printf("set %ld %f\n",l,psi->co[l][r][c]);

}
