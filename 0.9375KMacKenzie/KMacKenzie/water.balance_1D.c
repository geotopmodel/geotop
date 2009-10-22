
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
#include "water.balance_1D.h"
#include "pedo.funct.h"
#include "util_math.h"

extern T_INIT *UV;
extern STRINGBIN *files;
extern char *error_file_name;
extern long Nl, Nr, Nc;
extern double NoV;

/*--------------------------------------------*/
void water_balance_1D(TOPO *top, SOIL *sl, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, double time){

	double Dt, DPsiMax=0.0, te, tb, dt;
	long i, n, r, c, l, rref, cref, lref;
	FILE *f;// file pointer
	DOUBLETENSOR *P1, *P2;

	P1=new_doubletensor(Nl,Nr,Nc);
	P2=new_doubletensor(Nl,Nr,Nc);

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

				subflow(dt, land, top, sl, par, wat, &DPsiMax, &rref, &cref, &lref, P1, P2);

				if(DPsiMax>par->Dpsi) n++;	//decrease time step

			}while(dt/(double)n>par->dtmin && DPsiMax>par->Dpsi);

			//f=fopen(files->co[ferr]+1,"a");
			//fprintf(f,"\nWATER BALANCE: tb:%f dt:%f DPsiMax:%f maxadmitted:%f\n",tb,dt,DPsiMax,par->Dpsi);
			//fclose(f);

			if(DPsiMax>par->Dpsi){
				printf("NOT ABLE TO REDUCE MAX PSI INCREMENT AFTER SUBSURFACE FLOW BELOW MAX ADMITTED VALUE, with dt:%f s\n",dt);
				printf("MAX PSI increment after lateral subsurface flow: %f, max admitted:%f\n",DPsiMax,par->Dpsi);
				/*printf("r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref,cref,lref,sl->P->co[lref][rref][cref],P1->co[lref][rref][cref]);
				printf("r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref+1,cref,lref,sl->P->co[lref][rref+1][cref],P1->co[lref][rref+1][cref]);
				printf("r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref-1,cref,lref,sl->P->co[lref][rref-1][cref],P1->co[lref][rref-1][cref]);
				printf("r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref,cref+1,lref,sl->P->co[lref][rref][cref+1],P1->co[lref][rref][cref+1]);
				printf("r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref,cref-1,lref,sl->P->co[lref][rref][cref-1],P1->co[lref][rref][cref-1]);*/

				f=fopen(error_file_name,"a");
				fprintf(f,"NOT ABLE TO REDUCE MAX PSI INCREMENT AFTER SUBSURFACE FLOW BELOW MAX ADMITTED VALUE, with dt:%f s\n",dt);
				fprintf(f,"MAX PSI increment after lateral subsurface flow: %f, max admitted:%f\n",DPsiMax,par->Dpsi);
				/*fprintf(f,"r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref,cref,lref,sl->P->co[lref][rref][cref],P1->co[lref][rref][cref]);
				fprintf(f,"r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref+1,cref,lref,sl->P->co[lref][rref+1][cref],P1->co[lref][rref+1][cref]);
				fprintf(f,"r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref-1,cref,lref,sl->P->co[lref][rref-1][cref],P1->co[lref][rref-1][cref]);
				fprintf(f,"r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref,cref+1,lref,sl->P->co[lref][rref][cref+1],P1->co[lref][rref][cref+1]);
				fprintf(f,"r:%ld c:%ld l:%ld Pbefore:%f Pafter:%f\n",rref,cref-1,lref,sl->P->co[lref][rref][cref-1],P1->co[lref][rref][cref-1]);*/

				fclose(f);

				stop_execution();
			}

			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(land->LC->co[r][c]!=NoV){
						for(l=1;l<=Nl;l++){
							sl->P->co[l][r][c]=P1->co[l][r][c];

							if(sl->P->co[l][r][c]!=sl->P->co[l][r][c]){
								printf("2.NOvalue Psi %ld %ld %ld\n",l,r,c);
								stop_execution();
							}

						}
					}
				}
			}

			//subflow_channel(dt, Dt, cnet, sl, par);

			vertical_water_balance(dt, land->LC, sl, wat, time, par);

			supflow(dt, Dt, top, land, wat, cnet, par);

			te=tb+dt;

		}while(te<Dt);

		routing(cnet);

		output_waterbalance(Dt, wat, sl, par, land->LC);

	}

	free_doubletensor(P1);
	free_doubletensor(P2);

}


/*--------------------------------------------*/

void vertical_water_balance(double Dt, DOUBLEMATRIX *Z, SOIL *sl, WATER *wat, double time, PAR *par)

{

	long l;					/* the index of depth is l not d because d is used for sl-thickness vector*/
	long r,c;				/* counters of rows and columns*/
	long i;
	DOUBLEVECTOR *psi;
	short sy;
	double masserrorbasin=0.0, masserror;
	double Pnbasin=0.0, Infbasin=0.0;
	//double DW, th0, th1, theq, h, psisat, oversat, ice0, T0;
	FILE *f;                    /* file which contains the errors of simulation*/

	psi=new_doublevector(Nl);

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(Z->co[r][c]!=NoV){

				sy=sl->type->co[r][c];

				/*control of net precipitation:*/
				if (wat->Pn->co[r][c]<-0.01){
					if(par->n_error<par->max_error){
						par->n_error++;
						f=fopen(error_file_name,"a");
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
				Richards_1D(Dt, r, c, sl, psi, &(wat->h_sup->co[r][c]), wat->Pn->co[r][c], time, par, &masserror, wat->q_sub);

 				if(sl->Jinf->co[r][c]!=sl->Jinf->co[r][c]){
					printf("ERROR ON INFILTRATION r:%ld c:%ld\n",r,c);
					stop_execution();
				}

				wat->error->co[r][c]+=masserror;
				masserrorbasin+=masserror/(double)par->total_pixel;
				Pnbasin+=wat->Pn->co[r][c]*Dt/(double)par->total_pixel;
				Infbasin+=sl->Jinf->co[r][c]*Dt/(double)par->total_pixel;

				/*write q_out for specified pixels*/
				for(i=1;i<=par->chkpt->nrh;i++){
					wat->out1->co[28][i]=masserror*3600.0/Dt;
					if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){
						wat->out1->co[15][i]+=sl->J->co[1][r][c];
					}
				}


				/*update P*/
				for(l=1;l<=Nl;l++){
					/*if(sl->T->co[l][r][c]<=Tfreezing && par->en_balance==1){

						psisat=psi_saturation(sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
							1-1/sl->pa->co[sy][jns][l]);

						th1=teta_psi(Fmin(psi->co[l],psisat),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
							sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil);

						theq=teta_psi(Psif(sl->T->co[l][r][c]),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],
							1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil);

						if(fabs(th1-theq)>1.E-8){

							th0=teta_psi(Fmin(sl->P->co[l][r][c],psisat),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil);

							oversat=teta_psi(psi->co[l],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil)
								-th1;

							DW=(th1-th0)*(sl->pa->co[sy][jdz][l]);
							h=internal_energy_soil(th0, sl->thice->co[l][r][c], sl->T->co[l][r][c], sl->pa->co[sy][jdz][l], sl->pa->co[sy][jct][l], sl->pa->co[sy][jsat][l]);

							ice0=sl->thice->co[l][r][c];
							T0=sl->T->co[l][r][c];

							from_internal_soil_energy(r, c, l, h+Lf*DW, &th1, &(sl->thice->co[l][r][c]), &(sl->T->co[l][r][c]), sl->pa->co[sy], par->psimin);

							sl->P->co[l][r][c]=psi_teta(th1+oversat,sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil);

						}else{
							sl->P->co[l][r][c]=psi->co[l];
						}

					}else{*/
						sl->P->co[l][r][c]=psi->co[l];
					//}
				}

				//error check
				if (wat->h_sup->co[r][c]<0.0){
					if (wat->h_sup->co[r][c]<-1E-8){
						if(par->n_error<par->max_error){
							par->n_error++;
							f=fopen(error_file_name,"a");
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
						f=fopen(error_file_name,"a");
						fprintf(f,"h_sup[%ld,%ld], at the end of vertical_water_balance subroutine, is not a number!\n",r,c);
						fclose(f);
					}
				}

			}
		}
	}

	free_doublevector(psi);

	wat->out2->co[8]+=masserrorbasin*3600.0/Dt;

	f=fopen(error_file_name,"a");
	fprintf(f,"\nERROR MASS BALANCE: %20.18fmm/h - Pnet: %20.18fmm/h - Infiltration: %20.18fmm/h\n",masserrorbasin*3600.0/Dt,Pnbasin*3600.0/Dt,Infbasin*3600.0/Dt);
	fclose(f);

}


/*-----------------------------------------------------------------------------------------------------------------------------*/
void supflow(double Dt, double Dtmax, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par)
{
	long r,c,R,C,ch,s;                                    /* counters*/
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; /* differential of number-pixel for rows and*/
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; /* columns, depending on Drainage Directions*/
	double dx,dy;                                         /* the two dimensions of a pixel*/
	double Ks;											  /* the Strickler's coefficent calculated with a standard deviation*/
	//double q_sup_max;                                     /* maximum superficial flow (in the first part) or the incoming flow in a channel pixel (in the second part)*/
	double b[10];                                         /* area perpendicular to the superficial flow divided by h_sup*/
	double i;											  /*hydraulic gradient*/
	double q,tb,te,dt;
	short lu;


if(par->point_sim==0){	//distributed simulations

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
				if(wat->h_sup->co[r][c]>0 && top->DD->co[r][c]>=1 && top->DD->co[r][c]<=8){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
					//i=top->i_DD->co[r][c];
					if(i<0) i=0;  //correct false slopes
					q=b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;	//mm/s
					//q=wat->h_sup->co[r][c]/dt;
					if(q>0){
						if(wat->h_sup->co[r][c]/q<dt) dt=wat->h_sup->co[r][c]/q;
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
		initialize_doublematrix(wat->q_sup,0.0);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(wat->h_sup->co[r][c]>0 && top->DD->co[r][c]>=1 && top->DD->co[r][c]<=8){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
					//i=top->i_DD->co[r][c];
					if(i<0) i=0.0;
					wat->q_sup->co[r][c]=b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;
					//if(wat->h_sup->co[r][c]>0)printf("%ld %ld h:%f q:%f\n",r,c,wat->h_sup->co[r][c],wat->q_sup->co[r][c]);
					//wat->q_sup->co[r][c]=wat->h_sup->co[r][c]/dt;
					if(wat->q_sup->co[r][c]!=wat->q_sup->co[r][c]){
						printf("NO VALUE SUP:%ld %ld %ld %ld %f %f\n",r,c,r+r_DD[top->DD->co[r][c]],c+c_DD[top->DD->co[r][c]],wat->h_sup->co[r][c],wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]]);
						printf("i:%f Dh:%f Ks:%f pow:%f iF:%f\n",i,wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]],Ks,pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m),top->i_DD->co[r][c]);
					}

				}else if(top->DD->co[r][c]==10){
					wat->q_sup->co[r][c]=wat->h_sup->co[r][c]/dt; /*[mm/s]*/
					//if(wat->h_sup->co[r][c]>0)printf("ch. %ld %ld h:%f q:%f\n",r,c,wat->h_sup->co[r][c],wat->q_sup->co[r][c]);
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
							if(wat->q_sup->co[r][c]!=wat->q_sup->co[r][c]) printf("A\n");
						}
						//if(wat->q_sup->co[r][c]>0)printf("corr %ld %ld qdt:%f h:%f dt:%f\n",r,c,wat->q_sup->co[r][c]*dt,wat->h_sup->co[r][c],dt);
					}

					/*the superficial flow is added to the land pixels (code 0):*/
					if (top->pixel_type->co[R][C]==0) wat->h_sup->co[R][C]+=wat->q_sup->co[r][c]*dt;

					/*the superficial flow is added to flow which flows into the channel pixels (code 10):*/
					if (top->pixel_type->co[R][C]==10){
						for(ch=1;ch<=cnet->r->nh;ch++){
							if(R==cnet->r->co[ch] && C==cnet->c->co[ch]){
								cnet->Q->co[ch]+=wat->q_sup->co[r][c]*dt/par->Dt;
								if(cnet->Q->co[ch]!=cnet->Q->co[ch]){
									printf("qsup no value: r:%ld c:%ld ch:%ld R:%ld C:%ld qsup:%f hsup:%f\n",r,c,ch,R,C,wat->q_sup->co[r][c],wat->h_sup->co[r][c]);
								}
								//if(wat->q_sup->co[r][c]>0) printf("channel r:%ld c:%ld ch:%ld q:%f %f dt:%f\n",R,C,ch,wat->q_sup->co[r][c],cnet->Q->co[ch],dt);
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
			wat->h_sup->co[r][c]=Fmax(wat->h_sup->co[r][c], 0.0);
			//q_sup_max=Fmin(2.0*b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(top->i_ch->co[r][c])*1000.0/dx/dy,cnet->Q->co[ch]);
			//q_sup_max=cnet->Q->co[ch];
			//if (top->DD->co[r][c]==10) q_sup_max=cnet->Q->co[ch];	//closure section
			//if(q_sup_max!=q_sup_max) printf("q_sup_max no value: qsup:%f r:%ld c:%ld ch:%ld h:%f ich:%f\n",cnet->Q->co[ch],r,c,ch,wat->h_sup->co[r][c],top->i_ch->co[r][c]);

			/*wat of channel pixel which becames lateral flow and doesn't go into the channel;
			the variable cnet->q_sup is use in this case improperly:*/
			//cnet->Q->co[ch]-=q_sup_max;

			/*the water on surface in channel pixel is updated:*/
			//wat->h_sup->co[r][c]-=wat->q_sup->co[r][c]*dt;
			//wat->h_sup->co[r][c]+=cnet->Q->co[ch]*dt;

			/*now is assingned to cnet->q_sup is right value:*/
			//cnet->Q->co[ch]=q_sup_max;

			/*The new water arrived in the channel-network is spread immediately along a
			virtual coordinate s (=distance from the outlet) so that it can travel to the
			outlet, after simple traslation, compling with the convolution of the De
			Saint Venant equations. This is made either for superficial flow or for the
			subsuperficial flow:*/
			for(s=1;s<=cnet->Q_sup_spread->nh;s++){
				cnet->Q_sup_spread->co[s]+=(cnet->Q->co[ch])*(cnet->fraction_spread->co[ch][s])*0.001*dx*dy; /*in mc/s*/
				//if(cnet->Q_sup_spread->co[s]>0)printf("s:%ld Q:%f\n",s,cnet->Q_sup_spread->co[s]);
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
void subflow(double Dt, LAND *land, TOPO *top, SOIL *sl, PAR *par, WATER *wat, double *DeltaPsiMax, long *rref, long *cref, long *lref,
	DOUBLETENSOR *P1, DOUBLETENSOR *P2){

	long l;  /*the index of depth is l not d because d is used for the vector of the sl-thickness*/
	long r,c,R,C;/*counters of rows and columns*/
	short sy, sy1, flux; //sl type
	double dx,dy; /*the two dimensions of a pixel*/
	double ds,dn;
	double q;    /*auxiliar variable to calculate the subflow for unit of land-surface*/
	double i; //lateral gradient of pressure [-]
	double psisat, psisat1, Kh, Kh1;

	//debug
	//long ri=0, ci=0;

	*DeltaPsiMax=0.0;

	dx=UV->U->co[1];	//m
	dy=UV->U->co[2];	//m

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				for(l=1;l<=Nl;l++){

					if(sl->P->co[l][r][c]!=sl->P->co[l][r][c]){
						printf("1.NOvalue Psi %ld %ld %ld\n",l,r,c);
						stop_execution();
					}

				}
			}
		}
	}

	initialize_doubletensor(wat->q_sub, 0.0);

	if(par->point_sim==0){	//2D simulations

		for(r=1;r<=Nr;r++){

			for(c=1;c<=Nc-1;c++){

				R=r;
				C=c+1;
				ds=dx;
				dn=dy;

				if((top->pixel_type->co[r][c]==0 || top->pixel_type->co[r][c]==10) && (top->pixel_type->co[R][C]==0 || top->pixel_type->co[R][C]==10)){

					sy=sl->type->co[r][c];
					sy1=sl->type->co[R][C];

					for(l=1;l<=Nl;l++){

						flux=0;
						if( sl->pa->co[sy][jlatfl][l]==1 || sl->pa->co[sy1][jlatfl][l]==1 ) flux=1;
						if( sl->pa->co[sy][jlatfl][l]==2 || sl->pa->co[sy1][jlatfl][l]==2 ) flux=2;

						psisat=-20.0;
						psisat1=-20.0;

						if( flux==2 || (flux==1 && (sl->P->co[l][r][c]>=psisat || sl->P->co[l][R][C]>=psisat1)) ){

							Kh=K(sl->P->co[l][r][c], sl->pa->co[sy][jKh][l], par->imp, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l],
									sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], sl->pa->co[sy][jv][l],
									sl->pa->co[sy][jpsimin][l], sl->T->co[l][r][c]);
							Kh1=K(sl->P->co[l][R][C], sl->pa->co[sy1][jKh][l], par->imp, sl->thice->co[l][R][C], sl->pa->co[sy1][jsat][l],
									sl->pa->co[sy1][jres][l], sl->pa->co[sy1][ja][l], sl->pa->co[sy1][jns][l], 1-1/sl->pa->co[sy1][jns][l], sl->pa->co[sy1][jv][l],
									sl->pa->co[sy1][jpsimin][l], sl->T->co[l][R][C]);

							i=top->dz_dx->co[r][c]+0.001*(sl->P->co[l][r][c]-sl->P->co[l][R][C])/ds;

							q=Harmonic_Mean(ds, ds, Kh, Kh1)*i;	//[mm/s]

							check_q_2(sl->P, P1, &q, sl, Dt, l, r, c, R, C, par->psimin, par->Esoil);

							wat->q_sub->co[l][r][c] += q;
							wat->q_sub->co[l][R][C] -= q;

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

				if((top->pixel_type->co[r][c]==0 || top->pixel_type->co[r][c]==10) && (top->pixel_type->co[R][C]==0 || top->pixel_type->co[R][C]==10)){

					sy=sl->type->co[r][c];
					sy1=sl->type->co[R][C];

					for(l=1;l<=Nl;l++){

						flux=0;
						if( sl->pa->co[sy][jlatfl][l]==1 || sl->pa->co[sy1][jlatfl][l]==1 ) flux=1;
						if( sl->pa->co[sy][jlatfl][l]==2 || sl->pa->co[sy1][jlatfl][l]==2 ) flux=2;

						psisat=-20.0;
						psisat1=-20.0;

						if( flux==2 || (flux==1 && (sl->P->co[l][r][c]>=psisat || sl->P->co[l][R][C]>=psisat1)) ){

							Kh=K(sl->P->co[l][r][c], sl->pa->co[sy][jKh][l], par->imp, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l],
									sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], sl->pa->co[sy][jv][l],
									sl->pa->co[sy][jpsimin][l], sl->T->co[l][r][c]);

							Kh1=K(sl->P->co[l][R][C], sl->pa->co[sy1][jKh][l], par->imp, sl->thice->co[l][R][C], sl->pa->co[sy1][jsat][l],
									sl->pa->co[sy1][jres][l], sl->pa->co[sy1][ja][l], sl->pa->co[sy1][jns][l], 1-1/sl->pa->co[sy1][jns][l], sl->pa->co[sy1][jv][l],
									sl->pa->co[sy1][jpsimin][l], sl->T->co[l][R][C]);

							i=top->dz_dy->co[r][c]+0.001*(sl->P->co[l][r][c]-sl->P->co[l][R][C])/ds;

							q=Harmonic_Mean(ds, ds, Kh, Kh1)*i;	//[mm/s]

							check_q_2(P1, P2, &q, sl, Dt, l, r, c, R, C, par->psimin, par->Esoil);

							wat->q_sub->co[l][r][c] += q;
							wat->q_sub->co[l][R][C] -= q;
						}
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
								sl->pa->co[sy][jpsimin][l], sl->T->co[l][r][c]);
						i=top->dz_dx->co[r][c]+top->dz_dy->co[r][c];
						q=Kh*i;	//[mm/s]

						//if(l==1)printf("->%ld q:%e i:%e Kh:%e",l,q,i,Kh);

						check_q_1(sl->P, P2, &q, sl, Dt, l, r, c, par->psimin, par->Esoil, par->Dpsi);

						wat->q_sub->co[l][r][c] += q;

						//if(l==1)printf("->%ld q:%e qsub:%e",l,q,wat->q_sub->co[l][r][c]);

					}
				}
			}
		}
	}

	//calculates max increment in Psi, so that it is larger than a prefixed value Dt is reduced
	for(c=1;c<=Nc;c++){
		for(r=1;r<=Nr;r++){
			if(land->LC->co[r][c]!=NoV){
				for(l=1;l<=Nl;l++){

					set_psi(P1, wat->q_sub, sl, Dt, l, r, c, par->psimin, par->Esoil);

					if(*DeltaPsiMax<(P1->co[l][r][c]-sl->P->co[l][r][c])){
						*DeltaPsiMax=P1->co[l][r][c]-sl->P->co[l][r][c];
						*rref=r;
						*cref=c;
						*lref=l;
					}
				}
			}
		}
	}

}


/*--------------------------------------------*/
void subflow_channel(double Dt, double Dtmax, CHANNEL *cnet, SOIL *sl, PAR *par){

//	short sy; /* commented by Emanuele Cordano on 24/9/9 */
	long ch,s; //l,s,r,c;
	//double q, Kh, dP, dpixel=0.5*(UV->U->co[1]+UV->U->co[2]);

	/* Note: cnet->q_sup can be negative, this means that the channel cells are recharging the
	surrounding pixels, but a channel cell cannot recharge more than the amount of water contained in it:
	!no control for this! */
	initialize_doublevector(cnet->Q, 0.0);

	/*for(ch=1;ch<=cnet->r->nh;ch++){

		for(l=1;l<=Nl;l++){

			r=cnet->r->co[ch];
			c=cnet->c->co[ch];

			sy=sl->type->co[r][c];

			Kh=K(sl->P->co[l][r][c], sl->pa->co[sy][jKh][l], par->imp, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
				sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], sl->pa->co[sy][jv][l], sl->pa->co[sy][jpsimin][l], sl->T->co[l][r][c]);

			dP=Fmax(sl->P->co[l][r][c] - psi_saturation(sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
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
	}*/

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

void check_q_1(DOUBLETENSOR *Pbeg, DOUBLETENSOR *Pend, double *Q, SOIL *sl, double dt, long l, long r, long c, double psimin, double Esoil, double max_dpsi){

	double thetamin, q=*Q, theta0, theta1, e=1.E-3;
	short sy; //, sgnchg=0; /*commented by Emanuele Cordano on 24/9/09 */

	sy=sl->type->co[r][c];

	theta0=teta_psi(Pbeg->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	thetamin=teta_psi(sl->pa->co[sy][jpsimin][l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	//if(l==1) printf("th:%f thmin:%f ",theta0,thetamin);

	if(theta0 <= thetamin){
		q = 0.0;
	}else if(theta0 - q*dt/sl->pa->co[sy][jdz][l] < thetamin){
		q = (theta0-thetamin)*sl->pa->co[sy][jdz][l]/dt;
	}

	theta1 = theta0 - q*dt/sl->pa->co[sy][jdz][l];

	Pend->co[l][r][c]=psi_teta(theta1, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	if(Pend->co[l][r][c]-Pbeg->co[l][r][c] > max_dpsi - e){
		Pend->co[l][r][c] = Pbeg->co[l][r][c] + max_dpsi - e;
		theta1 = teta_psi(Pend->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
			sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);
		q = (theta0 - theta1)*sl->pa->co[sy][jdz][l]/dt;
	}

	//if(l==1) printf("q:%e theta:%f Pend:%f\n",q,theta1,Pend->co[l][r][c]);

	*Q = q;

}

/*--------------------------------------------*/

void check_q_2(DOUBLETENSOR *Pbeg, DOUBLETENSOR *Pend, double *Q, SOIL *sl, double dt, long l, long r1, long c1, long r2, long c2,
	double psimin, double Esoil){

	long rout, cout, rin, cin, r, c;
	double thetamin, q=*Q, theta;
	short sy, sgnchg=0;

//	long ri=0, ci=0;

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

	theta=teta_psi(Pbeg->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	thetamin=teta_psi(sl->pa->co[sy][jpsimin][l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	if(theta <= thetamin){
		q = 0.0;
	}else if(theta - q*dt/sl->pa->co[sy][jdz][l] < thetamin){
		q = (theta-thetamin)*sl->pa->co[sy][jdz][l]/dt;
	}

	if(q<0) q=0.0;

	theta -= q*dt/sl->pa->co[sy][jdz][l];

	Pend->co[l][r][c]=psi_teta(theta, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);


	r=rin;
	c=cin;
	sy=sl->type->co[r][c];

	theta=teta_psi(Pbeg->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	theta += q*dt/sl->pa->co[sy][jdz][l];

	Pend->co[l][r][c]=psi_teta(theta, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
		sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	if(sgnchg==0){
		*Q=q;
	}else{
		*Q=-q;
	}

}

/*--------------------------------------------*/

void set_psi(DOUBLETENSOR *psi, DOUBLETENSOR *q, SOIL *sl, double dt, long l, long r, long c, double psimin, double Esoil){

	double theta;
	short sy;
	FILE *f;

	double p;

	sy=sl->type->co[r][c];

	theta=teta_psi(sl->P->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	p=sl->P->co[l][r][c];

	//if(l==1)printf("...l:%ld P:%f theta:%f ",l,sl->P->co[l][r][c],theta);

	theta -= q->co[l][r][c]*dt/sl->pa->co[sy][jdz][l];

	psi->co[l][r][c]=psi_teta(theta, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], psimin, Esoil);

	//if(l==1)printf(";;;l:%ld P:%f theta:%f n:%f q:%e",l,psi->co[l][r][c],theta,sl->pa->co[sy][jsat][l]-sl->thice->co[l][r][c],q->co[l][r][c]);

	if(fabs(p-psi->co[l][r][c])>12000) {
		printf("ATTENTION!!!!!!! water.balance_1D.c: P_init - P_fin >12000");
		f=fopen(error_file_name,"a");
		fprintf(f,"ATTENTION!!!!!!! water.balance_1D.c: P_init - P_fin >12000");
		fclose(f);
		stop_execution();
	}

}


/*--------------------------------------------*/




void Richards_1D(double Dt, long r, long c, SOIL *sl, DOUBLEVECTOR *psi, double *h, double Pnet, double t, PAR *par, double *masslosscum, DOUBLETENSOR *q_sub)

{
	double Inf,massloss; //,K1,Kdw,K1dw,dz,ddw,dzdw; /* commented by Emanuele Cordano (unuseful variable)*/
	DOUBLEVECTOR *e1;
	long l; //, m;/* commented by Emanuele Cordano (unuseful variable)*/
//	short sy=sl->type->co[r][c]; /* commented by Emanuele Cordano (unuseful variable)*/

	*masslosscum=0.0;

	e1=new_doublevector(Nl);

	//initialization of the fluxes
	sl->Jinf->co[r][c]=0.0;
	for(l=1;l<=Nl;l++){
		sl->J->co[l][r][c]=0.0;
	}

	Inf=Pnet+(*h)/Dt;

	solve_Richards_1D(r, c, Dt, &Inf, *h, psi, e1, &massloss, sl, par, q_sub);

	*masslosscum = *masslosscum + massloss;

	sl->Jinf->co[r][c]+=Inf;

	*h+=(Pnet-Inf)*Dt; /*in mm*/
	if(*h<0 && *h>-1.E-3) *h=0.0;
	if(*h<0){
		printf("hsup negative:%e %ld %ld Pnet:%f Inf:%f\n",*h,r,c,Pnet,Inf);
		stop_execution();
		*h=0.0;
	}


	//fluxes update
	/*for(l=1;l<=Nl-1;l++){
		m=l;
		K1=K(e1->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
			sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
		dz=sl->pa->co[sy][jdz][m];

		m=l+1;
		Kdw=K(e1->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
			sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
		ddw=sl->pa->co[sy][jdz][m];

		dzdw=0.5*dz+0.5*ddw;
		K1dw=Harmonic_Mean(dz, ddw, K1, Kdw);

		sl->J->co[l][r][c]+=(K1dw*(e1->co[l]-e1->co[l+1])/dzdw + K1dw);
	}

	l=Nl;
	m=l;
	K1=K(e1->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
		sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);

	sl->J->co[Nl][r][c]+=par->f_bound_Richards*K1;*/

	//update psi
	for(l=1;l<=Nl;l++){
		psi->co[l]=e1->co[l];
	}

	/*Deallocate*/
	free_doublevector(e1);

}

/*--------------------------------------------*/


void find_coeff_Richards_1D(long r, long c, double dt, short *bc, double Inf, double h, DOUBLEVECTOR *psit, DOUBLEVECTOR *e0, DOUBLEVECTOR *adi,
	DOUBLEVECTOR *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, SOIL *sl, PAR *par, DOUBLETENSOR *q_sub){

	long l, m;
	short sy=sl->type->co[r][c];
	double dz, K1, Kdw, Kup, ddw, dup, dzdw, dzup, K1dw, K1up, C1, theta0, theta1;

	l=1;
	m=l;
	K1=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m], sl->pa->co[sy][ja][m],
		sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
	dz=sl->pa->co[sy][jdz][m];

	m=l+1;
	Kdw=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m], sl->pa->co[sy][ja][m],
		sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m],sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
	ddw=sl->pa->co[sy][jdz][m];

	dzdw=0.5*dz+0.5*ddw;

	K1dw=Harmonic_Mean(dz, ddw, K1, Kdw);


	C1=dteta_dpsi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);
	theta0=teta_psi(psit->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);
	theta1=teta_psi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);


	ad->co[l]= C1*dz/dt + K1dw/dzdw;
	ads->co[l]=-K1dw/dzdw;
	b->co[l]=(C1*e0->co[l]+theta0-theta1)*dz/dt - K1dw;

	if(Inf < K1 + K1*(h-e0->co[l])/(dz/2.)){
		*bc=0;	//Neumann
		b->co[l]+=Inf;
	}else{
		*bc=1;	//Dirichlet
		ad->co[l]+=K1/(dz/2.);
		b->co[l]+=(K1 + K1*h/(dz/2.));
	}

	//middle layers
	for(l=2;l<=Nl-1;l++){

		m=l;
		K1=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
			sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
		dz=sl->pa->co[sy][jdz][m];

		m=l+1;
		Kdw=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
			sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
		ddw=sl->pa->co[sy][jdz][m];

		m=l-1;
		Kup=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
			sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
		dup=sl->pa->co[sy][jdz][m];

		dzdw=0.5*dz+0.5*ddw;
		dzup=0.5*dz+0.5*dup;

		K1dw=Harmonic_Mean(dz, ddw, K1, Kdw);
		K1up=Harmonic_Mean(dz, dup, K1, Kup);


		C1=dteta_dpsi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
			sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);
		theta0=teta_psi(psit->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
			sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);
		theta1=teta_psi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
			sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);


		ad->co[l]=C1*dz/dt + (K1dw/dzdw+K1up/dzup);
		ads->co[l]=-K1dw/dzdw;
		adi->co[l-1]=-K1up/dzup;
		b->co[l]=(C1*e0->co[l]+theta0-theta1)*dz/dt + (K1up-K1dw);

	}


	//last layer
	l=Nl;

	m=l;
	K1=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
		sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
	dz=sl->pa->co[sy][jdz][m];

	m=l-1;
	Kup=K(e0->co[m], sl->pa->co[sy][jKv][m], par->imp, sl->thice->co[m][r][c], sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m],
		sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], sl->pa->co[sy][jv][m], sl->pa->co[sy][jpsimin][l], sl->T->co[m][r][c]);
	dup=sl->pa->co[sy][jdz][m];

	dzup=0.5*dz+0.5*dup;

	K1up=Harmonic_Mean(dz, dup, K1, Kup);


	C1=dteta_dpsi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);
	theta0=teta_psi(psit->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);
	theta1=teta_psi(e0->co[l], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin2, par->Esoil);

	ad->co[l]=C1*dz/dt + (K1up/dzup);
	adi->co[l-1]=-K1up/dzup;
	b->co[l]=(C1*e0->co[l]+theta0-theta1)*dz/dt + (K1up-K1*par->f_bound_Richards);

}


/*--------------------------------------------*/



void solve_Richards_1D(long r, long c, double dt, double *Inf, double h, DOUBLEVECTOR *psit, DOUBLEVECTOR *e1, double *massloss, SOIL *sl, PAR *par, DOUBLETENSOR *q_sub){

	double mass0, nw, massloss0, mass1, res, k, ActInf, q_lat;
	long cont, cont2, l;
	short sy=sl->type->co[r][c], bc=0, bc0=0, cond_out; /* modified bc and bc0 by Emanuele Cordano on 24/9/9 */
	DOUBLEVECTOR *e0, *adi, *ad, *ads, *b, *de;
	double max_error;
	ad=new_doublevector(Nl);
	adi=new_doublevector(Nl-1);
	ads=new_doublevector(Nl-1);
	b=new_doublevector(Nl);
	de=new_doublevector(Nl);
	e0=new_doublevector(Nl);

	*massloss=1.E99;
	//error=1.E99;
	if(par->superfast==1) max_error=1E-1;
	else max_error=1E-2;
	mass0=0.0;
	for(l=1;l<=Nl;l++){
		mass0+=sl->pa->co[sy][jdz][l]*teta_psi(psit->co[l],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],
		   sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2, par->Esoil);
	}

	cont=0;

	for(l=1;l<=Nl;l++){
		e1->co[l]=psit->co[l];
	}

	do{ //loop of Picard iteration

		for(l=1;l<=Nl;l++){
			e0->co[l]=e1->co[l];
		}

		find_coeff_Richards_1D(r, c, dt, &bc, *Inf, h, psit, e0, adi, ad, ads, b, sl, par, q_sub);
		tridiag(2,r,c,Nl,adi,ad,ads,b,e1);

		cont++;

		for(l=1;l<=Nl;l++){
			de->co[l]=e1->co[l]-e0->co[l];
		}

		cont2=0;
		nw=1.0;
		massloss0=(*massloss);
		//error0=error;

		do{
			for(l=1;l<=Nl;l++){
				e1->co[l]=e0->co[l]+nw*de->co[l];
			}

			k=K(e1->co[1], sl->pa->co[sy][jKv][1], par->imp, sl->thice->co[1][r][c], sl->pa->co[sy][jsat][1], sl->pa->co[sy][jres][1],
				sl->pa->co[sy][ja][1], sl->pa->co[sy][jns][1], 1-1/sl->pa->co[sy][jns][1], sl->pa->co[sy][jv][1], sl->pa->co[sy][jpsimin][1],
				sl->T->co[1][r][c]);

			ActInf=Fmin(*Inf, k + k*(h-e1->co[1])/(0.5*sl->pa->co[sy][jdz][1]));
			//printf("cont:%ld cont2:%ld %f %f %f\n",cont,cont2,ActInf,*Inf,k + k*(h-e1->co[1])/(0.5*sl->pa->co[sy][jdz][1]));

			mass1=0.0;
			q_lat=0.0;
			for(l=1;l<=Nl;l++){
				mass1+=sl->pa->co[sy][jdz][l]*teta_psi(e1->co[l],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],
					sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2, par->Esoil);
				//q_lat+=q_sub->co[l][r][c]*dt;
			}
			*massloss=mass1-mass0-ActInf*dt;

			if(ActInf>*Inf){
				printf("mass0:%f mass1:%f error:%e error0:%e\n\n\n",mass0,mass1,*massloss,massloss0);
				stop_execution();
			}

			cont2++;
			nw/=4.0;

		}while(fabs(*massloss)>fabs(massloss0) && cont2<5);

		res=0.0;
		for(l=1;l<=Nl;l++){
			if(e1->co[l]<par->psimin2) e1->co[l]=par->psimin2;
			res+=pow((e1->co[l]-e0->co[l]),2.0);
		}
		res=pow(res,0.5);

		cond_out=0;
		if(res<=par->TolVWb && fabs(*massloss)<=max_error) cond_out=1;

		}while(cont<par->MaxiterVWB && cond_out==0);

		if(fabs(*massloss)>max_error){
				printf("Error too high %ld %ld massloss:%f Inf:%f Infpot:%f lat:%f bc:%d bc0:%d psi1:%f e1:%f\n",r,c,*massloss,ActInf,*Inf,q_lat,bc,bc0,psit->co[1],e1->co[1]);
				//stop_execution();
			}


	*Inf=ActInf;

	free_doublevector(ad);
	free_doublevector(adi);
	free_doublevector(ads);
	free_doublevector(b);
	free_doublevector(de);
	free_doublevector(e0);

}

