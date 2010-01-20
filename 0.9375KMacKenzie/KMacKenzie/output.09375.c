
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Matteo Dall'Amico, Emanuele Cordano

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


#include "keywords_file.h"
#include "struct.geotop.09375.h"
#include "vegetation.h"
#include "pedo.funct.h"
#include "geo_statistic.09375.h"
#include "networks.h"
#include "rw_maps.h"
#include "constant.h"
#include "extensions.h"
#include "times.h"
#include "energy.balance.h"
#include "input.09375.h"
#include "output.09375.h"
#include "tabs.h"

#include "frost_table.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern char *error_file_name;
extern long Nl, Nr, Nc;
extern double NoV;


/****************************************************************************************************/
/* write_output: stampa i dati degli output variabili nel tempo                                     */
/****************************************************************************************************/
void write_output(TIMES *times, WATER *wat, CHANNEL *cnet, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac, METEO *met)

{
 /*internal auxiliary variables:*/
 long i,j,r,c=0,l,s; /*counters*/
 long n_file;      /*number of file of the type "TETAxySSSlZZ"(i.e. number of the basin-time-step)*/
 double t_i=0;         /*time of begin of an interval time for the output*/
 char SSSS[ ]={"SSSS"};
 char *name=NULL; /*modified by Emanuele Cordano on 24/9/9 */
 FILE *f;

 static double wt0_basin; /*mean intercepted precipitation [mm] in the previous output-basin Dt*/
 static double Ssup;      /*supercial Storage of water in all the basin [mm]*/
 static double Ssub;      /*subsuperficial Storage of water in all the basin [mm]*/
 static double Rout;      /*sum of the output flows from the last output-basin for unit of area[mm]*/
 static double R_G;
 static double S_ch0;     /*wat in the channel at the start of a basin-output step-time*/
 static double S_ch1;     /*wat in the channel at the end of a basin-output step-time*/
 static double Qsub_ch,Qsup_ch,Q_G; /*averaged output flows*/
 static double SWE_previous;
 static double GWE_previous;
 static double Smelt;   /*Snow melt [mm] during the time interval*/
 static double Ssubl;   /*Snow sublimation [mm] during the time interval*/
 static double Sevap;   /*Snow evaporation [mm] during the time interval*/
 static double Gmelt;   /*Glacier melt [mm] during the time interval*/
 static double Gsubl;   /*Glacier sublimation [mm] during the time interval*/
 static double Gevap;   /*Glacier evaporation [mm] during the time interval*/
 static double Smelt_previous;
 static double Ssubl_previous;
 static double Sevap_previous;
 static double Gmelt_previous;
 static double Gsubl_previous;
 static double Gevap_previous;
 static long isavings;

 /*internal variables to memorize input par:*/
 double DS_ch;            /*difference of water in the channel (S_ch0-S_ch1)*/
 double time_max;      /*time of all the simulation [s]*/
 double total_pixel;   /*total number of pixel which are not novalue*/
 double snowD, SWE, snowdensity, snowtemperature;
 double glacD=0.0, GWE=0.0, glacdensity=UV->V->co[2], glactemperature=UV->V->co[2];
 double D;
 double Q_GG;
 DOUBLEMATRIX *M;
 DOUBLETENSOR *Q;
 double JD,z;
 long d2,mo2,y2,h2,mi2;
 short sy=0; /* modified by Emanuele Cordano on 24 September 2009 */
 double fwet;
 double theta;
char *temp,*temp1; /* added by Emanuele Cordano on 21 September 2009 */

//double jd;
//long day,month,year,hour0,hour1,hour2,min0,min1,min2;

 /* Assignment to some internal variables of some input par:*/
time_max=times->TH*3600.0;
total_pixel=(double)par->total_pixel;

/* Print on the screen and into error-file the temporal coordinate of simulation: */
if (times->i_pixel==times->n_pixel){
	if(times->mm<=9){

		printf("%ld/%ld/%ld %ld:0%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));

		f=fopen(error_file_name,"a");
		fprintf(f,"%ld/%ld/%ld %ld:0%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));
		fclose(f);

	}else{

		printf("%ld/%ld/%ld %ld:%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));

		f=fopen(error_file_name,"a");
		fprintf(f,"%ld/%ld/%ld %ld:%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));
		fclose(f);

	}
}



//DISCHARGE
//****************************************************************************************************************
//****************************************************************************************************************
 if(times->time==0){
	Rout=0.0; Qsub_ch=0.0; Qsup_ch=0.0; R_G=0.0; Q_G=0.0;
 }

 Q_GG=0.0;
 for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
       if (land->LC->co[r][c]!=NoV){/*if the pixel is not a novalue*/
          Q_GG+=sl->J->co[Nl][r][c]*UV->U->co[1]*UV->U->co[2]*0.001; /*[mc/s]*/
       }
    }
 }


Qsub_ch+=cnet->Q_sub_s->co[1];
Qsup_ch+=cnet->Q_sup_s->co[1];
Q_G+=Q_GG;

//Calculation of the total outflow for unit of area [mm] for the current time step: Rout=Rsup_ch-Rsub_ch-Rsup_sea-Rsub_sea [mm]:
Rout+=(cnet->Q_sup_s->co[1]+cnet->Q_sub_s->co[1]+Q_GG)*par->Dt/(UV->U->co[1]*UV->U->co[2]*total_pixel)*1000.0;
R_G+=Q_GG*par->Dt/(UV->U->co[1]*UV->U->co[2]*total_pixel)*1000.0;

if (times->i_pixel==times->n_pixel){/*Print the outlet flows*/

	t_i=times->time-par->Dt*times->n_pixel;
	date_time(times->time+par->Dt, par->year0, par->JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
	temp=join_strings(files->co[fQ]+1,textfile);
	f=fopen(temp,"a");
	free(temp);
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
	fprintf(f,",%f,%f,%f,%f\n",(Qsub_ch+Qsup_ch+Q_G)/((double)times->n_pixel),Qsub_ch/((double)times->n_pixel),Qsup_ch/((double)times->n_pixel),Q_G/((double)times->n_pixel));
	fclose(f);

	Qsub_ch=0.0; Qsup_ch=0.0; Q_G=0.0;

	//melt fluxes
	if(par->output_balancesn!=0){

		for(i=1;i<=land->clax->nh;i++){

			write_suffix(SSSS, i, 0);

			temp=join_strings(files->co[fmeltlu]+1,"_snowcovered_");
			temp1=join_strings(temp,SSSS);
			name=join_strings(temp1,textfile);
			free(temp);
			free(temp1);
			f=fopen(name,"a");
			if(land->cont->co[i][1]>0){
				if(land->cont->co[i][2]>0){
					fprintf(f,"%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,
						(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
						wat->outfluxes->co[1][i]*times->n_pixel/(double)land->cont->co[i][2],
						wat->outfluxes->co[2][i]*times->n_pixel/(double)land->cont->co[i][2],
						wat->outfluxes->co[3][i]*times->n_pixel/(double)land->cont->co[i][2],
						land->cont->co[i][2]/times->n_pixel,
						land->cont->co[i][2]*UV->U->co[1]*UV->U->co[2]*1.0E-6/times->n_pixel,
						1.0E2*land->cont->co[i][2]/(land->cont->co[i][1]*times->n_pixel),
						1.0E2*land->cont->co[i][2]/(times->n_pixel*total_pixel));
				}else{
					fprintf(f,"%f,%f,%f,%f,%f,%f,%d,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
						0.0,0.0,0.0,0,0.0,0.0,0.0);
				}
			}else{
				fprintf(f,"%f,%f,%f,%f,%f,%f,%d,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
					0.0,0.0,0.0,0,0.0,0.0,0.0);
			}
			fclose(f);
			free(name);
			temp=join_strings(files->co[fmeltlu]+1,"_snowfree_");
			temp1=join_strings(temp,SSSS);

			name=join_strings(temp1,textfile);
			free(temp);
			free(temp1);

			f=fopen(name,"a");
			if(land->cont->co[i][1]>0){
				if(land->cont->co[i][2]<land->cont->co[i][1]*times->n_pixel){
					fprintf(f,"%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
						wat->outfluxes->co[4][i]*times->n_pixel/(double)(land->cont->co[i][1]*times->n_pixel-land->cont->co[i][2]),
						wat->outfluxes->co[5][i]*times->n_pixel/(double)(land->cont->co[i][1]*times->n_pixel-land->cont->co[i][2]),
						wat->outfluxes->co[6][i]*times->n_pixel/(double)(land->cont->co[i][1]*times->n_pixel-land->cont->co[i][2]),
						land->cont->co[i][1]-land->cont->co[i][2]/times->n_pixel,
						(land->cont->co[i][1]*times->n_pixel-land->cont->co[i][2])*UV->U->co[1]*UV->U->co[2]*1.0E-6/times->n_pixel,
						1.0E2*(land->cont->co[i][1]*times->n_pixel-land->cont->co[i][2])/(land->cont->co[i][1]*times->n_pixel),
						1.0E2*(land->cont->co[i][1]*times->n_pixel-land->cont->co[i][2])/(times->n_pixel*total_pixel));
				}else{
					fprintf(f,"%f,%f,%f,%f,%f,%f,%d,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
						0.0,0.0,0.0,0,0.0,0.0,0.0);
				}
			}else{
				fprintf(f,"%f,%f,%f,%f,%f,%f,%d,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
					0.0,0.0,0.0,0,0.0,0.0,0.0);
			}
			fclose(f);
			free(name);
			land->cont->co[i][2]=0;
			for(j=1;j<=6;j++){
				wat->outfluxes->co[j][i]=0.0;
			}
		}
	}
}


//DATA POINT
//****************************************************************************************************************
//****************************************************************************************************************

if(par->state_pixel==1){

	for(i=1;i<=par->chkpt->nrh;i++){

		write_suffix(SSSS, i, 0);
		r=par->rc->co[i][1];
		c=par->rc->co[i][2];

		sy=sl->type->co[r][c];

		// update of mean punctual values
		for(l=1;l<=Nl;l++){
			sl->Tmean->co[l][i]+=sl->T->co[l][r][c]/(double)times->n_pixel;
			//printf("time=%f,l=%ld,sl->T=%f,sl->Tmean=%f,times->i_pixel=%ld,time->n_pixel=%ld",times->time,l,sl->T->co[l][r][c],sl->Tmean->co[l],times->i_pixel,times->n_pixel);stop_execution();
			sl->thetai_mean->co[l][i]+=sl->thice->co[l][r][c]/(double)times->n_pixel;
			sl->thetaw_mean->co[l][i]+=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
					sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
					PSImin,par->Esoil)/(double)times->n_pixel;
			sl->psi_mean->co[l][i]+=sl->P->co[l][r][c]/(double)times->n_pixel;
		}
		/*Print of pixel-output every times->n_pixel time step */
		if (times->i_pixel==times->n_pixel){

			wat->out1->co[11][i]=wat->h_sup->co[r][c];
			for(l=1;l<=Nl;l++){
				wat->out1->co[12][i]+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
					sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
					PSImin,par->Esoil);
			}

			//wat->wcan_rain [mm] is the water (in both the liquid and solid form) that is currently on the leaves
			wat->out1->co[13][i]=-wat->out1->co[14][i];
			wat->out1->co[14][i]=wat->wcan_rain->co[r][c];
			wat->out1->co[13][i]+=wat->wcan_rain->co[r][c];

			//Further calculation
 			wat->out1->co[11][i]-=wat->out1->co[1][i]; /*wat store on sl-surface [mm]*/
			wat->out1->co[12][i]-=wat->out1->co[2][i]; /*wat store in subsoil [mm]*/
			wat->out1->co[15][i]*=3600.0*24.0/par->Dt/(double)times->n_pixel; /*from mm in the pixel-output-interval to mm/d*/

			//Snow
			if(snow->lnum->co[r][c]>0){
				snowD=0.0;
				SWE=0.0;
				snowtemperature=0.0;
				for(l=1;l<=snow->lnum->co[r][c];l++){
					snowD+=snow->Dzl->co[l][r][c];
					SWE+=1.0E+3*(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c])/rho_w;
					snowtemperature+=snow->T->co[l][r][c]*snow->Dzl->co[l][r][c];
				}
				snowdensity=SWE*rho_w/snowD;
				snowtemperature/=snowD;
			}else{
				snowD=0.0;
				SWE=0.0;
				snowdensity=0.0;
				snowtemperature=5.0;
			}

			//Glacier
			if(par->glaclayer_max>0){
				if(glac->lnum->co[r][c]>0){
					glacD=0.0;
					GWE=0.0;
					glactemperature=0.0;
					for(l=1;l<=glac->lnum->co[r][c];l++){
						glacD+=glac->Dzl->co[l][r][c];
						GWE+=1.0E+3*(glac->w_liq->co[l][r][c]+glac->w_ice->co[l][r][c])/rho_w;
						glactemperature+=glac->T->co[l][r][c]*glac->Dzl->co[l][r][c];
					}
					glacdensity=GWE*rho_w/glacD;
					glactemperature/=glacD;
				}else{
					glacD=0.0;
					GWE=0.0;
					glacdensity=0.0;
					glactemperature=5.0;
				}
			}


			if(egy->out1->co[56][i]>0){
				fwet=pow(wat->out1->co[14][i]/(0.1*egy->out1->co[56][i]),2./3.);
				if(fwet>1) fwet=1.0;
			}else{
				fwet=0.0;
			}

			/*update of the file "f" every times->n_pixel times: */
			t_i=times->time-par->Dt*times->n_pixel;
			temp=join_strings(files->co[fpoint]+1,SSSS);
			name=join_strings(temp,textfile);
			free(temp);
			f=fopen(name,"a");
			date_time(times->time+par->Dt, par->year0, par->JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
			write_date(f, d2, mo2, y2, h2, mi2);
			fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
			fprintf(f,",%f,%f,%f,",0.5*(t_i+par->Dt+times->time+par->Dt)/(double)86400,(t_i+par->Dt),(times->time+par->Dt));
			fprintf(f,"%14.3f,%14.9f,%f,%14.4f,%14.6f,%14.6f,%14.6f,%14.6f,%14.12f,%14.6f,%14.12f,",
					egy->out1->co[17][i]/*v*/,egy->out1->co[55][i]/*vdir*/,egy->out1->co[18][i]/*RH*/,egy->out1->co[19][i]/*P*/,egy->out1->co[20][i]/*Ta*/,egy->out1->co[64][i]/*Tsurface*/,
					egy->out1->co[46][i],egy->out1->co[47][i],egy->out1->co[48][i],egy->out1->co[49][i],egy->out1->co[50][i]);
			fprintf(f,"%14.5f,%14.5f,%14.5f,%14.5f,%14.12f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,",
					egy->out1->co[14][i]/*Rsw*/,egy->out1->co[53][i]/*SWbeam*/,egy->out1->co[54][i]/*SWdiff*/,egy->out1->co[45][i],
					(180/Pi)*egy->hsun,(180/Pi)*egy->dsun,(180/Pi)*acos(cos(top->slopes->co[r][c])*sin(egy->hsun)+sin(top->slopes->co[r][c])*cos(egy->hsun)*cos(-top->aspect->co[r][c]+egy->dsun)),
					egy->out1->co[12][i]/*Rlwin*/,egy->out1->co[13][i]/*Rlwout*/,egy->out1->co[11][i]/*Rnet*/,egy->out1->co[14][i]+egy->out1->co[45][i]/*SW*/,
					egy->out1->co[12][i]+egy->out1->co[13][i]/*LW*/,egy->out1->co[6][i]/*H*/,egy->out1->co[5][i]/*L*/,egy->out1->co[9][i]/*Qrain*/,egy->out1->co[8][i]/*G*/,egy->out1->co[7][i]/*surfEB*/);
			fprintf(f,"%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,",
					egy->out1->co[24][i]/*Rswin_cum*/,egy->out1->co[28][i]/*Rswout_cum*/,egy->out1->co[22][i]/*Rlwin_cum*/,egy->out1->co[23][i]/*Rlwout_cum*/,
					egy->out1->co[24][i]+egy->out1->co[28][i],egy->out1->co[22][i]+egy->out1->co[23][i],egy->out1->co[21][i]/*Rnet_cum*/,
					egy->out1->co[26][i]/*H_cum*/,egy->out1->co[25][i]/*ET_cum*/,egy->out1->co[27][i]/*G_cum*/);

			/*Eg[mm],Sg[mm],Etc[mm],Psnow[mm],Prain[mm],Psnow_c[mm],Prain_SOILc[mm],Prain_SNOWc[mm],Ptot_c[mm],Wtrain[mm],Wtsnow - 11
			Ptot_atm,Rain_atm,Snow_atm,Ptot_atm_cum,Prain_atm_cum,Psnow_atm_cum,Pn[mm],Runoff[mm] - 8
			q_sup[mm],q_sub[mm],DS_sup[mm],DS_sub[mm],q_G[mm] - 5*/
			fprintf(f,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
					 egy->out1->co[1][i]/*Eg*/,egy->out1->co[2][i]/*Sg*/,egy->out1->co[4][i]/*Etc*/,
					 wat->out1->co[3][i]/*Psnow*/,wat->out1->co[4][i]/*Prain*/,wat->out1->co[21][i]/*Psnow_cum*/,wat->out1->co[22][i]/*Prain_cum_soil*/,
					 wat->out1->co[27][i]/*rain_on_snow_cum*/,wat->out1->co[21][i]+wat->out1->co[22][i]/*Ptot_cum*/,
					 wat->wcan_rain->co[r][c]/*wcan_rain*/,wat->wcan_snow->co[r][c]/*wcan_snow*/,
					 wat->out1->co[16][i]+wat->out1->co[17][i]/*Ptot atm*/,wat->out1->co[17][i]/*rain*/,wat->out1->co[16][i]/*snow*/,
					 wat->out1->co[23][i]+wat->out1->co[24][i],wat->out1->co[24][i],wat->out1->co[23][i],
					 wat->out1->co[5][i]/*net liquid precipitation*/,wat->out1->co[6][i]/*runoff*/,
					 wat->out1->co[7][i]/*q_sup*/,wat->out1->co[10][i]/*q_sub*/,wat->out1->co[11][i]/*wat stored on sl*/,wat->out1->co[12][i]/*wat stored in subsoil*/,
					 wat->out1->co[8][i]/*Q_g*/);

			fprintf(f,"%14.8f,%14.8f,%14.8f,%14.8f,%f,%f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,",
					 snowD,SWE,snowdensity,snowtemperature,snow->out_bs->co[9][i],snow->out_bs->co[10][i],snow->melted->co[i],snow->subl->co[i],snow->evap->co[i],
					 glacD,GWE,glacdensity,glactemperature,glac->melted->co[i],glac->subl->co[i],glac->evap->co[i]);
			fprintf(f,"%14.9f,",wat->out1->co[15][i]/*q_v*/);

			fprintf(f,"%14.12f,%14.3f,%14.6f,", egy->out1->co[16][i],egy->out1->co[29][i],egy->out1->co[30][i]);
			fprintf(f,"%f,%f,%f,%f,%f,%20.18f,",egy->out1->co[51][i],egy->out1->co[52][i],egy->out1->co[56][i],egy->out1->co[57][i],egy->out1->co[58][i],wat->out1->co[28][i]/(double)times->n_pixel);
			fprintf(f,"%f,%f,%f,%f,%f,%f,%f,",egy->out1->co[64][i],egy->out1->co[10][i],egy->out1->co[63][i],egy->out1->co[59][i],egy->out1->co[60][i],egy->out1->co[61][i],egy->out1->co[62][i]);
			fprintf(f,"%f,%f,%f,%f,%f,%f,%f,",egy->out1->co[6][i]+egy->out1->co[61][i],egy->out1->co[5][i]+egy->out1->co[62][i],egy->out1->co[65][i],egy->out1->co[66][i],egy->out1->co[67][i],egy->out1->co[68][i],egy->out1->co[69][i]);
			fprintf(f,"%e,%e,%e,%e,%e,%e,%e,%f,",egy->out1->co[70][i],egy->out1->co[71][i],egy->out1->co[72][i],egy->out1->co[73][i],egy->out1->co[74][i],egy->out1->co[77][i],egy->out1->co[29][i],egy->out1->co[30][i]);
			fprintf(f,"%f,%f,%14.8f,%14.8f,%14.8f,%14.8f,%f,",egy->out1->co[75][i],egy->out1->co[76][i],egy->out1->co[78][i],egy->out1->co[79][i],egy->out1->co[80][i],egy->out1->co[81][i],egy->out1->co[82][i]);
			fprintf(f,"%f,%14.8e\n",egy->out1->co[83][i],egy->out1->co[84][i]);
			fclose(f);

			/*update of the file "f" every times->n_pixel times: */
			z=0.0;
			free(name);
			temp=join_strings(files->co[fsnz]+1,SSSS);
			name=join_strings(temp,textfile);
			free(temp);
			f=fopen(name,"a");
			write_date(f, d2, mo2, y2, h2, mi2);
			fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
			fprintf(f,",%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f,%f,%f,%f",(t_i+par->Dt),(times->time+par->Dt),snow->type->co[r][c],snowD,SWE,snowtemperature,snowdensity,
				snow->melted->co[i],snow->subl->co[i],snow->evap->co[i],snow->lnum->co[r][c],snow->out_bs->co[1][i],snow->out_bs->co[2][i],
				snow->out_bs->co[3][i],snow->out_bs->co[4][i],snow->out_bs->co[9][i],snow->out_bs->co[10][i]);

			for(l=1;l<=par->snowlayer_max;l++){
				fprintf(f,",%f",snow->Dzl->co[l][r][c]);
			}
			for(l=1;l<=par->snowlayer_max;l++){
				z+=snow->Dzl->co[l][r][c];
				fprintf(f,",%f",z);
			}
			for(l=1;l<=par->snowlayer_max;l++){
				fprintf(f,",%f",snow->T->co[l][r][c]);
			}
			for(l=1;l<=par->snowlayer_max;l++){
				if(snow->Dzl->co[l][r][c]>0){
					fprintf(f,",%f",( snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c] )/ ( snow->Dzl->co[l][r][c]*0.001 ) );
				}else{
					fprintf(f,",%f",0.0);
				}
			}
			for(l=1;l<=par->snowlayer_max;l++){
				if(snow->Dzl->co[l][r][c]>0){
					fprintf(f,",%f",1000.0*(snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c])/ rho_w );
				}else{
					fprintf(f,",%f",0.0);
				}
			}
			for(l=1;l<=par->snowlayer_max;l++){
				if(snow->Dzl->co[l][r][c]>0){
					fprintf(f,",%f",snow->w_ice->co[l][r][c]/(rho_i*0.001*snow->Dzl->co[l][r][c]));
				}else{
					fprintf(f,",%f",0.0);
				}
			}
			for(l=1;l<=par->snowlayer_max;l++){
				if(snow->Dzl->co[l][r][c]>0){
					fprintf(f,",%f",snow->w_liq->co[l][r][c]/(rho_w*0.001*snow->Dzl->co[l][r][c]));
				}else{
					fprintf(f,",%f",0.0);
				}
			}
			for(l=1;l<=par->snowlayer_max;l++){
				fprintf(f,",%f",snow->CR1m->co[l]);
			}
			for(l=1;l<=par->snowlayer_max;l++){
				fprintf(f,",%f",snow->CR2m->co[l]);
			}
			for(l=1;l<=par->snowlayer_max;l++){
				fprintf(f,",%f",snow->CR3m->co[l]);
			}
			fprintf(f,"\n");
			fclose(f);

			/*update of the file "f" every times->n_pixel times: */
			if(par->glaclayer_max>0){
				temp=join_strings(files->co[fglz]+1,SSSS);
				free(name);
				name=join_strings(temp,textfile);
				free(temp);
				f=fopen(name,"a");
				free(name); //added by Emanuele
				write_date(f, d2, mo2, y2, h2, mi2);
				fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
				fprintf(f,",%f,%f,%f,%f,%f,%ld,%f,%f,%f,%f,%f",(t_i+par->Dt),(times->time+par->Dt),glacD,GWE,glactemperature,glac->lnum->co[r][c],glacdensity,
					egy->out1->co[15][i]/*albedo*/,glac->melted->co[i],glac->subl->co[i],glac->evap->co[i]);
				for(l=1;l<=par->glaclayer_max;l++){
					fprintf(f,",%f",glac->T->co[l][r][c]);
				}
				for(l=1;l<=par->glaclayer_max;l++){
					fprintf(f,",%f",glac->Dzl->co[l][r][c]);
				}
				for(l=1;l<=par->glaclayer_max;l++){
					if(glac->Dzl->co[l][r][c]>0){
						fprintf(f,",%f",( glac->w_ice->co[l][r][c]+glac->w_liq->co[l][r][c] )/ ( glac->Dzl->co[l][r][c]*0.001 ) );
					}else{
						fprintf(f,",%f",UV->V->co[2]);
					}
				}
				for(l=1;l<=par->glaclayer_max;l++){
					fprintf(f,",%f",1000.0*( glac->w_ice->co[l][r][c]+glac->w_liq->co[l][r][c] )/ rho_w );
				}
				for(l=1;l<=par->glaclayer_max;l++){
					fprintf(f,",%f",glac->w_ice->co[l][r][c]);
				}
				for(l=1;l<=par->glaclayer_max;l++){
					fprintf(f,",%f",glac->w_liq->co[l][r][c]);
				}
				fprintf(f,"\n");
				fclose(f);
			}

			//sl output
			write_soil_output(times->n_pixel, i, times->time, par->Dt, par->year0, par->JD0, par->rc, sl, PSImin, par->Esoil);

		}//end(times->i_pixel==times->n_pixel)
	}//end for(i=1;i<=par->chkpt->nrh;i++)

	if(times->i_pixel==times->n_pixel){
		for(i=1;i<=par->chkpt->nrh;i++){
			for(j=1;j<=20;j++) { egy->out1->co[j][i]=0.0; }
			for(j=29;j<=30;j++) { egy->out1->co[j][i]=0.0; }
			for(j=35;j<=84;j++) { egy->out1->co[j][i]=0.0; }
			for(j=1;j<=12;j++) { wat->out1->co[j][i]=0.0; }
			for(j=15;j<=20;j++) { wat->out1->co[j][i]=0.0; }
			for(j=28;j<=28;j++) { wat->out1->co[j][i]=0.0; }
			for(j=1;j<=5;j++) { snow->out_bs->co[2*j-1][i]=0.0; }

			for(l=1;l<=Nl;l++){
				wat->out1->co[2][i]+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
								sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
								PSImin,par->Esoil);
				sl->Tmean->co[l][i]=0.0;
				sl->thetaw_mean->co[l][i]=0.0;
				sl->thetai_mean->co[l][i]=0.0;
				sl->psi_mean->co[l][i]=0.0;
				sl->Tmin->co[l][i]=99.0;
				sl->thetaw_min->co[l][i]=99.0;
				sl->thetai_min->co[l][i]=99.0;
				sl->Tmax->co[l][i]=-99.0;
				sl->thetaw_max->co[l][i]=-99.0;
				sl->thetai_max->co[l][i]=-99.0;
			}
			wat->out1->co[1][i]=wat->h_sup->co[r][c];

		}
	}
}// end if(par->state_pixel==1)


//BASIN DATA
//****************************************************************************************************************
//****************************************************************************************************************
if(times->time==0){
	/*NOTE: these are all static variables*/
	wt0_basin=0.0;		Ssup=0.0;				Ssub=0.0;
	SWE_previous=0.0;   Smelt_previous=0.0;		Ssubl_previous=0.0;		Sevap_previous=0.0;   Smelt=0.0;	Ssubl=0.0;	Sevap=0.0;
	GWE_previous=0.0;   Gmelt_previous=0.0;		Gsubl_previous=0.0;		Gevap_previous=0.0;   Gmelt=0.0;	Gsubl=0.0;	Gevap=0.0;

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if (land->LC->co[r][c]!=NoV){ /*if the pixel is not a novalue*/

				sy=sl->type->co[r][c];

				/*find the total water on the soil-surface (volume for unit of area):*/
				Ssup+=wat->h_sup->co[r][c];

				/*find the total water in the subsoil (volume for unit of area):*/
				for(l=1;l<=Nl;l++){
					Ssub+=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
						sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],PSImin,par->Esoil)*sl->pa->co[sy][jdz][l];
				}

				/*water on the leaves*/
				wt0_basin+=wat->wcan_rain->co[r][c];

				/*find the mean water equivalent snow and glacier height at the initial condition:*/
				for(l=1;l<=snow->lnum->co[r][c];l++){
					SWE_previous+=(1.0E+3*(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c])/rho_w)/total_pixel;
				}

				if(par->glaclayer_max>0){
					for(l=1;l<=glac->lnum->co[r][c];l++){
						GWE_previous+=(1.0E+3*(glac->w_liq->co[l][r][c]+glac->w_ice->co[l][r][c])/rho_w)/total_pixel;
					}
				}

				 S_ch0=0.0; S_ch1=0.0;

			}
		}
	}
}

/* Updating snow and glacier melting : */
Smelt+=snow->melted_basin/total_pixel;	//[mm]
Ssubl+=snow->subl_basin/total_pixel;	//[mm]
Sevap+=snow->evap_basin/total_pixel;	//[mm]
Gmelt+=glac->melted_basin/total_pixel;	//[mm]
Gsubl+=glac->subl_basin/total_pixel;	//[mm]
Gevap+=glac->evap_basin/total_pixel;	//[mm]

if (times->i_basin==times->n_basin){/*Print of the output of all the basin*/

	wat->out2->co[5]=-wt0_basin;
	wat->out2->co[6]=-Ssup;
	wat->out2->co[7]=-Ssub;
	wt0_basin=0.0;	Ssup=0.0;	Ssub=0.0;	SWE=0.0;	GWE=0.0;

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if (land->LC->co[r][c]!=NoV){ /*if the pixel is not a novalue(=9)*/

				/*Snow water equivalent [mm]*/
				for(l=1;l<=snow->lnum->co[r][c];l++){
					SWE+=(1.0E+3*(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c])/rho_w)/total_pixel;
				}

				/*Glacier water equivalent [mm]*/
				if(par->glaclayer_max>0){
					for(l=1;l<=glac->lnum->co[r][c];l++){
						GWE+=(1.0E+3*(glac->w_liq->co[l][r][c]+glac->w_ice->co[l][r][c])/rho_w)/total_pixel;
					}
				}

				/*Water on the leaves*/
				wat->out2->co[5]+=wat->wcan_rain->co[r][c];
				wt0_basin+=wat->wcan_rain->co[r][c];

				/*Total water storage on the sl-surface (volume for unit of area[mm]):*/
				wat->out2->co[6]+=wat->h_sup->co[r][c];
				Ssup+=wat->h_sup->co[r][c];

				/*Total water storage in the subsoil (volume for unit of area[mm]):*/
				for(l=1;l<=Nl;l++){
					wat->out2->co[7]+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
								sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
								PSImin,par->Esoil);
					Ssub+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
								sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
								PSImin,par->Esoil);
				}
			}
		}
	}

	/*computation of the total water present in the channels:*/
	for(s=1;s<=cnet->Q_sup_s->nh;s++){
		S_ch1+=cnet->Q_sup_s->co[s];
		S_ch1+=cnet->Q_sub_s->co[s];
	}
	DS_ch=(S_ch1-S_ch0)*par->Dt*1000.0/(UV->U->co[1]*UV->U->co[2]*total_pixel);/*[mm]*/
	S_ch0=S_ch1;
	S_ch1=0.0;

	/*computation of the final means to write in the file "O_5BASINname": for global output each dt pixel times->n_basin*/
	for(j=1;j<=7;j++){
		wat->out2->co[j]/=total_pixel;
	}

	for(j=1;j<=11;j++){
		egy->out2->co[j]/=total_pixel;
	}

	/*update of the file "O_5BASINname" every times->n_basin->co[1] times:*/
	temp=join_strings(files->co[fbas]+1,textfile);
	f=fopen(temp,"a");
	free(temp);
	t_i=times->time-par->Dt*times->n_basin;
	date_time(times->time+par->Dt, par->year0, par->JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
	fprintf(f,",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		(t_i+par->Dt),(times->time+par->Dt),egy->out2->co[3],egy->out2->co[4],egy->out2->co[5],egy->out2->co[8],egy->out2->co[9],egy->out2->co[7],egy->out2->co[6],//9
		egy->out2->co[12],egy->out2->co[13],egy->out2->co[11],egy->out2->co[10],egy->out2->co[1],egy->out2->co[2],egy->out2->co[15],//16
		egy->out2->co[14],wat->out2->co[1]/*Prain*/,wat->out2->co[2]/*Psnow*/,wat->out2->co[5]/*Dwt*/,wat->out2->co[3]/*Pn*/, wat->out2->co[4]/*runoff*/,
		wat->out2->co[6]/*Dssup*/,wat->out2->co[7]/*Dssub*/,DS_ch,R_G,Rout,Ssup/total_pixel,Ssub/total_pixel,wt0_basin/total_pixel,//30
		SWE,SWE-SWE_previous,Smelt,Smelt-Smelt_previous,Ssubl,Ssubl-Ssubl_previous,Sevap,Sevap-Sevap_previous,GWE,GWE-GWE_previous,Gmelt,Gmelt-Gmelt_previous,Gsubl,//43
		Gsubl-Gsubl_previous,Gevap,Gevap-Gevap_previous,wat->out2->co[8]/(double)times->n_basin);//47
	fclose(f);

	/*writing of some results on the screen:*/
	printf(" SW=%f W/m2  LW=%f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Prain=%6.2f mm  Psnow=%6.2f mm  Rout=%6.2f mm/h \n Max Error Richards=%20.18f mm/h\n\n",
		egy->out2->co[8],egy->out2->co[9],egy->out2->co[7],egy->out2->co[6],wat->out2->co[1],wat->out2->co[2],Rout*3600.0/(double)(par->Dt*times->n_basin),wat->out2->co[8]/(double)times->n_basin);

	Rout=0.0; R_G=0.0;
	SWE_previous=SWE;   Smelt_previous=Smelt;	Ssubl_previous=Ssubl;	Sevap_previous=Sevap;
	GWE_previous=GWE;   Gmelt_previous=Gmelt;	Gsubl_previous=Gsubl;	Gevap_previous=Gevap;

	for(j=1;j<=15;j++){
		egy->out2->co[j]=0.0;
	}

	for(j=1;j<=8;j++){
		wat->out2->co[j]=0.0;
	}
}

//ALTIMETRIC RANKS
//****************************************************************************************************************
//****************************************************************************************************************
for(i=1;i<=par->ES_num;i++){

	write_suffix(SSSS, i, 0);
	temp=join_strings(files->co[farank]+1,SSSS);

	name=join_strings(temp,textfile);
	free(temp);
	if (times->i_basin==times->n_basin){/*Print of the output of all the basin*/

		for(l=1;l<=snow->lnum->co[r][c];l++){
			egy->out3->co[19][i]+=(1.0E+3*(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c])/rho_w);
		}

		if(par->glaclayer_max>0){
			for(l=1;l<=glac->lnum->co[r][c];l++){
				egy->out3->co[23][i]+=(1.0E+3*(glac->w_liq->co[l][r][c]+glac->w_ice->co[l][r][c])/rho_w);
			}
		}

		for(j=1;j<=23;j++){
			egy->out3->co[j][i]/=(double)top->ES_pixel->co[i];
		}

		f=fopen(name,"a");
		free(name);
		fprintf(f,"%14.2f %f %f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f\n",
							(0.5*(t_i+par->Dt)+0.5*(times->time+par->Dt))/86400.0,(t_i+par->Dt),(times->time+par->Dt),egy->out3->co[1][i],egy->out3->co[2][i],egy->out3->co[3][i],
							egy->out3->co[4][i],egy->out3->co[5][i],egy->out3->co[6][i],egy->out3->co[11][i],egy->out3->co[7][i],
							egy->out3->co[8][i],egy->out3->co[9][i],egy->out3->co[10][i],egy->out3->co[12][i],egy->out3->co[13][i],
							egy->out3->co[14][i],egy->out3->co[15][i],egy->out3->co[19][i],egy->out3->co[16][i],egy->out3->co[18][i],
							egy->out3->co[17][i],egy->out3->co[23][i],egy->out3->co[20][i],egy->out3->co[22][i],egy->out3->co[21][i]);
		fclose(f);

		for(j=1;j<=23;j++){
			egy->out3->co[j][i]=0.0;
		}

	}
}


//DISTRIBUTED OUTPUTS
//****************************************************************************************************************
//****************************************************************************************************************
if(par->output_h_sup>0){
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV) wat->hsupav->co[r][c]+=wat->h_sup->co[r][c]/((par->output_h_sup*3600.0)/(par->Dt));
		}
	}
}

//TETA
if(par->output_TETAxy>0 && fmod(times->time+par->Dt,par->output_TETAxy*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_TETAxy*3600.0));
	Q=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(Q,UV->V->co[2]);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				sy=sl->type->co[r][c];
				for(l=1;l<=Nl;l++){
					Q->co[l][r][c]=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],
									sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],PSImin, par->Esoil);
				}
			}
		}
	}
	write_tensorseries2(n_file, files->co[fliq]+1, 0, par->format_out, Q, UV,times->time+par->Dt); //USE_NETCDF_MAP
	//modified by Emanuele Cordano
	free_doubletensor(Q);

#ifdef USE_NETCDF_MAP
	temp=join_strings(files->co[fliq]+1,"$_error");
	write_map(temp, 0, par->format_out, wat->error, UV,times->time+par->Dt,n_file);
#else
	write_suffix(SSSS, n_file, 0);
	temp=join_strings(files->co[fliq]+1,"error");
	name=join_strings(temp,SSSS);
	write_map(name, 0, par->format_out, wat->error, UV,0,0);
	initmatrix(0.0, wat->error, land->LC, NoV);
	if (name!=NULL) free(name);
#endif
	free(temp);
}



//T  EMANUELE CORDANO 27 nov 2009
if(par->output_Txy>0 && fmod(times->time+par->Dt,par->output_Txy*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Txy*3600.0));
	write_tensorseries2(n_file, files->co[fT]+1, 0, par->format_out, sl->T, UV,times->time+par->Dt); //USE_NETCDF_MAP

	//***************************************************************************
	//calculate active layer depth
	/*M=new_doublematrix(Nr,Nc);
	initialize_doublematrix(M,UV->V->co[2]);
	V=new_doublevector(Nl);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				sy=sl->type->co[r][c];
				for(l=1;l<=Nl;l++){
					V->co[l]=sl->T->co[l][r][c];
				}
				M->co[r][c]=find_activelayer(V, sl->pa->co[sy][jdz]);
			}
		}
	}
	free_doublevector(V);
	write_suffix(SSSS, n_file, 0);
	temp1=join_strings(WORKING_DIRECTORY,"tables/01thawed_soil_depth");
	temp2=join_strings(temp1,SSSS);
	write_map(temp2, 0, par->format_out, M, UV);
	free(temp1);
	free(temp2);*/
	//end of active layer depth calculation
	//***************************************************************************

}

//TETAICE
if(par->output_TETAICExy>0 && fmod(times->time+par->Dt,par->output_TETAICExy*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_TETAICExy*3600.0));
	write_tensorseries2(n_file, files->co[fice]+1, 0, par->format_out, sl->thice, UV,times->time+par->Dt); //USE_NETCDF_MAP
}

//PSI
if(par->output_PSIxy>0 && fmod(times->time+par->Dt,par->output_PSIxy*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_PSIxy*3600.0));
#ifdef USE_NETCDF_MAP
	//LIQ file correspond exactly to netcdf variable
	temp=join_strings(files->co[fpsi]+1,"");
	write_tensorseries2(n_file,temp, 0, par->format_out, sl->P, UV,times->time+par->Dt);
	free(temp);
#else
	temp=join_strings(files->co[fpsi]+1,"LIQ");
	write_tensorseries2(n_file,temp, 0, par->format_out, sl->P, UV,0);
	free(temp);
#endif
	Q=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(Q,UV->V->co[2]);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				for(l=1;l<=Nl;l++){
					theta=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],
						sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],PSImin, par->Esoil);
					Q->co[l][r][c]=psi_teta(theta + sl->thice->co[l][r][c], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
						sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PSImin, par->Esoil);
				}
			}
		}
	}
#ifdef USE_NETCDF_MAP
	//TOT file correspond netcdfvarname_tot variable
	temp=join_strings(files->co[fpsi]+1,"$_tot"); //$ is token field separator
	printf("%s\n",temp);
	write_tensorseries2(n_file,temp, 0, par->format_out, Q, UV,times->time+par->Dt);
	free(temp);
#else
	temp=join_strings(files->co[fpsi]+1,"TOT");
	write_tensorseries2(n_file,temp, 0, par->format_out, Q, UV,times->time+par->Dt);
	free(temp);
#endif
	//***************************************************************************
	//calculate saturation front depth
	/*K=new_doublematrix(Nr,Nc);
		initialize_doublematrix(K,UV->V->co[2]);
		V=new_doublevector(Nl);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					sy=sl->type->co[r][c];
					for(l=1;l<=Nl;l++){
						V->co[l]=Q->co[l][r][c];
					}
					K->co[r][c]=find_watertable(V, sl->pa->co[sy][jdz]);
				}
			}
		}
		free_doublevector(V);
		write_suffix(SSSS, n_file, 0);
		temp1=join_strings(WORKING_DIRECTORY,"tables/02sat_front_depth");
		temp2=join_strings(temp1,SSSS);
		write_map(temp2, 0, par->format_out, K, UV);
		free(temp1);
		free(temp2);*/
	//end of saturation front depth calculation
	//***************************************************************************
	free_doubletensor(Q);
}


//***************************************************************************
//calculate frost table depth
/*if( par->output_Txy>0 && fmod(times->time+par->Dt,par->output_Txy*3600.0)==0 && par->output_PSIxy>0 && fmod(times->time+par->Dt,par->output_PSIxy*3600.0)==0 ){
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				//if thawed soil depth < saturation front depth
				if(M->co[r][c] < K->co[r][c]) M->co[r][c]=K->co[r][c];
			}
		}
	}
	n_file=(long)((times->time+par->Dt)/(par->output_Txy*3600.0));
	write_suffix(SSSS, n_file, 0);
	temp1=join_strings(WORKING_DIRECTORY,"tables/03frost_table_depth");
	temp2=join_strings(temp1,SSSS);
	write_map(temp2, 0, par->format_out, M, UV);
	free(temp1);
	free(temp2);
}
if( par->output_Txy>0 && fmod(times->time+par->Dt,par->output_Txy*3600.0)==0 ) free_doublematrix(M);
if( par->output_PSIxy>0 && fmod(times->time+par->Dt,par->output_PSIxy*3600.0)==0 ) free_doublematrix(K);*/
//end of frost table depth calculation
//***************************************************************************

//matrix to handle further data
M=new_doublematrix(Nr,Nc);
initialize_doublematrix(M,UV->V->co[2]);

//ALBEDO
if(par->output_albedo>0 && fmod(times->time+par->Dt,par->output_albedo*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_albedo*3600.0));
#ifdef USE_NETCDF_MAP
	//write_map(temp,0, par->format_out, land->albedo, UV,times->time+par->Dt,n_file);
	write_map(files->co[falb]+1,0, par->format_out, land->albedo, UV,times->time+par->Dt,n_file); //091209
#else
	write_suffix(SSSS, n_file, 0);
	temp=join_strings(files->co[falb]+1,SSSS);
	write_map(temp,0, par->format_out, land->albedo, UV,0,0);
	free(temp);
#endif
	initmatrix(0.0, land->albedo, top->Z0, NoV);

}

//SNOW DEPTH
if(par->output_snow>0 && fmod(times->time+par->Dt,par->output_snow*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_snow*3600.0));
	write_suffix(SSSS, n_file, 0);
	for(r=1;r<=snow->Dzl->nrh;r++){
		for(c=1;c<=snow->Dzl->nch;c++){
			if(land->LC->co[r][c]!=NoV){
				M->co[r][c]=0.0;
				for(l=1;l<=snow->lnum->co[r][c];l++){
					M->co[r][c]+=snow->Dzl->co[l][r][c];
				}
			}
		}
	}
#ifdef USE_NETCDF_MAP
	temp=join_strings(files->co[fsn]+1,"$_dist"); //$ is token field separator
	write_map(temp, 0, par->format_out, M, UV,times->time+par->Dt,n_file);
#else
	temp=join_strings(files->co[fsn]+1,"Dist");
	temp1=join_strings(temp,SSSS);
	write_map(temp1, 0, par->format_out, M, UV,0,0);
	free(temp1);
#endif
	free(temp);

	//write_map(join_strings(join_strings(files->co[fsn]+1,"Dmax"),SSSS), 0, par->format_out, snow->max, UV);
	//initmatrix(0.0, snow->max, top->Z0, NoV);
	//write_map(join_strings(join_strings(files->co[fsn]+1,"Daverage"),SSSS), 0, par->format_out, snow->average, UV);
	//initmatrix(0.0, snow->average, top->Z0, NoV);
	if(par->blowing_snow==1){
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsubl"),SSSS), 0, par->format_out, snow->Wsubl_cum, UV);
		//initmatrix(0.0, snow->Wsubl_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsalt"),SSSS), 0, par->format_out, snow->Wtrans_cum, UV);
		//initmatrix(0.0, snow->Wtrans_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsusp"),SSSS), 0, par->format_out, snow->Wsusp_cum, UV);
		//initmatrix(0.0, snow->Wsusp_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsbgr"),SSSS), 0, par->format_out, snow->Wsubgrid_cum, UV);
		//initmatrix(0.0, snow->Wsubgrid_cum, top->Z0, NoV);
	#ifdef USE_NETCDF_MAP
		temp=join_strings(files->co[fsn]+1,"$_bstot"); //$ is token field separator
		write_map(temp, 0, par->format_out, snow->Wtot, UV,times->time+par->Dt,n_file);
	#else
		temp=join_strings(files->co[fsn]+1,"BStot");
		temp1=join_strings(temp,SSSS);
		write_map(temp1, 0, par->format_out, snow->Wtot, UV,0,0);
		free(temp1);
	#endif
		free(temp);
		initmatrix(0.0, snow->Wtot, top->Z0, NoV);

	}
}

//GLACIER DEPTH
if(par->glaclayer_max>0 && par->output_glac>0 && fmod(times->time+par->Dt,par->output_glac*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_glac*3600.0));
	write_suffix(SSSS, n_file, 0);
	for(r=1;r<=glac->Dzl->nrh;r++){
		for(c=1;c<=glac->Dzl->nch;c++){
			if(land->LC->co[r][c]!=NoV){
				M->co[r][c]=0.0;
				for(l=1;l<=glac->lnum->co[r][c];l++){
					M->co[r][c]+=(glac->w_liq->co[l][r][c]+glac->w_ice->co[l][r][c]);
				}
			}
		}
	}
#ifdef USE_NETCDF_MAP
	write_map(files->co[fgl]+1, 0, par->format_out, M, UV,times->time+par->Dt,n_file);
#else
	temp=join_strings(files->co[fgl]+1,SSSS);
	write_map(temp, 0, par->format_out, M, UV,0,0);
	free(temp);
#endif

}

//WATER OVER THE SURFACE
if(par->output_h_sup>0 && fmod(times->time+par->Dt,par->output_h_sup*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_h_sup*3600.0));
#ifdef USE_NETCDF_MAP
	write_map(files->co[fhsup]+1, 0, par->format_out, wat->hsupav, UV,times->time+par->Dt,n_file);
#else
	write_suffix(SSSS, n_file, 0);
	temp=join_strings(files->co[fhsup]+1,SSSS);
	write_map(temp, 0, par->format_out, wat->hsupav, UV,0,0);
	free(temp);
#endif
	initmatrix(0.0, wat->hsupav, top->Z0, NoV);
}

//RADIATION
if(par->output_Rn>0 && fmod(times->time+par->Dt,par->output_Rn*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Rn*3600.0));
#ifdef USE_NETCDF_MAP
	name=join_strings(files->co[fRn]+1,"$_nmean");
	write_map(name, 0, par->format_out, egy->Rn_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->Rn_mean, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fRn]+1,"$_lwin");
	write_map(name, 0, par->format_out, egy->LW_in, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->LW_in, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fRn]+1,"$lwout");
	write_map(name, 0, par->format_out, egy->LW_out, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->LW_out, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fRn]+1,"$_sw");
	write_map(name, 0, par->format_out, egy->SW, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->SW, top->Z0, NoV);
	free(name);

	if(par->distr_stat==1){
		name=join_strings(files->co[fRn]+1,"$_nmax");
		write_map(name, 0, par->format_out, egy->Rn_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->Rn_max, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fRn]+1,"$_nmin");
		write_map(name, 0, par->format_out, egy->Rn_min, UV,times->time+par->Dt,n_file);
		initmatrix(1.0E+9, egy->Rn_min, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fRn]+1,"$_lwmax");
		write_map(name, 0, par->format_out, egy->LW_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->LW_max, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fRn]+1,"$_lwmin");
		write_map(name, 0, par->format_out, egy->LW_min, UV,times->time+par->Dt,n_file);
		initmatrix(1.0E+9, egy->LW_min, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fRn]+1,"$_swmax");
		write_map(name, 0, par->format_out, egy->SW_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->SW_max, top->Z0, NoV);
		free(name);
	}
#else
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fRn]+1,"nmean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->Rn_mean, UV,0,0);
	free(temp);
	initmatrix(0.0, egy->Rn_mean, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fRn]+1,"LWin");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->LW_in, UV,0,0);
	free(temp);
	initmatrix(0.0, egy->LW_in, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fRn]+1,"LWout");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->LW_out, UV,0,0);
	free(temp);
	initmatrix(0.0, egy->LW_out, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fRn]+1,"SW");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->SW, UV,0,0);
	free(temp);
	initmatrix(0.0, egy->SW, top->Z0, NoV);
	free(name);

	if(par->distr_stat==1){
		name=join_strings(files->co[fRn]+1,"nmax");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->Rn_max, UV,0,0);
		free(temp);
		initmatrix(-1.0E+9, egy->Rn_max, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fRn]+1,"nmin");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->Rn_min, UV,0,0);
		free(temp);
		initmatrix(1.0E+9, egy->Rn_min, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fRn]+1,"LWmax");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->LW_max, UV,0,0);
		free(temp);
		initmatrix(-1.0E+9, egy->LW_max, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fRn]+1,"LWmin");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->LW_min, UV,0,0);
		free(temp);
		initmatrix(1.0E+9, egy->LW_min, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fRn]+1,"SWmax");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->SW_max, UV,0,0);
		free(temp);
		initmatrix(-1.0E+9, egy->SW_max, top->Z0, NoV);
		free(name);
	}
#endif
}

//GROUND HEAT FLUX
if(par->output_G>0 && fmod(times->time+par->Dt,par->output_G*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_G*3600.0));
#ifdef USE_NETCDF_MAP
	name=join_strings(files->co[fG]+1,"$_mean");
	write_map(name, 0, par->format_out, egy->G_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->G_mean, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fG]+1,"$_snowsoil");
	write_map(name, 0, par->format_out, egy->G_snowsoil, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->G_snowsoil, top->Z0, NoV);
	free(name);

	if(par->distr_stat==1){

		name=join_strings(files->co[fG]+1,"$_max");
		write_map(name, 0, par->format_out, egy->G_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->G_max, top->Z0, NoV);
		free(name);

		name=join_strings(files->co[fG]+1,"$_min");
		write_map(name, 0, par->format_out, egy->G_min, UV,times->time+par->Dt,n_file);
		initmatrix(1.0E+9, egy->G_min, top->Z0, NoV);
		free(name);

	}
#else
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fG]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->G_mean, UV,0,0);
	initmatrix(0.0, egy->G_mean, top->Z0, NoV);
	free(temp);
	free(name);

	name=join_strings(files->co[fG]+1,"snowsoil");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->G_snowsoil, UV,0,0);
	initmatrix(0.0, egy->G_snowsoil, top->Z0, NoV);
	free(temp);
	free(name);

	if(par->distr_stat==1){

		name=join_strings(files->co[fG]+1,"max");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->G_max, UV,0,0);
		initmatrix(-1.0E+9, egy->G_max, top->Z0, NoV);
		free(temp);
		free(name);

		name=join_strings(files->co[fG]+1,"min");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->G_min, UV,0,0);
		initmatrix(1.0E+9, egy->G_min, top->Z0, NoV);
		free(name);
		free(temp);

	}
#endif

}

//SENSIBLE HEAT FLUX
if(par->output_H>0 && fmod(times->time+par->Dt,par->output_H*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_H*3600.0));
#ifdef USE_NETCDF_MAP
	name=join_strings(files->co[fH]+1,"$_mean");
	write_map(name, 0, par->format_out, egy->H_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->H_mean, top->Z0, NoV);
	free(name);

	if(par->distr_stat==1){
		name=join_strings(files->co[fH]+1,"$_max");
		write_map(name, 0, par->format_out, egy->H_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->H_max, top->Z0, NoV);
		free(name);
		name=join_strings(files->co[fH]+1,"$_min");
		write_map(name, 0, par->format_out, egy->H_min, UV,times->time+par->Dt,n_file);
		initmatrix(1.0E+9, egy->H_min, top->Z0, NoV);
		free(name);
	}
#else
	write_suffix(SSSS, n_file, 0);

	//name=join_strings(files->co[fH]+1,"mean+");
	name=join_strings(files->co[fH]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->H_mean, UV,0,0);
	initmatrix(0.0, egy->H_mean, top->Z0, NoV);
	free(name);
	free(temp);
	//name=join_strings(files->co[fH]+1,"mean-");
	//write_map(join_strings(name,SSSS), 0, par->format_out, egy->H_mean2, UV);

	if(par->distr_stat==1){
		name=join_strings(files->co[fH]+1,"max");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->H_max, UV,0,0);
		initmatrix(-1.0E+9, egy->H_max, top->Z0, NoV);
		free(name);
		free(temp);
		name=join_strings(files->co[fH]+1,"min");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->H_min, UV,0,0);
		initmatrix(1.0E+9, egy->H_min, top->Z0, NoV);
		free(name);
		free(temp);
	}
#endif


}

//LATENT HEAT FLUX
if(par->output_ET>0 && fmod(times->time+par->Dt,par->output_ET*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_ET*3600.0));
#ifdef USE_NETCDF_MAP
	name=join_strings(files->co[fLE]+1,"$_mean");
	write_map(name, 0, par->format_out, egy->ET_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->ET_mean, top->Z0, NoV);
	free(name);

	if(par->distr_stat==1){
		name=join_strings(files->co[fLE]+1,"$_max");
		write_map(name, 0, par->format_out, egy->ET_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->ET_max, top->Z0, NoV);
		free(name);
		name=join_strings(files->co[fLE]+1,"$_min");
		write_map(name, 0, par->format_out, egy->ET_min, UV,times->time+par->Dt,n_file);
		initmatrix(1.0E+9, egy->ET_min, top->Z0, NoV);
		free(name);
	}
#else
	write_suffix(SSSS, n_file, 0);

	//name=join_strings(files->co[fLE]+1,"mean+");
	name=join_strings(files->co[fLE]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->ET_mean, UV,0,0);
	initmatrix(0.0, egy->ET_mean, top->Z0, NoV);
	//name=join_strings(files->co[fLE]+1,"mean-");
	//write_map(join_strings(name,SSSS), 0, par->format_out, egy->ET_mean2, UV);
	free(temp);
	free(name);

	if(par->distr_stat==1){
		name=join_strings(files->co[fLE]+1,"max");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->ET_max, UV,0,0);
		initmatrix(-1.0E+9, egy->ET_max, top->Z0, NoV);
		free(temp);
		free(name);
		name=join_strings(files->co[fLE]+1,"min");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->ET_min, UV,0,0);
		initmatrix(1.0E+9, egy->ET_min, top->Z0, NoV);
		free(name);
		free(temp);
	}
#endif
}

//SURFACE TEMPERATURE
if(par->output_Ts>0 && fmod(times->time+par->Dt,par->output_Ts*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Ts*3600.0));
#ifdef USE_NETCDF_MAP
	name=join_strings(files->co[fTs]+1,"$_mean");
	write_map(name, 0, par->format_out, egy->Ts_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->Ts_mean, top->Z0, NoV);
	free(name);
	if(par->distr_stat==1){
		name=join_strings(files->co[fTs]+1,"$_max");
		write_map(name, 0, par->format_out, egy->Ts_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->Ts_max, top->Z0, NoV);
		free(name);
		name=join_strings(files->co[fTs]+1,"$_min");
		write_map(name, 0, par->format_out, egy->Ts_min, UV,times->time+par->Dt,n_file);
		initmatrix(1.0E+9, egy->Ts_min, top->Z0, NoV);
		free(name);
	}
#else
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fTs]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->Ts_mean, UV,0,0);
	initmatrix(0.0, egy->Ts_mean, top->Z0, NoV);
	free(temp);
	free(name);
	if(par->distr_stat==1){
		name=join_strings(files->co[fTs]+1,"max");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->Ts_max, UV,0,0);
		initmatrix(-1.0E+9, egy->Ts_max, top->Z0, NoV);
		free(temp);
		free(name);
		name=join_strings(files->co[fTs]+1,"min");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->Ts_min, UV,0,0);
		initmatrix(1.0E+9, egy->Ts_min, top->Z0, NoV);
		free(name);
		free(temp);
	}
#endif

}

//PRECIPITATION
if(par->output_P>0 && fmod(times->time+par->Dt,par->output_P*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_P*3600.0));
#ifdef USE_NETCDF_MAP
	temp=join_strings(files->co[fprec]+1,"$_total");
	write_map(temp, 0, par->format_out, wat->PrTOT_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, wat->PrTOT_mean, top->Z0, NoV);
	free(temp);
	temp=join_strings(files->co[fprec]+1,"$_snow");
	write_map(temp, 0, par->format_out, wat->PrSNW_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, wat->PrSNW_mean, top->Z0, NoV);
	free(temp);
#else
	write_suffix(SSSS, n_file, 0);
	temp=join_strings(files->co[fprec]+1,"TOTAL");
	name=join_strings(temp,SSSS);
	write_map(name, 0, par->format_out, wat->PrTOT_mean, UV,0,0);
	initmatrix(0.0, wat->PrTOT_mean, top->Z0, NoV);
	free(name);
	free(temp);
	temp=join_strings(files->co[fprec]+1,"SNOW");
	name=join_strings(temp,SSSS);
	write_map(name, 0, par->format_out, wat->PrSNW_mean, UV,0,0);
	initmatrix(0.0, wat->PrSNW_mean, top->Z0, NoV);
	free(name);
	free(temp);
#endif
}

//INTERCEPTED PRECIPITATION
if(par->output_Wr>0 && fmod(times->time+par->Dt,par->output_Wr*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Wr*3600.0));
#ifdef USE_NETCDF_MAP
	write_map(files->co[fcint]+1, 0, par->format_out, wat->wcan_rain, UV,times->time+par->Dt,n_file);
#else
	write_suffix(SSSS, n_file, 0);
	name=join_strings(files->co[fcint]+1,SSSS);
	write_map(name, 0, par->format_out, wat->wcan_rain, UV,0,0);
	free(name);
#endif
}

//SNOW ENERGY BALANCE
if(par->output_balancesn>0 && fmod(times->time+par->Dt,par->output_balancesn*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_balancesn*3600.0));
#ifdef USE_NETCDF_MAP
	write_map(files->co[fmsn]+1, 0, par->format_out, snow->MELTED, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, snow->MELTED, top->Z0, NoV);
	write_map(files->co[fssn]+1, 0, par->format_out, snow->SUBL, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, snow->SUBL, top->Z0, NoV);
#else
	write_suffix(SSSS, n_file, 0);
	temp=join_strings(files->co[fmsn]+1,SSSS);
	write_map(temp, 0, par->format_out, snow->MELTED, UV,0,0);
	free(temp);
	initmatrix(0.0, snow->MELTED, top->Z0, NoV);
	temp=join_strings(files->co[fssn]+1,SSSS);
	write_map(temp, 0, par->format_out, snow->SUBL, UV,0,0);
	free(temp);
	initmatrix(0.0, snow->SUBL, top->Z0, NoV);
#endif

	for(r=1;r<=snow->Dzl->nrh;r++){
		for(c=1;c<=snow->Dzl->nch;c++){
			if(land->LC->co[r][c]!=NoV){
				M->co[r][c]=0.0;
				D=0.0;
				for(l=1;l<=snow->lnum->co[r][c];l++){
					M->co[r][c]+=(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c]);
					D+=0.001*snow->Dzl->co[l][r][c];
				}
				if(D>0) M->co[r][c]/=D;
			}
		}
	}
#ifdef USE_NETCDF_MAP
	write_map(files->co[fsnd]+1, 0, par->format_out, M, UV,times->time+par->Dt,n_file);
#else
	temp=join_strings(files->co[fsnd]+1,SSSS);
	write_map(temp, 0, par->format_out, M, UV,0,0);
	free(temp);
#endif
}

//GLACIER ENERGY BALANCE
if(par->glaclayer_max>0 && par->output_balancegl>0 && fmod(times->time+par->Dt,par->output_balancegl*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_balancegl*3600.0));
#ifdef USE_NETCDF_MAP
	write_map(files->co[fmgl]+1, 0, par->format_out, glac->MELTED, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, glac->MELTED, top->Z0, NoV);
	write_map(files->co[fsgl]+1, 0, par->format_out, glac->SUBL, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, glac->SUBL, top->Z0, NoV);
#else
	write_suffix(SSSS, n_file, 0);
	temp=join_strings(files->co[fmgl]+1,SSSS);
	write_map(temp, 0, par->format_out, glac->MELTED, UV,0,0);
	free(temp);
	initmatrix(0.0, glac->MELTED, top->Z0, NoV);
	temp=join_strings(files->co[fsgl]+1,SSSS);
	write_map(temp, 0, par->format_out, glac->SUBL, UV,0,0);
	free(temp);
	initmatrix(0.0, glac->SUBL, top->Z0, NoV);
#endif

	for(r=1;r<=glac->Dzl->nrh;r++){
		for(c=1;c<=glac->Dzl->nch;c++){
			if(land->LC->co[r][c]!=NoV){
				M->co[r][c]=0.0;
				D=0.0;
				for(l=1;l<=glac->lnum->co[r][c];l++){
					M->co[r][c]+=(glac->w_liq->co[l][r][c]+glac->w_ice->co[l][r][c]);
					D+=0.001*glac->Dzl->co[l][r][c];
				}
				if(D>0) M->co[r][c]/=D;
			}else{
				M->co[r][c]=UV->V->co[2];
			}
		}
	}
#ifdef USE_NETCDF_MAP
	write_map(files->co[fgld]+1, 0, par->format_out, M, UV,0,n_file);
#else
	temp=join_strings(files->co[fgld]+1,SSSS);
	write_map(temp, 0, par->format_out, M, UV,0,0);
	free(temp);
#endif

}

//SOLAR RADIATION
if(par->output_Rswdown>0 && fmod(times->time+par->Dt,par->output_Rswdown*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Rswdown*3600.0));
#ifdef USE_NETCDF_MAP
	name=join_strings(files->co[fSW]+1,"$_mean");
	write_map(name, 0, par->format_out, egy->Rswdown_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->Rswdown_mean, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fSW]+1,"$_beam_mean");
	write_map(name, 0, par->format_out, egy->Rswbeam, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->Rswbeam, top->Z0, NoV);
	free(name);

	if(par->distr_stat==1){
		name=join_strings(files->co[fSW]+1,"$_max");
		write_map(name, 0, par->format_out, egy->Rswdown_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->Rswdown_max, top->Z0, NoV);
		free(name);
	}

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				if(egy->nDt_sun->co[r][c]>0){
					M->co[r][c]=egy->nDt_shadow->co[r][c]/(double)(egy->nDt_sun->co[r][c]);
				}else{
					M->co[r][c]=0.0;
				}
			}
		}
	}
	name=join_strings(files->co[fSW]+1,"$_shadowfractime");
	write_map(name, 0, par->format_out, M, UV,times->time+par->Dt,n_file);
	initlongmatrix(0, egy->nDt_shadow, top->Z0, NoV);
	initlongmatrix(0, egy->nDt_sun, top->Z0, NoV);
	free(name);
#else
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fSW]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->Rswdown_mean, UV,0,0);
	initmatrix(0.0, egy->Rswdown_mean, top->Z0, NoV);
	free(name);
	free(temp);

	name=join_strings(files->co[fSW]+1,"beam_mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->Rswbeam, UV,0,0);
	initmatrix(0.0, egy->Rswbeam, top->Z0, NoV);
	free(name);
	free(temp);
	if(par->distr_stat==1){
		name=join_strings(files->co[fSW]+1,"max");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->Rswdown_max, UV,0,0);
		initmatrix(-1.0E+9, egy->Rswdown_max, top->Z0, NoV);
		free(name);
		free(temp);
	}

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				if(egy->nDt_sun->co[r][c]>0){
					M->co[r][c]=egy->nDt_shadow->co[r][c]/(double)(egy->nDt_sun->co[r][c]);
				}else{
					M->co[r][c]=0.0;
				}
			}
		}
	}
	name=join_strings(files->co[fSW]+1,"shadowfractime");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, M, UV,0,0);
	initlongmatrix(0, egy->nDt_shadow, top->Z0, NoV);
	initlongmatrix(0, egy->nDt_sun, top->Z0, NoV);
	free(name);
	free(temp);
#endif
}

//METEO
if(par->output_meteo>0 && fmod(times->time+par->Dt,par->output_meteo*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_meteo*3600.0));
#ifdef USE_NETCDF_MAP
	name=join_strings(files->co[fTa]+1,"$_mean");
	write_map(name, 0, par->format_out, egy->Ta_mean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, egy->Ta_mean, top->Z0, NoV);
	free(name);
	if(par->distr_stat==1){
		name=join_strings(files->co[fTa]+1,"$_max");
		write_map(name, 0, par->format_out, egy->Ta_max, UV,times->time+par->Dt,n_file);
		initmatrix(-1.0E+9, egy->Ta_max, top->Z0, NoV);
		free(name);
		name=join_strings(files->co[fTa]+1,"$_min");
		write_map(name, 0, par->format_out, egy->Ta_min, UV,times->time+par->Dt,n_file);
		initmatrix(1.0E+9, egy->Ta_min, top->Z0, NoV);
		free(name);
	}
#else
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fTa]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, egy->Ta_mean, UV,0,0);
	initmatrix(0.0, egy->Ta_mean, top->Z0, NoV);
	free(temp);
	free(name);
	if(par->distr_stat==1){
		name=join_strings(files->co[fTa]+1,"max");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->Ta_max, UV,0,0);
		initmatrix(-1.0E+9, egy->Ta_max, top->Z0, NoV);
		free(name);
		free(temp);
		name=join_strings(files->co[fTa]+1,"min");
		temp=join_strings(name,SSSS);
		write_map(temp, 0, par->format_out, egy->Ta_min, UV,0,0);
		initmatrix(1.0E+9, egy->Ta_min, top->Z0, NoV);
		free(name);
		free(temp);
	}
#endif

	if(par->micromet==1){
#ifdef USE_NETCDF_MAP
	name=join_strings(files->co[fwspd]+1,"$_mean");
	write_map(name, 0, par->format_out, met->Vspdmean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, met->Vspdmean, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[fwdir]+1,"$_mean");
	write_map(name, 0, par->format_out, met->Vdirmean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, met->Vdirmean, top->Z0, NoV);
	free(name);

	name=join_strings(files->co[frh]+1,"$_mean");
	write_map(name, 0, par->format_out, met->RHmean, UV,times->time+par->Dt,n_file);
	initmatrix(0.0, met->RHmean, top->Z0, NoV);
	free(name);
#else
	name=join_strings(files->co[fwspd]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, met->Vspdmean, UV,0,0);
	initmatrix(0.0, met->Vspdmean, top->Z0, NoV);
	free(temp);
	free(name);

	name=join_strings(files->co[fwdir]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, met->Vdirmean, UV,0,0);
	initmatrix(0.0, met->Vdirmean, top->Z0, NoV);
	free(temp);
	free(name);

	name=join_strings(files->co[frh]+1,"mean");
	temp=join_strings(name,SSSS);
	write_map(temp, 0, par->format_out, met->RHmean, UV,0,0);
	initmatrix(0.0, met->RHmean, top->Z0, NoV);
	free(temp);
	free(name);
#endif

	}
}

free_doublematrix(M);





/**********************************************************************************************************/
/**********************************************************************************************************/
//SPECIAL PLOTS AT SOME DAYS
/**********************************************************************************************************/
/**********************************************************************************************************/

if(times->i_plot==times->n_plot){

	printf("\nWriting output data for JD:%ld year:%ld file:%ld\n",times->d_plot, times->AAAA, times->nt_plot);
	f=fopen(error_file_name,"a");
	fprintf(f,"\nWriting output data for JD:%ld year:%ld file:%ld  ",times->d_plot, times->AAAA, times->nt_plot);
	write_date(f, times->DD, times->MM, times->AAAA, times->hh, times->mm);
	fprintf(f,"\n");
	fclose(f);

	n_file=(long)((times->time+par->Dt)/(times->n_plot*par->Dt));

	M=new_doublematrix(Nr,Nc);

	initialize_doublematrix(M,UV->V->co[2]);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				M->co[r][c]=egy->Hgplot->co[r][c] + egy->Hvplot->co[r][c];
			}
		}
	}
	plot(files->co[pH]+1, times->d_plot, times->AAAA, times->nt_plot, M, par->format_out);

	initialize_doublematrix(M,UV->V->co[2]);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				M->co[r][c]=egy->LEgplot->co[r][c] + egy->LEvplot->co[r][c];
			}
		}
	}
	plot(files->co[pLE]+1, times->d_plot, times->AAAA, times->nt_plot, M, par->format_out);

	plot(files->co[pHg]+1, times->d_plot, times->AAAA, times->nt_plot, egy->Hgplot, par->format_out);
	plot(files->co[pLEg]+1, times->d_plot, times->AAAA, times->nt_plot, egy->LEgplot, par->format_out);
	plot(files->co[pHv]+1, times->d_plot, times->AAAA, times->nt_plot, egy->Hvplot, par->format_out);
	plot(files->co[pLEv]+1, times->d_plot, times->AAAA, times->nt_plot, egy->LEvplot, par->format_out);
	plot(files->co[pSWin]+1, times->d_plot, times->AAAA, times->nt_plot, egy->SWinplot, par->format_out);
	plot(files->co[pSWg]+1, times->d_plot, times->AAAA, times->nt_plot, egy->SWgplot, par->format_out);
	plot(files->co[pSWv]+1, times->d_plot, times->AAAA, times->nt_plot, egy->SWvplot, par->format_out);
	plot(files->co[pLWin]+1, times->d_plot, times->AAAA, times->nt_plot, egy->LWinplot, par->format_out);
	plot(files->co[pLWg]+1, times->d_plot, times->AAAA, times->nt_plot, egy->LWgplot, par->format_out);
	plot(files->co[pLWv]+1, times->d_plot, times->AAAA, times->nt_plot, egy->LWvplot, par->format_out);
	plot(files->co[pTs]+1, times->d_plot, times->AAAA, times->nt_plot, egy->Tsplot, par->format_out);
	plot(files->co[pTg]+1, times->d_plot, times->AAAA, times->nt_plot, egy->Tgplot, par->format_out);
	plot(files->co[pTv]+1, times->d_plot, times->AAAA, times->nt_plot, egy->Tvplot, par->format_out);
	plot(files->co[pTa]+1, times->d_plot, times->AAAA, times->nt_plot, met->Taplot, par->format_out);
	plot(files->co[pD]+1, times->d_plot, times->AAAA, times->nt_plot, snow->Dplot, par->format_out);
	if(par->micromet==1){
		plot(files->co[pVspd]+1, times->d_plot, times->AAAA, times->nt_plot, met->Vspdplot, par->format_out);
		plot(files->co[pVdir]+1, times->d_plot, times->AAAA, times->nt_plot, met->Vdirplot, par->format_out);
		plot(files->co[pRH]+1, times->d_plot, times->AAAA, times->nt_plot, met->RHplot, par->format_out);
	}

	initialize_doublematrix(M,UV->V->co[2]);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				sy=sl->type->co[r][c];
				l=1;
				//printf("%ld %ld psi:%f\n",r,c,sl->P->co[l][r][c]);
				M->co[r][c]=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],
					sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],PSImin,par->Esoil);
			}
		}
	}
	plot(files->co[pth]+1, times->d_plot, times->AAAA, times->nt_plot, M, par->format_out);

	free_doublematrix(M);

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				egy->Hgplot->co[r][c]=0.0;
				egy->LEgplot->co[r][c]=0.0;
				egy->Hvplot->co[r][c]=0.0;
				egy->LEvplot->co[r][c]=0.0;
				egy->SWinplot->co[r][c]=0.0;
				egy->SWgplot->co[r][c]=0.0;
				egy->SWvplot->co[r][c]=0.0;
				egy->LWinplot->co[r][c]=0.0;
				egy->LWgplot->co[r][c]=0.0;
				egy->LWvplot->co[r][c]=0.0;
				egy->Tsplot->co[r][c]=0.0;
				egy->Tgplot->co[r][c]=0.0;
				egy->Tvplot->co[r][c]=0.0;

				snow->Dplot->co[r][c]=0.0;
				met->Taplot->co[r][c]=0.0;
				if(par->micromet==1){
					met->Vspdplot->co[r][c]=0.0;
					met->Vdirplot->co[r][c]=0.0;
					met->RHplot->co[r][c]=0.0;
				}
			}
		}
	}
}

/**********************************************************************************************************/
/**********************************************************************************************************/
//TRANSECTS
/**********************************************************************************************************/
/**********************************************************************************************************/

/*date_time(times->time, par->year0, par->JD0, 0.0, &jd, &day, &month, &year, &hour0, &min0);
date_time(times->time+par->Dt, par->year0, par->JD0, 0.0, &jd, &day, &month, &year, &hour1, &min1);
date_time(times->time+2.0*par->Dt, par->year0, par->JD0, 0.0, &jd, &day, &month, &year, &hour2, &min2);

M=new_doublematrix(Nr,Nc);
for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){

			if(snow->lnum->co[r][c]>=1){
				M->co[r][c]=snow->T->co[snow->lnum->co[r][c]][r][c];
			}else if(par->glaclayer_max>0){
				if(snow->lnum->co[r][c]==0 && glac->lnum->co[r][c]>=1){
					M->co[r][c]=glac->T->co[1][r][c];
				}else{
					M->co[r][c]=sl->T->co[1][r][c];
				}
			}else{
				M->co[r][c]=sl->T->co[1][r][c];
			}

		}else{
			M->co[r][c]=NoV;
		}
	}
}

for(j=0;j<2;j++){

	for(i=0;i<dim2(par->transect[j]);i++){

		if(par->transect[j][i][2]>=hour1+min1/60.0 && par->transect[j][i][2]<hour2+min2/60.0){

			par->vtrans[j][i]=interp_value(par->transect[j][i][0], par->transect[j][i][1], M, top->Z0)*(par->transect[j][i][2]-(hour1+min1/60.0))/((hour2+min2/60.0)-(hour1+min1/60.0));

			//printf("1. %ld %ld %f\n",i,j,par->vtrans[j][i]);

		}else if(par->transect[j][i][2]>hour0+min0/60.0 && par->transect[j][i][2]<hour1+min1/60.0){

			par->vtrans[j][i]+=interp_value(par->transect[j][i][0], par->transect[j][i][1], M, top->Z0)*((hour1+min1/60.0)-par->transect[j][i][2])/((hour1+min1/60.0)-(hour0+min0/60.0));

			//printf("2. %ld %ld %f %ld\n",i,j,par->vtrans[j][i],par->ibeg->co[j+1]);

			if(par->ibeg->co[j+1]==-1){
				par->ibeg->co[j+1]=0;
				par->cont_trans->co[j+1]+=1;
			}

			if(i!=0){
				//printf("%ld %ld %f\n",j,i,pow(pow(par->transect[j][i][0]-par->transect[j][i-1][0],2.0)+pow(par->transect[j][i][1]-par->transect[j][i-1][1],2.0),0.5));

				if(pow(pow(par->transect[j][i][0]-par->transect[j][i-1][0],2.0)+pow(par->transect[j][i][1]-par->transect[j][i-1][1],2.0),0.5)>40){
					par->ibeg->co[j+1]=i;
					par->cont_trans->co[j+1]+=1;
				}
			}

			if(fabs(par->vtrans[j][i]-NoV)>1.0E-2){
				f=fopen(namefile_i(join_strings(WORKING_DIRECTORY,"_transectMOD"),j+1),"a");
				fprintf(f,"%ld,%ld,%f,%f,%f,%f,%f,%f\n",par->cont_trans->co[j+1],i,pow(pow(par->transect[j][i][0]-par->transect[j][par->ibeg->co[j+1]][0],2.0)+pow(par->transect[j][i][1]-par->transect[j][par->ibeg->co[j+1]][1],2.0),0.5),
					par->transect[j][i][0],par->transect[j][i][1],par->transect[j][i][2],par->transect[j][i][3],par->vtrans[j][i]);
				fclose(f);
			}
		}
	}


}

free_doublematrix(M);*/



/**********************************************************************************************************/
/**********************************************************************************************************/
//SAVING POINTS
/**********************************************************************************************************/
/**********************************************************************************************************/

if(times->time==0) isavings=0;

if(isavings < par->saving_points->nh){

	if(par->saving_points->nh==1 && par->saving_points->co[1]==0.0){

		isavings=1;

	}else{

		if(times->time+par->Dt >= 86400.0*par->saving_points->co[isavings+1]){
			isavings+=1;

			write_suffix(SSSS, isavings, 0);

#ifdef USE_NETCDF_MAP
			//force output map in esriascii format
			short int tmp_format_output;
			if (par->format_out >= 4){
				tmp_format_output=par->format_out;
				par->format_out=3;
			}
#endif

			for(l=1;l<=Nl;l++){
				write_tensorseries(1, l, isavings, files->co[rpsi]+1, 0, par->format_out, sl->P, UV,0);
				write_tensorseries(1, l, isavings, files->co[riceg]+1, 0, par->format_out, sl->thice, UV,0);
				write_tensorseries(1, l, isavings, files->co[rTg]+1, 0, par->format_out, sl->T, UV,0);
			}

			temp=join_strings(files->co[rhsup]+1,SSSS);
			write_map(temp, 0, par->format_out, wat->h_sup, UV,0,0); //USE_NETCDF_MAP
			free(temp);
			temp=join_strings(files->co[rwcrn]+1,SSSS);
			write_map(temp, 0, par->format_out, wat->wcan_rain, UV,0,0); //USE_NETCDF_MAP
			free(temp);
			temp=join_strings(files->co[rwcsn]+1,SSSS);
			write_map(temp, 0, par->format_out, wat->wcan_snow, UV,0,0); //USE_NETCDF_MAP
			free(temp);
			temp=join_strings(files->co[rTv]+1,SSSS);
			write_map(temp, 0, par->format_out, sl->Tv, UV,0,0); //USE_NETCDF_MAP
			free(temp);




			for(l=1;l<=par->snowlayer_max;l++){
				write_tensorseries(1, l, isavings, files->co[rDzs]+1, 0, par->format_out, snow->Dzl, UV,0);
				write_tensorseries(1, l, isavings, files->co[rwls]+1, 0, par->format_out, snow->w_liq, UV,0);
				write_tensorseries(1, l, isavings, files->co[rwis]+1, 0, par->format_out, snow->w_ice, UV,0);
				write_tensorseries(1, l, isavings, files->co[rTs]+1, 0, par->format_out, snow->T, UV,0);

			}
			temp=join_strings(files->co[rsnag_adim]+1,SSSS);
			write_map(temp, 0, par->format_out, snow->nondimens_age, UV,0,0); //USE_NETCDF_MAP
			free(temp);
			temp=join_strings(files->co[rsnag_dim]+1,SSSS);
			write_map(temp, 0, par->format_out, snow->dimens_age, UV,0,0); //USE_NETCDF_MAP
			free(temp);
			M=copydouble_longmatrix(snow->lnum);
			temp=join_strings(files->co[rns]+1, SSSS);
			write_map(temp, 1, par->format_out, M, UV,0,0); //USE_NETCDF_MAP
			free(temp);
			free_doublematrix(M);


			if(par->glaclayer_max>0){
				for(l=1;l<=par->glaclayer_max;l++){
					write_tensorseries(1, l, isavings, files->co[rDzi]+1, 0, par->format_out, glac->Dzl, UV,0);
					write_tensorseries(1, l, isavings, files->co[rwli]+1, 0, par->format_out, glac->w_liq, UV,0);
					write_tensorseries(1, l, isavings, files->co[rwii]+1, 0, par->format_out, glac->w_ice, UV,0);
					write_tensorseries(1, l, isavings, files->co[rTi]+1, 0, par->format_out, glac->T, UV,0);

				}
				M=copydouble_longmatrix(glac->lnum);
				temp=join_strings(files->co[rni]+1,SSSS);
				write_map(temp, 1, par->format_out, M, UV,0,0); //USE_NETCDF_MAP
				free(temp);
				free_doublematrix(M);

			}

#ifdef USE_NETCDF_MAP
			//restore original output map parameter
			if (par->format_out >= 4){
				par->format_out=tmp_format_output;
			}
#endif
			temp=join_strings(files->co[rQch]+1,SSSS);
			name=join_strings(temp,textfile);
			f=t_fopen(name,"w");
			fprintf(f,"/**Channel discharges(m3/s) for each channel-pixel\n");
			fprintf(f," Q_sup_s		Q_sub_s*/\n");
			fprintf(f,"index{1}\n");
			fprintf(f,"1:double matrix Q_channel {%ld,2}\n",cnet->Q_sup_s->nh);
			for(j=1;j<=cnet->Q_sup_s->nh;j++){
				fprintf(f,"%20.16f  %20.16f\n",cnet->Q_sup_s->co[j],cnet->Q_sub_s->co[j]);
			}
			t_fclose(f);
			free(name);
			free(temp);
			for(i=1;i<=par->nLC;i++){
				temp=join_strings(files->co[rSFA]+1,SSSS);
				temp1=join_strings(temp,"L");
				name=namefile_i(temp1,i);
				f=t_fopen(name,"w");
				fprintf(f,"/**Fraction of snow free area and average sensible heat flux (W/m2) from snow free area\n");
				fprintf(f," VSFA		HSFA*/\n");
				fprintf(f,"index{1}\n");
				fprintf(f,"1:double vector SFA {%f,%f}\n",egy->VSFA->co[i],egy->HSFA->co[i]);
				t_fclose(f);
				free(name);
				free(temp);
				free(temp1);
			}

		}
	}
}
// file bugs
/*long cell1=1, cell2=16,cell3=26;
double T1,Thw1,Thi1,lambdat1,CT1,T2,Thw2,Thi2,lambdat2,CT2,T3,Thw3,Thi3,lambdat3,CT3;
r=1;c=1; // just for the 1D version, otherwise has to be changed
long ns=snow->lnum->co[r][c];
long ng=0;
long nsng=ns+ng;
char* namebugs=join_strings(files->element[fbugs]+1,textfile);
f=t_fopen(namebugs,"a");
date_time(times->time+par->Dt, par->year0, par->JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
write_date(f, d2, mo2, y2, h2, mi2);
fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
l=cell1;
	T1=sl->T->co[l][r][c];
	Thw1=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
			sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],PSImin, par->Esoil);
	Thi1=sl->thice->co[l][r][c];
	lambdat1=k_thermal_soil(Thw1,Thi1,sl->pa->co[sy][jsat][l],T1, sl->pa->co[sy][jkt][l]);
	CT1=sl->pa->co[sy][jct][l-nsng]*(1.-sl->pa->co[sy][jsat][l-nsng]) + 1000*c_ice*Thi1 + 1000*c_liq*Thw1;
l=cell2;
	T2=sl->T->co[l][r][c];
	Thw2=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
			sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],PSImin, par->Esoil);
	Thi2=sl->thice->co[l][r][c];
	lambdat2=k_thermal_soil(Thw2,Thi2,sl->pa->co[sy][jsat][l],T2, sl->pa->co[sy][jkt][l]);
	CT2=sl->pa->co[sy][jct][l-nsng]*(1.-sl->pa->co[sy][jsat][l-nsng]) + 1000*c_ice*Thi2 + 1000*c_liq*Thw2;
l=cell3;
	T3=sl->T->co[l][r][c];
	Thw3=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
			sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],PSImin, par->Esoil);
	Thi3=sl->thice->co[l][r][c];
	lambdat3=k_thermal_soil(Thw3,Thi3,sl->pa->co[sy][jsat][l],T3, sl->pa->co[sy][jkt][l]);
	CT3=sl->pa->co[sy][jct][l-nsng]*(1.-sl->pa->co[sy][jsat][l-nsng]) + 1000*c_ice*Thi3 + 1000*c_liq*Thw3;
fprintf(f,",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",T1,T2,T3,Thw1,Thw2,Thw3,Thi1,Thi2,Thi3,lambdat1,lambdat2,lambdat3,CT1,CT2,CT3);
fclose(f);
free(namebugs);*/
}

/*--------------------------------------------------------------------------------------------------*/













/****************************************************************************************************/
/* All the structs and substructs of the simulation are deallocated:                                */
/****************************************************************************************************/
void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met, TIMES *times)
{
long r,c,l; /* added by Emanuele Cordano for de-allocating double or triple pointers; */
 /* Deallocation of struct TOPO "top": */
 printf("Deallocating top\n");
 free_doublematrix(top->Z0);//mance z1 z0dp;
 //if (top->Z1!=NULL) free_doublematrix(top->Z1); /* added by Emanuele Cordano */
 free_doublematrix(top->sky);
 free_shortmatrix(top->pixel_type);
 free_shortmatrix(top->DD);
 free_doublematrix(top->i_DD);
 free_doublematrix(top->dz_dx);
 free_doublematrix(top->dz_dy);
 free_shortmatrix(top->curv);
 /*free_doublematrix(top->area); */
 free_doublematrix(top->aspect);
 free_doublematrix(top->slopes);
 free_doublematrix(top->i_ch);
 free_doublematrix(top->pixel_distance);
 if(par->ES_num>0){
	free_longvector(top->ES_pixel);
	free_doublevector(top->ES_aspect);
	free_doublevector(top->ES_slope);
 }
 if(par->point_sim==1) {
	 int num_points=top->Z0->nch;
	 for(r=1;r<=top->Z0->nrh;r++){// top->Z0->nrh=1 therefore there is only one line...I don't know why
		 for(c=1;c<=num_points;c++){// for every point of the simulation 1D
			 if(top->horizon_height[r-1][c-1]!=NULL) {
				printf("\ndeallocating matrix top->horizon_height[%ld][%ld]",r,c);
				free_alloc2(&top->horizon_height[r-1][c-1]);
			 }
			 /*if(top->horizon_height[r-1][c-1]!=NULL){
				 printf("\ndeallocating vector top->horizon_height[%ld][%ld]",r,c);
				 free(top->horizon_height[r-1][c-1]);
			 }*/
			 if(top->horizon_height[r-1]!=NULL){
				 printf("\nr=%ld, c=%ld, ciao3",r,c);
				 free(top->horizon_height[r-1]);
			}
		 }
	 }
	 if(top->horizon_height!=NULL)	 free(top->horizon_height);
	 /* added by Emanuele Cordano on 22 September 2009
	 r=0;
	 c=0;
	 l=0;
	 while (top->horizon_height[l]!=NULL) {
		r=0;
	 	while (top->horizon_height[l][r]!=NULL) {
			c=0;
	//		 while(top->horizon_height[l][r][c]!=NULL){
	//			 free(top->horizon_height[l][r][c]);
	//			 c++;
	//		 }
	  		 free(top->horizon_height[l][r]);
	  		 r++;
	  	 }
	 	free(top->horizon_height[l]);
	  	l++;
	 }*/
 }
 if(par->point_sim==1 && par->micromet==1) free_doublematrix(top->Z1);
 if(par->micromet==1){
	free_doublematrix(top->Zm);
	free_doublematrix(top->curv_m);
	free_doublematrix(top->slope_m);
	free_doublematrix(top->slopeaz_m);
 }
 free_longmatrix(top->lrc_cont);

/*deallocating top content data
  *
  * \author Emanuele Cordano
  * \date September 2009
 *
 *
 */


// long r,c,l;
r=0;
c=0;
l=1;
while (top->i_cont[l]!=NULL && l<=Nl) {
//		r=0;
//		while (top->i_cont[l][r]!=NULL) {
//			/*(while (top->i_cont[l][r][c]!=NULL) { CONTROLLARE QUI
//				free(top->i_cont[l][r][c]);
//				c++;
//			}*/
// 		 free(top->i_cont[l][r]);
// 		 r++;
// 	 }
	free(top->i_cont[l]);
 	l++;
}

free(top->i_cont);

free_doubletensor(top->Z);
free_doublematrix(top->slope_H);
free(top);

 /* Deallocation of struct SOIL "sl": */
 printf("Deallocating sl\n");
 free_doubletensor(sl->P);
 free_doubletensor(sl->T);
 free_doublematrix(sl->Tv);
 free_doublematrix(sl->Tmean);
 free_doublematrix(sl->psi_mean);
 free_doublematrix(sl->thetai_mean);
 free_doublematrix(sl->thetaw_mean);
 free_doublematrix(sl->Tmin);
 free_doublematrix(sl->thetai_min);
 free_doublematrix(sl->thetaw_min);
 free_doublematrix(sl->Tmax);
 free_doublematrix(sl->thetai_max);
 free_doublematrix(sl->thetaw_max);

 free_doubletensor(sl->thice);
 free_doublematrix(sl->Jinf);
 free_shortmatrix(sl->bc);
 free_doubletensor(sl->J);
 free_shortmatrix(sl->type);
 free_doubletensor(sl->pa);
 int i;
 if(par->superfast==1){
	 for(i=1;i<par->num_of_time-1;i++){
		 //if(sl->output[i]!=NULL)
			 free_alloc2(&sl->output[i]);
	 }
 }
 free(sl);

 /* Deallocation of struct LAND "land": */
 printf("Deallocating land\n");
 free_doublematrix(land->LC);
 if(par->output_albedo>0) free_doublematrix(land->albedo);
 free_shortmatrix(land->shadow);
 free_longvector(land->clax);
 free_longmatrix(land->cont);
 free_doublematrix(land->ty);
 if(files->index->nh>nfiles && par->point_sim==0){if(existing_file(files->co[nfiles+1]+1)>0) free_shortmatrix(land->LC2);}
 free(land);

 /* Deallocation of struct WATER "water": */
 printf("Deallocating water\n");
 free_doublematrix(wat->weights_Kriging);
 free_doublematrix(wat->h_sup);
 free_doublematrix(wat->q_sup);
 free_doubletensor(wat->q_sub);
 free_doublematrix(wat->total);
 free_doublematrix(wat->Pn);
 free_doublematrix(wat->wcan_rain);
 free_doublematrix(wat->wcan_snow);

 if (par->output_P>0){
	free_doublematrix(wat->PrTOT_mean);
	free_doublematrix(wat->PrSNW_mean);
 }
 free_doublematrix(wat->out1);
 free_doublevector(wat->out2);
 free_doublematrix(wat->outfluxes);
 if(par->output_h_sup>0) free_doublematrix(wat->hsupav);
 free(wat);

 /* Deallocation of struct CHANNEL "channel": */
 printf("Deallocating channel network\n");
 free_longvector(cnet->r);
 free_longvector(cnet->c);
 if(par->point_sim!=1) free_longmatrix(cnet->ch);
 free_doublevector(cnet->Q);
 free_doublevector(cnet->s0);
 free_doublematrix(cnet->fraction_spread);
 free_doublevector(cnet->Q_sup_s);
 free_doublevector(cnet->Q_sub_s);
 free_doublevector(cnet->Qsup_spread);
 free_doublevector(cnet->Qsub_spread);
 free_doublevector(cnet->Qsup);
 free_doublevector(cnet->Qsub);
 free(cnet);

 /* Deallocation of struct ENERGY "egy": */
 printf("Deallocating egy\n");
 free_doublematrix(egy->Hgrid);
 free_doublematrix(egy->Tsgrid);

 if(par->output_Rn>0){
	free_doublematrix(egy->Rn_mean);
	if(par->distr_stat==1)free_doublematrix(egy->Rn_max);
	if(par->distr_stat==1)free_doublematrix(egy->Rn_min);
	if(par->distr_stat==1)free_doublematrix(egy->LW_max);
	if(par->distr_stat==1)free_doublematrix(egy->LW_min);
	free_doublematrix(egy->LW_in);
	free_doublematrix(egy->LW_out);
	free_doublematrix(egy->SW);
	if(par->distr_stat==1)free_doublematrix(egy->SW_max);
 }

 if(par->output_ET>0){
	free_doublematrix(egy->ET_mean);
	//free_doublematrix(egy->ET_mean2);
	if(par->distr_stat==1)free_doublematrix(egy->ET_max);
	if(par->distr_stat==1)free_doublematrix(egy->ET_min);
 }

 if(par->output_H>0){
	free_doublematrix(egy->H_mean);
	//free_doublematrix(egy->H_mean2);
	if(par->distr_stat==1)free_doublematrix(egy->H_max);
	if(par->distr_stat==1)free_doublematrix(egy->H_min);
 }

 if(par->output_G>0){
	free_doublematrix(egy->G_mean);
	if(par->distr_stat==1)free_doublematrix(egy->G_max);
	if(par->distr_stat==1)free_doublematrix(egy->G_min);
	free_doublematrix(egy->G_snowsoil);
 }
 if(par->output_Ts>0){
	free_doublematrix(egy->Ts_mean);
	if(par->distr_stat==1)free_doublematrix(egy->Ts_max);
	if(par->distr_stat==1)free_doublematrix(egy->Ts_min);
 }
 if(par->output_Rswdown>0){
	free_doublematrix(egy->Rswdown_mean);
	if(par->distr_stat==1)free_doublematrix(egy->Rswdown_max);
	free_doublematrix(egy->Rswbeam);
 }

 free_longmatrix(egy->nDt_shadow);
 free_longmatrix(egy->nDt_sun);

 if(par->output_meteo>0){
	free_doublematrix(egy->Ta_mean);
	if(par->distr_stat==1)free_doublematrix(egy->Ta_max);
	if(par->distr_stat==1)free_doublematrix(egy->Ta_min);
 }

 free_doublematrix(egy->out1);
 free_doublevector(egy->out2);

 if(par->ES_num>0) free_doublematrix(egy->out3);

 if(par->JD_plots->co[1]>=0){
	free_doublematrix(egy->Hgplot);
	free_doublematrix(egy->LEgplot);
	free_doublematrix(egy->Hvplot);
	free_doublematrix(egy->LEvplot);
	free_doublematrix(egy->SWinplot);
	free_doublematrix(egy->SWgplot);
	free_doublematrix(egy->SWvplot);
	free_doublematrix(egy->LWinplot);
	free_doublematrix(egy->LWgplot);
	free_doublematrix(egy->LWvplot);
	free_doublematrix(egy->Tsplot);
	free_doublematrix(egy->Tgplot);
	free_doublematrix(egy->Tvplot);
 }

 free_doublevector(egy->HSFA);
 free_doublevector(egy->VSFA);

 free(egy);
 /* Deallocation of struct SNOW "snow": */
 printf("Deallocating snow\n");
 if(par->JD_plots->co[1]>=0) free_doublematrix(snow->Dplot);
 free_shortmatrix(snow->type);
 free_longmatrix(snow->lnum);
 free_doubletensor(snow->Dzl);
 free_doubletensor(snow->w_liq);
 free_doubletensor(snow->w_ice);
 free_doubletensor(snow->T);
 free_doublematrix(snow->rho_newsnow);
 free_doublematrix(snow->Psnow);
 if(par->blowing_snow==1){
	free_doublematrix(snow->Wtrans);
	free_doublematrix(snow->Qsub);
	free_doublematrix(snow->Qtrans);
	free_doublematrix(snow->Qtrans_x);
	free_doublematrix(snow->Qtrans_y);
	//free_doublematrix(snow->ListonSWE);
	//free_doublematrix(snow->softSWE);
	//free_doublematrix(snow->softSWE1);
	if(par->output_snow>0){
		free_doublematrix(snow->Wtot);
		/*free_doublematrix(snow->Wtrans_cum);
		free_doublematrix(snow->Wsusp_cum);
		free_doublematrix(snow->Wsubl_cum);
		free_doublematrix(snow->Wsubgrid_cum);*/
	}
 }
 free_doublematrix(snow->out_bs);
 free_doublematrix(snow->nondimens_age);
 free_doublematrix(snow->dimens_age);
 free_doublevector(snow->evap);
 free_doublevector(snow->subl);
 free_doublevector(snow->melted);
 if(par->output_snow>0){
	free_doublematrix(snow->max);
	free_doublematrix(snow->average);
 }
 if(par->output_balancesn>0){
	free_doublematrix(snow->MELTED);
	free_doublematrix(snow->SUBL);
 }
 free_doublevector(snow->CR1);
 free_doublevector(snow->CR2);
 free_doublevector(snow->CR3);
 free_doublevector(snow->CR1m);
 free_doublevector(snow->CR2m);
 free_doublevector(snow->CR3m);
 if(par->blowing_snow==1) free_longvector(snow->change_dir_wind);
 free(snow);

 printf("Deallocating glacier\n");
 free_doublevector(glac->evap);
 free_doublevector(glac->subl);
 free_doublevector(glac->melted);
 if(par->glaclayer_max>0){
	free_longmatrix(glac->lnum);
	free_doubletensor(glac->Dzl);
	free_doubletensor(glac->w_liq);
	free_doubletensor(glac->w_ice);
	free_doubletensor(glac->T);
	if(par->output_balancegl>0) free_doublematrix(glac->MELTED);
	if(par->output_balancegl>0)	free_doublematrix(glac->SUBL);
 }
 free(glac);

 printf("Deallocating met\n");
 free_doublematrix(met->Tgrid);
 free_doublematrix(met->Pgrid);
 if(par->micromet==1){
	free_doublematrix(met->Vgrid);
	free_doublematrix(met->Vdir);
	free_doublematrix(met->RHgrid);
	if(par->output_meteo>0){
		free_doublematrix(met->Vspdmean);
		free_doublematrix(met->Vdirmean);
		free_doublematrix(met->RHmean);
	}
 }
 if(par->JD_plots->co[1]>=0){
	free_doublematrix(met->Taplot);
	if(par->micromet==1){
		free_doublematrix(met->Vspdplot);
		free_doublematrix(met->Vdirplot);
		free_doublematrix(met->RHplot);
	}
 }

 /*! deallocating met content data
  *
  * \author Emanuele Cordano
  * \date September 2009
  * */
 //long r,c,l;
 r=0;
 c=0;
 l=0;
 /* while (met->data[l]!=NULL) {
	 r=0;
	 while (met->data[l][r]!=NULL) {
	 while (met->data[l][r][c]!=NULL) {
			 free(met->data[l][r][c]);
		 }*/
	 /*
		 free(met->data[l][r]);
		 r++;
	 }
	free(met->data[l]);
	l++;
 }

 r=0;
 c=0;
 l=0;
 while (met->horizon[l]!=NULL) {
	 r=0;
 	 while (met->horizon[l][r]!=NULL) {

 		 free(met->horizon[l][r]);
 		 r++;
 	 }
 	free(met->horizon[l]);
 	l++;
  }
 r=0;
 c=0;
 l=0;
 while (met->column[l]!=NULL) {
	 free(met->column[l]);
	 l++;
 }
 l=0;
  while (met->var[l]!=NULL) {
 	 free(met->var[l]);
 	 l++;
 }*/
 /* END deallocating met */

 int num_meteo_st=met->st->E->nh;
 for(i=0;i<num_meteo_st;i++){
	 //printf("\nfree data of meteo_st %d",i+1);stop_execution();
	 if(met->data[i]!=NULL)
		free_alloc2(&met->data[i]);
	 if(met->column[i]!=NULL)
		free(met->column[i]);// deallocates each dynamic vector
	 if(met->horizon[i]!=NULL)
		free_alloc2(&met->horizon[i]);
 }
 free(met->data);

 //if(met->var[i]!=NULL) free_alloc2(&met->var);
 //free(met->var);// was like this
 //free(met->horizon);// was like this
 free(met->column);



 /*l=0;
 while(met->LRs[l]!=NULL){
	 free(met->LRs[l]);
	 l++;
 }*/

 if(existing_file_text(files->co[fLRs]+1)==1){
	 if(met->LRs!=NULL)
	 free_alloc2(&met->LRs);// was free(met->LRs)
 }
 free(met->LRv);
 free_longvector(met->LRp);

 dealloc_meteostations(met->st);

 free(met);

 /* Deallocation of struct PAR "par": */
 printf("Deallocating par\n");
 free_doublevector(par->Dmin);
 free_doublevector(par->Dmax);
 if(par->glaclayer_max>0){
	free_doublevector(par->Dmin_glac);
	free_doublevector(par->Dmax_glac);
 }
 free_doublematrix(par->chkpt);
 if(par->state_pixel==1) free_longmatrix(par->rc);
 free_doublevector(par->saving_points);
 free_longvector(par->JD_plots); /* corrected by Emanuele Cordano on 21 September 2009 */
 if(par->point_sim==1){
	if(par->micromet==1){
		free_longvector(par->r_points);
		free_longvector(par->c_points);
	}
 }
 /*free(par->transect);
 free(par->vtrans);
 free_longvector(par->cont_trans);
 free_longvector(par->ibeg);*/
 free(par);

  /* Deallocation of struct FILENAMES "filenames": */
 printf("Deallocating files\n");
 free_stringbin(files);
 free(error_file_name);
 /* Deallocation of struct T_INIT "UV": */
 printf("Deallocating UV\n");
 free_doublevector(UV->U);
 free_doublevector(UV->V);
 free(UV);

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void dealloc_meteostations(METEO_STATIONS *st){

	free_doublevector(st->E);
	free_doublevector(st->N);
	free_doublevector(st->lat);
	free_doublevector(st->lon);
	free_doublevector(st->Z);
	free_doublevector(st->sky);
	free_doublevector(st->ST);
	free_doublevector(st->Vheight);
	free_doublevector(st->Theight);
	free_doublevector(st->JD0);
	free_longvector(st->Y0);
	free_doublevector(st->Dt);
	free_longvector(st->offset);
	free(st);

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/



void write_date(FILE *f, long day, long month, long year, long hour, long min)

{
	if(day<10){
		fprintf(f,"0%1ld/",day);
	}else{
		fprintf(f,"%2ld/",day);
	}

	if(month<10){
		fprintf(f,"0%1ld/",month);
	}else{
		fprintf(f,"%2ld/",month);
	}

	fprintf(f,"%4ld ",year);

	if(hour<10){
		fprintf(f,"0%1ld:",hour);
	}else{
		fprintf(f,"%2ld:",hour);
	}

	if(min<10){
		fprintf(f,"0%1ld",min);
	}else{
		fprintf(f,"%2ld",min);
	}
}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void plot(char *name, long JD, long y, long i, DOUBLEMATRIX *M, short format){
	// USE_NETCDF_MAP
	//force output map in esriascii format
	if (format >= 4) {
		printf("Warning: output map are writtenm in esriiasci format. Poject not built to support NetCDF files \n");//added by Emanuele Cordano
		format=3;
		}
	// USE_NETCDF_MAP

	char ADS[ ]={"aaaaddddLssss"};

	write_suffix(ADS, y, 0);
	write_suffix(ADS, JD, 4);
	write_suffix(ADS, i, 9);
	write_map(join_strings(name,ADS), 0, format, M, UV,0,0); //USE_NETCDF_MAP

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void write_init_condit(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac){

	/*internal auxiliary variables:*/
	long i,j,r,c,l;
	char *name,SSSS[ ]={"SSSS"};
	double total_pixel, z;
	short sy,lu;
	DOUBLEVECTOR *root_fraction;
	FILE *f;
	char *temp1,*temp; /* added by Emanuele Cordano on 21 september 2009 */
	//DISCHARGE
	temp=join_strings(files->co[fQ]+1,textfile);
	f=t_fopen(temp,"w");
	free(temp);
	fprintf(f,"DATE,JDfrom0,JD,Q_tot[mc/s],Qsub_ch[mc/s],Qsup_ch[mc/s],Q_G[mc/s]\n");
	t_fclose(f);

	root_fraction=new_doublevector(Nl);

	//MELT FLUXES
	if(par->output_balancesn!=0){
		for(i=1;i<=land->clax->nh;i++){
			write_suffix(SSSS, i, 0);
			temp=join_strings(files->co[fmeltlu]+1,"_snowcovered_");
			temp1=join_strings(temp,SSSS);
			name=join_strings(temp1,textfile);/* bug found by Xujun on 12/11/09*/
			free(temp);
			free(temp1);
			f=t_fopen(name,"w");
			fprintf(f,"/** LAND USE CLASS: %ld\n",land->clax->co[i]);
			fprintf(f,"SNOW COVERED AREA\n");
			total_pixel=-99; /* data unavailable note by Emanuele Cordano on 24/9/9 */
			fprintf(f,"Number of pixels for this land use: %ld/%ld\n",land->cont->co[i][1],(long)total_pixel);
			fprintf(f,"Surface of this land use (km2): %10.4f/%10.4f */\n",land->cont->co[i][1]*UV->U->co[1]*UV->U->co[2]*1.0E-6,
				total_pixel*UV->U->co[1]*UV->U->co[2]*1.0E-6);
			fprintf(f,"1)t0[s],2)t1[s],3)t[d], 4)rain[mm],5)meltsnow[mm],6)meltglac[mm],7)pixels,8)area[km2],9)frac.landuse,10)frac.basin\n");
			t_fclose(f);
			free(name);

			temp=join_strings(files->co[fmeltlu]+1,"_snowfree_");
			temp1=join_strings(temp,SSSS);
			name=join_strings(temp1,textfile);
			free(temp1);
			free(temp);
			f=t_fopen(name,"w");
			fprintf(f,"/** LAND USE CLASS: %ld\n",land->clax->co[i]);
			fprintf(f,"SNOW FREE AREA\n");
			fprintf(f,"Number of pixels for this land use: %ld/%ld\n",land->cont->co[i][1],(long)total_pixel);
			fprintf(f,"Surface of this land use (km2): %10.4f/%10.4f\n",land->cont->co[i][1]*UV->U->co[1]*UV->U->co[2]*1.0E-6,
				total_pixel*UV->U->co[1]*UV->U->co[2]*1.0E-6);
			fprintf(f,"1)t0[s],2)t1[s],3)t[d], 4)rain[mm],5)meltsnow[mm],6)meltglac[mm],7)pixels,8)area[km2],9)frac.landuse,10)frac.basin\n");
			t_fclose(f);
			free(name);

		}
	}




	//DATA POINTS
	for(i=1;i<=par->chkpt->nrh;i++){

		write_suffix(SSSS, i, 0);
		r=par->rc->co[i][1];
		c=par->rc->co[i][2];
		sy=sl->type->co[r][c];
		lu=(short)land->LC->co[r][c];

		for(l=1;l<=Nl;l++){
			wat->out1->co[2][i]+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
								sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
								PSImin,par->Esoil);
		}
		wat->out1->co[1][i]=wat->h_sup->co[r][c];
		wat->out1->co[14][i]=wat->wcan_rain->co[r][c];


		temp=join_strings(files->co[fpoint]+1,"_info_");
		temp1=join_strings(temp,SSSS);
		name=join_strings(temp1,textfile);
		free(temp);
		free(temp1);
		f=t_fopen(name,"w");
		fprintf(f,"GEOtop KMackenzie: Summary of the main properties of the simulation in the point\nEast [m],North [m]\n%.2f,%.2f\nrow identification number,col identification number,Number of layers:\n%ld,%ld,%ld\n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c,Nl);
		fprintf(f,"Integration time, Output time interval (sec)\n%.0f,%.0f\n",par->Dt,(double)(times->n_pixel*par->Dt));
		fprintf(f,"soil layer depth [mm]\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.0f,",sl->pa->co[sy][jdz][l]);
			}fprintf(f,"%.0f\n",sl->pa->co[sy][jdz][Nl]);
		fprintf(f,"z coordinate [mm] of the center of each soil layer\n");
		z=0.0;
		for(l=1;l<=Nl-1;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,"%.2f,",z-0.5*sl->pa->co[sy][jdz][l]);
			}
		fprintf(f,"%.2f\n",z+0.5*sl->pa->co[sy][jdz][Nl]);
		fprintf(f,"Residual water content[-] in each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%f,",sl->pa->co[sy][jres][l]);
			}
		fprintf(f,"%f\n",sl->pa->co[sy][jres][Nl]);
		fprintf(f,"Saturated water content[-] in each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.2f,",sl->pa->co[sy][jsat][l]);
			}
		fprintf(f,"%.2f\n",sl->pa->co[sy][jsat][Nl]);
		fprintf(f,"Alpha of Van Genuchten [mm^-1] for each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.2f,",sl->pa->co[sy][ja][l]);
			}
		fprintf(f,"%.2f\n",sl->pa->co[sy][ja][Nl]);
		fprintf(f,"n of Van Genuchten [-] for each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.2f,",sl->pa->co[sy][jns][l]);
			}
		fprintf(f,"%.2f\n",sl->pa->co[sy][jns][Nl]);
		fprintf(f,"m of Van Genuchten [mm^-1] for each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.2f,",1-1/sl->pa->co[sy][jns][l]);
			}
		fprintf(f,"%.2f\n",1-1/sl->pa->co[sy][jns][Nl]);
		fprintf(f,"v of Van Genuchten [mm^-1] for each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.2f,",sl->pa->co[sy][jv][l]);
			}
		fprintf(f,"%.2f\n",sl->pa->co[sy][jv][Nl]);
		fprintf(f,"Kv_sat [mm/s] for each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.3f,",sl->pa->co[sy][jKv][l]);
			}
		fprintf(f,"%.3f\n",sl->pa->co[sy][jKv][Nl]);
		fprintf(f,"Kh_sat [mm/s] for each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.3f,",sl->pa->co[sy][jKh][l]);
			}
		fprintf(f,"%.3f\n",sl->pa->co[sy][jKh][Nl]);
		fprintf(f,"Thermal capacity of the soil skeleton [J/(m3 K)] for each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.0f,",sl->pa->co[sy][jct][l]);
			}
		fprintf(f,"%.0f\n",sl->pa->co[sy][jct][Nl]);
		fprintf(f,"Thermal conductivity of the soil skeleton [W/(m K)] for each layer\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%.2f,",sl->pa->co[sy][jkt][l]);
		}
		fprintf(f,"%.2f\n",sl->pa->co[sy][jkt][Nl]);
		fprintf(f,"Water content of wilting point [-]:\n");fprintf(f,"%.3f\n",land->ty->co[lu][jtwp]);
		fprintf(f,"Water content of field capacity [-]:\n");fprintf(f,"%.3f\n",land->ty->co[lu][jtfc]);
		fprintf(f,"Elevation above sea level [m]:\n");fprintf(f,"%f m\n",top->Z0->co[r][c]);
		fprintf(f,"Gauckler-Strickler [m^1/3/s]:\n");fprintf(f,"%f\n",land->ty->co[lu][jcm]);
		fprintf(f,"Sky view factor [-]:\n");fprintf(f,"%.3f\n",top->sky->co[r][c]);
		if (top->pixel_type->co[r][c]==0){
			fprintf(f,"The pixel-type is land\n0\n");
		}else if (top->pixel_type->co[r][c]==9){
			fprintf(f,"The pixel-type is novalue\n9 \n");
		}else if (top->pixel_type->co[r][c]==10){
			fprintf(f,"The pixel-type is channel\n10 \n");
		}else if (top->pixel_type->co[r][c]==11){
			fprintf(f,"The pixel-type is lake\n11\n");
		}else if (top->pixel_type->co[r][c]==12){
			fprintf(f,"The pixel-type is sea\n12 \n");
		}
		fprintf(f,"Drainage Direction:\n%d \n",top->DD->co[r][c]);
		fprintf(f,"Slope along Drainage Direction [-]:\n%f \n",top->i_DD->co[r][c]);
		fprintf(f,"Slope along positive x direction [-]:\n%f \n",top->dz_dx->co[r][c]);
		fprintf(f,"Slope along negative y direction [-]:\n%f \n",top->dz_dy->co[r][c]);
		fprintf(f,"Topology of curvature (0-1) [-]:\n%d \n",top->curv->co[r][c]);
		/*fprintf(f,"Area considering the slope [m^2]:\n%f \n",top->area->co[r][c]);*/
		fprintf(f,"Aspect [deg] [0=Nord, clockwise]:\n%f \n",top->aspect->co[r][c]*180.0/Pi);
		fprintf(f,"Mean slope of the pixel [deg]:\n%f \n",top->slopes->co[r][c]*180.0/Pi);
		fprintf(f,"Slope to calculate the surface velocity of the channel incoming flow [-]:\n%f \n",top->i_ch->co[r][c]);
		fprintf(f,"Land use number:\n%d\n",(short)land->LC->co[r][c]);

		root(land->ty->co[lu][jroot], sl->pa->co[sy][jdz], root_fraction->co);
		fprintf(f,"Root fraction [-]:\n");
		for(l=1;l<=Nl-1;l++){
			fprintf(f,"%f,",root_fraction->co[l]);
		}fprintf(f,"%f\n",root_fraction->co[Nl]);

		//fprintf(f," Albedo without snow and alpha=0 [-]: %f \n",land->ty->co[lu][jalbedo]);
		fprintf(f," Surface fraction of land covered by vegetation [-]: %f \n",land->ty->co[lu][jcf]);
		fprintf(f," Momentum roughness length z0soil [m]: %f \n",land->ty->co[lu][jz0]);
		fprintf(f," Vegetation height [m]: %f \n",land->ty->co[lu][jHveg]);
		fprintf(f," KRIGING WEIGHTS=\n");
		for(j=1;j<=n;j++){
			fprintf(f," STATION %ld = %f\n",j,wat->weights_Kriging->co[(r-1)*Nc+c][j]);
		}
		fprintf(f," */ \n");
		t_fclose(f);
		free(name);
		temp=join_strings(files->co[fpoint]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		fprintf(f,"DATE,JDfrom0,JD,t[d],t_i[s],t_f[s],v[m/s],Vdir,RH[-],P[hPa],Tair[C],Tsurface[C],Tdew[C],eair[mbar],Qair[-],esurf[mbar]");//16
		fprintf(f,",Qsurf[-],SWin[W/m2],SWin_beam,SWin_diff,SWout[W/m2],alpha[deg],direction[deg],phi[deg],LWin[W/m2],LWout[W/m2],Rnet[W/m2],SW[W/m2],LW[W/m2],H[W/m2],LE[W/m2]");//32
		fprintf(f,",Qrain[W/m2],Gsoil[W/m2],SurfaceEB[W/m2],SWin_c[MJ],SWout_c[MJ],LWin_c[MJ],LWout_c[MJ],SW_cum[MJ],LWn_cum[MJ],Rnet_cm[MJ],H_cum[MJ],LE_cum[MJ]");//44
		fprintf(f,",G_cum[MJ],Eg[mm],Sg[mm],Etc[mm],Psnow[mm],Prain[mm],Psnow_c[mm],Prain_SOILc[mm],Prain_SNOWc[mm],Ptot_c[mm],Wtrain[mm],Wtsnow[mm]");//56
		fprintf(f,",Ptot_atm,Rain_atm,Snow_atm,Ptot_atm_cum,Prain_atm_cum,Psnow_atm_cum,Pn[mm],Runoff[mm]");//69
		fprintf(f,",q_sup[mm],q_sub[mm],DS_sup[mm],DS_sub[mm],q_G[mm],snowDEPTH[mm],SWE[mm],snowDENSITY[kg/m3],snowT[C],BStot[mm],BStot_cum[mm]");//80
		fprintf(f,",snowMELT[mm],snowSUBL[mm],snowEVAP[mm],glacierDEPTH[mm],GWE[mm],gDENSITY[kg/m3],glcT[C],glcMELT[mm],glcSUBL[mm],glcEVAP[mm],qv(mm/d)");//91
		fprintf(f,",LWemissiv,LObukhov[m],numb.iter.,LWmin[W/m2],LWmax[W/m2],LAI,z0[m],d0[m],ErrorRichards[mm/h]");//100
		fprintf(f,",Ts,Tg,Tv,SWv,LWv,Hv,LEv,Htot,LEtot,Hg0,LEg0,Hg1,LEg1,fc");
		fprintf(f,",Ch,Cv,Cb,Cc,Ch_ic,Cv_ic,L_Obukhov,n_iter");
		fprintf(f,",Hv,LEv,Qv,Qg,Qa,Qs,u_top,decay,Locc\n");//121
		t_fclose(f);
		free(name);
		/*name=join_strings(files->co[fveg]+1,SSSS);
		name=join_strings(temp,textfile);
		f=t_fopen(name,"w");
		fprintf(f,"DATE,JDfrom0,JD,v[m/s],Tair,Tcan,Tg,Tsur,snowD,SWE,Lobukhov,Lobukhov_wc,a,LAI,fc,z0[m],d0[m],SWv,LWv,Hv,LEv,Ev[mm],Etrans,");
		fprintf(f,"Melt_over_canopy,Prain,Psnow,Wcrain,Wcrain_max,Wcsnow,Wcsnow_max,Hv,LEv,Hg_uc,LEg_ug,Hg_bg,LEg_bg,Hg,LEg,Qv,Qg,Qa,Qs,Ch,Cv,Cb,Cc,Ch_ic,Cv_ic\n");
		t_fclose(f);*/

		temp=join_strings(files->co[fsnz]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/**Snow profile for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */\n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,t_i[s],t_f[s],Snowtype,Snowdepth[mm],SWE[mm],Temp_aver[C],Density_aver[kg/m3],melting[mm],sublimation[mm],evaporation[mm],nlayer,BStrans[mm],BStrans_cum[mm],BSsubl[mm],BSsubl_cum[mm],BStot[mm],BStot_cum[mm]");
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",Dlayer[mm]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",Zlayer[mm]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",Temp[C]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",Density[kg/m3]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",SWE[mm]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",ice[kg/m2]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",water[kg/m2]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",CR_destr_met[1/h]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",CR_overburden[1/h]_%1ld",l);
		}
		for(l=1;l<=par->snowlayer_max;l++){
			fprintf(f,",CR_melting[1/h]_%1ld",l);
		}
		fprintf(f,"\n");
		t_fclose(f);
		free(name);
		if(par->glaclayer_max>0){
			temp=join_strings(files->co[fglz]+1,SSSS);
			name=join_strings(temp,textfile);
			free(temp);
			f=t_fopen(name,"w");
			//fprintf(f,"/** Ice profile for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: \n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
			fprintf(f,"DATE,JDfrom0,JD,t_i[s],t_f[s],Glacdepth[mm],GWE[mm],Temp_aver[C],Num.glac_layer,Rhoav[kg/m3],albedo,Gl.melt[mm],Gl.subl[mm],Gl.evap[mm]");
			for(l=1;l<=par->glaclayer_max;l++){
				fprintf(f,",Temp[C]_%1ld",l);
			}
			for(l=1;l<=par->glaclayer_max;l++){
				fprintf(f,",Depth[mm]_%1ld",l);
			}
			for(l=1;l<=par->glaclayer_max;l++){
				fprintf(f,",Rho[kg/m3]_%1ld",l);
			}
			for(l=1;l<=par->glaclayer_max;l++){
				fprintf(f,",GWE[mm]_%1ld",l);
			}
			for(l=1;l<=par->glaclayer_max;l++){
				fprintf(f,",th_ice_%1ld",l);
			}
			for(l=1;l<=par->glaclayer_max;l++){
				fprintf(f,",th_water_%1ld",l);
			}
			fprintf(f," \n");
			t_fclose(f);
			free(name);
		}

		/*creation of the file "Tz.txt": */
		temp=join_strings(files->co[fTz]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of sl temperature for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,time");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);
		free(name);
		/*creation of the file "Tz_MEAN.txt": */
		temp=join_strings(files->co[fTz_mean]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of soil temperature average for the time interval for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",parameters->checkpoints->element[i][1],parameters->checkpoints->element[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,time");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);
		free(name);
		/*creation of the file "PSIz.txt": */
		temp=join_strings(files->co[fpsiz]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of water pressure (mm) for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,time");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);
		free(name);
		/*creation of the file "PSI_MEANz.txt": */
		temp=join_strings(files->co[fpsiz_mean]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of water pressure (mm) average in the time interval for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,time");
		z=0.0;
		for(l=1;l<=Nl;l++){
		z+=sl->pa->co[sy][jdz][l];
		fprintf(f,",%.0f",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);
		free(name);
		/*creation of the file "TETAz.txt": */
		temp=join_strings(files->co[fliqz]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of water content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,time");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);
		free(name);
		/*creation of the file "TETAz_MEAN.txt": */
		temp=join_strings(files->co[fliqz_mean]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of water content average in the time interval for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,time");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);
		free(name);
		/*creation of the file "TETAICEz.txt": */
		temp=join_strings(files->co[ficez]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of ice content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld:\n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,time");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);
		free(name);
		/*creation of the file "TETAICEz_MEAN.txt": */
		temp=join_strings(files->co[ficez_mean]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of ice content average in the time interval for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,time");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);
		free(name);

	}

	//DATA BASIN
	temp=join_strings(files->co[fbas]+1,textfile);
	f=t_fopen(temp,"w");
	free(temp);
	fprintf(f,"DATE,JDfrom0,JD");
	fprintf(f,",t_i[s],t_f[s],Esoil[mm],Ecv[mm],Etc[mm],SW[W/m2],LW[W/m2],H[W/m2],LE[W/m2],SWv[W/m2],LWv[W/m2],Hv[W/m2],LEv[W/m2]");
	fprintf(f,",Ta[C],Tg[C],Tv[C],SWin[C],Prain[mm],Psnow[mm],Dwt[mm]");
	fprintf(f,",Pn[mm],Infiltration[mm],DS_sup[mm],DS_sub[mm],DS_ch[mm],R_G[mm],R_tot[mm],SSup[mm],SSub[mm],wt[mm]");
	fprintf(f,",SWE[mm],D_SWE[mm],MeltSnow[mm],D_MeltSnow[mm],SublSnow[mm],D_SublSnow[mm],EvapSnow[mm]");
	fprintf(f,",D_EvapSnow[mm],GWE[mm],D_GWE[mm],MeltGlac[mm],D_MeltGlac[mm],SublGlac[mm],D_SublGlac[mm],EvapGlac[mm]");
	fprintf(f,",D_EvapGlac[mm],ErrorRichards[mm/h]\n");
	t_fclose(f);


	//ALTIMETRIC RANKS
	for(i=1;i<=par->ES_num;i++){

		write_suffix(SSSS, i, 0);
		temp=join_strings(files->co[farank]+1,SSSS);
		name=join_strings(temp,textfile);
		free(temp);
		f=t_fopen(name,"w");
		fprintf(f,"/** File with calculated mean quantities for altimetric stripes.\n");
		fprintf(f," Minimum elevation: %10.3f m a.s.l.\n",top->Zmin+(i-1)*(top->Zmax-top->Zmin)/(double)par->ES_num);
		fprintf(f," Maximum elevation: %10.3f m a.s.l.\n",top->Zmin+i*(top->Zmax-top->Zmin)/(double)par->ES_num);
		fprintf(f," Note: t_i and t_f are the time from the start of simulation of the begin and of the end of the time interval in which the means are been calculated.\n");
		fprintf(f," Number of not NO_VALUE pixel %4ld  \n",top->ES_pixel->co[i]);
		fprintf(f," Area %10.4f km^2\n",top->ES_pixel->co[i]*UV->U->co[1]*UV->U->co[2]/1.0E6);
		fprintf(f," mean aspect (deg): %15.5f\n",top->ES_aspect->co[i]);
		fprintf(f," mean slope  (deg): %15.5f\n",top->ES_slope->co[i]);
		fprintf(f," t[d]         t_i[s]         t_f[s]         Eg[mm]         Sg[mm]         Ecv[mm]        Etc[mm]        ET[W/m2]       ");
		fprintf(f,"H[W/m2]        Rnet[W/m2]	  G[W/mq]        Ecanopy[W/m2]  Qrain[W/m2]    Ts[C]          Prain[mm]      Psnow[mm]		");
		fprintf(f,"Pn[mm]         Runoff[mm]     SWE[mm]        MeltSnow[mm]   SublSnow[mm]   EvapSnow[mm]   GWE[mm]        MeltGlac[mm]   SublGlac[mm]   EvapGlac[mm] */\n");
		t_fclose(f);
		free(name);
	}

	free_doublevector(root_fraction);


}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void write_soil_output(long n, long i, double t, double dt, long y0, double JD0, LONGMATRIX *rc, SOIL *sl, double psimin, double Esoil){
	/* Author:  Stefano Endrizzi  Year:
	* function writes the soil parameters of specified checkpoints pixels. It is called when the current time equals the time of output
	* Input:
	* 			n:	number of Dt after which the output of a pixel is printed (Es. se Dt=900sec e time output=1h, allora n=4)
	* 			i:	number of pixel (among the checkpoints) that is currently written (see chkpt->nrh)
	* 			t:	current time
	* 			dt: integration interval
	* 			y0: year of the beginning of the simulation
	* 			rc: matrix containing the number of row and column of the checkpoint i
	* 			sl: soil stucture
	* 		psimin:	Absolute minimum admitted suction potential
	* 		 Esoil:Soil comprimibility if oversaturated water is added
	* Output:
	* 		soil parameters are written
	* comment: Matteo Dall'Amico, May 2009 memory leakage corrections: Emanuele Cordano, September 2009 */
	char *name,*temp,SSSS[ ]={"SSSS"};
	double t_i, JD;
	long d2, mo2, y2, h2, mi2, l, r=rc->co[i][1], c=rc->co[i][2];
	FILE *f;
	short sy=sl->type->co[r][c];

	write_suffix(SSSS, i, 0);
	t_i=t-dt*n;
	//date_time(0.5*(t_i+dt+t+dt), y0, JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
	date_time(t+dt, y0, JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);

	/*update of the sl profile temperature at the end of the output interval in the control pixel:*/
	temp=join_strings(files->co[fTz]+1,SSSS);
	name=join_strings(temp,textfile);
	free(temp);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f,%.0f",JD+(double)(daysfrom0(y2)),JD,t+dt);
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",sl->T->co[l][r][c]);
	}
	fprintf(f," \n");
	fclose(f);
	free(name);
	/*writing psi at the end of the output interval: */
	temp=join_strings(files->co[fpsiz]+1,SSSS);
	name=join_strings(temp,textfile);
	free(temp);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f,%.0f",JD+(double)(daysfrom0(y2)),JD,t+dt);
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",sl->P->co[l][r][c]);
	}
	fprintf(f," \n");
	fclose(f);
	free(name);
	/*writing the water content at the end of the output interval: */
	temp=join_strings(files->co[fliqz]+1,SSSS);
	name=join_strings(temp,textfile);
	free(temp);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f,%.0f",JD+(double)(daysfrom0(y2)),JD,t+dt);
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
					sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],psimin,Esoil));
		//printf("alpha:%f\n",sl->pa->co[sy][ja][l]);
	}
	fprintf(f," \n");
	fclose(f);
	free(name);
	/*writing the ice content at the end of the output interval: */
	temp=join_strings(files->co[ficez]+1,SSSS);
	name=join_strings(temp,textfile);
	free(temp);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f,%.0f",JD+(double)(daysfrom0(y2)),JD,t+dt);
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",sl->thice->co[l][r][c]);
	}
	fprintf(f," \n");
	fclose(f);
	free(name);
	// find the date at the middle of the interval
	t_i=t-dt*n;
	date_time(0.5*(t_i+t)+dt, y0, JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);

	/*write the soil profile temperature as an AVERAGE of the output interval in the control pixel:*/
	temp=join_strings(files->co[fTz_mean]+1,SSSS);
	name=join_strings(temp,textfile);
	free(temp);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f,%.0f",JD+(double)(daysfrom0(y2)),JD,0.5*(t_i+t)+dt);
	for(l=1;l<=Nl;l++){
		if (t==0){
			fprintf(f,",%f",sl->T->co[l][r][c]);}
		else {
			fprintf(f,",%f",sl->Tmean->co[l][i]);
		}
	}
	fprintf(f," \n");
	fclose(f);
	free(name);
	/*write the  soil profile suction as an AVERAGE of the output interval in the control pixel:*/
	temp=join_strings(files->element[fpsiz_mean]+1,SSSS);
	name=join_strings(temp,textfile);
	free(temp);
	f=fopen(name,"a");
	//printf("\ntimes+dt=%f, t_i=%f, n=%ld, t_i+dt=%f, 0.5*(t_i+dt+t+dt)=%f",t+dt,t_i,n,t_i+dt,0.5*(t_i+dt+t+dt)); stop_execution();
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);fprintf(f,",%.0f", 0.5*(t_i+t)+dt);
	for(l=1;l<=Nl;l++){
		if(t==0){
			fprintf(f,",%f",sl->P->co[l][r][c]);
		}else{
			fprintf(f,",%f",sl->psi_mean->co[l][i]);
		}
	}
	fprintf(f," \n");
	fclose(f);
	free(name);
	/*write the  soil profile water content as an AVERAGE of the output interval in the control pixel:*/
	temp=join_strings(files->element[fliqz_mean]+1,SSSS);
	name=join_strings(temp,textfile);
	free(temp);
	f=fopen(name,"a");
	//printf("\ntimes+dt=%f, t_i=%f, n=%ld, t_i+dt=%f, 0.5*(t_i+dt+t+dt)=%f",t+dt,t_i,n,t_i+dt,0.5*(t_i+dt+t+dt)); stop_execution();
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);fprintf(f,",%.0f", 0.5*(t_i+t)+dt);
	for(l=1;l<=Nl;l++){
		if(t==0){
			fprintf(f,",%f",teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],psimin,Esoil));
		}else{
			fprintf(f,",%f",sl->thetaw_mean->co[l][i]);
		}
	}
	fprintf(f," \n");
	fclose(f);
	free(name);
	/*write the  soil profile ice content as an AVERAGE of the output interval in the control pixel:*/
	temp=join_strings(files->element[ficez_mean]+1,SSSS);
	name=join_strings(temp,textfile);
	free(temp);
	f=fopen(name,"a");
	//printf("\ntimes+dt=%f, t_i=%f, n=%ld, t_i+dt=%f, 0.5*(t_i+dt+t+dt)=%f",t+dt,t_i,n,t_i+dt,0.5*(t_i+dt+t+dt)); stop_execution();
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);fprintf(f,",%.0f", 0.5*(t_i+t)+dt);
	for(l=1;l<=Nl;l++){
		if(t==0){
			fprintf(f,",%f",sl->thice->co[l][r][c]);
		}else{
			fprintf(f,",%f",sl->thetai_mean->co[l][i]);
		}
	}
	fprintf(f," \n");
	fclose(f);
	free(name);


}
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void initmatrix(double val, DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue){

	long r,c;
	for(r=1;r<=destination->nrh;r++){
		for(c=1;c<=destination->nch;c++){
			if(origin->co[r][c]!=novalue) destination->co[r][c]=val;
		}
	}
}

void initlongmatrix(long val, LONGMATRIX *destination, DOUBLEMATRIX *origin, double novalue){

	long r,c;
	for(r=1;r<=destination->nrh;r++){
		for(c=1;c<=destination->nch;c++){
			if(origin->co[r][c]!=novalue) destination->co[r][c]=val;
		}
	}
}

void inittensor(double val, DOUBLETENSOR *destination, DOUBLEMATRIX *origin, double novalue){

	long l,r,c;
	for(r=1;r<=destination->nrh;r++){
		for(c=1;c<=destination->nch;c++){
			if(origin->co[r][c]!=novalue){
				for(l=1;l<=destination->ndh;l++){
					destination->co[l][r][c]=val;
				}
			}
		}
	}
}
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

double interp_value(double E, double N, DOUBLEMATRIX *M, DOUBLEMATRIX *Z){

	double  dN, dE, N0, E0, DN, DE, w1, V;
	long r, c;

	r=row1(N, Nr, 0, UV);
	c=col1(E, Nc, 0, UV);

	if(r==0 || c==0 || Z->co[r][c]==NoV){

		V=NoV;

	}else{

		dN=UV->U->co[1];
		dE=UV->U->co[2];

		N0=UV->U->co[3] + (Nr-r+0.5)*dN;
		E0=UV->U->co[4] + (c-0.5)*dE;

		DN=(N-N0)/(0.5*dN);
		DE=(E-E0)/(0.5*dE);

		//printf("No:%f E0:%f N:%f E:%f Dn:%f De:%f\n",N0,E0,N,E,DN,DE);

		if(DN<=DE && DN>=-DE){
			//printf("%ld %ld %f\n",r,c,M->co[r][c]);

			if(Z->co[r][c+1]!=NoV){
				w1=(E0+dE-E)/dE;
				V=w1*M->co[r][c]+(1.0-w1)*M->co[r][c+1];
			}else{
				V=M->co[r][c];
			}

		}else if(DN>=DE && DN>=-DE){
			//printf("2\n");
			//printf("%ld %ld %f\n",r,c,M->co[r][c]);
			if(Z->co[r-1][c]!=NoV){
				w1=(N0+dN-N)/dN;
				V=w1*M->co[r][c]+(1.0-w1)*M->co[r-1][c];
			}else{
				V=M->co[r][c];
			}

		}else if(DN>=DE && DN<=-DE){
			//printf("3\n");
			//printf("%ld %ld %f\n",r,c,M->co[r][c]);
			if(Z->co[r][c-1]!=NoV){
				w1=(E-(E0-dE))/dE;
				V=w1*M->co[r][c]+(1.0-w1)*M->co[r][c-1];
			}else{
				V=M->co[r][c];
			}

		}else{
			//printf("4\n");
			//printf("%ld %ld %f\n",r,c,M->co[r][c]);
			if(Z->co[r+1][c]!=NoV){
				w1=(N-(N0-dN))/dN;
				V=w1*M->co[r][c]+(1.0-w1)*M->co[r+1][c];
			}else{
				V=M->co[r][c];
			}

		}
	}

	return(V);

}

long row1(double N, long nrows, long i, T_INIT *UV)
{

	long cont=0;

	if(N>=UV->U->co[3] && N<=UV->U->co[3]+nrows*UV->U->co[1]){
		do{
			cont+=1;
		}while(UV->U->co[3]+(nrows-cont)*UV->U->co[1]>N);
	}

	return(cont);
}

long col1(double E, long ncols, long i, T_INIT *UV)
{

	long cont=0;

	if(E>=UV->U->co[4] && E<=UV->U->co[4]+ncols*UV->U->co[2]){
		do{
			cont+=1;
		}while(UV->U->co[4]+cont*UV->U->co[2]<E);
	}

	return(cont);
}


/****************************************************************************************************/
/* write_output_superfast: memorizza i dati degli output in un apposito tensore e li stampa su file alla fine della simulazione                                     */
/****************************************************************************************************/
void write_output_superfast(TIMES *times, WATER *wat, CHANNEL *cnet, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac, METEO *met)
/*
 * Author: Matteo Dall'Amico, on 06/10/2009 in Trento
 * */
{
 /*internal auxiliary variables:*/
 long i,r=0,c=0,l; /*counters*/
 char *name=NULL; /*modified by Emanuele Cordano on 24/9/9 */
 FILE *f;

 /*internal variables to memorize input par:*/
 double time_max;      /*time of all the simulation [s]*/
 double total_pixel;   /*total number of pixel which are not novalue*/
 short sy=0; /* modified by Emanuele Cordano on 24 September 2009 */

 /* Assignment to some internal variables of some input par:*/
time_max=times->TH*3600.0;
total_pixel=(double)par->total_pixel;

/* Print on the screen and into error-file the temporal coordinate of simulation: */
if (times->i_pixel==times->n_pixel){
	if(times->mm<=9){

		printf("%ld/%ld/%ld %ld:0%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));

		f=fopen(error_file_name,"a");
		fprintf(f,"%ld/%ld/%ld %ld:0%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));
		fclose(f);

	}else{

		printf("%ld/%ld/%ld %ld:%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));

		f=fopen(error_file_name,"a");
		fprintf(f,"%ld/%ld/%ld %ld:%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));
		fclose(f);

	}
}

//****************************************************************************************************************
//DATA POINT
//****************************************************************************************************************

if(par->state_pixel==1){

	for(i=1;i<=par->chkpt->nrh;i++){

		//write_suffix(SSSS, i, 0);//Mat
		r=par->rc->co[i][1];
		c=par->rc->co[i][2];
		sy=sl->type->co[r][c];

		// update of min, max and mean punctual values
		for(l=1;l<=Nl;l++){
			sl->Tmean->co[l][i]+=sl->T->co[l][r][c]/(double)times->n_pixel;
			//printf("time=%f,l=%ld,sl->T=%f,sl->Tmean=%f,times->i_pixel=%ld,time->n_pixel=%ld",times->time,l,sl->T->co[l][r][c],sl->Tmean->co[l],times->i_pixel,times->n_pixel);stop_execution();
			sl->thetai_mean->co[l][i]+=sl->thice->co[l][r][c]/(double)times->n_pixel;
			sl->thetaw_mean->co[l][i]+=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
					sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
					PSImin,par->Esoil)/(double)times->n_pixel;

			sl->Tmin->co[l][i]=Fmin(sl->T->co[l][r][c], sl->Tmin->co[l][i]);
			sl->thetai_min->co[l][i]=Fmin(sl->thice->co[l][r][c], sl->thetai_min->co[l][i]);
			sl->thetaw_min->co[l][i]=Fmin(teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
					sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
					PSImin,par->Esoil),sl->thetaw_min->co[l][i]);

			sl->Tmax->co[l][i]=Fmax(sl->T->co[l][r][c], sl->Tmax->co[l][i]);
			sl->thetai_max->co[l][i]=Fmax(sl->thice->co[l][r][c], sl->thetai_max->co[l][i]);
			sl->thetaw_max->co[l][i]=Fmax(teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
					sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
					PSImin,par->Esoil),sl->thetaw_max->co[l][i]);
		}

		/*Print of pixel-output every times->n_pixel time step */
		if (times->i_pixel==times->n_pixel){
			/* matrix 9xNl that contains the value of the output variable for each layer
			 * row0: Tmin
			 * row1: Tmax
			 * row2: Tmean
			 * row3: ThetaWmin
			 * row4: ThetaWmax
			 * row5: ThetaWmean
			 * row6: ThetaImin
			 * row7: ThetaImax
			 * row8: ThetaImean
			*/
			for(l=1;l<=Nl;l++){
				/* record the variables in the supertensor */
				sl->output[times->count_super_printed+1][0][l-1]=sl->Tmin->co[l][i];
				sl->output[times->count_super_printed+1][1][l-1]=sl->Tmax->co[l][i];
				sl->output[times->count_super_printed+1][2][l-1]=sl->Tmean->co[l][i];
				sl->output[times->count_super_printed+1][3][l-1]=sl->thetaw_min->co[l][i];
				sl->output[times->count_super_printed+1][4][l-1]=sl->thetaw_max->co[l][i];
				sl->output[times->count_super_printed+1][5][l-1]=sl->thetaw_mean->co[l][i];
				sl->output[times->count_super_printed+1][6][l-1]=sl->thetai_min->co[l][i];
				sl->output[times->count_super_printed+1][7][l-1]=sl->thetai_max->co[l][i];
				sl->output[times->count_super_printed+1][8][l-1]=sl->thetai_mean->co[l][i];
				/* update the variables */
				sl->thetaw_mean->co[l][i]=0.0;
				sl->Tmean->co[l][i]=0.0;
				sl->thetai_mean->co[l][i]=0.0;
				sl->Tmin->co[l][i]=99.0;
				sl->thetaw_min->co[l][i]=99.0;
				sl->thetai_min->co[l][i]=99.0;
				sl->Tmax->co[l][i]=-99.0;
				sl->thetaw_max->co[l][i]=-99.0;
				sl->thetai_max->co[l][i]=-99.0;
				}
			}//end(times->i_pixel==times->n_pixel)
		}//end for(i=1;i<=par->chkpt->nrh;i++)

		// now I have to print the output "just" at the end of the simulation
		if(times->time>=time_max){
			double Dt_output=times->n_pixel*par->Dt; // [sec]
			name=join_strings(files->co[supfs]+1,textfile);
			write_supertensor(r,c,sl->type->co[r][c], Dt_output, par->num_of_time, times->time, par->Dt, par->year0, par->JD0, name, sl->pa->co[sy][jdz], times, sl->output);
		}
		free(name);
	}// end if(par->state_pixel==1)
}// end function


/****************************************************************************************************/
/* write_supertensor                    */
/****************************************************************************************************/
void write_supertensor(long r, long c,short sy,double Dt_output, short n,double t, double dt, long y0, double JD0, char* file, double* Dz, TIMES *times, double *** out){
	/*Author: Matteo Dall'Amico, 6 october 2009 in Caltrano
	 * writes the super_output_tensor in a file
	 * r: row
	 * c: column
	 * sy: soil type
	 * Dt_output: print output
	 * n: number of times the output must be written (num of times)
	 * t: time
	 * dt: integration interval
	 * y0: year of beginning
	 * JD0: julian day of beginning
	 * file: file name
	 * Dz: vector of layer depth
	 * times: TIMES structure
	 * tensor of ouput results */

	//char* name=NULL;
	double z;
	int i,j,l;
	double JD;
	long d2, mo2, y2, h2, mi2;
	FILE *f;
	//double Dt_output=n*dt;

	f=fopen(file,"w");
	// write the header
	fprintf(f,"DATE,time");
	// write Tmin
	for(l=1;l<=Nl;l++){
		fprintf(f,",Tmin ");
	}
	// write Tmax
	for(l=1;l<=Nl;l++){
		fprintf(f,",Tmax ");
	}
	// write Tmean
	for(l=1;l<=Nl;l++){
		fprintf(f,",Tmean ");
	}
	// write Theta_w min
	for(l=1;l<=Nl;l++){
		fprintf(f,",ThetaW min ");
	}
	// write Theta_w max
	for(l=1;l<=Nl;l++){
		fprintf(f,",ThetaW max ");
	}
	// write Theta_w mean
	for(l=1;l<=Nl;l++){
		fprintf(f,",ThetaW mean ");
	}
	// write Theta_I min
	for(l=1;l<=Nl;l++){
		fprintf(f,",ThetaI min ");
	}
	// write Theta_I max
	for(l=1;l<=Nl;l++){
		fprintf(f,",ThetaI max ");
	}
	// write Theta_I mean
	for(l=1;l<=Nl;l++){
		fprintf(f,",ThetaI mean ");
		}
	fprintf(f," \n");
	// write the header
	fprintf(f,"DATE,time");
	for(j=1;j<=9;j++){
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=Dz[l];
			fprintf(f,",%.0f ",z-0.5*Dz[l]);
		}
	}
	fprintf(f," \n");
	// ORA DEVI SCRIVERE LA Tmin, Tmax, Tmean... ecc.
	// Tmin
	for(i=0;i<n;i++){// for every superfast output time

		date_time(i*Dt_output, y0, JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
		write_date(f, d2, mo2, y2, h2, mi2);
		fprintf(f,",%.0f",i*Dt_output);
		for(j=1;j<=9;j++){// for every variable
			for(l=1;l<=Nl;l++){// for every layer
				fprintf(f,",%f",out[i][j-1][l-1]);
				}
			}
		fprintf(f," \n");
		}
	fclose(f);
}
