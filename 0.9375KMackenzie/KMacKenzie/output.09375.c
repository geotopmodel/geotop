
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


#include "keywords_file.h"
#include "struct.geotop.09375.h"
#include "liston.h"
#include "output.09375.h"
#include "pedo.funct.h"
#include "write_dem.h"
#include "t_datamanipulation.h"
#include "geo_statistic.09375.h"
#include "t_random.h"
#include "networks.h"
#include "t_utilities.h"
#include "rw_maps.h"
#include "constant.h"
#include "extensions.h"
#include "times.h"
#include "energy.balance.h"
#include "input.09375.h"



extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;


/****************************************************************************************************/
/* write_output: stampa i dati degli output variabili nel tempo                                     */
/****************************************************************************************************/
void write_output(TIMES *times, WATER *wat, CHANNEL *cnet, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac, METEO *met)

{
 /*internal auxiliary variables:*/
 long i,j,r,c,l,s; /*counters*/
 long n_file;      /*number of file of the type "TETAxySSSlZZ"(i.e. number of the basin-time-step)*/
 double t_i;         /*time of begin of an interval time for the output*/
 char SSSS[ ]={"SSSS"};
 char *name;
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
 short sy;
 double fwet;


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

		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"%ld/%ld/%ld %ld:0%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));
		fclose(f);

	}else{

		printf("%ld/%ld/%ld %ld:%ld JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
	       times->DD,times->MM,times->AAAA,times->hh,times->mm,times->JD,(long)(floor(times->time/86400))+1,
	       (100.0*(double)times->time/time_max));

		f=fopen(files->co[ferr]+1,"a");
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
	f=fopen(join_strings(files->co[fQ]+1,textfile),"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
	fprintf(f,",%f,%f,%f,%f\n",(Qsub_ch+Qsup_ch+Q_G)/((double)times->n_pixel),Qsub_ch/((double)times->n_pixel),Qsup_ch/((double)times->n_pixel),Q_G/((double)times->n_pixel));
	fclose(f);

	Qsub_ch=0.0; Qsup_ch=0.0; Q_G=0.0;

	//melt fluxes
	if(par->output_balancesn!=0){

		for(i=1;i<=land->clax->nh;i++){

			write_suffix(SSSS, i, 0);

			name=join_strings(files->co[fmeltlu]+1,"_snowcovered_");
			name=join_strings(name,SSSS);
			name=join_strings(name,textfile);
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
					fprintf(f,"%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
						0.0,0.0,0.0,0,0.0,0.0,0.0);
				}
			}else{
				fprintf(f,"%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
					0.0,0.0,0.0,0,0.0,0.0,0.0);
			}
			fclose(f);

			name=join_strings(files->co[fmeltlu]+1,"_snowfree_");
			name=join_strings(name,SSSS);
			name=join_strings(name,textfile);
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
					fprintf(f,"%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
						0.0,0.0,0.0,0,0.0,0.0,0.0);
				}
			}else{
				fprintf(f,"%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f\n",times->time+par->Dt*(1-times->n_pixel),times->time+par->Dt,(0.5*(times->time+par->Dt*(1-times->n_pixel))+0.5*(times->time+par->Dt))/(double)86400,
					0.0,0.0,0.0,0,0.0,0.0,0.0);
			}
			fclose(f);

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

		/*Print of pixel-output every times->n_pixel time step */
		if (times->i_pixel==times->n_pixel){

			wat->out1->co[11][i]=wat->h_sup->co[r][c];
			for(l=1;l<=Nl;l++){
				wat->out1->co[12][i]+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
					sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
					fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);
			}

			//wat->wt [mm] is the water (in both the liquid and solid form) that is currently on the leaves
			wat->out1->co[13][i]=-wat->out1->co[14][i];
			wat->out1->co[14][i]=wat->wt->co[r][c];
			wat->out1->co[13][i]+=wat->wt->co[r][c];

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
			name=join_strings(files->co[fpoint]+1,SSSS);
			name=join_strings(name,textfile);
			f=fopen(name,"a");
			date_time(times->time+par->Dt, par->year0, par->JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
			write_date(f, d2, mo2, y2, h2, mi2);
			fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
			fprintf(f,",%f,%f,%f,",0.5*(t_i+par->Dt+times->time+par->Dt)/(double)86400,(t_i+par->Dt),(times->time+par->Dt));
			fprintf(f,"%14.3f,%14.9f,%f,%14.4f,%14.6f,%14.6f,%14.6f,%14.6f,%14.12f,%14.6f,%14.12f,",
					egy->out1->co[17][i]/*v*/,egy->out1->co[55][i]/*vdir*/,egy->out1->co[18][i]/*RH*/,egy->out1->co[19][i]/*P*/,egy->out1->co[20][i]/*Ta*/,egy->out1->co[64][i]/*Tsurface*/,
					egy->out1->co[46][i],egy->out1->co[47][i],egy->out1->co[48][i],egy->out1->co[49][i],egy->out1->co[50][i]);
			fprintf(f,"%14.5f,%14.5f,%14.5f,%14.5f,%14.12f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,%14.5f,",
					egy->out1->co[14][i]/*Rsw*/,egy->out1->co[53][i]/*SWbeam*/,egy->out1->co[54][i]/*SWdiff*/,egy->out1->co[45][i],egy->out1->co[15][i]/*albedo*/,
					(180/Pi)*egy->hsun,(180/Pi)*egy->dsun,(180/Pi)*acos(cos(top->slopes->co[r][c])*sin(egy->hsun)+sin(top->slopes->co[r][c])*cos(egy->hsun)*cos(-top->aspect->co[r][c]+egy->dsun)),
					egy->out1->co[12][i]/*Rlwin*/,egy->out1->co[13][i]/*Rlwout*/,egy->out1->co[11][i]/*Rnet*/,egy->out1->co[14][i]+egy->out1->co[45][i]/*SW*/,
					egy->out1->co[12][i]+egy->out1->co[13][i]/*LW*/,egy->out1->co[6][i]/*H*/,egy->out1->co[5][i]/*L*/,egy->out1->co[9][i]/*Qrain*/,egy->out1->co[8][i]/*G*/,egy->out1->co[7][i]/*surfEB*/);
			fprintf(f,"%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,",
					egy->out1->co[24][i]/*Rswin_cum*/,egy->out1->co[28][i]/*Rswout_cum*/,egy->out1->co[22][i]/*Rlwin_cum*/,egy->out1->co[23][i]/*Rlwout_cum*/,
					egy->out1->co[24][i]+egy->out1->co[28][i],egy->out1->co[22][i]+egy->out1->co[23][i],egy->out1->co[21][i]/*Rnet_cum*/,
					egy->out1->co[26][i]/*H_cum*/,egy->out1->co[25][i]/*ET_cum*/,egy->out1->co[27][i]/*G_cum*/);
			fprintf(f,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
					 egy->out1->co[1][i]/*Eg*/,egy->out1->co[2][i]/*Sg*/,egy->out1->co[4][i]/*Etc*/,
					 wat->out1->co[3][i]/*Psnow*/,wat->out1->co[4][i]/*Prain*/,wat->out1->co[21][i]/*Psnow_cum*/,wat->out1->co[22][i]/*Prain_cum_soil*/,wat->out1->co[27][i]/*rain_on_snow_cum*/,
					 wat->out1->co[21][i]+wat->out1->co[22][i],wat->out1->co[14][i]/*wt*/,wat->out1->co[18][i]/*maxstorage*/,wat->out1->co[13][i]/*Delta Water intercepted by leaves*/,
					 wat->out1->co[16][i]+wat->out1->co[17][i]/*Ptot atm*/,wat->out1->co[17][i]/*rain*/,wat->out1->co[16][i]/*snow*/,wat->out1->co[19][i]/*evap*/,
					 wat->out1->co[20][i]/*drip*/,wat->out1->co[23][i]+wat->out1->co[24][i],wat->out1->co[24][i],wat->out1->co[23][i],wat->out1->co[25][i],wat->out1->co[26][i],
					 wat->out1->co[5][i]/*net liquid precipitation*/,wat->out1->co[6][i]/*runoff*/,
					 wat->out1->co[7][i]/*q_sup*/,wat->out1->co[10][i]/*q_sub*/,wat->out1->co[11][i]/*wat stored on sl*/,wat->out1->co[12][i]/*wat stored in subsoil*/,
					 wat->out1->co[8][i]/*Q_g*/);
			fprintf(f,"%14.8f,%14.8f,%14.8f,%14.8f,%f,%f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,",
					 snowD,SWE,snowdensity,snowtemperature,snow->out_bs->co[9][i],snow->out_bs->co[10][i],snow->melted->co[i],snow->subl->co[i],snow->evap->co[i],
					 glacD,GWE,glacdensity,glactemperature,glac->melted->co[i],glac->subl->co[i],glac->evap->co[i]);
			fprintf(f,"%14.9f,",wat->out1->co[15][i]/*q_v*/);
			/*fprintf(f,"%14.12f,%14.3f,%14.6f,%14.11f,%14.11f,%14.9f,%14.11f,%14.11f,%14.11f,%14.11f,%14.11f,%14.11f,%14.11f,",
					 egy->out1->co[16][i],egy->out1->co[29][i],egy->out1->co[30][i],egy->out1->co[35][i],egy->out1->co[36][i],egy->out1->co[37][i],egy->out1->co[38][i],
					 egy->out1->co[39][i],egy->out1->co[40][i],egy->out1->co[41][i],egy->out1->co[42][i],egy->out1->co[43][i],egy->out1->co[44][i]);*/
			fprintf(f,"%14.12f,%14.3f,%14.6f,", egy->out1->co[16][i],egy->out1->co[29][i],egy->out1->co[30][i]);
			fprintf(f,"%f,%f,%f,%f,%f,%20.18f,",egy->out1->co[51][i],egy->out1->co[52][i],egy->out1->co[56][i],egy->out1->co[57][i],egy->out1->co[58][i],wat->out1->co[28][i]/(double)times->n_pixel);
			fprintf(f,"%f,%f,%f,%f,%f,%f,%f,%f,",egy->out1->co[64][i],egy->out1->co[10][i],egy->out1->co[63][i],egy->out1->co[59][i],egy->out1->co[60][i],egy->out1->co[61][i],egy->out1->co[62][i],fwet);
			fprintf(f,"%f,%f,%f,%f,%f,%f,%f,",egy->out1->co[6][i]+egy->out1->co[61][i],egy->out1->co[5][i]+egy->out1->co[62][i],egy->out1->co[65][i],egy->out1->co[66][i],egy->out1->co[67][i],egy->out1->co[68][i],egy->out1->co[69][i]);
			fprintf(f,"%f,%f,%f,%f,%f,%f,%f,%f,%e,%e,%e,%e\n",egy->out1->co[70][i],egy->out1->co[71][i],egy->out1->co[72][i],egy->out1->co[73][i],egy->out1->co[74][i],egy->out1->co[77][i],egy->out1->co[75][i],egy->out1->co[76][i],egy->out1->co[78][i],egy->out1->co[79][i],egy->out1->co[80][i],egy->out1->co[81][i]);
			fclose(f);

			/*update of the file "f" every times->n_pixel times: */
			z=0.0;
			name=join_strings(files->co[fsnz]+1,SSSS);
			name=join_strings(name,textfile);
			f=fopen(name,"a");
			write_date(f, d2, mo2, y2, h2, mi2);
			fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
			fprintf(f,",%f,%f,%ld,%f,%f,%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",(t_i+par->Dt),(times->time+par->Dt),snow->type->co[r][c],snowD,SWE,snowtemperature,snowdensity,egy->out1->co[15][i]/*albedo*/,
				snow->melted->co[i],snow->subl->co[i],snow->evap->co[i],snow->lnum->co[r][c],snow->out_bs->co[1][i],snow->out_bs->co[2][i],
				snow->out_bs->co[3][i],snow->out_bs->co[4][i],snow->out_bs->co[5][i],snow->out_bs->co[6][i],snow->out_bs->co[7][i],
				snow->out_bs->co[8][i],snow->out_bs->co[9][i],snow->out_bs->co[10][i]);
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
				name=join_strings(files->co[fglz]+1,SSSS);
				name=join_strings(name,textfile);
				f=fopen(name,"a");
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
			write_soil_output(times->n_pixel, i, times->time, par->Dt, par->year0, par->JD0, par->rc, sl, par->psimin, par->Esoil);

		}//end(times->i_pixel==times->n_pixel)
	}//end(i point cycle)

	if(times->i_pixel==times->n_pixel){
		for(i=1;i<=par->chkpt->nrh;i++){
			for(j=1;j<=20;j++) { egy->out1->co[j][i]=0.0; }
			for(j=29;j<=30;j++) { egy->out1->co[j][i]=0.0; }
			for(j=35;j<=81;j++) { egy->out1->co[j][i]=0.0; }
			for(j=1;j<=12;j++) { wat->out1->co[j][i]=0.0; }
			for(j=15;j<=20;j++) { wat->out1->co[j][i]=0.0; }
			for(j=28;j<=28;j++) { wat->out1->co[j][i]=0.0; }
			for(j=1;j<=5;j++) { snow->out_bs->co[2*j-1][i]=0.0; }

			for(l=1;l<=Nl;l++){
				wat->out1->co[2][i]+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
								sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
								fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);
			}
			wat->out1->co[1][i]=wat->h_sup->co[r][c];

		}
	}
}


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
						sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil)*sl->pa->co[sy][jdz][l];
				}

				/*water on the leaves*/
				wt0_basin+=wat->wt->co[r][c];

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
				wat->out2->co[5]+=wat->wt->co[r][c];
				wt0_basin+=wat->wt->co[r][c];

				/*Total water storage on the sl-surface (volume for unit of area[mm]):*/
				wat->out2->co[6]+=wat->h_sup->co[r][c];
				Ssup+=wat->h_sup->co[r][c];

				/*Total water storage in the subsoil (volume for unit of area[mm]):*/
				for(l=1;l<=Nl;l++){
					wat->out2->co[7]+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
								sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
								fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);
					Ssub+=sl->pa->co[sy][jdz][l]*teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],
								sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],
								fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);
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
	f=fopen(join_strings(files->co[fbas]+1,textfile),"a");
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
	printf(" SW=%f W/m2  LW:%f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Prain=%6.2f mm  Psnow=%6.2f mm  Rout=%6.2f mm/h \n Max Error Richards=%20.18f mm/h\n\n",
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
	name=join_strings(files->co[farank]+1,SSSS);
	name=join_strings(name,textfile);

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
if(par->output_TETAxy>0){
	Q=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(Q,UV->V->co[2]);

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				sy=sl->type->co[r][c];
				for(l=1;l<=Nl;l++){
					Q->co[l][r][c]=teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],
									sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l],Psif(sl->T->co[l][r][c],par->psimin)),
									par->Esoil);
					sl->thwav->co[l][r][c]+=(Q->co[l][r][c]/((par->output_TETAxy*3600.0)/(par->Dt)));
				}
			}
		}
	}

	if( fmod(times->time+par->Dt,par->output_TETAxy*3600.0)==0 ){
		n_file=(long)((times->time+par->Dt)/(par->output_TETAxy*3600.0));
		write_tensorseries2(n_file, join_strings(files->co[fliq]+1, "_IST_"), 0, par->format_out, Q, UV);
		write_tensorseries2(n_file, join_strings(files->co[fliq]+1, "_MEAN_"), 0, par->format_out, sl->thwav, UV);

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					for(l=1;l<=Nl;l++){
						sl->thwav->co[l][r][c]=0.0;
					}
				}
			}
		}
	}

	free_doubletensor(Q);
}



//T
if(par->output_Txy>0){

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				sy=sl->type->co[r][c];
				for(l=1;l<=Nl;l++){
					sl->Tav->co[l][r][c]+=(sl->T->co[l][r][c]/((par->output_Txy*3600.0)/(par->Dt)));
				}
			}
		}
	}

	if( fmod(times->time+par->Dt,par->output_Txy*3600.0)==0 ){
		n_file=(long)((times->time+par->Dt)/(par->output_Txy*3600.0));
		write_tensorseries2(n_file, join_strings(files->co[fT]+1, "IST"), 0, par->format_out, sl->T, UV);
		write_tensorseries2(n_file, join_strings(files->co[fT]+1, "MEAN"), 0, par->format_out, sl->Tav, UV);

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					for(l=1;l<=Nl;l++){
						sl->Tav->co[l][r][c]=0.0;
					}
				}
			}
		}
	}
}

//TETAICE
if(par->output_TETAICExy>0){

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				sy=sl->type->co[r][c];
				for(l=1;l<=Nl;l++){
					sl->thiav->co[l][r][c]+=(sl->thice->co[l][r][c]/((par->output_TETAICExy*3600.0)/(par->Dt)));
				}
			}
		}
	}

	if( fmod(times->time+par->Dt,par->output_TETAICExy*3600.0)==0 ){
		n_file=(long)((times->time+par->Dt)/(par->output_TETAICExy*3600.0));
		write_tensorseries2(n_file, join_strings(files->co[fice]+1, "IST"), 0, par->format_out, sl->thice, UV);
		write_tensorseries2(n_file, join_strings(files->co[fice]+1, "MEAN"), 0, par->format_out, sl->thiav, UV);

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					for(l=1;l<=Nl;l++){
						sl->thiav->co[l][r][c]=0.0;
					}
				}
			}
		}
	}
}

//PSI
if(par->output_PSIxy>0 && fmod(times->time+par->Dt,par->output_PSIxy*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_PSIxy*3600.0));
	write_tensorseries2(n_file, files->co[fpsi]+1, 0, par->format_out, sl->P, UV);
}

//matrix to handle data
M=new_doublematrix(Nr,Nc);
initialize_doublematrix(M,UV->V->co[2]);

//ALBEDO
if(par->output_albedo>0 && fmod(times->time+par->Dt,par->output_albedo*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_albedo*3600.0));
	write_suffix(SSSS, n_file, 0);
	write_map(join_strings(files->co[falb]+1,SSSS), 0, par->format_out, land->albedo, UV);
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
	write_map(join_strings(join_strings(files->co[fsn]+1,"Dist"),SSSS), 0, par->format_out, M, UV);

	//write_map(join_strings(join_strings(files->co[fsn]+1,"Dmax"),SSSS), 0, par->format_out, snow->max, UV);
	//initmatrix(0.0, snow->max, top->Z0, NoV);
	//write_map(join_strings(join_strings(files->co[fsn]+1,"Daverage"),SSSS), 0, par->format_out, snow->average, UV);
	//initmatrix(0.0, snow->average, top->Z0, NoV);
	if(par->snowtrans==1){
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsubl"),SSSS), 0, par->format_out, snow->Wsubl_cum, UV);
		//initmatrix(0.0, snow->Wsubl_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsalt"),SSSS), 0, par->format_out, snow->Wsalt_cum, UV);
		//initmatrix(0.0, snow->Wsalt_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsusp"),SSSS), 0, par->format_out, snow->Wsusp_cum, UV);
		//initmatrix(0.0, snow->Wsusp_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsbgr"),SSSS), 0, par->format_out, snow->Wsubgrid_cum, UV);
		//initmatrix(0.0, snow->Wsubgrid_cum, top->Z0, NoV);
		write_map(join_strings(join_strings(files->co[fsn]+1,"BStot"),SSSS), 0, par->format_out, snow->Wtot_cum, UV);
		initmatrix(0.0, snow->Wtot_cum, top->Z0, NoV);
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
	write_map(join_strings(files->co[fgl]+1,SSSS), 0, par->format_out, M, UV);
}

//WATER OVER THE SURFACE
if(par->output_h_sup>0 && fmod(times->time+par->Dt,par->output_h_sup*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_h_sup*3600.0));
	write_suffix(SSSS, n_file, 0);
	write_map(join_strings(files->co[fhsup]+1,SSSS), 0, par->format_out, wat->hsupav, UV);
	initmatrix(0.0, wat->hsupav, top->Z0, NoV);
}

//RADIATION
if(par->output_Rn>0 && fmod(times->time+par->Dt,par->output_Rn*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Rn*3600.0));
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fRn]+1,"nmean");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->Rn_mean, UV);
	initmatrix(0.0, egy->Rn_mean, top->Z0, NoV);

	name=join_strings(files->co[fRn]+1,"LWin");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->LW_in, UV);
	initmatrix(0.0, egy->LW_in, top->Z0, NoV);

	name=join_strings(files->co[fRn]+1,"LWout");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->LW_out, UV);
	initmatrix(0.0, egy->LW_out, top->Z0, NoV);

	name=join_strings(files->co[fRn]+1,"SW");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->SW, UV);
	initmatrix(0.0, egy->SW, top->Z0, NoV);

	if(par->distr_stat==1){
		name=join_strings(files->co[fRn]+1,"nmax");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->Rn_max, UV);
		initmatrix(-1.0E+9, egy->Rn_max, top->Z0, NoV);

		name=join_strings(files->co[fRn]+1,"nmin");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->Rn_min, UV);
		initmatrix(1.0E+9, egy->Rn_min, top->Z0, NoV);

		name=join_strings(files->co[fRn]+1,"LWmax");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->LW_max, UV);
		initmatrix(-1.0E+9, egy->LW_max, top->Z0, NoV);

		name=join_strings(files->co[fRn]+1,"LWmin");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->LW_min, UV);
		initmatrix(1.0E+9, egy->LW_min, top->Z0, NoV);

		name=join_strings(files->co[fRn]+1,"SWmax");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->SW_max, UV);
		initmatrix(-1.0E+9, egy->SW_max, top->Z0, NoV);
	}
}

//GROUND HEAT FLUX
if(par->output_G>0 && fmod(times->time+par->Dt,par->output_G*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_G*3600.0));
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fG]+1,"mean");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->G_mean, UV);
	initmatrix(0.0, egy->G_mean, top->Z0, NoV);

	name=join_strings(files->co[fG]+1,"snowsoil");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->G_snowsoil, UV);
	initmatrix(0.0, egy->G_snowsoil, top->Z0, NoV);

	if(par->distr_stat==1){
		name=join_strings(files->co[fG]+1,"max");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->G_max, UV);
		initmatrix(-1.0E+9, egy->G_max, top->Z0, NoV);

		name=join_strings(files->co[fG]+1,"min");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->G_min, UV);
		initmatrix(1.0E+9, egy->G_min, top->Z0, NoV);
	}
}

//SENSIBLE HEAT FLUX
if(par->output_H>0 && fmod(times->time+par->Dt,par->output_H*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_H*3600.0));
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fH]+1,"mean");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->H_mean, UV);
	initmatrix(0.0, egy->H_mean, top->Z0, NoV);

	if(par->distr_stat==1){
		name=join_strings(files->co[fH]+1,"max");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->H_max, UV);
		initmatrix(-1.0E+9, egy->H_max, top->Z0, NoV);

		name=join_strings(files->co[fH]+1,"min");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->H_min, UV);
		initmatrix(1.0E+9, egy->H_min, top->Z0, NoV);
	}
}

//LATENT HEAT FLUX
if(par->output_ET>0 && fmod(times->time+par->Dt,par->output_ET*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_ET*3600.0));
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fLE]+1,"mean");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->ET_mean, UV);
	initmatrix(0.0, egy->ET_mean, top->Z0, NoV);

	if(par->distr_stat==1){
		name=join_strings(files->co[fLE]+1,"max");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->ET_max, UV);
		initmatrix(-1.0E+9, egy->ET_max, top->Z0, NoV);

		name=join_strings(files->co[fLE]+1,"min");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->ET_min, UV);
		initmatrix(1.0E+9, egy->ET_min, top->Z0, NoV);
	}
}

//SURFACE TEMPERATURE
if(par->output_Ts>0 && fmod(times->time+par->Dt,par->output_Ts*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Ts*3600.0));
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fTs]+1,"mean");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->Ts_mean, UV);
	initmatrix(0.0, egy->Ts_mean, top->Z0, NoV);

	if(par->distr_stat==1){
		name=join_strings(files->co[fTs]+1,"max");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->Ts_max, UV);
		initmatrix(-1.0E+9, egy->Ts_max, top->Z0, NoV);

		name=join_strings(files->co[fTs]+1,"min");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->Ts_min, UV);
		initmatrix(1.0E+9, egy->Ts_min, top->Z0, NoV);
	}
}

//PRECIPITATION
if(par->output_P>0 && fmod(times->time+par->Dt,par->output_P*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_P*3600.0));
	write_suffix(SSSS, n_file, 0);

	write_map(join_strings(join_strings(files->co[fprec]+1,"TOTAL"),SSSS), 0, par->format_out, wat->PrTOT_mean, UV);
	initmatrix(0.0, wat->PrTOT_mean, top->Z0, NoV);

	write_map(join_strings(join_strings(files->co[fprec]+1,"SNOW"),SSSS), 0, par->format_out, wat->PrSNW_mean, UV);
	initmatrix(0.0, wat->PrSNW_mean, top->Z0, NoV);
}

//INTERCEPTED PRECIPITATION
if(par->output_Wr>0 && fmod(times->time+par->Dt,par->output_Wr*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Wr*3600.0));
	write_suffix(SSSS, n_file, 0);

	write_map(join_strings(files->co[fcint]+1,SSSS), 0, par->format_out, wat->wt, UV);
}

//SNOW ENERGY BALANCE
if(par->output_balancesn>0 && fmod(times->time+par->Dt,par->output_balancesn*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_balancesn*3600.0));
	write_suffix(SSSS, n_file, 0);

	write_map(join_strings(files->co[fmsn]+1,SSSS), 0, par->format_out, snow->MELTED, UV);
	initmatrix(0.0, snow->MELTED, top->Z0, NoV);
	write_map(join_strings(files->co[fssn]+1,SSSS), 0, par->format_out, snow->SUBL, UV);
	initmatrix(0.0, snow->SUBL, top->Z0, NoV);

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
	write_map(join_strings(files->co[fsnd]+1,SSSS), 0, par->format_out, M, UV);
}

//GLACIER ENERGY BALANCE
if(par->glaclayer_max>0 && par->output_balancegl>0 && fmod(times->time+par->Dt,par->output_balancegl*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_balancegl*3600.0));
	write_suffix(SSSS, n_file, 0);

	write_map(join_strings(files->co[fmgl]+1,SSSS), 0, par->format_out, glac->MELTED, UV);
	initmatrix(0.0, glac->MELTED, top->Z0, NoV);
	write_map(join_strings(files->co[fsgl]+1,SSSS), 0, par->format_out, glac->SUBL, UV);
	initmatrix(0.0, glac->SUBL, top->Z0, NoV);

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
	write_map(join_strings(files->co[fgld]+1,SSSS), 0, par->format_out, M, UV);
}

//SOLAR RADIATION
if(par->output_Rswdown>0 && fmod(times->time+par->Dt,par->output_Rswdown*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_Rswdown*3600.0));
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fSW]+1,"mean");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->Rswdown_mean, UV);
	initmatrix(0.0, egy->Rswdown_mean, top->Z0, NoV);

	name=join_strings(files->co[fSW]+1,"beam_mean");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->Rswbeam, UV);
	initmatrix(0.0, egy->Rswbeam, top->Z0, NoV);

	if(par->distr_stat==1){
		name=join_strings(files->co[fSW]+1,"max");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->Rswdown_max, UV);
		initmatrix(-1.0E+9, egy->Rswdown_max, top->Z0, NoV);
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
	write_map(join_strings(name,SSSS), 0, par->format_out, M, UV);
	initlongmatrix(0, egy->nDt_shadow, top->Z0, NoV);
	initlongmatrix(0, egy->nDt_sun, top->Z0, NoV);
}

//METEO
if(par->output_meteo>0 && fmod(times->time+par->Dt,par->output_meteo*3600.0)==0){
	n_file=(long)((times->time+par->Dt)/(par->output_meteo*3600.0));
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fTa]+1,"mean");
	write_map(join_strings(name,SSSS), 0, par->format_out, egy->Ta_mean, UV);
	initmatrix(0.0, egy->Ta_mean, top->Z0, NoV);

	if(par->distr_stat==1){
		name=join_strings(files->co[fTa]+1,"max");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->Ta_max, UV);
		initmatrix(-1.0E+9, egy->Ta_max, top->Z0, NoV);

		name=join_strings(files->co[fTa]+1,"min");
		write_map(join_strings(name,SSSS), 0, par->format_out, egy->Ta_min, UV);
		initmatrix(1.0E+9, egy->Ta_min, top->Z0, NoV);
	}

	if(par->micromet1==1){
		name=join_strings(files->co[fwspd]+1,"mean");
		write_map(join_strings(name,SSSS), 0, par->format_out, met->Vspdmean, UV);
		initmatrix(0.0, met->Vspdmean, top->Z0, NoV);

		name=join_strings(files->co[fwdir]+1,"mean");
		write_map(join_strings(name,SSSS), 0, par->format_out, met->Vdirmean, UV);
		initmatrix(0.0, met->Vdirmean, top->Z0, NoV);

		name=join_strings(files->co[frh]+1,"mean");
		write_map(join_strings(name,SSSS), 0, par->format_out, met->RHmean, UV);
		initmatrix(0.0, met->RHmean, top->Z0, NoV);
	}
}

free_doublematrix(M);



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


			for(l=1;l<=Nl;l++){
				write_tensorseries(1, l, isavings, files->co[rpsi]+1, 0, par->format_out, sl->P, UV);
				write_tensorseries(1, l, isavings, files->co[riceg]+1, 0, par->format_out, sl->thice, UV);
				write_tensorseries(1, l, isavings, files->co[rTg]+1, 0, par->format_out, sl->T, UV);
			}
			write_map(join_strings(files->co[rhsup]+1,SSSS), 0, par->format_out, wat->h_sup, UV);


			for(l=1;l<=par->snowlayer_max;l++){
				write_tensorseries(1, l, isavings, files->co[rDzs]+1, 0, par->format_out, snow->Dzl, UV);
				write_tensorseries(1, l, isavings, files->co[rwls]+1, 0, par->format_out, snow->w_liq, UV);
				write_tensorseries(1, l, isavings, files->co[rwis]+1, 0, par->format_out, snow->w_ice, UV);
				write_tensorseries(1, l, isavings, files->co[rTs]+1, 0, par->format_out, snow->T, UV);
			}
			write_map(join_strings(files->co[rsnag]+1,SSSS), 0, par->format_out, snow->age, UV);
			M=copydouble_longmatrix(snow->lnum);
			write_map(join_strings(files->co[rns]+1, SSSS), 1, par->format_out, M, UV);
			free_doublematrix(M);


			if(par->glaclayer_max>0){
				for(l=1;l<=par->glaclayer_max;l++){
					write_tensorseries(1, l, isavings, files->co[rDzi]+1, 0, par->format_out, glac->Dzl, UV);
					write_tensorseries(1, l, isavings, files->co[rwli]+1, 0, par->format_out, glac->w_liq, UV);
					write_tensorseries(1, l, isavings, files->co[rwii]+1, 0, par->format_out, glac->w_ice, UV);
					write_tensorseries(1, l, isavings, files->co[rTi]+1, 0, par->format_out, glac->T, UV);
				}
				M=copydouble_longmatrix(glac->lnum);
				write_map(join_strings(files->co[rni]+1,SSSS), 1, par->format_out, M, UV);
				free_doublematrix(M);
			}


			name=join_strings(files->co[rQch]+1,SSSS);
			name=join_strings(name,textfile);
			f=t_fopen(name,"w");
			fprintf(f,"/**Channel discharges(m3/s) for each channel-pixel\n");
			fprintf(f," Q_sup_s		Q_sub_s*/\n");
			fprintf(f,"index{1}\n");
			fprintf(f,"1:double matrix Q_channel {%ld,2}\n",cnet->Q_sup_s->nh);
			for(j=1;j<=cnet->Q_sup_s->nh;j++){
				fprintf(f,"%20.16f  %20.16f\n",cnet->Q_sup_s->co[j],cnet->Q_sub_s->co[j]);
			}
			t_fclose(f);


		}
	}
}

}

/*--------------------------------------------------------------------------------------------------*/













/****************************************************************************************************/
/* All the structs and substructs of the simulation are deallocated:                                */
/****************************************************************************************************/
void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met, LISTON *liston)
{

 /* Deallocation of struct TOPO "top": */
 printf("Deallocating top\n");
 free_doublematrix(top->Z0);
 free_doublematrix(top->sky);
 free_shortmatrix(top->pixel_type);
 free_shortmatrix(top->DD);
 free_doublematrix(top->i_DD);
 free_doublematrix(top->dz_dx);
 free_doublematrix(top->dz_dy);
 free_shortmatrix(top->curv);
 free_doublematrix(top->area);
 free_doublematrix(top->aspect);
 free_doublematrix(top->slopes);
 free_doublematrix(top->i_ch);
 free_doublematrix(top->pixel_distance);
 if(par->ES_num>0){
	free_longvector(top->ES_pixel);
	free_doublevector(top->ES_aspect);
	free_doublevector(top->ES_slope);
 }
 if(par->point_sim==1) free(top->horizon_height);
 free(top);

 /* Deallocation of struct SOIL "sl": */
 printf("Deallocating sl\n");
 free_doubletensor(sl->P);
 free_doubletensor(sl->T);
 free_doublematrix(sl->Tv);
 free_doubletensor(sl->thice);
 free_doublematrix(sl->Jinf);
 free_doubletensor(sl->J);
 free_shortmatrix(sl->type);
 free_doubletensor(sl->pa);
 if(par->output_TETAICExy>0) free_doubletensor(sl->thiav);
 if(par->output_TETAxy>0) free_doubletensor(sl->thwav);
 if(par->output_Txy>0) free_doubletensor(sl->Tav);
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
 free_doublematrix(wat->total);
 free_doublematrix(wat->Psnow);
 free_doublematrix(wat->Pn);
 free_doublematrix(wat->wt);
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
 free_doublevector(cnet->Q);
 free_doublevector(cnet->s0);
 free_doublematrix(cnet->fraction_spread);
 free_doublevector(cnet->Q_sup_s);
 free_doublevector(cnet->Q_sub_s);
 free_doublevector(cnet->Q_sup_spread);
 free_doublevector(cnet->Q_sub_spread);
 free(cnet);

 /* Deallocation of struct FILENAMES "filenames": */
 printf("Deallocating files\n");
 free(files);

 /* Deallocation of struct T_INIT "UV": */
 printf("Deallocating UV\n");
 free_doublevector(UV->U);
 free_doublevector(UV->V);
 free(UV);

 /* Deallocation of struct ENERGY "egy": */
 printf("Deallocating egy\n");
 free_doublematrix(egy->Hgrid);
 free_doublematrix(egy->Tsgrid);

 //if(par->micromet2==1) free_doublematrix(egy->SWin);
 if(par->micromet3==1) free_doublematrix(egy->LWin);

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
	if(par->distr_stat==1)free_doublematrix(egy->ET_max);
	if(par->distr_stat==1)free_doublematrix(egy->ET_min);
 }

 if(par->output_H>0){
	free_doublematrix(egy->H_mean);
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

 if(par->JD_plots->co[1]!=0){
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
 if(par->JD_plots->co[1]!=0) free_doublematrix(snow->Dplot);
 free_shortmatrix(snow->type);
 free_longmatrix(snow->lnum);
 free_doubletensor(snow->Dzl);
 free_doubletensor(snow->w_liq);
 free_doubletensor(snow->w_ice);
 free_doubletensor(snow->T);
 free_doublematrix(snow->rho_newsnow);
 if(par->snowtrans==1){
	free_doublematrix(snow->Wsalt);
	free_doublematrix(snow->Wsusp);
	free_doublematrix(snow->Wsubl);
	free_doublematrix(snow->Wsubgrid);
	free_doublematrix(snow->ListonSWE);
	free_doublematrix(snow->softSWE);
	free_doublematrix(snow->softSWE1);
	if(par->output_snow>0){
		free_doublematrix(snow->Wtot_cum);
		free_doublematrix(snow->Wsalt_cum);
		free_doublematrix(snow->Wsusp_cum);
		free_doublematrix(snow->Wsubl_cum);
		free_doublematrix(snow->Wsubgrid_cum);
	}
 }
 free_doublematrix(snow->out_bs);
 free_doublematrix(snow->age);
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
 if(par->micromet1==1){
	free_doublematrix(met->Vgrid);
	free_doublematrix(met->Vdir);
	free_doublematrix(met->RHgrid);
	if(par->output_meteo>0){
		free_doublematrix(met->Vspdmean);
		free_doublematrix(met->Vdirmean);
		free_doublematrix(met->RHmean);
	}
 }
 if(par->micromet2==1) free_doublematrix(met->CFgrid);
 if(par->JD_plots->co[1]!=0){
	free_doublematrix(met->Taplot);
	if(par->micromet1==1){
		free_doublematrix(met->Vspdplot);
		free_doublematrix(met->Vdirplot);
		free_doublematrix(met->RHplot);
	}
 }

 free(met->LT);
 free(met->Lrh);
 free(met->Lws);
 free(met->Lwd);
 free(met->LP);

 free(met->data);
 free(met->column);
 free(met->horizon);
 free(met->var);

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
 free_longvector(par->JD_plots);
 if(par->point_sim==1){
	if(par->micromet1==1 || par->micromet2==1 || par->micromet3==1){
		free_longvector(par->r_points);
		free_longvector(par->c_points);
	}
 }
 /*free(par->transect);
 free(par->vtrans);
 free_longvector(par->cont_trans);
 free_longvector(par->ibeg);*/
 free(par);

 printf("Deallocating liston\n");
 deallocate_liston(liston);
 free(liston);

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

	char ADS[ ]={"aaaaddddLssss"};

	write_suffix(ADS, y, 0);
	write_suffix(ADS, JD, 4);
	write_suffix(ADS, i, 9);
	write_map(join_strings(name,ADS), 0, format, M, UV);
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

	//DISCHARGE
	f=t_fopen(join_strings(files->co[fQ]+1,textfile),"w");
	fprintf(f,"DATE,JDfrom0,JD,Q_tot[mc/s],Qsub_ch[mc/s],Qsup_ch[mc/s],Q_G[mc/s]\n");
	t_fclose(f);

	root_fraction=new_doublevector(Nl);

	//MELT FLUXES
	if(par->output_balancesn!=0){
		for(i=1;i<=land->clax->nh;i++){
			write_suffix(SSSS, i, 0);
			name=join_strings(files->co[fmeltlu]+1,"_snowcovered_");
			name=join_strings(name,SSSS);
			name=join_strings(name,textfile);
			f=t_fopen(name,"w");
			fprintf(f,"/** LAND USE CLASS: %ld\n",land->clax->co[i]);
			fprintf(f,"SNOW COVERED AREA\n");
			fprintf(f,"Number of pixels for this land use: %ld/%ld\n",land->cont->co[i][1],(long)total_pixel);
			fprintf(f,"Surface of this land use (km2): %10.4f/%10.4f */\n",land->cont->co[i][1]*UV->U->co[1]*UV->U->co[2]*1.0E-6,
				total_pixel*UV->U->co[1]*UV->U->co[2]*1.0E-6);
			fprintf(f,"1)t0[s],2)t1[s],3)t[d], 4)rain[mm],5)meltsnow[mm],6)meltglac[mm],7)pixels,8)area[km2],9)frac.landuse,10)frac.basin\n");
			t_fclose(f);

			name=join_strings(files->co[fmeltlu]+1,"_snowfree_");
			name=join_strings(name,SSSS);
			name=join_strings(name,textfile);
			f=t_fopen(name,"w");
			fprintf(f,"/** LAND USE CLASS: %ld\n",land->clax->co[i]);
			fprintf(f,"SNOW FREE AREA\n");
			fprintf(f,"Number of pixels for this land use: %ld/%ld\n",land->cont->co[i][1],(long)total_pixel);
			fprintf(f,"Surface of this land use (km2): %10.4f/%10.4f\n",land->cont->co[i][1]*UV->U->co[1]*UV->U->co[2]*1.0E-6,
				total_pixel*UV->U->co[1]*UV->U->co[2]*1.0E-6);
			fprintf(f,"1)t0[s],2)t1[s],3)t[d], 4)rain[mm],5)meltsnow[mm],6)meltglac[mm],7)pixels,8)area[km2],9)frac.landuse,10)frac.basin\n");
			t_fclose(f);
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
								fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)),par->Esoil);
		}
		wat->out1->co[1][i]=wat->h_sup->co[r][c];
		wat->out1->co[14][i]=wat->wt->co[r][c];


		name=join_strings(files->co[fpoint]+1,"_info_");
		name=join_strings(name,SSSS);
		name=join_strings(name,textfile);
		f=t_fopen(name,"w");
		fprintf(f,"/** The main properties of the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld are:\n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f," Elevation above sea level: %10.3f m\n",top->Z0->co[r][c]);
		fprintf(f," Gauckler-Strickler [m^1/3/s]: %f\n",land->ty->co[lu][jcm]);
		for(l=1;l<=Nl;l++){
			fprintf(f," Residual water content[-] of the layer %ld: %f\n",l,sl->pa->co[sy][jres][l]);
		}
		for(l=1;l<=Nl;l++){
			fprintf(f," Saturated water content[-] of the layer %ld: %f\n",l,sl->pa->co[sy][jsat][l]);
		}
		for(l=1;l<=Nl;l++){
			fprintf(f," Alpha of van Genuchten[mm^-1] of the layer %ld: %f\n",l,sl->pa->co[sy][ja][l]);
		}
		for(l=1;l<=Nl;l++){
			fprintf(f," n of van Genuchten[-] of the layer %ld: %f\n",l,sl->pa->co[sy][jns][l]);
		}
		for(l=1;l<=Nl;l++){
			fprintf(f," m of van Genuchten[-] of the layer %ld: %f\n",l,1-1/sl->pa->co[sy][jns][l]);
		}
		for(l=1;l<=Nl;l++){
			fprintf(f," v of van Genuchten[-] of the layer %ld: %f\n",l,sl->pa->co[sy][jv][l]);
		}
		fprintf(f," Water content of wilting point [-]: %f\n",land->ty->co[lu][jtwp]);
		fprintf(f," Water content of field capacity [-]: %f\n",land->ty->co[lu][jtfc]);
		for(l=1;l<=Nl;l++){
			fprintf(f," Kv_sat of layer %ld [mm/s]: %f\n",l,sl->pa->co[sy][jKv][l]);
		}
		for(l=1;l<=Nl;l++){
			fprintf(f," Kh_sat of layer %ld [mm/s]: %f\n",l,sl->pa->co[sy][jKh][l]);
		}

		fprintf(f," Terrain elevation [m]: %f\n",top->Z0->co[r][c]);
		fprintf(f," Sky view factor [-]: %f\n",top->sky->co[r][c]);
		if (top->pixel_type->co[r][c]==0){
			fprintf(f," The pixel-type is land (0) \n");
		}else if (top->pixel_type->co[r][c]==9){
			fprintf(f," The pixel-type is novalue (9) \n");
		}else if (top->pixel_type->co[r][c]==10){
			fprintf(f," The pixel-type is channel (10) \n");
		}else if (top->pixel_type->co[r][c]==11){
			fprintf(f," The pixel-type is lake (11) \n");
		}else if (top->pixel_type->co[r][c]==12){
			fprintf(f," The pixel-type is sea (12) \n");
		}
		fprintf(f," Drainage Direction is %d \n",top->DD->co[r][c]);
		fprintf(f," Slope along Drainage Direction [-]: %f \n",top->i_DD->co[r][c]);
		fprintf(f," Slope along positive x direction [-]: %f \n",top->dz_dx->co[r][c]);
		fprintf(f," Slope along negative y direction [-]: %f \n",top->dz_dy->co[r][c]);
		fprintf(f," Topology of curvature (0-1) [-]: %d \n",top->curv->co[r][c]);
		fprintf(f," Area considering the slope [m^2]: %f \n",top->area->co[r][c]);
		fprintf(f," Aspect [deg] [0=Nord, clockwise]: %f \n",top->aspect->co[r][c]*180.0/Pi);
		fprintf(f," Mean slope of the pixel [deg]: %f \n",top->slopes->co[r][c]*180.0/Pi);
		fprintf(f," Slope to calculate the surface velocity of the channel incoming flow [-]: %f \n",top->i_ch->co[r][c]);
		fprintf(f," Land use number is %d \n",(short)land->LC->co[r][c]);

		root(land->ty->co[lu][jroot], sl->pa->co[sy][jdz], root_fraction->co);
		for(l=1;l<=Nl;l++){
			fprintf(f," The root fraction [-] of layer %ld: %f\n",l,root_fraction->co[l]);
		}

		//fprintf(f," Albedo without snow and alpha=0 [-]: %f \n",land->ty->co[lu][jalbedo]);
		fprintf(f," Surface fraction of land covered by vegetation [-]: %f \n",land->ty->co[lu][jcf]);
		fprintf(f," Momentum roughness length z0soil [m]: %f \n",land->ty->co[lu][jz0soil]);
		fprintf(f," Momentum roughness length z0veg [m]: %f \n",land->ty->co[lu][jz0veg]);
		fprintf(f," KRIGING WEIGHTS=\n");
		for(j=1;j<=n;j++){
			fprintf(f," STATION %ld = %f\n",j,wat->weights_Kriging->co[(r-1)*Nc+c][j]);
		}
		fprintf(f," */ \n");
		t_fclose(f);

		name=join_strings(files->co[fpoint]+1,SSSS);
		name=join_strings(name,textfile);
		f=t_fopen(name,"w");
		fprintf(f,"DATE,JDfrom0,JD,t[d],t_i[s],t_f[s],v[m/s],Vdir,RH[-],P[hPa],Tair[C],Tsurface[C],Tdew[C],eair[mbar],Qair[-],esurf[mbar]");//16
		fprintf(f,",Qsurf[-],SWin[W/m2],SWin_beam,SWin_diff,SWout[W/m2],albedo[-],alpha[deg],direction[deg],phi[deg],LWin[W/m2],LWout[W/m2],Rnet[W/m2],SW[W/m2],LW[W/m2],H[W/m2],LE[W/m2]");//32
		fprintf(f,",Qrain[W/m2],Gsoil[W/m2],SurfaceEB[W/m2],SWin_c[MJ],SWout_c[MJ],LWin_c[MJ],LWout_c[MJ],SW_cum[MJ],LWn_cum[MJ],Rnet_cm[MJ],H_cum[MJ],LE_cum[MJ]");//44
		fprintf(f,",G_cum[MJ],Eg[mm],Sg[mm],Etc[mm],Psnow[mm],Prain[mm],Psnow_c[mm],Prain_SOILc[mm],Prain_SNOWc[mm],Ptot_c[mm],Wt[mm],maxStor");//56
		fprintf(f,",DWt[mm],Ptot_atm,Rain_atm,Snow_atm,Evap_can,Drip_can,Ptot_atm_cum,Prain_atm_cum,Psnow_atm_cum,Evap_can_cum,Drip_can_cum,Pn[mm],Runoff[mm]");//69
		fprintf(f,",q_sup[mm],q_sub[mm],DS_sup[mm],DS_sub[mm],q_G[mm],snowDEPTH[mm],SWE[mm],snowDENSITY[kg/m3],snowT[C],BStot[mm],BStot_cum[mm]");//80
		fprintf(f,",snowMELT[mm],snowSUBL[mm],snowEVAP[mm],glacierDEPTH[mm],GWE[mm],gDENSITY[kg/m3],glcT[C],glcMELT[mm],glcSUBL[mm],glcEVAP[mm],qv(mm/d)");//91
		fprintf(f,",LWemissiv,LObukhov[m],numb.iter.,LWmin[W/m2],LWmax[W/m2],LAI,z0[m],d0[m],ErrorRichards[mm/h]");//100
		fprintf(f,",Ts,Tg,Tv,SWv,LWv,Hv,LEv,fwet,Htot,LEtot,Hg0,LEg0,Hg1,LEg1,fc,Ch,Cv,Cb,Cc,Ch_ic,Cv_ic,Hv,LEv,Qv,Qg,Qa,Qs\n");//120
		t_fclose(f);

		name=join_strings(files->co[fsnz]+1,SSSS);
		name=join_strings(name,textfile);
		f=t_fopen(name,"w");
		//fprintf(f,"/**Snow profile for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */\n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD,t_i[s],t_f[s],Snowtype,Snowdepth[mm],SWE[mm],Temp_aver[C],Density_aver[kg/m3],albedo,melting[mm],sublimation[mm],evaporation[mm],nlayer,BSsalt[mm],BSsalt_cum[mm],BSsusp[mm],BSsusp_cum[mm],BSsubl[mm],BSsubl_cum[mm],BSsbgr[mm],BSsbgr_cum[mm],BStot[mm],BStot_cum[mm]");
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

		if(par->glaclayer_max>0){
			name=join_strings(files->co[fglz]+1,SSSS);
			name=join_strings(name,textfile);
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
		}

		/*creation of the file "Tz.txt": */
		name=join_strings(files->co[fTz]+1,SSSS);
		name=join_strings(name,textfile);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of sl temperature for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);

		/*creation of the file "PSIz.txt": */
		name=join_strings(files->co[fpsiz]+1,SSSS);
		name=join_strings(name,textfile);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of water pressure (mm) for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);

		/*creation of the file "TETAz.txt": */
		name=join_strings(files->co[fliqz]+1,SSSS);
		name=join_strings(name,textfile);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of water content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);

		/*creation of the file "TETAICEz.txt": */
		name=join_strings(files->co[ficez]+1,SSSS);
		name=join_strings(name,textfile);
		f=t_fopen(name,"w");
		//fprintf(f,"/** Profiles of ice content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld:\n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
		fprintf(f,"DATE,JDfrom0,JD");
		z=0.0;
		for(l=1;l<=Nl;l++){
			z+=sl->pa->co[sy][jdz][l];
			fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
		}
		fprintf(f," \n");
		t_fclose(f);

	}

	//DATA BASIN
	f=t_fopen(join_strings(files->co[fbas]+1,textfile),"w");
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
		name=join_strings(files->co[farank]+1,SSSS);
		name=join_strings(name,textfile);

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
	}

	free_doublevector(root_fraction);


}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void write_soil_output(long n, long i, double t, double dt, long y0, double JD0, LONGMATRIX *rc, SOIL *sl, double psimin, double Esoil){
	/* Author:  Stefano Endrizzi  Year:
	* function writes the soil parameters of specified checkpoints pixels
	* Input:
	* 			n:	number of Dt after which the output of a pixel is printed
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
	* comment: Matteo Dall'Amico, May 2009 */
	char *name,SSSS[ ]={"SSSS"};
	double t_i, JD;
	long d2, mo2, y2, h2, mi2, l, r=rc->co[i][1], c=rc->co[i][2];
	FILE *f;
	short sy=sl->type->co[r][c];

	write_suffix(SSSS, i, 0);
	t_i=t-dt*n;
	//date_time(0.5*(t_i+dt+t+dt), y0, JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
	date_time(t+dt, y0, JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);

	/*update of the sl profile temperature in the control pixel:*/
	name=join_strings(files->co[fTz]+1,SSSS);
	name=join_strings(name,textfile);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",sl->T->co[l][r][c]);
	}
	fprintf(f," \n");
	fclose(f);

	/*writing psi as weighted mean: */
	name=join_strings(files->co[fpsiz]+1,SSSS);
	name=join_strings(name,textfile);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",sl->P->co[l][r][c]);
	}
	fprintf(f," \n");
	fclose(f);

	/*writing the water content as weighted mean: */
	name=join_strings(files->co[fliqz]+1,SSSS);
	name=join_strings(name,textfile);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",teta_psi(sl->P->co[l][r][c],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],fmin(sl->pa->co[sy][jpsimin][l],Psif(sl->T->co[l][r][c],psimin)),Esoil));
		//printf("alpha:%f\n",sl->pa->co[sy][ja][l]);
	}
	fprintf(f," \n");
	fclose(f);

	/*writing the ice content: */
	name=join_strings(files->co[ficez]+1,SSSS);
	name=join_strings(name,textfile);
	f=fopen(name,"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f",JD+(double)(daysfrom0(y2)),JD);
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",sl->thice->co[l][r][c]);
	}
	fprintf(f," \n");
	fclose(f);
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
			//printf("1\n");
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

