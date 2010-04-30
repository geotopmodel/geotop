
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */
    
    
    
#include "struct.geotop.h"
#include "output.h"
#include "pedo.funct.h"
#include "geo_statistic.h"
#include "networks.h"
#include "rw_maps.h"
#include "constant.h"
#include "keywords_file.h"
#include "extensions.h"
#include "times.h"
#include "energy.balance.h"
#include "input.h"
#include "tabs.h"
#include "vegetation.h"
#include "tables.h"

#include <time.h>
 
extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;
extern DOUBLEMATRIX *outdata_point;
extern DOUBLEVECTOR *outdata_basin;
extern char *MISSING_FILE;

/****************************************************************************************************/
/* write_output: stampa i dati degli output variabili nel tempo                                     */
/****************************************************************************************************/
void write_output(TIMES *times, WATER *wat, CHANNEL *cnet, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac, METEO *met)
				  
{
 /*internal auxiliary variables:*/
 long i,j,r,c,l; /*counters*/
 long n_file;      /*number of file of the type "TETAxySSSlZZ"(i.e. number of the basin-time-step)*/
 double t_i;         /*time of begin of an interval time for the output*/ 
 char SSSS[ ]={"SSSS"};
 char *name;	
 char *temp1,*temp2;  
 FILE *f;
 //time variables
 static time_t start_time, stop_time;
 double percent_done, elapsed_time, remaining_time, total_time;
	
// static double Qsub_ch, Qsup_ch;
 static double Qtot; /*averaged output flows*/
 static long isavings;
 static double mass_error_tot;
 
 /*internal variables to memorize input par:*/
 double time_max;      /*time of all the simulation [s]*/
 double snowD, SWE, snowdensity, snowtemperature;
 double glacD=0.0, GWE=0.0, glacdensity=UV->V->co[2], glactemperature=UV->V->co[2];
 DOUBLEMATRIX *M, *K;
 double JD,z;
 long d2,mo2,y2,h2,mi2;
 long sy;
 
 //initialize static variables
 if(times->time==0){
	 mass_error_tot=0.0;
	 time( &start_time );
 }
 
 //Assignment to some internal variables of some input par
 time_max=times->TH*3600.0;//seconds
 
 //times
 date_time(times->time+par->Dt, par->year0, par->JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);	

//DISCHARGE
//****************************************************************************************************************
//****************************************************************************************************************

if(par->state_discharge==1 && strcmp(files->co[fQ]+1 , MISSING_FILE) != 0){

	if(times->time==0) Qtot=0.0;
	Qtot += cnet->Q_out;
	
	if (times->i_discharge==times->n_discharge){/*Print the outlet flows*/

		t_i=times->time-par->Dt*times->n_discharge;
		name=join_strings(files->co[fQ]+1,textfile);
		f=fopen(name,"a");
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(times->time+par->Dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		fprintf(f,",%f\n",Qtot/((double)times->n_discharge));
		fclose(f);
		free(name);

		Qtot = 0.0;

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

		for(l=1;l<=Nl;l++){
			sl->T_av->co[l][i] += sl->T->co[l][r][c]/(double)times->n_pixel;
			sl->th_av->co[l][i] += sl->th->co[l][r][c]/(double)times->n_pixel;
			sl->thice_av->co[l][i] += sl->thice->co[l][r][c]/(double)times->n_pixel;
		}
		
		outdata_point->co[othawed][i]+=find_activelayerdepth(r, c, sl)/(double)times->n_pixel;					
		outdata_point->co[owtable][i]+=find_watertabledepth(r, c, sl)/(double)times->n_pixel;					
						
		//Print of pixel-output every times->n_pixel time step
		if (times->i_pixel==times->n_pixel){	
			
			//Point data
			if(strcmp(files->co[fpoint]+1 , MISSING_FILE) != 0){
				t_i=times->time-par->Dt*times->n_pixel;
				name=join_strings(files->co[fpoint]+1,SSSS);
				name=join_strings(name,textfile);
				f=fopen(name,"a");
				fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
				fprintf(f,",%f,%f,%f",(times->time+par->Dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   				
				for(j=1;j<=otot;j++){
					fprintf(f,",%f",outdata_point->co[j][i]);
				}
				fprintf(f, "\n");
				fclose(f);
			}
			
			//Snow
			if(strcmp(files->co[fsnz]+1 , MISSING_FILE) != 0){
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
						
				z=0.0;	   
				name=join_strings(files->co[fsnz]+1,SSSS);
				temp1=join_strings(name,textfile);
				f=fopen(temp1,"a");
				fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
				fprintf(f,",%f,%f,%f",(times->time+par->Dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
				fprintf(f,",%d,%f,%f,%f,%f,%f,%f,%f,%ld,%f,%f,%f,%f,%f,%f",snow->type->co[r][c],snowD,
					SWE,snowtemperature,snowdensity,outdata_point->co[omrsnow][i],outdata_point->co[osrsnow][i],outdata_point->co[oersnow][i],
					snow->lnum->co[r][c], snow->out_bs->co[1][i],snow->out_bs->co[2][i],snow->out_bs->co[3][i],snow->out_bs->co[4][i],
					snow->out_bs->co[9][i], snow->out_bs->co[10][i]);	   
						
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
							
				fprintf(f,"\n");	  
				fclose(f);	 
				free(temp1);
				free(name);
			}
			
			//Glacier
			if(par->glaclayer_max>0  && strcmp(files->co[fglz]+1 , MISSING_FILE) != 0){
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

				name=join_strings(files->co[fglz]+1,SSSS);
				temp1=join_strings(name,textfile);
				f=fopen(temp1,"a");			
				fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
				fprintf(f,",%f,%f,%f",(times->time+par->Dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
				fprintf(f,",%f,%f,%f,%ld,%f,%f,%f,%f,%f",glacD,GWE,glactemperature,glac->lnum->co[r][c],glacdensity,
					outdata_point->co[15][i]/*albedo*/,outdata_point->co[omrglac][i],outdata_point->co[osrglac][i],outdata_point->co[oerglac][i]);	   
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
				free(temp1);
				free(name);			
			}
			
			//sl output
			write_soil_output(times->n_pixel, i, times->time, par->Dt, par->year0, par->JD0, par->rc, sl, PsiMin);
							
		}//end(times->i_pixel==times->n_pixel)			
	}//end(i point cycle)
	
	if(times->i_pixel==times->n_pixel){
		for(i=1;i<=par->chkpt->nrh;i++){
			for(j=1;j<=otot;j++) { outdata_point->co[j][i]=0.0; }
			for(j=1;j<=5;j++) { snow->out_bs->co[2*j-1][i]=0.0; }
		}	

		percent_done = 100.0*(double)times->time/time_max;
		time( &stop_time );
		elapsed_time = difftime( stop_time, start_time );
		if( percent_done > 1.0e-6 ){
			total_time = elapsed_time * 100.0 / percent_done;
		}else{
			total_time = -999.9;
		}
		remaining_time = (total_time - elapsed_time);
	
	   	
		printf("%ld/%ld/%ld %2.0f:%02.0f %.2f%% - Time elapsed (hours:min:sec) %2.0f:%02.0f:%02.0f Time remaining (hours:min) %2.0f:%02.0f  \n",
			   d2,mo2,y2,(float)h2,(float)mi2,(100.0*(double)times->time/time_max), 
			   floor(elapsed_time / 3600.0), floor(((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60.),
			   floor((((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. - floor( ((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. ))*60.),
			   floor(remaining_time / 3600.0), floor(((remaining_time/3600)-floor(remaining_time / 3600.0))*60.) );
	
		
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f %.2f%% - Time elapsed (hours:min:sec) %2.0f:%02.0f:%02.0f Time remaining (hours:min) %2.0f:%02.0f  \n",
				d2,mo2,y2,(float)h2,(float)mi2,(100.0*(double)times->time/time_max), 
			floor(elapsed_time / 3600.0), floor(((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60.),
		   floor((((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. - floor( ((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. ))*60.),
		   floor(remaining_time / 3600.0), floor(((remaining_time/3600)-floor(remaining_time / 3600.0))*60.) );
		fclose(f);
	}
		
}else{
	
	percent_done = 100.0*(double)times->time/time_max;
	time( &stop_time );
	elapsed_time = difftime( stop_time, start_time );
	if( percent_done > 1.0e-6 ){
		total_time = elapsed_time * 100.0 / percent_done;
	}else{
		total_time = -999.9;
	}
	remaining_time = (total_time - elapsed_time);
	
	
	printf("%ld/%ld/%ld %2.0f:%02.0f %.2f%% - Time elapsed (hours:min:sec) %2.0f:%02.0f:%02.0f Time remaining (hours:min) %2.0f:%02.0f  \n",
		   d2,mo2,y2,(float)h2,(float)mi2,(100.0*(double)(times->time+par->Dt)/time_max), 
		   floor(elapsed_time / 3600.0), floor(((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60.),
		   floor((((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. - floor( ((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. ))*60.),
		   floor(remaining_time / 3600.0), floor(((remaining_time/3600)-floor(remaining_time / 3600.0))*60.) );
	
	
	f=fopen(files->co[ferr]+1,"a");
	fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f %.2f%% - Time elapsed (hours:min:sec) %2.0f:%02.0f:%02.0f Time remaining (hours:min) %2.0f:%02.0f  \n",
			d2,mo2,y2,(float)h2,(float)mi2,(100.0*(double)(times->time+par->Dt)/time_max), 
			floor(elapsed_time / 3600.0), floor(((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60.),
			floor((((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. - floor( ((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. ))*60.),
			floor(remaining_time / 3600.0), floor(((remaining_time/3600)-floor(remaining_time / 3600.0))*60.) );
	fclose(f);
	
}

//BASIN DATA
//****************************************************************************************************************
//****************************************************************************************************************

if(par->state_basin==1){
	
	if (times->i_basin==times->n_basin){

		if(strcmp(files->co[fbas]+1 , MISSING_FILE) != 0){
			name=join_strings(files->co[fbas]+1,textfile);
			f=fopen(name,"a");
			t_i=times->time-par->Dt*times->n_basin;
			fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
			fprintf(f,",%f,%f,%f",(times->time+par->Dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
			for(j=1;j<=ootot;j++){
				fprintf(f,",%f",outdata_basin->co[j]);
			}
			fprintf(f, "\n");
			fclose(f);
			free(name);
		}
		
		
		printf("\n%ld/%ld/%ld %2.0f:%02.0f JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
			   d2,mo2,y2,(float)h2,(float)mi2,times->JD,(long)(floor(times->time/86400))+1,
				(100.0*(double)(times->time+par->Dt)/time_max));

		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"\n%ld/%ld/%ld %2.0f:%02.0f JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
			   d2,mo2,y2,(float)h2,(float)mi2,times->JD,(long)(floor(times->time/86400))+1,
				(100.0*(double)(times->time+par->Dt)/time_max));
		fclose(f);
		
	
		mass_error_tot += outdata_basin->co[oomasserror];
		printf(" SW=%f W/m2  LW:%f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm\n\n",
			   outdata_basin->co[ooSW],outdata_basin->co[ooLW],outdata_basin->co[ooH],outdata_basin->co[ooLE],outdata_basin->co[oorainover],
			   outdata_basin->co[oosnowover],outdata_basin->co[oomasserror]*(3600.0/par->Dt)/(double)times->n_basin,mass_error_tot);
	
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f," SW=%f W/m2  LW:%f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm\n\n",
				outdata_basin->co[ooSW],outdata_basin->co[ooLW],outdata_basin->co[ooH],outdata_basin->co[ooLE],outdata_basin->co[oorainover],
				outdata_basin->co[oosnowover],outdata_basin->co[oomasserror]*(3600.0/par->Dt)/(double)times->n_basin,mass_error_tot);
		fclose(f);		
		

		for(j=1;j<=ootot;j++){ 
			outdata_basin->co[j]=0.0; 
		}
		
	}
}


//DISTRIBUTED OUTPUTS
//****************************************************************************************************************
//****************************************************************************************************************

//TETA
if(par->output_TETAxy>0 && fmod(times->time+par->Dt,par->output_TETAxy*3600.0)<1.E-5){
	n_file=(long)((times->time+par->Dt)/(par->output_TETAxy*3600.0));
	write_suffix(SSSS, n_file, 0);
	
	if(strcmp(files->co[fliq]+1 , MISSING_FILE) != 0) write_tensorseries2(n_file, files->co[fliq]+1, 0, par->format_out, sl->th, UV);	
	
	if(strcmp(files->co[fliqsup]+1 , MISSING_FILE) != 0){
		M=new_doublematrix(Nr, Nc);
		initialize_doublematrix(M,UV->V->co[2]);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					M->co[r][c] = sl->th->co[1][r][c];
				}
			}
		}
		temp1=join_strings(files->co[fliqsup]+1,SSSS);
		write_map(temp1, 0, par->format_out, M, UV);	
		free(temp1);
		free_doublematrix(M);
	}
}

	
//T
if(par->output_Txy>0 && fmod(times->time+par->Dt,par->output_Txy*3600.0)<1.E-5){
	n_file=(long)((times->time+par->Dt)/(par->output_Txy*3600.0));	
	write_suffix(SSSS, n_file, 0);
	
	//write T tensor
	if(strcmp(files->co[fT]+1 , MISSING_FILE) != 0) write_tensorseries2(n_file, files->co[fT]+1, 0, par->format_out, sl->T, UV);
	
	//calculate active layer depth		   
	if( strcmp(files->co[fthawed]+1 , MISSING_FILE) != 0 ){
		
		M=new_doublematrix(Nr,Nc);
		initialize_doublematrix(M,UV->V->co[2]);	
	
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					M->co[r][c]=find_activelayerdepth(r, c, sl);					
				}
			}
		}
		temp1=join_strings(files->co[fthawed]+1,SSSS);
		write_map(temp1, 0, par->format_out, M, UV);	
		free(temp1);
	}

}

//TETAICE
if(par->output_TETAICExy>0 && fmod(times->time+par->Dt,par->output_TETAICExy*3600.0)<1.E-5 && strcmp(files->co[fice]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_TETAICExy*3600.0));	
	write_tensorseries2(n_file, files->co[fice]+1, 0, par->format_out, sl->thice, UV);
}

//PSI
if(par->output_PSIxy>0 && fmod(times->time+par->Dt,par->output_PSIxy*3600.0)<1.E-5){
	n_file=(long)((times->time+par->Dt)/(par->output_PSIxy*3600.0));	
	write_suffix(SSSS, n_file, 0);
	
	//write psi tensors
	if(strcmp(files->co[fpsi]+1 , MISSING_FILE) != 0){
		temp1=join_strings(files->co[fpsi]+1,"TOT");
		write_tensorseries2(n_file, temp1, 0, par->format_out, sl->Ptot, UV);
		free(temp1);
		
		temp1=join_strings(files->co[fpsi]+1,"LIQ");
		write_tensorseries2(n_file, temp1, 0, par->format_out, sl->P, UV);
		free(temp1);		
	}
	
	//calculate saturation front depth
	if( strcmp(files->co[fwtable]+1 , MISSING_FILE) != 0 ){
		K=new_doublematrix(Nr,Nc);
		initialize_doublematrix(K,UV->V->co[2]);	
	
		initialize_doublematrix(K,UV->V->co[2]);				
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					K->co[r][c]=find_watertabledepth(r, c, sl);				
				}
			}
		}
		temp1=join_strings(files->co[fwtable]+1,SSSS);
		write_map(temp1, 0, par->format_out, K, UV);	
		free(temp1);
	}
	
}	


//***************************************************************************
//calculate frost table depth
if( par->output_Txy>0 && fmod(times->time+par->Dt,par->output_Txy*3600.0)<1.E-5 && par->output_PSIxy>0 && fmod(times->time+par->Dt,par->output_PSIxy*3600.0)<1.E-5 ){
	if( strcmp(files->co[fthawed]+1 , MISSING_FILE) != 0 && strcmp(files->co[fwtable]+1 , MISSING_FILE) != 0 && strcmp(files->co[fftable]+1 , MISSING_FILE) != 0){
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
		temp1=join_strings(files->co[fftable]+1,SSSS);
		write_map(temp1, 0, par->format_out, M, UV);	
		free(temp1);
	}
}
if( par->output_Txy>0 && fmod(times->time+par->Dt,par->output_Txy*3600.0)<1.E-5 && strcmp(files->co[fthawed]+1 , MISSING_FILE) != 0  ) free_doublematrix(M);
if( par->output_PSIxy>0 && fmod(times->time+par->Dt,par->output_PSIxy*3600.0)<1.E-5 && strcmp(files->co[fwtable]+1 , MISSING_FILE) != 0  ) free_doublematrix(K);
//end of frost table depth calculation
//***************************************************************************

//matrix to handle further data
M=new_doublematrix(Nr,Nc);
initialize_doublematrix(M,UV->V->co[2]);					

//SNOW DEPTH
if(par->output_snow>0 && fmod(times->time+par->Dt,par->output_snow*3600.0)<1.E-5 && strcmp(files->co[fsn]+1 , MISSING_FILE) != 0){
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
	
	//write_map(join_strings(join_strings(files->co[fsn]+1,"Dist"),SSSS), 0, par->format_out, M, UV);
	temp1=join_strings(files->co[fsn]+1,SSSS);
	write_map(temp1, 0, par->format_out, M, UV);
	free(temp1);
				
	//write_map(join_strings(join_strings(files->co[fsn]+1,"Dmax"),SSSS), 0, par->format_out, snow->max, UV);
	//initmatrix(0.0, snow->max, top->Z0, NoV);
	//write_map(join_strings(join_strings(files->co[fsn]+1,"Daverage"),SSSS), 0, par->format_out, snow->average, UV);
	//initmatrix(0.0, snow->average, top->Z0, NoV);
	//if(par->blowing_snow==1){
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsubl"),SSSS), 0, par->format_out, snow->Wsubl_cum, UV);
		//initmatrix(0.0, snow->Wsubl_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsalt"),SSSS), 0, par->format_out, snow->Wtrans_cum, UV);
		//initmatrix(0.0, snow->Wtrans_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsusp"),SSSS), 0, par->format_out, snow->Wsusp_cum, UV);
		//initmatrix(0.0, snow->Wsusp_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BSsbgr"),SSSS), 0, par->format_out, snow->Wsubgrid_cum, UV);
		//initmatrix(0.0, snow->Wsubgrid_cum, top->Z0, NoV);
		//write_map(join_strings(join_strings(files->co[fsn]+1,"BStot"),SSSS), 0, par->format_out, snow->Wtot, UV);
		//initmatrix(0.0, snow->Wtot, top->Z0, NoV);
	//}
}

//GLACIER DEPTH
if(par->glaclayer_max>0 && par->output_glac>0 && fmod(times->time+par->Dt,par->output_glac*3600.0)<1.E-5 && strcmp(files->co[fgl]+1 , MISSING_FILE) != 0){
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
	temp1=join_strings(files->co[fgl]+1,SSSS);
	write_map(temp1, 0, par->format_out, M, UV);
	free(temp1);
}

//WATER OVER THE SURFACE
if(par->output_h_sup>0 && strcmp(files->co[fhsup]+1 , MISSING_FILE) != 0){
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				wat->hsupav->co[r][c] += Fmax(0.0, sl->P->co[0][r][c])/((par->output_h_sup*3600.0)/(par->Dt));
				if(cnet->ch->co[r][c]>0) cnet->hsupav->co[r][c] += cnet->h_sup->co[ cnet->ch->co[r][c] ]/((par->output_h_sup*3600.0)/(par->Dt));
			}
		}
	}

	if(fmod(times->time+par->Dt,par->output_h_sup*3600.0)<1.E-5){
		n_file=(long)((times->time+par->Dt)/(par->output_h_sup*3600.0));	
		write_suffix(SSSS, n_file, 0);
	
		temp1 = join_strings(files->co[fhsup]+1,"LAND");
		temp2 = join_strings(temp1, SSSS);
		write_map(temp2, 0, par->format_out, wat->hsupav, UV);
		free(temp1);
		free(temp2);
	
		temp1 = join_strings(files->co[fhsup]+1,"CHANNEL");
		temp2 = join_strings(temp1, SSSS);
		write_map(temp2, 0, par->format_out, cnet->hsupav, UV);
		free(temp1);
		free(temp2);
	
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					wat->hsupav->co[r][c] = 0.;
					if(cnet->ch->co[r][c]>0) cnet->hsupav->co[r][c] = 0.;
				}
			}
		}
	}
}

//RADIATION
if(par->output_Rn>0 && fmod(times->time+par->Dt,par->output_Rn*3600.0)<1.E-5 && strcmp(files->co[fRn]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_Rn*3600.0));	
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fRn]+1,"nmean");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->Rn_mean, UV);
	initmatrix(0.0, egy->Rn_mean, land->LC, NoV);
	free(temp1);
	free(name);
			
	name=join_strings(files->co[fRn]+1,"LWin");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->LW_in, UV);
	initmatrix(0.0, egy->LW_in, land->LC, NoV);
	free(temp1);
	free(name);
			
	name=join_strings(files->co[fRn]+1,"LW");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->LW, UV);
	initmatrix(0.0, egy->LW, land->LC, NoV);
	free(temp1);
	free(name);

	name=join_strings(files->co[fRn]+1,"SW");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->SW, UV);
	initmatrix(0.0, egy->SW, land->LC, NoV);
	free(temp1);
	free(name);
			
	if(par->distr_stat==1){
		name=join_strings(files->co[fRn]+1,"nmax");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->Rn_max, UV);
		initmatrix(-1.0E+9, egy->Rn_max, land->LC, NoV);
		free(temp1);
		free(name);

		name=join_strings(files->co[fRn]+1,"nmin");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->Rn_min, UV);
		initmatrix(1.0E+9, egy->Rn_min, land->LC, NoV);
		free(temp1);
		free(name);

		name=join_strings(files->co[fRn]+1,"LWmax");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->LW_max, UV);
		initmatrix(-1.0E+9, egy->LW_max, land->LC, NoV);
		free(temp1);
		free(name);
	
		name=join_strings(files->co[fRn]+1,"LWmin");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->LW_min, UV);
		initmatrix(1.0E+9, egy->LW_min, land->LC, NoV);
		free(temp1);
		free(name);

		name=join_strings(files->co[fRn]+1,"SWmax");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->SW_max, UV);
		initmatrix(-1.0E+9, egy->SW_max, land->LC, NoV);	
		free(temp1);
		free(name);

	}
}

//GROUND HEAT FLUX
if(par->output_G>0 && fmod(times->time+par->Dt,par->output_G*3600.0)<1.E-5 && strcmp(files->co[fG]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_G*3600.0));	
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fG]+1,"mean");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->SEB_mean, UV);
	initmatrix(0.0, egy->SEB_mean, land->LC, NoV);
	free(temp1);
	free(name);
	
	name=join_strings(files->co[fG]+1,"snowsoil");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->G_snowsoil, UV);
	initmatrix(0.0, egy->G_snowsoil, land->LC, NoV);
	free(temp1);
	free(name);

	if(par->distr_stat==1){
		name=join_strings(files->co[fG]+1,"max");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->G_max, UV);
		initmatrix(-1.0E+9, egy->G_max, land->LC, NoV);
		free(temp1);
		free(name);

		name=join_strings(files->co[fG]+1,"min");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->G_min, UV);
		initmatrix(1.0E+9, egy->G_min, land->LC, NoV);	
		free(temp1);
		free(name);
	}
}

//SENSIBLE HEAT FLUX
if(par->output_H>0 && fmod(times->time+par->Dt,par->output_H*3600.0)<1.E-5 && strcmp(files->co[fH]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_H*3600.0));	
	write_suffix(SSSS, n_file, 0);

	//name=join_strings(files->co[fH]+1,"mean+");
	name=join_strings(files->co[fH]+1,"mean");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->H_mean, UV);
	initmatrix(0.0, egy->H_mean, land->LC, NoV);	
	//name=join_strings(files->co[fH]+1,"mean-");
	//write_map(join_strings(name,SSSS), 0, par->format_out, egy->H_mean2, UV);
	free(temp1);
	free(name);
				
	if(par->distr_stat==1){
		name=join_strings(files->co[fH]+1,"max");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->H_max, UV);
		initmatrix(-1.0E+9, egy->H_max, land->LC, NoV);
		free(temp1);
		free(name);
		
		name=join_strings(files->co[fH]+1,"min");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->H_min, UV);
		initmatrix(1.0E+9, egy->H_min, land->LC, NoV);
		free(temp1);
		free(name);
	}
}

//LATENT HEAT FLUX
if(par->output_ET>0 && fmod(times->time+par->Dt,par->output_ET*3600.0)<1.E-5 && strcmp(files->co[fLE]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_ET*3600.0));	
	write_suffix(SSSS, n_file, 0);

	//name=join_strings(files->co[fLE]+1,"mean+");
	name=join_strings(files->co[fLE]+1,"mean");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->ET_mean, UV);
	initmatrix(0.0, egy->ET_mean, land->LC, NoV);	
	//name=join_strings(files->co[fLE]+1,"mean-");
	//write_map(join_strings(name,SSSS), 0, par->format_out, egy->ET_mean2, UV);
	free(temp1);
	free(name);
		
	if(par->distr_stat==1){			
		name=join_strings(files->co[fLE]+1,"max");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->ET_max, UV);
		initmatrix(-1.0E+9, egy->ET_max, land->LC, NoV);			
		free(temp1);
		free(name);
		
		name=join_strings(files->co[fLE]+1,"min");	
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->ET_min, UV);	
		initmatrix(1.0E+9, egy->ET_min, land->LC, NoV);	
		free(temp1);
		free(name);
	}
}

//SURFACE TEMPERATURE
if(par->output_Ts>0 && fmod(times->time+par->Dt,par->output_Ts*3600.0)<1.E-5 && strcmp(files->co[fTs]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_Ts*3600.0));	
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fTs]+1,"mean");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, egy->Ts_mean, UV);		
	initmatrix(0.0, egy->Ts_mean, land->LC, NoV);
	free(temp1);
	free(name);
	
	if(par->distr_stat==1){
		name=join_strings(files->co[fTs]+1,"max");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->Ts_max, UV);	
		initmatrix(-1.0E+9, egy->Ts_max, land->LC, NoV);
		free(temp1);
		free(name);
		
		name=join_strings(files->co[fTs]+1,"min");
		temp1=join_strings(name, SSSS);
		write_map(temp1, 0, par->format_out, egy->Ts_min, UV);
		initmatrix(1.0E+9, egy->Ts_min, top->Z0, NoV);
		free(temp1);
		free(name);
	}
}

//PRECIPITATION
if(par->output_P>0 && fmod(times->time+par->Dt,par->output_P*3600.0)<1.E-5 && strcmp(files->co[fprec]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_P*3600.0));	
	write_suffix(SSSS, n_file, 0);
	
	name=join_strings(files->co[fprec]+1,"TOTAL");
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, wat->PrTOT_mean, UV);	
	initmatrix(0.0, wat->PrTOT_mean, top->Z0, NoV);
	free(temp1);
	free(name);
	
	name=join_strings(files->co[fprec]+1,"SNOW");	
	temp1=join_strings(name, SSSS);
	write_map(temp1, 0, par->format_out, wat->PrSNW_mean, UV);	
	initmatrix(0.0, wat->PrSNW_mean, land->LC, NoV);
	free(temp1);
	free(name);
}

//INTERCEPTED PRECIPITATION
if(par->output_Wr>0 && fmod(times->time+par->Dt,par->output_Wr*3600.0)<1.E-5 && strcmp(files->co[fcint]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_Wr*3600.0));	
	write_suffix(SSSS, n_file, 0);

	temp1=join_strings(files->co[fcint]+1,"water");
	temp2=join_strings(temp1,SSSS);
	write_map(temp2, 0, par->format_out, wat->wcan_rain, UV);
	free(temp1);
	free(temp2);
	
	temp1=join_strings(files->co[fcint]+1,"snow");
	temp2=join_strings(temp1,SSSS);
	write_map(temp2, 0, par->format_out, wat->wcan_snow, UV);
	free(temp1);
	free(temp2);

}
	
//SNOW ENERGY BALANCE	
if(par->output_balancesn>0 && fmod(times->time+par->Dt,par->output_balancesn*3600.0)<1.E-5){
	n_file=(long)((times->time+par->Dt)/(par->output_balancesn*3600.0));	
	write_suffix(SSSS, n_file, 0);

	if(strcmp(files->co[fmsn]+1 , MISSING_FILE) != 0){
		temp1=join_strings(files->co[fmsn]+1,SSSS);
		write_map(temp1, 0, par->format_out, snow->MELTED, UV);		
		initmatrix(0.0, snow->MELTED, land->LC, NoV);	
		free(temp1);
	}
	
	if(strcmp(files->co[fssn]+1 , MISSING_FILE) != 0){
		temp1=join_strings(files->co[fssn]+1,SSSS);
		write_map(temp1, 0, par->format_out, snow->SUBL, UV);	
		initmatrix(0.0, snow->SUBL, land->LC, NoV);	
		free(temp1);
	}
	
	if(strcmp(files->co[fswe]+1 , MISSING_FILE) != 0){
		for(r=1;r<=snow->Dzl->nrh;r++){
			for(c=1;c<=snow->Dzl->nch;c++){
				if(land->LC->co[r][c]!=NoV){
					M->co[r][c]=0.0;
					//D=0.0;
					for(l=1;l<=snow->lnum->co[r][c];l++){
						M->co[r][c]+=(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c]);
						//D+=0.001*snow->Dzl->co[l][r][c];
					}
					//if(D>0) M->co[r][c]/=D;
				}
			}
		}
		temp1=join_strings(files->co[fswe]+1,SSSS);
		write_map(temp1, 0, par->format_out, M, UV);
		free(temp1);
	}
	
	if(strcmp(files->co[fsndur]+1 , MISSING_FILE) != 0){
		temp1=join_strings(files->co[fsndur]+1,SSSS);
		write_map(temp1, 0, par->format_out, snow->t_snow, UV);	
		initmatrix(0.0, snow->t_snow, land->LC, NoV);			
		free(temp1);
	}
	
}		

//GLACIER ENERGY BALANCE	
if(par->glaclayer_max>0 && par->output_balancegl>0 && fmod(times->time+par->Dt,par->output_balancegl*3600.0)<1.E-5){
	n_file=(long)((times->time+par->Dt)/(par->output_balancegl*3600.0));	
	write_suffix(SSSS, n_file, 0);

	if(strcmp(files->co[fmgl]+1 , MISSING_FILE) != 0){
		temp1=join_strings(files->co[fmgl]+1,SSSS);
		write_map(temp1, 0, par->format_out, glac->MELTED, UV);		
		initmatrix(0.0, glac->MELTED, land->LC, NoV);	
		free(temp1);
	}
	
	if(strcmp(files->co[fsgl]+1 , MISSING_FILE) != 0){
		temp1=join_strings(files->co[fsgl]+1,SSSS);
		write_map(temp1, 0, par->format_out, glac->SUBL, UV);
		initmatrix(0.0, glac->SUBL, land->LC, NoV);				
		free(temp1);
	}
	
	if(strcmp(files->co[fgwe]+1 , MISSING_FILE) != 0){
		for(r=1;r<=glac->Dzl->nrh;r++){
			for(c=1;c<=glac->Dzl->nch;c++){
				if(land->LC->co[r][c]!=NoV){
					M->co[r][c]=0.0;
					//D=0.0;
					for(l=1;l<=glac->lnum->co[r][c];l++){
						M->co[r][c]+=(glac->w_liq->co[l][r][c]+glac->w_ice->co[l][r][c]);
						//D+=0.001*glac->Dzl->co[l][r][c];
					}
					//if(D>0) M->co[r][c]/=D;
				}else{
					M->co[r][c]=UV->V->co[2];
				}
			}
		}
		temp1=join_strings(files->co[fgwe]+1,SSSS);
		write_map(temp1, 0, par->format_out, M, UV);	
		free(temp1);
	}
}	
		
//SOLAR RADIATION
if(par->output_Rswdown>0 && fmod(times->time+par->Dt,par->output_Rswdown*3600.0)<1.E-5 && strcmp(files->co[fSW]+1 , MISSING_FILE) != 0){
	n_file=(long)((times->time+par->Dt)/(par->output_Rswdown*3600.0));	
	write_suffix(SSSS, n_file, 0);

	name=join_strings(files->co[fSW]+1,"mean");
	temp1=join_strings(name,SSSS);
	write_map(temp1, 0, par->format_out, egy->Rswdown_mean, UV);
	initmatrix(0.0, egy->Rswdown_mean, land->LC, NoV);
	free(temp1);
	free(name);
		
	name=join_strings(files->co[fSW]+1,"beam_mean");
	temp1=join_strings(name,SSSS);
	write_map(temp1, 0, par->format_out, egy->Rswbeam, UV);	
	initmatrix(0.0, egy->Rswbeam, land->LC, NoV);			
	free(temp1);
	free(name);

	if(par->distr_stat==1){					
		name=join_strings(files->co[fSW]+1,"max");
		temp1=join_strings(name,SSSS);
		write_map(temp1, 0, par->format_out, egy->Rswdown_max, UV);
		initmatrix(-1.0E+9, egy->Rswdown_max, land->LC, NoV);	
		free(temp1);
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
	name=join_strings(files->co[fSW]+1,"shadowfractime");
	temp1=join_strings(name,SSSS);
	write_map(temp1, 0, par->format_out, M, UV);
	initlongmatrix(0, egy->nDt_shadow, land->LC, NoV);	
	initlongmatrix(0, egy->nDt_sun, land->LC, NoV);		
	free(temp1);
	free(name);
}

//METEO
if(par->output_meteo>0 && fmod(times->time+par->Dt,par->output_meteo*3600.0)<1.E-5){
	n_file=(long)((times->time+par->Dt)/(par->output_meteo*3600.0));	
	write_suffix(SSSS, n_file, 0);
		
	if(strcmp(files->co[fTa]+1 , MISSING_FILE) != 0){
		name=join_strings(files->co[fTa]+1,"mean");
		temp1=join_strings(name,SSSS);
		write_map(temp1, 0, par->format_out, egy->Ta_mean, UV);
		initmatrix(0.0, egy->Ta_mean, land->LC, NoV);	
		free(temp1);
		free(name);
	
		if(par->distr_stat==1){
			name=join_strings(files->co[fTa]+1,"max");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, egy->Ta_max, UV);
			initmatrix(-1.0E+9, egy->Ta_max, land->LC, NoV);
			free(temp1);
			free(name);

			name=join_strings(files->co[fTa]+1,"min");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, egy->Ta_min, UV);	
			initmatrix(1.0E+9, egy->Ta_min, land->LC, NoV);	
			free(temp1);
			free(name);
		}
	}
	
	if(par->micromet==1){
		if(strcmp(files->co[fwspd]+1 , MISSING_FILE) != 0){
			name=join_strings(files->co[fwspd]+1,"mean");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, met->Vspdmean, UV);
			initmatrix(0.0, met->Vspdmean, land->LC, NoV);	
			free(temp1);
			free(name);
		}
		
		if(strcmp(files->co[fwdir]+1 , MISSING_FILE) != 0){
			name=join_strings(files->co[fwdir]+1,"mean");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, met->Vdirmean, UV);	
			initmatrix(0.0, met->Vdirmean, land->LC, NoV);	
			free(temp1);
			free(name);
		}
		
		if(strcmp(files->co[frh]+1 , MISSING_FILE) != 0){
			name=join_strings(files->co[frh]+1,"mean");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, met->RHmean, UV);	
			initmatrix(0.0, met->RHmean, land->LC, NoV);
			free(temp1);
			free(name);
		}
	}
}	
	
free_doublematrix(M);

				



/**********************************************************************************************************/
/**********************************************************************************************************/
//SPECIAL PLOTS AT SOME DAYS
/**********************************************************************************************************/
/**********************************************************************************************************/

if(times->i_plot==times->n_plot){
	
	printf("\nWriting output data for JD:%ld year:%ld file:%ld\n",times->d_plot, times->year, times->nt_plot);
	f=fopen(files->co[ferr]+1,"a");
	fprintf(f,"\nWriting output data for JD:%ld year:%ld file:%ld  ",times->d_plot, times->year, times->nt_plot);
	fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",times->day,times->month,times->year,(float)times->hour,(float)times->min);

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
	plot(files->co[pH]+1, times->d_plot, times->year, times->nt_plot, M, par->format_out);

	initialize_doublematrix(M,UV->V->co[2]);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				M->co[r][c]=egy->LEgplot->co[r][c] + egy->LEvplot->co[r][c];
			}
		}
	}
	plot(files->co[pLE]+1, times->d_plot, times->year, times->nt_plot, M, par->format_out);
					
	plot(files->co[pHg]+1, times->d_plot, times->year, times->nt_plot, egy->Hgplot, par->format_out);
	plot(files->co[pLEg]+1, times->d_plot, times->year, times->nt_plot, egy->LEgplot, par->format_out);
	plot(files->co[pHv]+1, times->d_plot, times->year, times->nt_plot, egy->Hvplot, par->format_out);
	plot(files->co[pLEv]+1, times->d_plot, times->year, times->nt_plot, egy->LEvplot, par->format_out);
	plot(files->co[pSWin]+1, times->d_plot, times->year, times->nt_plot, egy->SWinplot, par->format_out);
	plot(files->co[pSWg]+1, times->d_plot, times->year, times->nt_plot, egy->SWgplot, par->format_out);
	plot(files->co[pSWv]+1, times->d_plot, times->year, times->nt_plot, egy->SWvplot, par->format_out);
	plot(files->co[pLWin]+1, times->d_plot, times->year, times->nt_plot, egy->LWinplot, par->format_out);
	plot(files->co[pLWg]+1, times->d_plot, times->year, times->nt_plot, egy->LWgplot, par->format_out);
	plot(files->co[pLWv]+1, times->d_plot, times->year, times->nt_plot, egy->LWvplot, par->format_out);
	plot(files->co[pTs]+1, times->d_plot, times->year, times->nt_plot, egy->Tsplot, par->format_out);
	plot(files->co[pTg]+1, times->d_plot, times->year, times->nt_plot, egy->Tgplot, par->format_out);
	plot(files->co[pTv]+1, times->d_plot, times->year, times->nt_plot, egy->Tvplot, par->format_out);
	plot(files->co[pTa]+1, times->d_plot, times->year, times->nt_plot, met->Taplot, par->format_out);
	plot(files->co[pD]+1, times->d_plot, times->year, times->nt_plot, snow->Dplot, par->format_out);
	if(par->micromet==1){
		plot(files->co[pVspd]+1, times->d_plot, times->year, times->nt_plot, met->Vspdplot, par->format_out);
		plot(files->co[pVdir]+1, times->d_plot, times->year, times->nt_plot, met->Vdirplot, par->format_out);
		plot(files->co[pRH]+1, times->d_plot, times->year, times->nt_plot, met->RHplot, par->format_out);
	}
	
	initialize_doublematrix(M,UV->V->co[2]);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				sy=sl->type->co[r][c];
				l=1;
				M->co[r][c]=sl->th->co[l][r][c];
			}
		}
	}
	plot(files->co[pth]+1, times->d_plot, times->year, times->nt_plot, M, par->format_out);
	
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

			for(l=0;l<=Nl;l++){
				write_tensorseries(1, l, isavings, files->co[rpsi]+1, 0, par->format_out, sl->P, UV);
			}
			
			for(l=1;l<=Nl;l++){
				write_tensorseries(1, l, isavings, files->co[riceg]+1, 0, par->format_out, sl->thice, UV);
				write_tensorseries(1, l, isavings, files->co[rTg]+1, 0, par->format_out, sl->T, UV);
			}
			write_map(join_strings(files->co[rwcrn]+1,SSSS), 0, par->format_out, wat->wcan_rain, UV);			
			write_map(join_strings(files->co[rwcsn]+1,SSSS), 0, par->format_out, wat->wcan_snow, UV);			

			write_map(join_strings(files->co[rTv]+1,SSSS), 0, par->format_out, sl->Tv, UV);			

			
			for(l=1;l<=par->snowlayer_max;l++){
				write_tensorseries(1, l, isavings, files->co[rDzs]+1, 0, par->format_out, snow->Dzl, UV);
				write_tensorseries(1, l, isavings, files->co[rwls]+1, 0, par->format_out, snow->w_liq, UV);
				write_tensorseries(1, l, isavings, files->co[rwis]+1, 0, par->format_out, snow->w_ice, UV);
				write_tensorseries(1, l, isavings, files->co[rTs]+1, 0, par->format_out, snow->T, UV);
			}
			write_map(join_strings(files->co[rsnag_adim]+1,SSSS), 0, par->format_out, snow->nondimens_age, UV);			
			write_map(join_strings(files->co[rsnag_dim]+1,SSSS), 0, par->format_out, snow->dimens_age, UV);			
			
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

			M=new_doublematrix(Nr, Nc);
			initialize_doublematrix(M, NoV);
			for(i=1;i<=cnet->r->nh;i++){
				if(cnet->r->co[i]>0) M->co[cnet->r->co[i]][cnet->c->co[i]] = cnet->h_sup->co[i];
			}
			write_map(join_strings(files->co[rQch]+1,SSSS), 0, par->format_out, M, UV);			
			free_doublematrix(M);
			
			/*name=join_strings(files->co[rQch]+1,SSSS);
			name=join_strings(name,textfile);
			f=t_fopen(name,"w");	
			fprintf(f,"index{1}\n");
			fprintf(f,"1:double matrix Q_channel {%ld,2}\n",cnet->Q_sup_s->nh);
			for(j=1;j<=cnet->Q_sup_s->nh;j++){
				fprintf(f,"%20.16f  %20.16f\n",cnet->Q_sup_s->co[j],cnet->Q_sub_s->co[j]);
			}
			t_fclose(f);*/

			
		}
	}
}
 
}

/*--------------------------------------------------------------------------------------------------*/













/****************************************************************************************************/
/* All the structs and substructs of the simulation are deallocated:                                */
/****************************************************************************************************/
void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met)
{

 long i,j,r,l;

 /* Deallocation of struct TOPO "top": */
 printf("Deallocating top\n");
 free_doublematrix(top->Z0);
 free_doublematrix(top->sky);
 free_shortmatrix(top->pixel_type);
 free_shortmatrix(top->DD);
 free_doublematrix(top->i_DD);
 //free_doublematrix(top->dz_dx);
 //free_doublematrix(top->dz_dy);
 free_doublematrix(top->aspect);
 free_doublematrix(top->slopes);
 if(par->point_sim==1) free(top->horizon_height);
 if(par->point_sim==1 && par->micromet==1) free_doublematrix(top->Z1);
 if(par->micromet==1){
	free_doublematrix(top->Zm);
	free_doublematrix(top->curv_m);
	free_doublematrix(top->slope_m);
	free_doublematrix(top->slopeaz_m);
 }
  
 free_longmatrix(top->lrc_cont);

 for(l=1;l<=Nl;l++){
	for(r=1;r<=Nr;r++){
		free(top->i_cont[l][r]);
	}
	free(top->i_cont[l]);
 }
 free(top->i_cont);
 
 free_longmatrix(top->rc_cont);
 
 for(r=1;r<=Nr;r++){
	free(top->j_cont[r]);
 }
 free(top->j_cont); 
 
 free_doubletensor(top->Z);
 free_doublematrix(top->Z0dp);
	 
 free_longvector(top->Lp);
 free_longvector(top->Li);
// free_longvector(top->Up);
// free_longvector(top->Ui);
  
 free(top);
 
 /* Deallocation of struct SOIL "sl": */
 printf("Deallocating sl\n");
 free_doubletensor(sl->P);
 free_doubletensor(sl->Ptot);
 free_doubletensor(sl->T);
 free_doublematrix(sl->Tv);
 free_doubletensor(sl->thice); 
 free_doubletensor(sl->th); 
 free_longmatrix(sl->type);
 free_doubletensor(sl->pa);
 free_doubletensor(sl->ET);	
 if(par->state_pixel==1){
	 free_doublematrix(sl->T_av);
	 free_doublematrix(sl->th_av);
	 free_doublematrix(sl->thice_av);	
 }
 free(sl);

 /* Deallocation of struct LAND "land": */
 printf("Deallocating land\n");
 free_doublematrix(land->LC);
 free_shortmatrix(land->shadow);
 free_doublematrix(land->ty);
	
 free_longmatrix(land->vegparp);
 free_doublevector(land->vegpar);
 for(i=1;i<=par->n_landuses;i++){
	free(land->vegparv[i]);
	if(par->vegflag->co[i]==1){ 
		for(j=0;j<dim2(land->vegpars[i]);j++){
			free(land->vegpars[i][j]);
		}
		free(land->vegpars[i]);
	}
 }
 free(land->vegpars);	
 free(land->vegparv);	
 
 free(land);
	
 /* Deallocation of struct WATER "water": */
 printf("Deallocating water\n"); 
 free_doublematrix(wat->weights_Kriging);
 free_doublematrix(wat->PrecTot);
 free_doublematrix(wat->Pnet);
 free_doublematrix(wat->wcan_rain);
 free_doublematrix(wat->wcan_snow);

 if (par->output_P>0){
	free_doublematrix(wat->PrTOT_mean);
	free_doublematrix(wat->PrSNW_mean);
 }

 if(par->output_h_sup>0) free_doublematrix(wat->hsupav);
	free_doublematrix(wat->dh_sup);
 //free_UMFPACK_REAL_MATRIX(wat->Jmatrix);
 //free_UMFPACK_REAL_TRIPLET(wat->Jtriplet);
 free_doublevector(wat->Lx);
 //free_doublevector(wat->Ux);
	free_doublevector(wat->H0);
	free_doublevector(wat->H1);
	free_doublevector(wat->dH);
	free_doublevector(wat->B);
	free_doublevector(wat->f);
	free_doublevector(wat->df);
 free(wat);
 
 /* Deallocation of struct CHANNEL "channel": */
 printf("Deallocating channel network\n"); 
 free_longvector(cnet->r);
 free_longvector(cnet->c);
 free_longmatrix(cnet->ch);
 free_doublevector(cnet->Qsup);
 free_doublevector(cnet->Qsub);
 /*free_doublevector(cnet->s0);
 free_doublematrix(cnet->fraction_spread);
 free_doublevector(cnet->Q_sup_s);
 free_doublevector(cnet->Q_sub_s);
 free_doublevector(cnet->Qsup_spread);
 free_doublevector(cnet->Qsub_spread);*/
	free_doublevector(cnet->h_sup);
	free_doublevector(cnet->dh_sup);
	free_doublevector(cnet->length);
if(par->output_h_sup>0) free_doublematrix(cnet->hsupav);
 free(cnet);
 
 /* Deallocation of struct FILENAMES "filenames": */
 printf("Deallocating files\n"); 
 free_stringbin(files);
 free(MISSING_FILE);
 
 /* Deallocation of struct T_INIT "UV": */
 printf("Deallocating UV\n"); 
 free_doublevector(UV->U);
 free_doublevector(UV->V);
 free(UV);
  
 /* Deallocation of struct ENERGY "egy": */
 printf("Deallocating egy\n");  
 if(par->output_Rn>0){
	free_doublematrix(egy->Rn_mean);
	if(par->distr_stat==1)free_doublematrix(egy->Rn_max);	
	if(par->distr_stat==1)free_doublematrix(egy->Rn_min);
	if(par->distr_stat==1)free_doublematrix(egy->LW_max);
	if(par->distr_stat==1)free_doublematrix(egy->LW_min);
	free_doublematrix(egy->LW_in);
	free_doublematrix(egy->LW);
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
	free_doublematrix(egy->SEB_mean);
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
 
 free_doublevector(egy->Dlay);
 free_doublevector(egy->wliq);
 free_doublevector(egy->wice);
 free_doublevector(egy->Temp); 
 free_doublevector(egy->deltaw);
 free_doublevector(egy->SWlayer);
 free_doublevector(egy->soil_transp_layer);
 free_doublevector(egy->dFenergy);
 free_doublevector(egy->Kth0);
 free_doublevector(egy->Kth1);
 free_doublevector(egy->Fenergy);
 free_doublevector(egy->Newton_dir);
 free_doublevector(egy->T0);
 free_doublevector(egy->T1);
 free_doublevector(egy->Tstar);
 free_doublevector(egy->THETA);
 free_doublevector(egy->soil_evap_layer_bare);
 free_doublevector(egy->soil_evap_layer_veg);
 
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
 if(par->blowing_snow==1) free_longvector(snow->change_dir_wind);
 free(snow);

 printf("Deallocating glacier\n");
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
 
 free(met->data);
 free(met->column);
 free(met->horizon);
 free(met->var);
 
 free(met->LRs);
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
 free_doublevector(par->JD_plots); 
 if(par->point_sim==1){
	if(par->micromet==1){
		free_longvector(par->r_points);
		free_longvector(par->c_points);
	}
 } 
 free_shortvector(par->vegflag);
 /*free(par->transect);
 free(par->vtrans);
 free_longvector(par->cont_trans);
 free_longvector(par->ibeg);*/
 free(par);
	
 free_doublematrix(outdata_point);
 free_doublevector(outdata_basin);


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

void write_output_headers(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac){

	/*internal auxiliary variables:*/
	long i,r,c,l,j;
	char *name,*temp,*temp2,SSSS[ ]={"SSSS"};
	double z;
	long sy;
	short lu;
	DOUBLEVECTOR *root_fraction;
	FILE *f; 
	char *header_point[otot+1];
	char *header_basin[ootot+1];
		
	//output matrix and vectors
	outdata_point=new_doublematrix(otot,par->chkpt->nrh);
	initialize_doublematrix(outdata_point,0.0);
	
	outdata_basin=new_doublevector(ootot);
	initialize_doublevector(outdata_basin,0.0);
	
	//DISCHARGE
	if (par->state_discharge==1 && strcmp(files->co[fQ]+1 , MISSING_FILE) != 0){
		temp=join_strings(files->co[fQ]+1,textfile);
		f=t_fopen(temp,"w");
		//fprintf(f,"DATE[day/month/year hour:min],JDfrom0,JD,Q_tot[mc/s],Qsub_ch[mc/s],Qsup_ch[mc/s]\n");
		fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,Q_tot[m3/s]\n");
		t_fclose(f);
		free(temp);
	}
	
	
	if(par->state_pixel==1){	
						
		header_point[osnowover] =  "Psnow_over_canopy[mm]"; //       prec_snow_atm;
		header_point[orainover] =  "Prain_over_canopy[mm]"; //        prec_rain_atm;
		header_point[oprecsnow] =  "Psnow_under_canopy[mm]"; //        prec_snow;
		header_point[oprecrain] =  "Prain_under_canopy[mm]"; //         (prec_rain_on_soil+prec_rain_on_snow);
		header_point[orainonsnow] ="Prain_rain_on_snow[mm]"; //        prec_rain_on_snow;			
		header_point[oV] =  "Wind speed[m/s]"; //          Vpoint/(double)n;
		header_point[oVdir] =  "Wind direction[deg]";  //     (met->Vdir[r][c])/(double)n;
		header_point[oRH] =  "Relative Humidity[-]"; //    RHpoint/(double)n;
		header_point[oPa] =  "Pressure[kPa]";  //   Ppoint/(double)n;
		header_point[oTa] =  "Tair[C]"; //    Tpoint/(double)n;
		header_point[oTdew] = "Tdew[C]";//      Tdew/(double)n;
		header_point[oTg] =  "Tsurface[C]"; //    Tg/(double)n; //Ts[C]
		header_point[oTv] =  "Tvegetation[C]"; //    Tv/(double)n;
		header_point[oTs] =  "Tcanopyair[C]";//    Ts/(double)n;
		header_point[oEB] =  "Surface Energy balance[W/m2]"; //    surfEB/(double)n; 
		header_point[oG] =  "Soil heat flux[W/m2]"; //   G/(double)n; 
		header_point[oSWin] =  "SWin[W/m2]"; //      SWin/(double)n;
		header_point[oSWb] =  "SWbeam[W/m2]";//     (SWbeam/(double)n);
		header_point[oSWd] =  "SWdiff[W/m2]"; //     (SWdiff/(double)n);
		header_point[oLWin] =  "LWin[W/m2]";//      LWin/(double)n;	
		header_point[ominLWin] =  "LWin_min[W/m2]"; //      (epsa_min*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
		header_point[omaxLWin] =  "LWin_max[W/m2]";//     (epsa_max*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
		header_point[oSW] = "SWnet[W/m2]"; //     SW/(double)n;
		header_point[oLW] =  "LWnet[W/m2]"; //     LW/(double)n;
		header_point[oH] = "H[W/m2]"; //    H/(double)n; //H[W/m^2]
		header_point[oLE] =  "LE[W/m2]"; //     LE/(double)n; //ET[W/m^2]
		header_point[ofc] =  "Canopy fraction[-]"; //     fc/(double)n;
		header_point[oLSAI] =  "LSAI[m2/m2]"; //      LSAI/(double)n;
		header_point[oz0v] =  "z0veg[m]";//      z0/(double)n;
		header_point[od0v] =  "d0veg[m]";//      d0/(double)n;	
		header_point[oEcan] =  "Estored_canopy[W/m2]"; //       (SWv+LWv-Hv-LEv)/(double)n;
		header_point[oSWv] =  "SWv[W/m2]"; //      SWv/(double)n;
		header_point[oLWv] =  "LWv[W/m2]"; //      LWv/(double)n;
		header_point[oHv] =  "Hv[W/m2]"; //     Hv/(double)n;
		header_point[oLEv] =  "LEv[W/m2]";//      LEv/(double)n;
		header_point[oHg0] =  "Hg_unveg[W/m2]"; //      Hg0/(double)n;
		header_point[oLEg0] =  "LEg_unveg[W/m2]"; //       Levap(Tg)*Eg0/(double)n;
		header_point[oHg1] =  "Hg_veg[W/m2]"; //      Hg1/(double)n;
		header_point[oLEg1] =  "LEg_veg[W/m2]"; //       Levap(Tg)*Eg1/(double)n;
		header_point[oevapsur] =  "Evap_surface[mm]"; //        Er_soil*par->Dt;	//Eg[mm]
		header_point[otrasp] =  "Trasp_canopy[mm]"; //      Evt*(1000.0/rho_w)*par->Dt;	//Etc[mm]	
		header_point[owcan_rain] =  "Water on canopy[mm]";  //        wat->wcan_rain[r][c]/(double)n;
		header_point[owcan_snow] =  "Snow on canopy[mm]";   //    wat->wcan_snow[r][c]/(double)n;	
		header_point[oQv] =  "Qvegetation[-]"; //       (Qv)/(double)n;
		header_point[oQg] =   "Qsurface[-]";  //       (Qg)/(double)n;
		header_point[oQa] =   "Qair[-]";  //       (Qa)/(double)n;
		header_point[oQs] =   "Qcanopyair[-]";  //       (Qs)/(double)n;
		header_point[oLobuk] =   "LObukhov[m]";  //        turbulence[2]/(double)n;
		//header_point[oiterLo] =   "Number of iter."; //        turbulence[1]/(double)n;
		header_point[oLobukcan] =   "LObukhovcanopy[m]";//        (Locc)/(double)n;			
		header_point[outop] =   "Wind speed top canopy[m/s]";//     (u_top)/(double)n;
		header_point[odecay] =  "Decay of K in canopy[-]"; //     (decay)/(double)n;
		header_point[oSWup] =   "SWup[W/m2]"; //       SWup_above_v/(double)n;
		header_point[oLWup] =   "LWup[W/m2]"; //     LWup_above_v/(double)n;
		header_point[oHup] =  "Hup[W/m2]";  //      (H+fc*Hv)/(double)n;
		header_point[oLEup] =  "LEup[W/m2]";  //       (LE+fc*LEv)/(double)n;
		header_point[omrsnow] = "snow melted[mm]"; //       Mr_snow*par->Dt;	//[mm]
		header_point[oersnow] =  "snow evap[mm]"; //       Er_snow*par->Dt;	//[mm]
		header_point[osrsnow] =  "snow subl[mm]"; //       Sr_snow*par->Dt;	//[mm]
		header_point[omrglac] =  "glac melted[mm]"; //       Mr_glac*par->Dt;	//[mm]
		header_point[oerglac] =  "glac evap[mm]"; //       Er_glac*par->Dt;	//[mm]
		header_point[osrglac] =  "glac subl[mm]"; //       Sr_glac*par->Dt;	//[mm]*/
		header_point[othawed] =  "thawed soil depth[mm]";
		header_point[owtable] =  "water table depth[mm]";

		root_fraction=new_doublevector(Nl);
		
   
		//DATA POINTS
		for(i=1;i<=par->chkpt->nrh;i++){
						
			write_suffix(SSSS, i, 0);
			r=par->rc->co[i][1];
			c=par->rc->co[i][2];				
			sy=sl->type->co[r][c];
			lu=(short)land->LC->co[r][c];
			
			if(strcmp(files->co[fpoint]+1 , MISSING_FILE) != 0){			
				name=join_strings(files->co[fpoint]+1,"_info_");
				temp=join_strings(name,SSSS);
				temp2=join_strings(temp,textfile);
				f=t_fopen(temp2,"w");
				
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
				for(l=1;l<=Nl;l++){
					fprintf(f," Water content of wilting point [-] of the layer %ld: %f\n",l,sl->pa->co[sy][jwp][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Water content of field capacity [-] of the layer %ld: %f\n",l,sl->pa->co[sy][jfc][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Kv_sat of layer %ld [mm/s]: %f\n",l,sl->pa->co[sy][jKv][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Kh_sat of layer %ld [mm/s]: %f\n",l,sl->pa->co[sy][jKh][l]);
				}

				fprintf(f," Terrain elevation [m]: %f\n",top->Z0->co[r][c]);
				fprintf(f," Sky view factor [-]: %f\n",top->sky->co[r][c]);
				fprintf(f," The pixel-type is %d \n",top->pixel_type->co[r][c]);
				fprintf(f," Drainage Direction is %d \n",top->DD->co[r][c]);
				fprintf(f," Slope along Drainage Direction [-]: %f \n",top->i_DD->co[r][c]);
				fprintf(f," Aspect [deg] [0=Nord, clockwise]: %f \n",top->aspect->co[r][c]*180.0/Pi);
				fprintf(f," Mean slope of the pixel [deg]: %f \n",top->slopes->co[r][c]*180.0/Pi);
				fprintf(f," Land use number is %d \n",(short)land->LC->co[r][c]);

				for(l=1;l<=Nl;l++){
					fprintf(f," The root fraction [-] of layer %ld: %f\n",l,land->root_fraction->co[lu][l]);
				}
		
				fprintf(f," Surface fraction of land covered by vegetation [-]: %f \n",land->ty->co[lu][jcf]);
				fprintf(f," Leaf and Stem Area Index [-]: %f \n",land->ty->co[lu][jLSAI]);
				fprintf(f," Momentum roughness length z0soil [m]: %f \n",land->ty->co[lu][jz0]);
				fprintf(f," Vegetation height [m]: %f \n",land->ty->co[lu][jHveg]);

				fprintf(f," */ \n");
				t_fclose(f);
				free(temp2);
				free(temp);
				free(name);
			
				
				temp=join_strings(files->co[fpoint]+1,SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				for(j=1;j<=otot;j++){
					fprintf(f,",%s",header_point[j]);
				}
				fprintf(f,"\n");
				t_fclose(f);
				free(name);
				free(temp);
			}
			
			if(strcmp(files->co[fsnz]+1 , MISSING_FILE) != 0){
				temp=join_strings(files->co[fsnz]+1,SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");
				//fprintf(f,"/**Snow profile for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */\n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,Snowtype,Snowdepth[mm],SWE[mm],Temp_aver[C],Density_aver[kg/m3],melting[mm],sublimation[mm],evaporation[mm],nlayer,BStrans[mm],BStrans_cum[mm],BSsubl[mm],BSsubl_cum[mm],BStot[mm],BStot_cum[mm]");
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
				fprintf(f,"\n");	  
				t_fclose(f);
				free(name);
				free(temp);
			}
			
			if(par->glaclayer_max>0 && strcmp(files->co[fglz]+1 , MISSING_FILE) != 0){
				temp=join_strings(files->co[fglz]+1,SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");
				//fprintf(f,"/** Ice profile for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: \n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);	   
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,Glacdepth[mm],GWE[mm],Temp_aver[C],Num.glac_layer,Rhoav[kg/m3],albedo,Gl.melt[mm],Gl.subl[mm],Gl.evap[mm]");
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
				free(temp);
			}
		
			if(strcmp(files->co[fTz]+1 , MISSING_FILE) != 0){
				/*creation of the file "Tz.txt": */
				temp=join_strings(files->co[fTz]+1,SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");			
				//fprintf(f,"/** Profiles of sl temperature for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				z=0.0;
				for(l=1;l<=Nl;l++){
					z+=sl->pa->co[sy][jdz][l];
					fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
				}
				fprintf(f," \n");
				t_fclose(f);	
				free(name);
				free(temp);
			
				/*creation of the file "Tz.txt": averaged data*/
				temp=join_strings(files->co[fTz]+1,"MEAN");
				temp2=join_strings(temp,SSSS);
				name=join_strings(temp2,textfile);
				f=t_fopen(name,"w");		
				//fprintf(f,"/** Profiles of water content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				z=0.0;
				for(l=1;l<=Nl;l++){
					z+=sl->pa->co[sy][jdz][l];
					fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
				}	
				fprintf(f," \n");
				t_fclose(f);	
				free(name);
				free(temp2);
				free(temp);
			}
		
			if(strcmp(files->co[fpsiztot]+1 , MISSING_FILE) != 0){
				/*creation of the file "PSIz.txt": */
				temp=join_strings(files->co[fpsiztot]+1,SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				//fprintf(f,"/** Profiles of water pressure (mm) for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				z=0.0;
				for(l=1;l<=Nl;l++){
					z+=sl->pa->co[sy][jdz][l];
					fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
				}
				fprintf(f," \n");
				t_fclose(f);
				free(name);
				free(temp);
			}
			
			if(strcmp(files->co[fpsiz]+1 , MISSING_FILE) != 0){
				temp=join_strings(files->co[fpsiz]+1,SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				//fprintf(f,"/** Profiles of water pressure (mm) for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				z=0.0;
				for(l=1;l<=Nl;l++){
					z+=sl->pa->co[sy][jdz][l];
					fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
				}	
				fprintf(f," \n");
				t_fclose(f);
				free(name);
				free(temp);
			}
		
			if(strcmp(files->co[fliqz]+1 , MISSING_FILE) != 0){
				/*creation of the file "TETAz.txt": */
				temp=join_strings(files->co[fliqz]+1,SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				//fprintf(f,"/** Profiles of water content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				z=0.0;
				for(l=1;l<=Nl;l++){
					z+=sl->pa->co[sy][jdz][l];
					fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
				}	
				fprintf(f," \n");
				t_fclose(f);
				free(name);
				free(temp);			
			
				/*creation of the file "TETAz.txt": averaged data*/
				temp=join_strings(files->co[fliqz]+1,"MEAN");
				temp2=join_strings(temp,SSSS);
				name=join_strings(temp2,textfile);
				f=t_fopen(name,"w");		
				//fprintf(f,"/** Profiles of water content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				z=0.0;
				for(l=1;l<=Nl;l++){
					z+=sl->pa->co[sy][jdz][l];
					fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
				}	
				fprintf(f," \n");
				t_fclose(f);			
				free(name);
				free(temp2);
				free(temp);
			}
			
			if(strcmp(files->co[ficez]+1 , MISSING_FILE) != 0){
				/*creation of the file "TETAICEz.txt": */
				temp=join_strings(files->co[ficez]+1,SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");				
				//fprintf(f,"/** Profiles of ice content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld:\n",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				z=0.0;
				for(l=1;l<=Nl;l++){
					z+=sl->pa->co[sy][jdz][l];
					fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
				}
				fprintf(f," \n");
				t_fclose(f);	
				free(name);
				free(temp);
				
				/*creation of the file "TETAICEz.txt": averaged data*/
				temp=join_strings(files->co[ficez]+1,"MEAN");
				temp2=join_strings(temp,SSSS);
				name=join_strings(temp2,textfile);
				f=t_fopen(name,"w");		
				//fprintf(f,"/** Profiles of water content for the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld: */",par->chkpt->co[i][1],par->chkpt->co[i][2],r,c);
				fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
				z=0.0;
				for(l=1;l<=Nl;l++){
					z+=sl->pa->co[sy][jdz][l];
					fprintf(f,",%.0f ",z-0.5*sl->pa->co[sy][jdz][l]);
				}	
				fprintf(f," \n");
				t_fclose(f);	
				free(name);
				free(temp2);
				free(temp);
			}
			
		}
		
		free_doublevector(root_fraction);
		
	}
	
	if(par->state_basin==1){
	
		header_basin[ooprecrain] =  "Prain_below_canopy[mm]";
		header_basin[ooprecsnow] =  "Psnow_below_canopy[mm]";
		header_basin[oorainover] =  "Prain_above_canopy[mm]";
		header_basin[oosnowover] =  "Prain_above_canopy[mm]";
		header_basin[ooTa] =  "Tair[C]";
		header_basin[ooTg] =  "Tsurface[C]";
		header_basin[ooTv] =  "Tvegetation[C]";
		header_basin[ooevapsur] =  "Evap_surface[mm]";
		header_basin[ootrasp] =  "Transpiration_canopy[mm]";
		header_basin[ooLE] =  "LE[W/m2]";
		header_basin[ooH] =  "H[W/m2]";
		header_basin[ooSW] =  "SW[W/m2]";
		header_basin[ooLW] =  "LW[W/m2]";
		header_basin[ooLEv] =  "LEv[W/m2]";
		header_basin[ooHv] =  "Hv[W/m2]";
		header_basin[ooSWv] =  "SWv[W/m2]";
		header_basin[ooLWv] =  "LWv[W/m2]";
		header_basin[ooSWin] =  "SWin[W/m2]";
		header_basin[ooLWin] =  "LWin[W/m2]";
		header_basin[oomasserror] = "Mass_balance_error[mm]";
		
		//DATA BASIN
		if(strcmp(files->co[fbas]+1 , MISSING_FILE) != 0){
			temp=join_strings(files->co[fbas]+1,textfile);
			f=t_fopen(temp,"w");
			fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD");
			for(j=1;j<=ootot;j++){
				fprintf(f,",%s",header_basin[j]);
			}
			fprintf(f,"\n");
			t_fclose(f);
			free(temp);
		}
				
	}

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void write_soil_output(long n, long i, double t, double dt, long y0, double JD0, LONGMATRIX *rc, SOIL *sl, double psimin){

	char *name,*temp,*temp2,SSSS[ ]={"SSSS"};
	double t_i, JD;
	long d2, mo2, y2, h2, mi2, l, r=rc->co[i][1], c=rc->co[i][2];
	FILE *f;
	
	write_suffix(SSSS, i, 0);	
	t_i=t-dt*n;
	date_time(t+dt, y0, JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
	
	if(strcmp(files->co[fTz]+1 , MISSING_FILE) != 0){
		/*update of the sl profile temperature in the control pixel:*/
		temp=join_strings(files->co[fTz]+1,SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(t+dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",sl->T->co[l][r][c]);
		}
		fprintf(f," \n");
		fclose(f);
		free(name);
		free(temp);

		temp=join_strings(files->co[fTz]+1,"MEAN");
		temp2=join_strings(temp,SSSS);
		name=join_strings(temp2,textfile);
		f=fopen(name,"a");		
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(t+dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",sl->T_av->co[l][i]);
		}
		fprintf(f," \n");
		fclose(f);
		free(name);
		free(temp2);
		free(temp);
	}
	
	if(strcmp(files->co[fpsiztot]+1 , MISSING_FILE) != 0){
		temp=join_strings(files->co[fpsiztot]+1,SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(t+dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",sl->Ptot->co[l][r][c]);
		}
		fprintf(f," \n");
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files->co[fpsiz]+1 , MISSING_FILE) != 0){
		temp=join_strings(files->co[fpsiz]+1,SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(t+dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",sl->P->co[l][r][c]);
		}
		fprintf(f," \n");
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files->co[fliqz]+1 , MISSING_FILE) != 0){
		temp=join_strings(files->co[fliqz]+1,SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(t+dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",sl->th->co[l][r][c]);
		}
		fprintf(f," \n");
		fclose(f);
		free(name);
		free(temp);


		temp=join_strings(files->co[fliqz]+1,"MEAN");
		temp2=join_strings(temp,SSSS);
		name=join_strings(temp2,textfile);
		f=fopen(name,"a");		
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(t+dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",sl->th_av->co[l][i]);
		}
		fprintf(f," \n");
		fclose(f);
		free(name);
		free(temp2);
		free(temp);
	}
	
	if(strcmp(files->co[ficez]+1 , MISSING_FILE) != 0){
		temp=join_strings(files->co[ficez]+1,SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(t+dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",sl->thice->co[l][r][c]);
		}
		fprintf(f," \n");
		fclose(f);
		free(name);
		free(temp);
	
	
		temp=join_strings(files->co[ficez]+1,"MEAN");
		temp2=join_strings(temp,SSSS);
		name=join_strings(temp2,textfile);
		f=fopen(name,"a");		
		fprintf(f,"%ld/%ld/%ld %2.0f:%02.0f",d2,mo2,y2,(float)h2,(float)mi2);
		fprintf(f,",%f,%f,%f",(t+dt)/secinday,JD+(double)(daysfrom0(y2)),JD);   
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",sl->thice_av->co[l][i]);
		}
		fprintf(f," \n");
		fclose(f);
		free(name);
		free(temp2);
		free(temp);
	}
	
	for(l=1;l<=Nl;l++){
		sl->T_av->co[l][i] = 0.;
		sl->th_av->co[l][i] = 0.;
		sl->thice_av->co[l][i] = 0.;
	}
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
	
