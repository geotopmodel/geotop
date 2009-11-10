
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
#include "get_filenames.h"
#include "geomorphology.0875.h"
#include "pedo.funct.h"
#include "geo_statistic.09375.h"
#include "networks.h"
#include "constant.h"
#include "dtm_resolution.h"
#include "rw_maps.h"
#include "extensions.h"
#include "tabs.h"
#include "snow.09375.h"
#include "micromet.h"
#include "input.09375.h"
extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern char *error_file_name;
extern long Nl, Nr, Nc;
extern double NoV;

/*--------------------------------------------------------------------------------------------------*/
//! Subroutine which reads input data, performs  geomporphological analisys and allocates data
void get_all_input(int argc, char *argv[], TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet,
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times)

{

 /*counters of layers(l),rows(r),columns(c) and other things(i) (internal variables):*/
 long l,r,c,i,j,n_soiltypes,index,ncols,iland;
 /*auxiliary internal variables:*/
// static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0};/*differential of number-pixel for rows and*/
// static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0};/*columns, depending on Drainage Directions*/
 FILE *f;

 DOUBLEMATRIX *SNOW0;/*auxiliary internal vectors-variable*/
 DOUBLEMATRIX *GLACIER0=NULL;/*auxiliary internal vectors-variable*/ /* modified by Emanuele Cordano on 24 Sept 2009 */
 DOUBLEMATRIX *M;/*auxiliary internal vectors-variable*/
 //DOUBLEVECTOR *v;
 //Checking variables
 short sy;/*soil  type*/
 short a;
 double glac0, theta;
 INIT_TOOLS *IT;
 char *temp;

 IT=(INIT_TOOLS *)malloc(sizeof(INIT_TOOLS));
 if(!IT) t_error("IT was not allocated");

/****************************************************************************************************/
/*! Reading of program control par from geo_top.inpts file									*/
/*  (if it is in WORKING_DIRECTORY), otherwise from prompt        									*/
/****************************************************************************************************/
printf("\nGEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE\n");
printf("\nVersion 0.9375 Mackenzie\n");
printf("\nCopyright, 2008 Stefano Endrizzi, Riccardo Rigon, Matteo Dall'Amico and others following the footsteps of the FREE SOFTWARE FOUNDATION\n");

if(!argv[1]){
	WORKING_DIRECTORY=get_workingdirectory();
}else{
	WORKING_DIRECTORY=argv[1];
}

printf("\nWORKING DIRECTORY: %s",WORKING_DIRECTORY);

//READ INPUT PAR
//read_inpts_par(par, times, "__geotop", ".inpts", "2:");
//READ FILES
//files=read_filenames(WORKING_DIRECTORY, "__geotop", ".inpts", "1:");// old Mackenzie version
files=get_filenames_from_keys(WORKING_DIRECTORY,PROGRAM_NAME,PRINT); /* reads the keywords in __geotop.init and reads the filenames from __geotop.inpts */
//printf("files...=%s",files->element[I_CONTROL_PARAMETERS]+1); stop_execution();
read_inpts_par(par, times,files->element[I_CONTROL_PARAMETERS]+1,textfile, "1:"); /* reads the parameters in __control_parameters */

if(files->index->nh<nfiles){
	for(i=1;i<=files->index->nh;i++){
		printf("i:%ld file:%s\n",i,files->co[i]+1);
	}
	printf("files given:%ld files required%d\n",files->index->nh,nfiles);
	t_error("Wrong number of files in .inpts file");
}else if(files->index->nh>nfiles){
	printf("Warning: Excess files given:%ld files required%d\n",files->index->nh,nfiles);
}
for(i=1;i<=files->index->nh;i++){
	printf("File n:%ld %s\n",i,files->co[i]+1);
}
/* modidified by Emanuele Cordano on 23 September 2009 */
error_file_name=join_strings(files->co[ferr]+1,textfile);


f=t_fopen(error_file_name,"w");
fprintf(f,"Files of the errors in the simulation\n");
t_fclose(f);

/****************************************************************************************************/
/*! Reading of the Input files:                                                                      */
/****************************************************************************************************/
read_parameterfile(files->co[fpar]+1, par, IT);	/* reads the file __parameters.txt */
read_soil_parameters(files->co[fspar]+1, &n_soiltypes, sl);	/* reads the file _soil.txt and defines NL (number of soil layers) */

if(par->point_sim!=1){ /* distributed simulation */
	read_inputmaps(top, land, sl, par);
	read_optionsfile_distr(files->co[fopt]+1, par, times, land->LC);
}else{
	read_optionsfile_point(files->co[fopt]+1, par, top, land, sl, IT, times);
}


Nr=top->Z0->nrh;
Nc=top->Z0->nch;
NoV=UV->V->co[2];

//set INIZIAL TIME
times->time=0.0;

/****************************************************************************************************/
//Reading of RAIN data file,	METEO data file, 	and CLOUD data file

//METEO DATA
//meteo parameters

temp=join_strings(files->co[fmet]+1,textfile);/* reads the file meteo where all the parameters of the meteo stations are*/
f=t_fopen(temp,"r");
free(temp);
index=read_index(f,PRINT);
IT->met=read_doublematrix(f,"a",PRINT);
if(IT->met->nch!=13) t_error("Error in meteo parameters");
IT->met_col_names=read_stringarray(f,PRINT);
t_fclose(f);

//meteo_stations
met->st=(METEO_STATIONS *)malloc(sizeof(METEO_STATIONS));
if(!met->st) t_error("meteo_stations was not allocated");
init_meteo_stations(IT->met, met->st);// initializes the structure met->st with the elevation, aspect, longitude...of all the meteo stations

//meteo data
//met->column=alloc_long2(IT->met->nrh);/* allocate a matrix (n X 1) where n is the number of meteo stations */
met->data=(double ***)malloc(IT->met->nrh*sizeof(double**));
met->horizon=(double ***)malloc(IT->met->nrh*sizeof(double**));
met->var=alloc2(met->st->Z->nh,nmet);	/* allocates a matrix (n X nmet) where n is the # of meteo station and nmet is the number of meteorological variables as defined in constant.h */
met->column=(long **)malloc( (IT->met->nrh+1)*sizeof(long*) );
//meteo variables for the current instant
int numlines=0,ch_header;char ch1;
for(i=1;i<=IT->met->nrh;i++){
	temp=namefile_i(files->co[fmet]+1, i);// name of the meteo file: adds 0001
	f=t_fopen(temp,"r"); /* open the meteo file of each station */
	// find the number of lines of the meteo file
	ch1='\0';ch_header=0;
	while(ch1!='\n') {// to get rid of the header. In this case I am considering 1 line only of header and no comment (Matteo 28/9/09)
		ch1=fgetc(f);
		ch_header++;
	}
	numlines=count_meteo_lines(f);
	//printf("#station=%ld, name=%s, #of char in header=%d, numlines=%d",i,temp,ch_header-1,numlines); stop_execution();
	t_fclose(f);
	met->column[i-1]=alloc_long1(nmet); /* allocates a vector of "nmet" values of long to each meteo station */
	f=t_fopen(temp,"r"); /* open the meteo file of each station */
	ReadMeteoHeader(f, IT->met_col_names, met->st->offset->co[i], &ncols, met->column[i-1]);/* reads the header of the meteo file */
	met->data[i-1]=read_datameteo(f, met->st->offset->co[i], ncols, UV->V->co[2],numlines);
	met->horizon[i-1]=read_horizon(files->co[fhor]+1, i);
	free(temp);
	t_fclose(f);
}
//printf("met->column[0][iTsup]=%ld",met->column[0][iTsup]);stop_execution();
if (par->superfast==1 && met->column[0][iTsup]==-1) {
	t_error("The top Dirichlet boundary condition is not present, therefore you can't run the superfast version");
}
// superfast version
if(par->superfast==1){
	par->num_of_time=(long)(1+floor(3600.*times->TH/(times->n_pixel*par->Dt))); // number of times that the results have to be written in the super tensor. I have set 10 days
	sl->output=(double ***)malloc((par->num_of_time+1)*sizeof(double**));// one for the initial condition
	for(i=0;i<=par->num_of_time;i++){
		sl->output[i]=alloc2(9,Nl);
		/* matrix 9xNl that contains the value of the output variable for each layer
		 * row0: Tmin
		 * row1: Tmax
		 * row2: Tmean
		 * row3: ThetaWmin
		 * row4: ThetaWmax
		 * row5: ThetaWmean
		 * row6: ThetaImin
		 * row7: ThetaImax
		 * row8: ThetaImean */
	}
}
//read LAPSE RATES FILE  //0.JD 1.Year 2.Tair 3.Tdew 4.precipitation
met->LRp=new_longvector(3);	//column referring to the LRs
//year column 0 - JD column 1
met->LRp->co[1]=2;	//Column 2: Lapse rate for Tair (dT/dz in [C/m])
met->LRp->co[2]=3;	//Column 3: Lapse rate for Tdew (dT/dz in [C/m])
met->LRp->co[3]=4;	//Column 4: Lapse rate for prec (dP/dz in [mm/m])
if(existing_file_text(files->co[fLRs]+1)==1){// se esiste il file Lapse Rates
	temp=join_strings(files->co[fLRs]+1,textfile);
	f=t_fopen(temp,"r");
	met->LRs=read_datameteo(f, 1, 5, UV->V->co[2],count_meteo_lines(f));
	free(temp);
	t_fclose(f);
	par->LRflag=1;
}else{
	par->LRflag=0;
	printf("WARNING: LAPSE RATE FILE IS MISSING!!!\n");
	printf("BE AWARE OF THAT!!!! The default values will be used\n");
}
met->LRv=alloc1(5);// e' il vettore che contiene le variabili dei lapse rates interpolate per l'instante Dt di riferimento
met->LRv[met->LRp->co[1]]=NoV;
met->LRv[met->LRp->co[2]]=NoV;
met->LRv[met->LRp->co[3]]=NoV;


//FIND A STATION WITH SHORTWAVE RADIATION DATA
met->nstsrad=0;
do{
	met->nstsrad++;
	a=0;
	if( met->column[met->nstsrad-1][iSW]!=-1 || (met->column[met->nstsrad-1][iSWb]!=-1 && met->column[met->nstsrad-1][iSWd]!=-1) ) a=1;
}while(met->nstsrad<met->st->Z->nh && a==0);
if(a==0){
	printf("WARNING: NO shortwave radiation measurements available\n");
}else{
	printf("Shortwave radiation measurements from station %ld\n",met->nstsrad);
}

//FIND A STATION WITH CLOUD DATA
met->nstcloud=0;
do{
	met->nstcloud++;
	a=0;
	if( met->column[met->nstcloud-1][iC]!=-1 || met->column[met->nstcloud-1][itauC]!=-1) a=1;
}while(met->nstcloud<met->st->Z->nh && a==0);
if(a==0){
	printf("WARNING: NO cloudiness measurements available\n");
}else{
	printf("Cloudiness measurements from station %ld\n",met->nstcloud);
}

//FIND A STATION WITH LONGWAVE RADIATION DATA
met->nstlrad=0;
do{
	met->nstlrad++;
	a=0;
	if( met->column[met->nstlrad-1][iLWi]!=-1) a=1;
}while(met->nstlrad<met->st->Z->nh && a==0);
if(a==0){
	printf("WARNING: NO longwave radiation measurements available\n");
}else{
	printf("Longwave radiation measurements from station %ld\n",met->nstlrad);
}

/*if( (met->column[met->nstsrad-1][iSW]!=-1 || (met->column[met->nstsrad-1][iSWb]!=-1 && met->column[met->nstsrad-1][iSWd]!=-1))  ){
	if(met->column[met->nstsrad-1][iRh]==-1) printf("WARNING: No RH data for the radiation station\n");
	if(met->column[met->nstsrad-1][iT]==-1) printf("WARNING: No T data for the radiation station\n");
	if(met->column[met->nstsrad-1][iPt]==-1) printf("WARNING: No precipitation data for the radiation station\n");
}*/


/****************************************************************************************************/
/*! Completing the several time-indipendent input variables with the data of input-files:           */
/****************************************************************************************************/
/****************************************************************************************************/
// Completing of "land" (of the type LAND):

//Initialize matrix of shadow
land->shadow=new_shortmatrix(Nr,Nc);
initialize_shortmatrix(land->shadow,0);/* initialized as if it was always NOT in shadow*/

//Check that there aren't cell with an undefined land use value
for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if (land->LC->co[r][c]!=NoV){
			iland=0;
			for(i=1;i<=IT->land_classes->nch;i++){ if ((short)land->LC->co[r][c]==IT->land_classes->co[1][i]) iland++; }
			if(iland==0){
				printf("Cell r=%ld c=%ld has an undefined land class value (=%d)\n",r,c,(short)land->LC->co[r][c]);
				t_error("Fatal Error");
			}
		}
	}
}


par->n_landuses=0;// the maximum digit in the first row of block 1 in parameters file
for(i=1;i<=IT->land_classes->nch;i++){
	if(par->n_landuses<(long)IT->land_classes->co[1][i]) par->n_landuses=(long)IT->land_classes->co[1][i];
}
printf("\nNumber of land use types: %ld\n",par->n_landuses);

//properties for each land use
land->ty=new_doublematrix(par->n_landuses,nlandprop); /* land type (classes) */
initialize_doublematrix(land->ty,0.0);

for(i=1;i<=par->n_landuses;i++){

	for(j=1;j<=IT->land_classes->nch;j++){

		if(i==(long)IT->land_classes->co[1][j]){

			for(l=1;l<=nlandprop;l++){
				land->ty->co[i][l]=IT->land_classes->co[l+1][j];
			}

			//z0 (convert in m)
			land->ty->co[i][jz0]*=0.001;
			land->ty->co[i][jHveg]*=0.001;
			//land->ty->co[i][jz0veg]*=0.001;
			//land->ty->co[i][jd0]*=0.001;
			land->ty->co[i][jzb]*=0.001;

			for(l=1;l<=met->st->Z->nh;l++){
				if(land->ty->co[i][jHveg]>met->st->Vheight->co[l] || land->ty->co[i][jHveg]>met->st->Theight->co[l]){
					printf("hc:%f zmu:%f zmt:%f - set hc lower than measurement height - land cover %ld, meteo station %ld\n",
						land->ty->co[i][jHveg],met->st->Vheight->co[l],met->st->Theight->co[l],i,l);
					stop_execution();
					t_error("ERROR 1");
				}
			}
		}
	}
}


if(par->output_albedo>0){
	land->albedo=new_doublematrix(Nr,Nc);
	initialize_doublematrix(land->albedo,0.0);
}

land->clax=new_longvector(IT->land_classes->nch);
land->cont=new_longmatrix(IT->land_classes->nch,2);
initialize_longmatrix(land->cont,0);
for(i=1;i<=IT->land_classes->nch;i++){
	land->clax->co[i]=IT->land_classes->co[1][i];
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if (land->LC->co[r][c]!=UV->V->co[2]){
				if((short)land->LC->co[r][c]==land->clax->co[i]) land->cont->co[i][1]+=1;
			}
		}
	}
}


/****************************************************************************************************/
// Completing of "top" (of the type TOPO):

if(par->point_sim!=1){
	top->dz_dx=new_doublematrix(Nr,Nc);
	top->dz_dy=new_doublematrix(Nr,Nc);
	initialize_doublematrix(top->dz_dx,0.0);
	initialize_doublematrix(top->dz_dy,0.0);

	if(par->wat_balance==1){
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc-1;c++){
				if(land->LC->co[r][c]!=NoV){
					if(top->Z0->co[r][c+1]!=UV->V->co[2]){
						top->dz_dx->co[r][c]=(top->Z0dp->co[r][c]-top->Z0dp->co[r][c+1])/UV->U->co[1];
					}else{
						top->dz_dx->co[r][c]=0.0;
					}
				}
			}
		}
		for(c=1;c<=Nc;c++){
			for(r=1;r<=Nr-1;r++){
				if (land->LC->co[r][c]!=NoV){
					if (top->Z0->co[r+1][c]!=UV->V->co[2]){
						top->dz_dy->co[r][c]=(top->Z0dp->co[r][c]-top->Z0dp->co[r+1][c])/UV->U->co[2];
					}else{
						top->dz_dy->co[r][c]=0.0;
					}
				}
			}
		}
	}
}

/*It is found the total number of pixel which are not novalues:*/
par->total_pixel=0;
for(r=1;r<=Nr;r++){
    for(c=1;c<=Nc;c++){
       if (land->LC->co[r][c]!=NoV) par->total_pixel++;
    }
}
f=fopen(error_file_name,"a");
fprintf(f,"Valid pixels: %ld\n",par->total_pixel);
fprintf(f,"Novalue pixels: %ld\n",(Nr*Nc-par->total_pixel));
fprintf(f,"Basin area: %12.2f km2\n",(double)par->total_pixel*UV->U->co[1]*UV->U->co[2]/1.E6);

fclose(f);

//Cont for Richards 3D
top->i_cont=(long ***)malloc((Nl+1)*sizeof(long**));
for(l=1;l<=Nl;l++){
	top->i_cont[l]=(long **)malloc((Nr+1)*sizeof(long*));
	for(r=1;r<=Nr;r++){
		top->i_cont[l][r]=(long *)malloc((Nc+1)*sizeof(long));
	}
}
top->lrc_cont=new_longmatrix(Nl*par->total_pixel,3);
i_lrc_cont(land->LC, top->i_cont, top->lrc_cont);

/*! Modification of Dranaige Direction "top->DD" pixels lake (11), sea (12) or novalue (9) are put to the same value (0) */
for(r=1;r<=Nr;r++){
   for(c=1;c<=Nc;c++){
      if (land->LC->co[r][c]!=NoV){
         if (top->DD->co[r][c]>=9) top->DD->co[r][c]=9; /*novalue or errors*/
      }else{
         top->DD->co[r][c]=9;
      }
   }
}
/*! Channels pixels are put to the value (10) and other changes */
for(r=1;r<=Nr;r++){
   for(c=1;c<=Nc;c++){
      if (top->pixel_distance->co[r][c]>0.0 && top->pixel_distance->co[r][c]<UV->U->co[2]){
         top->pixel_distance->co[r][c]=UV->U->co[2]/2.0;
		 /*if(top->pixel_type->co[r][c]==0){
			top->DD->co[r][c]=0;
		 }else{
			top->DD->co[r][c]=10;
		 }*/
	  }
   }
}


/* Creating of the matrix of slope to calculate the surface velocity of the channel incoming flow ("top->i_ch"):*/
if(par->point_sim!=1 && par->wat_balance==1){
	top->i_ch=new_doublematrix(Nr,Nc);
	channel_lateral_slope(top->Z0dp,top->DD,UV,top->i_ch);
}else{
	top->i_ch=new_doublematrix(Nr,Nc);
	initialize_doublematrix(top->i_ch,0.0);
}


/****************************************************************************************************/
/*! Filling up of the struct "channel" (of the type CHANNEL):                                        */

/*The number of channel-pixel are counted:*/
i=0;
for(r=1;r<=Nr;r++){
   for(c=1;c<=Nc;c++){
      if (top->pixel_type->co[r][c]==10) i++;
   }
}
f=fopen(error_file_name,"a");
fprintf(f,"Channel pixels: %ld\n",i);
fclose(f);

if(par->point_sim==1){// point simulation: the distance from the outlet doesn't exist
	cnet->r=new_longvector(1);
	cnet->r->co[1]=1;
	cnet->c=new_longvector(1);
	cnet->c->co[1]=1;
	cnet->Q=new_doublevector(1);
	initialize_doublevector(cnet->Q,0.0);
	cnet->s0=new_doublevector(1);
	initialize_doublevector(cnet->s0,0.0);
	cnet->fraction_spread=new_doublematrix(1,1);
	cnet->fraction_spread->co[1][1]=1;
	cnet->Q_sup_s=new_doublevector(1);
	initialize_doublevector(cnet->Q_sup_s,0.0);
	cnet->Q_sub_s=new_doublevector(1);
	initialize_doublevector(cnet->Q_sub_s,0.0);
	cnet->Q_sup_spread=new_doublevector(1);
	initialize_doublevector(cnet->Q_sup_spread,0.0);
	cnet->Q_sub_spread=new_doublevector(1);
	initialize_doublevector(cnet->Q_sub_spread,0.0);
} else{// distributed simulation
/* Creation of the vectors with the position of the channel-pixels ("cnet->r" and "cnet->c"):*/
cnet->r=new_longvector(i);
cnet->c=new_longvector(i);
i=0;
for(r=1;r<=Nr;r++){
   for(c=1;c<=Nc;c++){
      if (top->pixel_type->co[r][c]==10){
         i++;
         cnet->r->co[i]=r;
         cnet->c->co[i]=c;
      }
   }
}

/* Initialization of the matrix with the sub/superficial flows which goes in a channel-cell ("cnet->q_sup"):*/
cnet->Q=new_doublevector(cnet->r->nh);

/* Extraction of the vector of channel-pixels distances ("cnet->s0") from the matrix of the distances of
   each pixel from the outlet ("top->pixel_distance")*/
cnet->s0=new_doublevector(cnet->r->nh);
for(i=1;i<=cnet->r->nh;i++){
   cnet->s0->co[i]=top->pixel_distance->co[cnet->r->co[i]][cnet->c->co[i]];
   if (cnet->s0->co[i]<=0.0){
		f=fopen(error_file_name,"a");
		fprintf(f,"Warning: negative distance from outlet of the channel pixel %ld, %ld\n",cnet->r->co[i],cnet->c->co[i]);
		fclose(f);
   }
}
if(par->print==1){
   f=fopen(error_file_name,"a");
   for(i=1;i<=cnet->r->nh;i++){
      fprintf(f,"cnet->s0->co[i] = %f con i = %ld r = %ld c = %ld\n",cnet->s0->co[i],i,cnet->r->co[i],cnet->c->co[i]);
   }
   fclose(f);
}

/* Creation of the matrix with the coefficient to spread the channel-flow for each channel-pixel ("cnet->fraction_spread"):*/

cnet->fraction_spread=De_Saint_Venant(cnet->s0,IT->u0,IT->D,par->Dt/(double)par->nDt_water);
if(par->print==1){
	temp=join_strings(WORKING_DIRECTORY,"ii_fraction_spread.txt");
   f=fopen(temp,"a");
   free(temp);
   fprintf(f,"/* Fraction spread calculated with De Saint Venant subroutine */ \n \n");
   for(r=1;r<=cnet->fraction_spread->nrh;r++){
      for(c=1;c<=cnet->fraction_spread->nch;c++){
         fprintf(f,"r=%4ld c=%4ld fraction_spread=%f ",r,c,cnet->fraction_spread->co[r][c]);
      }
      fprintf(f,"\n");
   }
   fclose(f);
}

/* Initialization of the vector with channel-flow (derived from q_sup) for each virtual channel-pixel with the
   same distance from outlet ("cnet->Q_sup_s"); note: vector's dimension is the number of virtual stretches
   of channel (starting from outlet):*/
cnet->Q_sup_s=new_doublevector(cnet->fraction_spread->nch);
initialize_doublevector(cnet->Q_sup_s,0.0);

cnet->Q_sup_spread=new_doublevector(cnet->fraction_spread->nch);

/* Initialization of the vector with channel-flow (derived from q_sub) for each virtual channel-pixel with the
   same distance from outlet ("cnet->Q_sub_s"); note: vector's dimension is the number of virtual stretches
   of channel (starting from outlet):*/
cnet->Q_sub_s=new_doublevector(cnet->fraction_spread->nch);
initialize_doublevector(cnet->Q_sub_s,0.0);

cnet->Q_sub_spread=new_doublevector(cnet->fraction_spread->nch);
}

if(par->recover==1){
	if(existing_file_text(files->co[rQch]+1)==1){
		temp=join_strings(files->co[rQch]+1,textfile);
		f=t_fopen(temp,"r");
		free(temp);
		index=read_index(f,PRINT);
		M=read_doublematrix(f,"a",PRINT);
		t_fclose(f);
		if(M->nrh!=cnet->fraction_spread->nch) t_error("The number of channel-pixels is not the same as the previous simulation");
		for(i=1;i<=M->nrh;i++){
			cnet->Q_sup_s->co[i]=M->co[i][1];
			cnet->Q_sub_s->co[i]=M->co[i][2];
		}
		free_doublematrix(M);
	}
}

/*M=new_doublematrix(Nr,Nc);
//distance_from_channel(M, top->DD, top->pixel_type);
distance_from_channel2(M, top->pixel_type, cnet->r, cnet->c);
for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			M->co[r][c]*=UV->U->co[1];
		}
	}
}
temp=join_strings(WORKING_DIRECTORY,"0dist_from_channel");
write_map(temp, 0, par->format_out, M, UV);
free(temp);
free_doublematrix(M);*/


/****************************************************************************************************/
/*! Completing of the initialization of SOIL structure                               */
/****************************************************************************************************/

/****************************************************************************************************/
/*! Initialization of the sl temperature tensor "T" [C]: */

/****************************************************************************************************/
/*! Completing of "sl" (of the type SOIL):                                                    */


sl->thice=new_doubletensor(Nl,Nr,Nc); /* soil theta_ice */
initialize_doubletensor(sl->thice,0.0);

sl->Jinf=new_doublematrix(Nr,Nc);/* infiltration */
initialize_doublematrix(sl->Jinf,0.0);

sl->J=new_doubletensor(Nl,Nr,Nc); /* water flux outgoing from one layer to the lower one */
initialize_doubletensor(sl->J,0.0);

sl->P=new_doubletensor(Nl,Nr,Nc); /* soil pressure */
initialize_doubletensor(sl->P,0.0);

sl->T=new_doubletensor(Nl,Nr,Nc); /* soil temperature */
initialize_doubletensor(sl->T,0.0);

sl->Tv=new_doublematrix(Nr,Nc);
initialize_doublematrix(sl->Tv,0.0);

sl->Tmean=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->Tmean,0.0);; // mean temperature in a layer for particular pixel in a Dt_output

sl->thetaw_mean=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->thetaw_mean,0.0);; // mean water content in a layer for particular pixel in a Dt_output

sl->thetai_mean=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->thetai_mean,0.0);; // mean ice content in a layer for particular pixel in a Dt_output

sl->psi_mean=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->psi_mean,0.0);; // mean water suction in a layer for particular pixel in a Dt_output

sl->Tmin=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->Tmin,99.0);; // min temperature in a layer for particular pixel in a Dt_output

sl->thetaw_min=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->thetaw_min,1.0);; // min water content in a layer for particular pixel in a Dt_output

sl->thetai_min=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->thetai_min,1.0);; // min ice content in a layer for particular pixel in a Dt_output

sl->Tmax=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->Tmax,-99.0);; // max temperature in a layer for particular pixel in a Dt_output

sl->thetaw_max=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->thetaw_max,-1.0);; // max water content in a layer for particular pixel in a Dt_output

sl->thetai_max=new_doublematrix(Nl,par->chkpt->nrh);
initialize_doublematrix(sl->thetai_max,-1.0);; // max ice content in a layer for particular pixel in a Dt_output



for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			for(i=1;i<=n_soiltypes;i++){
				if(sl->type->co[r][c]==i){
					sy=sl->type->co[r][c];

					//sl->Tv->co[r][c]=sl->pa->co[sy][jT][1];

					for(l=1;l<=Nl;l++){
						sl->P->co[l][r][c]=sl->pa->co[sy][jpsi][l];
						sl->T->co[l][r][c]=sl->pa->co[sy][jT][l];

						if(sl->T->co[l][r][c]<=Tfreezing){
							theta=teta_psi(sl->P->co[l][r][c], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
									1-1/sl->pa->co[sy][jns][l], par->psimin, par->Esoil);
							//printf("%ld P:%f thtot:%f T:%f\n",l,sl->P->co[l][r][c],theta,sl->T->co[l][r][c]);
							//Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
							sl->thice->co[l][r][c]=Fmin(theta, sl->pa->co[sy][jsat][l]) - teta_psi(Psif(sl->T->co[l][r][c]), 0.0, sl->pa->co[sy][jsat][l],
									sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin, par->Esoil);
							//printf("%ld thice:%f Psiunf:%f Thpsi:%f\n",l,sl->thice->co[l][r][c],Psif(sl->T->co[l][r][c]),teta_psi(Psif(sl->T->co[l][r][c]), 0.0, sl->pa->co[sy][jsat][l],
							//	sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin, par->Esoil));
							//if Theta(without freezing)<Theta_unfrozen(in equilibrium with T) Theta_ice is set at 0
							if(sl->thice->co[l][r][c]<0) sl->thice->co[l][r][c]=0.0;
							//Psi is updated taking into account the freezing
							theta-=sl->thice->co[l][r][c];
							sl->P->co[l][r][c]=psi_teta(theta, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
									sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], par->psimin, par->Esoil);
							//printf("%ld th:%f P:%f\n",l,theta,sl->P->co[l][r][c]);
						}
					}
				}
			}
		}else{

			for(l=1;l<=Nl;l++){
				sl->P->co[l][r][c]=NoV;
				sl->T->co[l][r][c]=NoV;
				sl->thice->co[l][r][c]=NoV;
			}

		}
	}
}


if(par->recover==1){
	assign_recovered(files->co[rTv]+1, sl->Tv->co, par, land->LC, IT->LU);
	for(l=1;l<=Nl;l++){
		temp=namefile_i_we(files->co[riceg]+1, l);
		assign_recovered(temp, sl->thice->co[l], par, land->LC, IT->LU);
		free(temp);
		temp=namefile_i_we(files->co[rpsi]+1, l);
		assign_recovered(temp, sl->P->co[l], par, land->LC, IT->LU);
		free(temp);
		temp=namefile_i_we(files->co[rTg]+1, l);
		assign_recovered(temp, sl->T->co[l], par, land->LC, IT->LU);
		free(temp);
	}
}

/****************************************************************************************************/
/*! Initialization of the struct "egy" (of the type ENERGY):*/

 egy->VSFA=new_doublevector(par->nLC);
 egy->HSFA=new_doublevector(par->nLC);

 egy->Hgrid=new_doublematrix(Nr,Nc);
 initialize_doublematrix(egy->Hgrid,0.0);
 egy->Tsgrid=new_doublematrix(Nr,Nc);
 initialize_doublematrix(egy->Tsgrid,0.0);
 if(par->output_Rn>0){
	egy->Rn_mean=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Rn_min=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Rn_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->LW_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->LW_min=new_doublematrix(Nr,Nc);
	egy->LW_in=new_doublematrix(Nr,Nc);
	egy->LW_out=new_doublematrix(Nr,Nc);
	egy->SW=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->SW_max=new_doublematrix(Nr,Nc);
 }
 if(par->output_ET>0){
	egy->ET_mean=new_doublematrix(Nr,Nc);
	//egy->ET_mean2=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->ET_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->ET_min=new_doublematrix(Nr,Nc);
 }
 if(par->output_H>0){
	egy->H_mean=new_doublematrix(Nr,Nc);
	//egy->H_mean2=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->H_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->H_min=new_doublematrix(Nr,Nc);
 }
 if(par->output_G>0){
	egy->G_mean=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->G_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->G_min=new_doublematrix(Nr,Nc);
	egy->G_snowsoil=new_doublematrix(Nr,Nc);
 }
 if(par->output_Ts>0){
	egy->Ts_mean=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Ts_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Ts_min=new_doublematrix(Nr,Nc);
 }
 if(par->output_Rswdown>0){
	egy->Rswdown_mean=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Rswdown_max=new_doublematrix(Nr,Nc);
	egy->Rswbeam=new_doublematrix(Nr,Nc);
 }

 egy->nDt_shadow=new_longmatrix(Nr,Nc);
 egy->nDt_sun=new_longmatrix(Nr,Nc);
 initialize_longmatrix(egy->nDt_shadow,0);
 initialize_longmatrix(egy->nDt_sun,0);

 if(par->output_meteo>0){
	egy->Ta_mean=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Ta_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Ta_min=new_doublematrix(Nr,Nc);
 }

 for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){

			if(par->output_Rn>0){
				egy->Rn_mean->co[r][c]=0.0;
				if(par->distr_stat==1)egy->Rn_max->co[r][c]=-9.0E99;
				if(par->distr_stat==1)egy->Rn_min->co[r][c]=+9.0E99;
				if(par->distr_stat==1)egy->LW_min->co[r][c]=+9.0E99;
				if(par->distr_stat==1)egy->LW_max->co[r][c]=-9.0E99;
				egy->LW_in->co[r][c]=0.0;
				egy->LW_out->co[r][c]=0.0;
				egy->SW->co[r][c]=0.0;
				if(par->distr_stat==1)egy->SW_max->co[r][c]=-9.0E99;
			}
			if(par->output_ET>0){
				egy->ET_mean->co[r][c]=0.0;
				//egy->ET_mean2->co[r][c]=0.0;
				if(par->distr_stat==1)egy->ET_max->co[r][c]=-9.0E99;
				if(par->distr_stat==1)egy->ET_min->co[r][c]=+9.0E99;
			}
			if(par->output_H>0){
				egy->H_mean->co[r][c]=0.0;
				//egy->H_mean2->co[r][c]=0.0;
				if(par->distr_stat==1)egy->H_max->co[r][c]=-9.0E99;
				if(par->distr_stat==1)egy->H_min->co[r][c]=+9.0E99;
			}
			if(par->output_G>0){
				egy->G_mean->co[r][c]=0.0;
				if(par->distr_stat==1)egy->G_max->co[r][c]=-9.0E99;
				if(par->distr_stat==1)egy->G_min->co[r][c]=+9.0E99;
				egy->G_snowsoil->co[r][c]=0.0;
			}
			if(par->output_Ts>0){
				egy->Ts_mean->co[r][c]=0.0;
				if(par->distr_stat==1)egy->Ts_max->co[r][c]=-9.0E99;
				if(par->distr_stat==1)egy->Ts_min->co[r][c]=+9.0E99;
			}
			if(par->output_Rswdown>0){
				egy->Rswdown_mean->co[r][c]=0.0;
				if(par->distr_stat==1)egy->Rswdown_max->co[r][c]=-9.0E99;
				egy->Rswbeam->co[r][c]=0.0;
			}
			if(par->output_meteo>0){
				egy->Ta_mean->co[r][c]=0.0;
				if(par->distr_stat==1)egy->Ta_max->co[r][c]=-9.0E99;
				if(par->distr_stat==1)egy->Ta_min->co[r][c]=+9.0E99;
			}

		}else{

			if(par->output_Rn>0){
				egy->Rn_mean->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->Rn_max->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->Rn_min->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->LW_min->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->LW_max->co[r][c]=UV->V->co[2];
				egy->LW_in->co[r][c]=UV->V->co[2];
				egy->LW_out->co[r][c]=UV->V->co[2];
				egy->SW->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->SW_max->co[r][c]=UV->V->co[2];
			}
			if(par->output_ET>0){
				egy->ET_mean->co[r][c]=UV->V->co[2];
				//egy->ET_mean2->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->ET_max->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->ET_min->co[r][c]=UV->V->co[2];
			}
			if(par->output_H>0){
				egy->H_mean->co[r][c]=UV->V->co[2];
				//egy->H_mean2->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->H_max->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->H_min->co[r][c]=UV->V->co[2];
			}
			if(par->output_G>0){
				egy->G_mean->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->G_max->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->G_min->co[r][c]=UV->V->co[2];
				egy->G_snowsoil->co[r][c]=UV->V->co[2];
			}
			if(par->output_Ts>0){
				egy->Ts_mean->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->Ts_max->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->Ts_min->co[r][c]=UV->V->co[2];
			}
			if(par->output_Rswdown>0){
				egy->Rswdown_mean->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->Rswdown_max->co[r][c]=UV->V->co[2];
				egy->Rswbeam->co[r][c]=UV->V->co[2];
			}
			if(par->output_meteo>0){
				egy->Ta_mean->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->Ta_max->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->Ta_min->co[r][c]=UV->V->co[2];
			}
		}
	}
 }

 egy->out1=new_doublematrix(84,par->chkpt->nrh);
 initialize_doublematrix(egy->out1,0.0);

 egy->out2=new_doublevector(15);
 initialize_doublevector(egy->out2,0.0);

 if(par->ES_num>0){
	egy->out3=new_doublematrix(23,par->ES_num);
	initialize_doublematrix(egy->out3,0.0);
 }

 if(par->JD_plots->co[1]>=0){
	egy->Hgplot=new_doublematrix(Nr,Nc);
	egy->LEgplot=new_doublematrix(Nr,Nc);
	egy->Hvplot=new_doublematrix(Nr,Nc);
	egy->LEvplot=new_doublematrix(Nr,Nc);

	egy->SWinplot=new_doublematrix(Nr,Nc);
	egy->SWgplot=new_doublematrix(Nr,Nc);
	egy->SWvplot=new_doublematrix(Nr,Nc);

	egy->LWinplot=new_doublematrix(Nr,Nc);
	egy->LWgplot=new_doublematrix(Nr,Nc);
	egy->LWvplot=new_doublematrix(Nr,Nc);

	egy->Tsplot=new_doublematrix(Nr,Nc);
	egy->Tgplot=new_doublematrix(Nr,Nc);
	egy->Tvplot=new_doublematrix(Nr,Nc);

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(top->Z0->co[r][c]==UV->V->co[2]){
				egy->Hgplot->co[r][c]=UV->V->co[2];
				egy->LEgplot->co[r][c]=UV->V->co[2];
				egy->Hvplot->co[r][c]=UV->V->co[2];
				egy->LEvplot->co[r][c]=UV->V->co[2];
				egy->SWinplot->co[r][c]=UV->V->co[2];
				egy->SWgplot->co[r][c]=UV->V->co[2];
				egy->SWvplot->co[r][c]=UV->V->co[2];
				egy->LWinplot->co[r][c]=UV->V->co[2];
				egy->LWgplot->co[r][c]=UV->V->co[2];
				egy->LWvplot->co[r][c]=UV->V->co[2];
				egy->Tsplot->co[r][c]=UV->V->co[2];
				egy->Tgplot->co[r][c]=UV->V->co[2];
				egy->Tvplot->co[r][c]=UV->V->co[2];
			}else{
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
			}
		}
	}
}

for(i=1;i<=par->nLC;i++){
	egy->VSFA->co[i]=1.0;
	egy->HSFA->co[i]=0.0;
	/*if(par->recover==1){
		f=t_fopen(namefile_i(join_strings(files->co[rSFA]+1,"L"),i),"r");
		//f=t_fopen(join_strings(files->co[rSFA]+1,textfile),"r");
		index=read_index(f,PRINT);
		v=read_doublearray(f,PRINT);
		t_fclose(f);
		egy->VSFA->co[i]=v->co[1];
		egy->HSFA->co[i]=v->co[2];
		free_doublevector(v);
	}*/
}

/****************************************************************************************************/
/*! Completing of the struct "water" (of the type WATER) with the initializations of the remanent
    matrices (wat->rain, wat->Pn, wat->wcan_rain):        */

/* Initialization of wat->Pn (liquid precipitation that reaches the sl surface in mm):*/
wat->Pn=new_doublematrix(Nr,Nc);
initialize_doublematrix(wat->Pn,0.0);

/* Initialization of wat->wcan_rain: (liquid precipitation intercepted by vegetation in mm):*/
wat->wcan_rain=new_doublematrix(Nr,Nc);
wat->wcan_snow=new_doublematrix(Nr,Nc);

/* Initialization of wat->h_sup: (height of water over the sl-surface not infiltrated in mm):*/
wat->h_sup=new_doublematrix(Nr,Nc);
wat->q_sup=new_doublematrix(Nr,Nc);
wat->q_sub=new_doubletensor(Nl,Nr,Nc);

wat->error=new_doublematrix(Nr,Nc);


/* Initialization of wat->total (total precipitation (rain+snow) precipitation):*/
if(par->point_sim==1 && par->micromet==1){
	wat->total=new_doublematrix(top->Z1->nrh,top->Z1->nch);
}else{
	wat->total=new_doublematrix(Nr,Nc);
}
initialize_doublematrix(wat->total,0.0);

/* Initialization of the matrices with the output of total precipitation and interception:*/
if(par->output_P>0){
	wat->PrTOT_mean=new_doublematrix(Nr,Nc);
	wat->PrSNW_mean=new_doublematrix(Nr,Nc);
}

if(par->output_h_sup>0) wat->hsupav=new_doublematrix(Nr,Nc);

for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			wat->error->co[r][c]=0.0;
			wat->wcan_rain->co[r][c]=0.0;
			wat->wcan_snow->co[r][c]=0.0;
			wat->h_sup->co[r][c]=0.0;
			if(par->output_P>0){
				wat->PrTOT_mean->co[r][c]=0.0;
				wat->PrSNW_mean->co[r][c]=0.0;
			}
			if(par->output_h_sup>0) wat->hsupav->co[r][c]=0.0;
		}else{
			wat->error->co[r][c]=NoV;
			wat->wcan_rain->co[r][c]=UV->V->co[2];
			wat->wcan_snow->co[r][c]=UV->V->co[2];
			wat->h_sup->co[r][c]=UV->V->co[2];
			if(par->output_P>0){
				wat->PrTOT_mean->co[r][c]=UV->V->co[2];
				wat->PrSNW_mean->co[r][c]=UV->V->co[2];
			}
			if(par->output_h_sup>0) wat->hsupav->co[r][c]=UV->V->co[2];
		}
	}
}

if(par->recover==1){
	assign_recovered(files->co[rhsup]+1, wat->h_sup->co, par, land->LC, IT->LU);
	assign_recovered(files->co[rwcrn]+1, wat->wcan_rain->co, par, land->LC, IT->LU);
	assign_recovered(files->co[rwcsn]+1, wat->wcan_rain->co, par, land->LC, IT->LU);
}


/* Initialization of vector "wat->output" with basin's means and pixel's values for the output:*/
wat->out1=new_doublematrix(28,par->chkpt->nrh);
initialize_doublematrix(wat->out1,0.0);
wat->out2=new_doublevector(8);
initialize_doublevector(wat->out2,0.0);
wat->outfluxes=new_doublematrix(6,IT->land_classes->nch);
initialize_doublematrix(wat->outfluxes,0.0);


/* Creation of Kriging weights to model the rainfall distribution ("wat->weights_Kriging"):*/
wat->weights_Kriging=new_doublematrix((Nr*Nc), met->st->Z->nh);
initialize_doublematrix(wat->weights_Kriging, 0.999999);
/* Call of the function ordi_kriging in geo_statistic.c:
	Output: 1) matrix with kriging weights (DOUBLEMATRIX [NR*NC,number of station] - wat->weights_Kriging)
	Input:  2) matrix with the coordinates of rain stations (par->rain_stations)
			3) matrix with elevations (top->Z0)
			4) vector with pixel dimension (UV->U)
			5) vector with turtle file header (UV->V)
			6) integral spatial scale (scala_integr)
			7) spatial variance (varianza) */
ordi_kriging2(wat->weights_Kriging, met->st->E, met->st->N, top->Z0, UV, par->integr_scale_rain, par->variance_rain);
//doublematrix_dem3(wat->weights_Kriging,UV->U,UV->V,"ii_kriging_weights.txt","Doublematrix of krigingweights",par->print);




/****************************************************************************************************/
/*! Initialization of the struct "snow" (of the type SNOW):*/

/***************************************************************************************************/
/*! Optional reading of initial real snow thickness map in the whole basin ("SNOW0"):    */
if(existing_file(files->co[fsn0]+1)>0){
	printf("Snow initial condition from file %s\n",files->co[fsn0]+1);
	SNOW0=read_map(2, files->co[fsn0]+1, land->LC, UV);
}else{
	SNOW0=copydoublematrix_const(IT->Dsnow0, land->LC, UV->V->co[2]);
}

/*! Optional reading of snow age in the whole basin     */
if(existing_file(files->co[fsnag0]+1)>0){
	printf("Snow age initial condition from file %s\n",files->co[fsnag0]+1);
	snow->dimens_age=read_map(2, files->co[fsnag0]+1, land->LC, UV);
}else{
    snow->dimens_age=copydoublematrix_const(IT->agesnow0, land->LC, UV->V->co[2]);
}

if(par->JD_plots->co[1]>=0){
	snow->Dplot=new_doublematrix(Nr,Nc);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]==NoV){
				snow->Dplot->co[r][c]=UV->V->co[2];
			}else{
				snow->Dplot->co[r][c]=0.0;
			}
		}
	}
}

snow->Psnow=new_doublematrix(Nr,Nc);

snow->CR1=new_doublevector(par->snowlayer_max);
snow->CR2=new_doublevector(par->snowlayer_max);
snow->CR3=new_doublevector(par->snowlayer_max);
snow->CR1m=new_doublevector(par->snowlayer_max);
snow->CR2m=new_doublevector(par->snowlayer_max);
snow->CR3m=new_doublevector(par->snowlayer_max);
initialize_doublevector(snow->CR1m,0.0);
initialize_doublevector(snow->CR2m,0.0);
initialize_doublevector(snow->CR3m,0.0);

snow->type=new_shortmatrix(Nr,Nc);
initialize_shortmatrix(snow->type,0);

snow->lnum=new_longmatrix(Nr,Nc);
initialize_longmatrix(snow->lnum,0);

snow->Dzl=new_doubletensor(par->snowlayer_max,Nr,Nc);   /*in mm*/
initialize_doubletensor(snow->Dzl,0.0);

snow->w_liq=new_doubletensor(par->snowlayer_max,Nr,Nc);
initialize_doubletensor(snow->w_liq,0.0);

snow->w_ice=new_doubletensor(par->snowlayer_max,Nr,Nc);
initialize_doubletensor(snow->w_ice,0.0);

snow->T=new_doubletensor(par->snowlayer_max,Nr,Nc);
initialize_doubletensor(snow->T,-99.0);

snow->rho_newsnow=new_doublematrix(Nr,Nc);
initialize_doublematrix(snow->rho_newsnow,UV->V->co[2]);

snow->nondimens_age=new_doublematrix(Nr,Nc);

if(par->blowing_snow==1){
	snow->Wtrans=new_doublematrix(Nr,Nc);
	snow->Qsub=new_doublematrix(Nr,Nc);
	snow->Qtrans=new_doublematrix(Nr,Nc);
	snow->Qtrans_x=new_doublematrix(Nr,Nc);
	snow->Qtrans_y=new_doublematrix(Nr,Nc);
	//snow->ListonSWE=new_doublematrix(Nr,Nc);
	//snow->softSWE=new_doublematrix(Nr,Nc);
	//initialize_doublematrix(snow->softSWE,UV->V->co[2]);
	//snow->softSWE1=new_doublematrix(Nr,Nc);
	if(par->output_snow>0){
		snow->Wtot=new_doublematrix(Nr,Nc);
		/*snow->Wtrans_cum=new_doublematrix(Nr,Nc);
		snow->Wsusp_cum=new_doublematrix(Nr,Nc);
		snow->Wsubl_cum=new_doublematrix(Nr,Nc);
		snow->Wsubgrid_cum=new_doublematrix(Nr,Nc);*/
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]==UV->V->co[2]){
					snow->Wtot->co[r][c]=UV->V->co[2];
					/*snow->Wtrans_cum->co[r][c]=UV->V->co[2];
					snow->Wsusp_cum->co[r][c]=UV->V->co[2];
					snow->Wsubl_cum->co[r][c]=UV->V->co[2];
					snow->Wsubgrid_cum->co[r][c]=UV->V->co[2];*/
				}else{
					snow->Wtot->co[r][c]=0.0;
					/*snow->Wtrans_cum->co[r][c]=0.0;
					snow->Wsusp_cum->co[r][c]=0.0;
					snow->Wsubl_cum->co[r][c]=0.0;
					snow->Wsubgrid_cum->co[r][c]=0.0;*/
				}
			}
		}
	}
}
snow->out_bs=new_doublematrix(10,par->chkpt->nrh);
initialize_doublematrix(snow->out_bs,0.0);
snow->evap=new_doublevector(par->chkpt->nrh);
initialize_doublevector(snow->evap,0.0);
snow->subl=new_doublevector(par->chkpt->nrh);
initialize_doublevector(snow->subl,0.0);
snow->melted=new_doublevector(par->chkpt->nrh);
initialize_doublevector(snow->melted,0.0);

if(par->output_balancesn>0){
	snow->MELTED=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->MELTED,NoV);
	snow->SUBL=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->SUBL,NoV);
	snow->t_snow=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->t_snow,NoV);
	snow->totav_snow=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->totav_snow,NoV);
}

if(par->output_snow>0){
	snow->max=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->max,NoV);
	snow->average=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->average,NoV);
}

if(par->blowing_snow==1){
	if(Nr>Nc){
		snow->change_dir_wind=new_longvector(Nr);
	}else{
		snow->change_dir_wind=new_longvector(Nc);
	}
}

for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){

		if(land->LC->co[r][c]!=NoV){

			if(SNOW0->co[r][c]<0){
				t_error("Error: negative snow data");
			}else if(SNOW0->co[r][c]<0.00001){
				SNOW0->co[r][c]=0.0;
				snow->type->co[r][c]=0;
				snow->dimens_age->co[r][c]=0.0;

			}else if(SNOW0->co[r][c]<par->Dmin->co[1]){
				snow->lnum->co[r][c]=1;
				snow->Dzl->co[1][r][c]=SNOW0->co[r][c];
				snow->type->co[r][c]=1;
				snow->dimens_age->co[r][c]*=86400.0;	//now in [s]

			}else{
				//just a first guess
				snow->type->co[r][c]=2;
				snow->lnum->co[r][c]=par->snowlayer_max;
				for(l=1;l<=par->snowlayer_max;l++){
					snow->Dzl->co[l][r][c]=SNOW0->co[r][c]/par->snowlayer_max;
				}
				snow->dimens_age->co[r][c]*=86400.0;	//now in [s]
			}

			if(snow->lnum->co[r][c]>0){
				for(l=1;l<=snow->lnum->co[r][c];l++){
					snow->w_ice->co[l][r][c]=IT->rhosnow0*0.001*snow->Dzl->co[l][r][c];
					snow->T->co[l][r][c]=-15.0;
					snow->w_liq->co[l][r][c]=0.0;
				}
			}

			if(par->output_snow>0){
				snow->max->co[r][c]=0.0;
				snow->average->co[r][c]=0.0;
			}

			if(par->output_balancesn>0){
				snow->MELTED->co[r][c]=0.0;
				snow->SUBL->co[r][c]=0.0;
				snow->t_snow->co[r][c]=0.0;
				snow->totav_snow->co[r][c]=0.0;
			}

			snow->nondimens_age->co[r][c]=snow->dimens_age->co[r][c];
			non_dimensionalize_snowage(&(snow->nondimens_age->co[r][c]), Tfreezing);

		}
	}
}


if(par->recover==1){

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			snow->type->co[r][c]=2;
		}
	}

	assign_recovered_long(files->co[rns]+1, snow->lnum->co, par, land->LC, IT->LU);
	assign_recovered(files->co[rsnag_adim]+1, snow->nondimens_age->co, par, land->LC, IT->LU);
	assign_recovered(files->co[rsnag_dim]+1, snow->dimens_age->co, par, land->LC, IT->LU);

    for(l=1;l<=par->snowlayer_max;l++){
    	temp=namefile_i_we(files->co[rDzs]+1, l);
		assign_recovered(temp, snow->Dzl->co[l], par, land->LC, IT->LU);
		free(temp);
		temp=namefile_i_we(files->co[rwls]+1, l);
		assign_recovered(temp, snow->w_liq->co[l], par, land->LC, IT->LU);
		free(temp);
		temp=namefile_i_we(files->co[rwis]+1, l);
		assign_recovered(temp, snow->w_ice->co[l], par, land->LC, IT->LU);
		free(temp);
		temp=namefile_i_we(files->co[rTs]+1, l);
		assign_recovered(temp, snow->T->co[l], par, land->LC, IT->LU);
		free(temp);
	}

}




/****************************************************************************************************/
/*! Initialization of the struct "glac" (of the type GLACIER):*/

/***************************************************************************************************/
/*! Optional reading of glacier depth in the whole basin ("GLACIER0"):    */
if(par->point_sim!=1 && existing_file(files->co[fgl0]+1)>0){
	if(par->glaclayer_max==0){
		printf("Warning: Glacier map present, but glacier represented with 0 layers\n");
		stop_execution();
	}
}

if(par->glaclayer_max>0){
	if(par->point_sim!=1 && existing_file(files->co[fgl0]+1)>0){
		printf("Glacier initial condition from file %s\n",files->co[fgl0]+1);
		GLACIER0=read_map(2, files->co[fgl0]+1, land->LC, UV);
	}else{
		GLACIER0=copydoublematrix_const(IT->Dglac0, land->LC, UV->V->co[2]);
	}

}else{
	//check
	if(IT->Dglac0>0){
		f=fopen(error_file_name,"a");
		printf("\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
		fprintf(f,"\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
		fclose(f);
	}
}

glac->evap=new_doublevector(par->chkpt->nrh);
glac->subl=new_doublevector(par->chkpt->nrh);
glac->melted=new_doublevector(par->chkpt->nrh);
initialize_doublevector(glac->evap,0.0);
initialize_doublevector(glac->subl,0.0);
initialize_doublevector(glac->melted,0.0);

//If the max number of glacier layers is greater than 1, the matrices (or tensors) lnum, Dzl. w_liq, w_ice, T and print matrices are defined, according to the respective flags
if(par->glaclayer_max>0){

	glac->lnum=new_longmatrix(Nr,Nc);
	initialize_longmatrix(glac->lnum,0);

	glac->Dzl=new_doubletensor(par->glaclayer_max,Nr,Nc);   /*in mm*/
	initialize_doubletensor(glac->Dzl,0.0);

	glac->w_liq=new_doubletensor(par->glaclayer_max,Nr,Nc);
	initialize_doubletensor(glac->w_liq,0.0);

	glac->w_ice=new_doubletensor(par->glaclayer_max,Nr,Nc);
	initialize_doubletensor(glac->w_ice,0.0);

	glac->T=new_doubletensor(par->glaclayer_max,Nr,Nc);
	initialize_doubletensor(glac->T,-99.0);

	if(par->output_balancegl>0){
		glac->MELTED=new_doublematrix(Nr,Nc);
		initialize_doublematrix(glac->MELTED,NoV);
		glac->SUBL=new_doublematrix(Nr,Nc);
		initialize_doublematrix(glac->SUBL,NoV);
	}

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){

				if(GLACIER0->co[r][c]<0){
					t_error("Error: negative glacier data");

				}else if(GLACIER0->co[r][c]<10){
					GLACIER0->co[r][c]=0.0;
					glac->lnum->co[r][c]=0;
					for(l=1;l<=par->glaclayer_max;l++){
						glac->Dzl->co[l][r][c]=0.0;
						glac->w_ice->co[l][r][c]=0.0;
						GLACIER0->co[r][c]=0.0;
					}

				}else{
					//just a first guess
					glac->lnum->co[r][c]=par->glaclayer_max;
					for(l=1;l<=par->glaclayer_max;l++){
						glac->Dzl->co[l][r][c]=GLACIER0->co[r][c]/par->glaclayer_max;
						glac->w_ice->co[l][r][c]=IT->rhoglac0*0.001*glac->Dzl->co[l][r][c];
						glac->T->co[l][r][c]=IT->Tglac0;
						glac->w_liq->co[l][r][c]=0.0;
					}
				}

				if(par->output_balancegl>0){
					glac->MELTED->co[r][c]=0.0;
					glac->SUBL->co[r][c]=0.0;
				}
			}
		}
	}

	/*if(par->recover==1){

		assign_recovered_long(files->co[rni]+1, glac->lnum->co, par, land->LC, IT->LU);

		for(l=1;l<=par->glaclayer_max;l++){

			assign_recovered(namefile_i_we(files->co[rDzi]+1, l),  glac->Dzl->co[l], par, land->LC, IT->LU);
			assign_recovered(namefile_i_we(files->co[rwli]+1, l),  glac->w_liq->co[l], par, land->LC, IT->LU);
			assign_recovered(namefile_i_we(files->co[rwii]+1, l),  glac->w_ice->co[l], par, land->LC, IT->LU);
			assign_recovered(namefile_i_we(files->co[rTi]+1, l),  glac->T->co[l], par, land->LC, IT->LU);

		}
	}*/
}



//*ALTIMETRIC RANKS******************************************************
if(par->ES_num>0){
	top->ES_pixel=new_longvector(par->ES_num);
	initialize_longvector(top->ES_pixel,0);
	top->ES_aspect=new_doublevector(par->ES_num);
	initialize_doublevector(top->ES_aspect,0.0);
	top->ES_slope=new_doublevector(par->ES_num);
	initialize_doublevector(top->ES_slope,0.0);
	top->Zmin=10000.0;
	top->Zmax=-1000.0;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				if(par->glaclayer_max>0){
					glac0=GLACIER0->co[r][c];
				}else{
					glac0=0.0;
				}
				if(glac0>=par->glac_thr){
					if(top->Z0->co[r][c]<top->Zmin) top->Zmin=top->Z0->co[r][c];
					if(top->Z0->co[r][c]>top->Zmax) top->Zmax=top->Z0->co[r][c];
				}
			}
		}
	}

	for(i=1;i<=par->ES_num;i++){
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					if(par->glaclayer_max>0){
						glac0=GLACIER0->co[r][c];
					}else{
						glac0=0.0;
					}
					if(glac0>=par->glac_thr){
						if( (top->Z0->co[r][c]>=top->Zmin+(i-1)*(top->Zmax-top->Zmin)/(double)par->ES_num &&
							top->Z0->co[r][c]<top->Zmin+i*(top->Zmax-top->Zmin)/(double)par->ES_num) ||
							(top->Z0->co[r][c]==top->Zmin+i*(top->Zmax-top->Zmin)/(double)par->ES_num &&
							i==par->ES_num) ){
							top->ES_pixel->co[i]+=1;
							top->ES_aspect->co[i]+=top->aspect->co[r][c]*180/Pi;
							top->ES_slope->co[i]+=top->slopes->co[r][c]*180/Pi;
						}
					}
				}
			}
		}
	}
	for(i=1;i<=par->ES_num;i++){
		top->ES_aspect->co[i]/=(double)top->ES_pixel->co[i];
		top->ES_slope->co[i]/=(double)top->ES_pixel->co[i];
	}
}



//***************************************************************************************************
// Filling up of the struct "met" (of the type METEO):

if(par->point_sim==1 && par->micromet==1){
	met->Tgrid=new_doublematrix(top->Z1->nrh,top->Z1->nch);
	met->Pgrid=new_doublematrix(top->Z1->nrh,top->Z1->nch);
}else{
	met->Tgrid=new_doublematrix(Nr,Nc);
	met->Pgrid=new_doublematrix(Nr,Nc);
}

initialize_doublematrix(met->Tgrid, 0.0);

if(par->micromet==1){
	if(par->point_sim==1){
		met->Vgrid=new_doublematrix(top->Z1->nrh,top->Z1->nch);
		met->Vdir=new_doublematrix(top->Z1->nrh,top->Z1->nch);
		met->RHgrid=new_doublematrix(top->Z1->nrh,top->Z1->nch);
	}else{
		met->Vgrid=new_doublematrix(Nr,Nc);
		met->Vdir=new_doublematrix(Nr,Nc);
		met->RHgrid=new_doublematrix(Nr,Nc);
	}
	if(par->output_meteo>0){
		met->Vspdmean=new_doublematrix(Nr,Nc);
		met->Vdirmean=new_doublematrix(Nr,Nc);
		met->RHmean=new_doublematrix(Nr,Nc);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]==UV->V->co[2]){
					met->Vspdmean->co[r][c]=UV->V->co[2];
					met->Vdirmean->co[r][c]=UV->V->co[2];
					met->RHmean->co[r][c]=UV->V->co[2];
				}else{
					met->Vspdmean->co[r][c]=0.0;
					met->Vdirmean->co[r][c]=0.0;
					met->RHmean->co[r][c]=0.0;
				}
			}
		}
	}
}

//plot output
if(par->JD_plots->co[1]>=0){
	met->Taplot=new_doublematrix(Nr,Nc);
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]==UV->V->co[2]){
				met->Taplot->co[r][c]=UV->V->co[2];
			}else{
				met->Taplot->co[r][c]=0.0;
			}
		}
	}
	if(par->micromet==1){
		met->Vspdplot=new_doublematrix(Nr,Nc);
		met->Vdirplot=new_doublematrix(Nr,Nc);
		met->RHplot=new_doublematrix(Nr,Nc);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]==UV->V->co[2]){
					met->Vspdplot->co[r][c]=UV->V->co[2];
					met->Vdirplot->co[r][c]=UV->V->co[2];
					met->RHplot->co[r][c]=UV->V->co[2];
				}else{
					met->Vspdplot->co[r][c]=0.0;
					met->Vdirplot->co[r][c]=0.0;
					met->RHplot->co[r][c]=0.0;
				}
			}
		}
	}
}


/****************************************************************************************************/
//LISTON's submodels
if(par->micromet==1){
	if(par->point_sim!=1){
		top->Zm=new_doublematrix(Nr,Nc);
		top->curv_m=new_doublematrix(Nr,Nc);
		top->slope_m=new_doublematrix(Nr,Nc);
		top->slopeaz_m=new_doublematrix(Nr,Nc);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				top->Zm->co[r][c]=top->Z0->co[r][c];
			}
		}
		topo_data(UV->U->co[1], UV->U->co[2], top->Zm, top->curv_m, top->slope_m, top->slopeaz_m, par->curve_len_scale, NoV);

	}else{
		top->Zm=new_doublematrix(top->Z1->nrh,top->Z1->nch);
		top->curv_m=new_doublematrix(top->Z1->nrh,top->Z1->nch);
		top->slope_m=new_doublematrix(top->Z1->nrh,top->Z1->nch);
		top->slopeaz_m=new_doublematrix(top->Z1->nrh,top->Z1->nch);
		for(r=1;r<=top->Z1->nrh;r++){
			for(c=1;c<=top->Z1->nch;c++){
				top->Zm->co[r][c]=top->Z1->co[r][c];
			}
		}
		topo_data(UV->U->co[1], UV->U->co[2], top->Zm, top->curv_m, top->slope_m, top->slopeaz_m, par->curve_len_scale, NoV);
	}
}


if(par->point_sim==1){
	if(par->micromet==1 || par->recover==1) free_doublematrix(IT->LU);
}

//SNOW COVERED AREA STATISTICS
if(par->point_sim==0){
	for(i=1;i<=par->nLC;i++){
		temp=namefile_i(files->co[fHpatch]+1,i);
		f=t_fopen(temp,"w");
		free(temp);
		fprintf(f,"DATE,JD,snowDav,perc.SFA,perc.SCA,H_SCA,H_SFA,MAX_H_SFA,MIN_H_SFA,RSM_H_SFA,TS_MEAN_SFA,TS_MIN_SFA,TS_MAX_SFA,TA_MEAN,Frac(Tsoil>Ta),Hvmean,Hvmin,Hvmax,Hvrsm,Tsvmean\n");
		t_fclose(f);
	}

	temp=join_strings(files->co[fHpatch]+1,textfile);
	f=t_fopen(temp,"w");
	free(temp);
	fprintf(f,"DATE,JD,snowDav,SWEav,Tav,Tsav,perc.SFA,perc.SCA\n");
	t_fclose(f);
}


/****************************************************************************************************/

/*Free the struct allocated in this subroutine:*/
free_doublematrix(SNOW0);
if(par->glaclayer_max>0) free_doublematrix(GLACIER0);
if(par->point_sim==0 && par->wat_balance==1) free_doublematrix(top->Z0dp);

free_doublematrix(IT->land_classes);
free_doublematrix(IT->met);
free_stringbin(IT->met_col_names);
free(IT);

/****************************************************************************************************/

/**********************************************************************************************************/
/**********************************************************************************************************/
//TRANSECTS
/**********************************************************************************************************/
/**********************************************************************************************************/
/*par->transect=(double***)malloc(2*sizeof(double**));
par->vtrans=(double**)malloc(2*sizeof(double*));

par->cont_trans=new_longvector(2);
initialize_longvector(par->cont_trans,0);

par->ibeg=new_longvector(2);
initialize_longvector(par->ibeg,-1);

for(j=0;j<2;j++){

	f=t_fopen(namefile_i(join_strings(WORKING_DIRECTORY,"_Tsrf_time_F14_"),j+1),"r");
	par->transect[j]=read_datameteo(f, 0, 4, NoV);
	t_fclose(f);

	f=t_fopen(namefile_i(join_strings(WORKING_DIRECTORY,"_transectMOD"),j+1),"w");
	fprintf(f,"flight,point,distance,E,N,time,Tmeas,Tmod\n");
	t_fclose(f);

	par->vtrans[j]=(double*)malloc(dim2(par->transect[j])*sizeof(double));

}*/

// file bugs
/*char* namebugs=join_strings(files->element[fbugs]+1,textfile);
f=t_fopen(namebugs,"w");
// need to plot the thermal and hydraulic parameters at various depth (5, 200, 400mm)
fprintf(f,"date,JDfrom0,JD,T#1,T#2,T#3,ThetaW#1,ThetaW#2,ThetaW#3,ThetaI#1,ThetaI#2,ThetaI#3,lambdaT#1,lambdaT#2,lambdaT#3,C_T#1,C_T#2,C_T#3\n");
t_fclose(f);
free(namebugs);*/

}
/*--------------------------------------------------------------------------------------------------*/







/****************************************************************************************************/
/* Subroutine to find the the coefficents to diffuse the channel-flows ("fraction_spread" matrix):  */
/****************************************************************************************************/
DOUBLEMATRIX *De_Saint_Venant(DOUBLEVECTOR *s0,double u0,double D,double Dt)
{
 DOUBLEMATRIX *fraction_spread;    /*the returned output matrix*/
 double cumulative_fraction_spread;/*the cumulative sum of the fractions for a virtual distance*/
 double s0max;                     /*the distance of the farthest channel-pixel [m]*/
 long smax;   /*the maximum number of Dx=u0*Dt from the farthest virtual channel-stretch to the outlet*/
 long ch;     /*counter of the channel-pixels*/
 long s;      /*counter of the virtual channel-stretch; s=1 is the channel-stretch next to outlet*/
 FILE *f;

 /* Finding out the distance of the farthest channel-pixel: */
 s0max=0.0;
 for(ch=1;ch<=s0->nh;ch++){
    if (s0max<s0->co[ch]) s0max=s0->co[ch];
 }

 f=fopen(error_file_name,"a");

 /* Finding out the distance of the farthest virtual channel-stretch; it has to be enough to allow an
    adeguate diffusion (to do this it is fixed a minimum fraction_spread=0.00005):*/
 smax=(long)(s0max/(u0*Dt))+1;
 do{ smax++; }while(((s0max*Dt)/pow((4.0*Pi*D*pow(smax*Dt,3.0)),0.5)*exp(-u0*pow(s0max-smax*u0*Dt,2.0)
                                                                          /(4.0*D*smax*u0*Dt)))>0.00001);

 fprintf(f,"s0max(farthest channel distance)= %f",s0max);
 fprintf(f,"\nsmax(number channel-pixels)=%ld\n",smax);

 /* Initialization of fraction_spread matrix:*/
 fraction_spread=new_doublematrix(s0->nh,smax);

 /* Filling up of fraction_spread matrix:*/
 for(ch=1;ch<=fraction_spread->nrh;ch++){
    cumulative_fraction_spread=0.0;
    for(s=1;s<=fraction_spread->nch;s++){
       /*now is used the following formula:
         DV(x0,x)=(x0*V*Dt)/(4*3.14*D*(x/u)^3)^0.5*exp(-(u*(x0-x)^2)/(4*D*x))
         where:
         DV = water volume fraction = fraction_spread
         V  = unitay volume (dirac) = 1
         x0 = distance from outlet of the letting of water into the channel = s0 [m]
         x  = distance from outlet at which DV is calculated = s*u0*Dt [m]
         D  = mean coefficent of dispersion = D [m^2/s]
         u  = mean speed of water in the channel = u0 [m/s] */
		if (s0->co[ch]<=0.0) t_error ("Negative distance of a channel pixel from outlet!");
		fraction_spread->co[ch][s]=(s0->co[ch]*Dt)/pow((4.0*Pi*D*pow(s*Dt,3.0)),0.5)*exp(-u0
                                        *pow(s0->co[ch]-s*u0*Dt,2.0)/(4.0*D*s*u0*Dt));


		if (fraction_spread->co[ch][s]<0.0) fraction_spread->co[ch][s]=0.0;
		if (fraction_spread->co[ch][s]>0.9) fraction_spread->co[ch][s]=0.9;

		cumulative_fraction_spread+=fraction_spread->co[ch][s];
	}

	if(cumulative_fraction_spread==0){
		fprintf(f,"\nWARNING: channel flow description not correct, reduce Dt or u0\n");
		printf("\nWARNING: channel flow description not correct, reduce Dt or u0\n");
	}

	for(s=1;s<=fraction_spread->nch;s++){
		fraction_spread->co[ch][s]*=(1.0/cumulative_fraction_spread);		//correction
	}
 }

 fclose(f);

 return fraction_spread;
}
/*--------------------------------------------------------------------------------------------------*/











/****************************************************************************************************/
/* year_avg_temp: calcola la temprarura media e l'escursione annuale dell'aria                   */
/* Input:	data_meteo: matrice dati meteorolgici                                                   */
/*  		ndays [gg]: intervallo in cui viene mediato il dato di temperatura per  */
/*	    	calcolare l'escursione media annuale                                                    */
/* Output:	T: temperatura media annuale                                                 */
/*		    DT: escursione annuale di temperatura                                       */
/****************************************************************************************************/
void year_avg_temp(DOUBLEVECTOR *Tdata, double ndays, double *T, double *DT, double Dt)

{
long n,i,j,k;
double temp=0;
double temp_max=0; /* added by Emanuele Cordano on 24/9/9 */
double temp_min=0;

/* calcolo il numero di elementi in ndays */
n=floor(ndays*86400/Dt);

/* calcolo DT */
for(i=n/2+1; i<Tdata->nh-n/2; i++){
	temp=0;
	k=0;
	for(j=i-n/2;j<=i+n/2;j++){
		k+=1;
		temp+=Tdata->co[j];
	}
	temp/=(double)k;
	if(i==n/2+1){
		temp_min=temp;
		temp_max=temp;
	}
	if(temp>temp_max) temp_max=temp;
	if(temp<temp_min) temp_min=temp;
}
*DT=temp_max-temp_min;

/* calcolo T */
*T=0;
for(i=1;i<=Tdata->nh;i++){
	*T+=Tdata->co[i];
}
*T/=(double)Tdata->nh;

}

/*--------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------*/

long row(double N, long nrows, long i, T_INIT *UV)
{

	long cont=0;

	if(N<UV->U->co[3] || N>UV->U->co[3]+nrows*UV->U->co[1]){
		printf("North coordinate %f in point #%ld out of region",N,i);
		t_error("Fatal error");
	}
	do{
		cont+=1;
	}while(UV->U->co[3]+(nrows-cont)*UV->U->co[1]>N);

	return(cont);
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

long col(double E, long ncols, long i, T_INIT *UV)
{

	long cont=0;

	if(E<UV->U->co[4] || E>UV->U->co[4]+ncols*UV->U->co[2]){
		printf("East coordinate %f in point #%ld out of region",E,i);
		t_error("Fatal error");
	}
	do{
		cont+=1;
	}while(UV->U->co[4]+cont*UV->U->co[2]<E);

	return(cont);
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------


void read_inputmaps(TOPO *top, LAND *land, SOIL *sl, PAR *par){

long r, c, i, rtot, ctot;//iLC;
DOUBLEMATRIX *M, *Q, *curv;
SHORTMATRIX *P;
LONGMATRIX *ca;/* map of contributing area */
T_INIT *UV2;
FILE *f;

//reading LAND USE AND TOPOGRAPHY
M=new_doublematrix(1,1);
top->Z0=read_map(0, files->co[fdem]+1, M, UV); //topography
free_doublematrix(M);
write_map(files->co[fdem]+1, 0, par->format_out, top->Z0, UV);

if(existing_file(files->co[flu]+1)>0){ //landuse (reading or write default values)
	land->LC=read_map(1, files->co[flu]+1, top->Z0, UV);
	//Check borders
	for(r=1;r<=land->LC->nrh;r++){
		for(i=1;i<=2;i++){
			land->LC->co[r][i]=UV->V->co[2];
			land->LC->co[r][land->LC->nch+1-i]=UV->V->co[2];
		}
	}
	for(c=1;c<=land->LC->nch;c++){
		for(i=1;i<=1;i++){
			land->LC->co[i][c]=UV->V->co[2];
			land->LC->co[land->LC->nrh+1-i][c]=UV->V->co[2];
		}
	}

	//Land use is the official mask
	for(r=1;r<=land->LC->nrh;r++){
		for(c=1;c<=land->LC->nch;c++){
			if(land->LC->co[r][c]!=UV->V->co[2]){
				if(top->Z0->co[r][c]==UV->V->co[2]){
					printf("ERROR Land use mask include DTM novalue pixels");
					printf("\nr:%ld c:%ld Z:%f landuse:%f\n",r,c,top->Z0->co[r][c],land->LC->co[r][c]);
					land->LC->co[r][c]=UV->V->co[2];
					printf("LANDUSE set at novalue where DTM is not available\n");
					stop_execution();
					//t_error("Land use mask include DTM novalue pixels");
				}
			}
		}
	}
}else{  //writes default value (1)
	//Write land->LC (land cover)
	land->LC=copydoublematrix_const(1.0, top->Z0, UV->V->co[2]);
	for(r=1;r<=land->LC->nrh;r++){
		land->LC->co[r][1]=UV->V->co[2];
		land->LC->co[r][land->LC->nch]=UV->V->co[2];
	}
	for(c=1;c<=land->LC->nch;c++){
		land->LC->co[1][c]=UV->V->co[2];
		land->LC->co[land->LC->nrh][c]=UV->V->co[2];
	}

}
write_map(files->co[flu]+1, 1, par->format_out, land->LC, UV);

/****************************************************************************************************/
//reading SKY VIEW FACTOR
if(existing_file(files->co[fsky]+1)>0){
	top->sky=read_map(2, files->co[fsky]+1, land->LC, UV);
}else{/*The sky view factor file "top->sky" must be calculated*/
	UV2=(T_INIT *)malloc(sizeof(T_INIT));
	if(!UV2) t_error("UV2 was not allocated");
	UV2->U=new_doublevector(4);
	UV2->V=new_doublevector(2);
	rtot=(long)(top->Z0->nrh/par->nsky);
	ctot=(long)(top->Z0->nch/par->nsky);
	if(rtot*par->nsky<top->Z0->nrh) rtot+=1;
	if(ctot*par->nsky<top->Z0->nch) ctot+=1;
	M=new_doublematrix(rtot,ctot);
	reduce_resolution(par->nsky, top->Z0, M, UV, UV2);
	P=new_shortmatrix(rtot,ctot);
	initialize_shortmatrix(P, (short)UV2->V->co[2]);
	nablaquadro_mask(M, P, UV2->U, UV2->V);
	Q=new_doublematrix(rtot,ctot);
	sky_view_factor(Q, 36, UV2, M, P);
	rtot=top->Z0->nrh-par->nsky*M->nrh;
	ctot=top->Z0->nch-par->nsky*M->nch;
	if(rtot<0) rtot=top->Z0->nrh-par->nsky*(M->nrh-1);
	if(ctot<0) ctot=top->Z0->nch-par->nsky*(M->nch-1);
	top->sky=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	amplify_resolution(par->nsky, rtot, ctot, Q, top->sky, UV2, UV);
	free_doublematrix(M);
	free_shortmatrix(P);
	free_doublematrix(Q);
	free_doublevector(UV2->U);
	free_doublevector(UV2->V);
	free(UV2);
}
write_map(files->co[fsky]+1, 0, par->format_out, top->sky, UV);

/****************************************************************************************************/
//reading DRAINAGE DIRECTIONS
if(existing_file(files->co[fdd]+1)>0){
	M=read_map(2, files->co[fdd]+1, land->LC, UV);
	top->DD=copyshort_doublematrix(M);
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			if(land->LC->co[r][c]==UV->V->co[2]) top->DD->co[r][c]=9;
			if(top->DD->co[r][c]>11) top->DD->co[r][c]=0;
		}
	}
	write_map(files->co[fdd]+1, 1, par->format_out, M, UV);
	free_doublematrix(M);
	top->Z0dp=depitted(top->DD, top->Z0);
}else{
	t_error("You have give to give the drainage directions as an input maps");
	/*//novalue (9)
	top->DD=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			top->DD->co[r][c]=9;
		}
	}
	if(par->wat_balance==1){
		par->wat_balance=0.0;
		f=fopen(error_file_name,"a");
		fprintf(f,"Warning: DRAINAGE DIRECTION MAPS NOT PRESENT, cannot resolve water balance\n");
		fclose(f);
		printf("Warning: DRAINAGE DIRECTION MAPS NOT PRESENT, cannot resolve water balance\n");
	}*/
}


/****************************************************************************************************/
//reading SOIL MAP
if(existing_file(files->co[fsoil]+1)>0){
	M=read_map(2, files->co[fsoil]+1, land->LC, UV);
	sl->type=copyshort_doublematrix(M);

}else{//default value (99)
	M=copydoublematrix_const(1.0, land->LC, UV->V->co[2]);
	sl->type=copyshort_doublematrix(M);
}
write_map(files->co[fsoil]+1, 1, par->format_out, M, UV);
free_doublematrix(M);

/****************************************************************************************************/
//reading CURVATURE
if(existing_file(files->co[fcurv]+1)>0){
	M=read_map(2, files->co[fcurv]+1, land->LC, UV);
	top->curv=copyshort_doublematrix(M);

}else{//calculating
	top->curv=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
	initialize_shortmatrix(top->curv,(short)UV->V->co[2]);
	nablaquadro_mask(top->Z0,top->curv,UV->U,UV->V);
	M=copydouble_shortmatrix(top->curv);
}
write_map(files->co[fcurv]+1, 1, par->format_out, M, UV);
free_doublematrix(M);

/****************************************************************************************************/
//SLOPE
top->slopes=new_doublematrix(top->Z0->nrh,top->Z0->nch);
if(existing_file(files->co[fslp]+1)>0){
	M=read_map(2, files->co[fslp]+1, land->LC, UV);		//reads in degrees
	fmultiplydoublematrix(top->slopes, M, Pi/180.0, UV->V->co[2]);	//top->slopes in radiants
}else{
	slopes0875(top->Z0, UV->U, UV->V, top->slopes);		//calculates in radiants
	M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	fmultiplydoublematrix(M, top->slopes, 180.0/Pi, UV->V->co[2]);	//prints in degrees
}
write_map(files->co[fslp]+1, 0, par->format_out, M, UV);
free_doublematrix(M);

/****************************************************************************************************/
//ASPECT
top->aspect=new_doublematrix(top->Z0->nrh,top->Z0->nch);
if(existing_file(files->co[fasp]+1)>0){
	M=read_map(2, files->co[fasp]+1, land->LC, UV);
	fmultiplydoublematrix(top->aspect, M, Pi/180.0, UV->V->co[2]);
}else{
	aspect0875(top->Z0, UV->U, UV->V, top->aspect);
	M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	fmultiplydoublematrix(M, top->aspect, 180.0/Pi, UV->V->co[2]);
}
write_map(files->co[fasp]+1, 0, par->format_out, M, UV);
free_doublematrix(M);

/****************************************************************************************************/
//GRADIENTS along drainage directions
if(par->wat_balance==1){
	if(existing_file(files->co[fgrad]+1)>0){
		top->i_DD=read_map(2, files->co[fgrad]+1, land->LC, UV);
	}else{/*The matrix the slope along the Drainage Direction "top->i_DD" must be calculated*/
		top->i_DD=new_doublematrix(top->Z0->nrh,top->Z0->nch);
		initialize_doublematrix(top->i_DD,UV->V->co[2]);
		gradients(top->Z0dp,top->DD,top->i_DD,UV);
	}
	write_map(files->co[fgrad]+1, 0, par->format_out, top->i_DD, UV);
}else{
	top->i_DD=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	initialize_doublematrix(top->i_DD,0.0);
}

/****************************************************************************************************/
//Effective area of each pixel
/*if(existing_file(files->co[farea]+1)>0){
	top->area=read_map(2, files->co[farea]+1, land->LC, UV);
}else{
	top->area=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	area0875(top->Z0,UV->U,UV->V,top->area);
}
write_map(files->co[farea]+1, 0, par->format_out, top->area, UV);
*/
/****************************************************************************************************/
//Channel network (in top->pixel_type)
if(existing_file(files->co[fnet]+1)>0){
	M=read_map(2, files->co[fnet]+1, land->LC, UV);
	top->pixel_type=copyshort_doublematrix(M);
	free_doublematrix(M);
}else{
	//reading or calculating TCA
	if(existing_file(files->co[ftca]+1)>0){
		M=read_map(2, files->co[ftca]+1, land->LC, UV);
		ca=copylong_doublematrix(M);
	}else{
		ca=new_longmatrix(top->Z0->nrh,top->Z0->nch);
		initialize_longmatrix(ca,0.0);
		tca(top->DD,ca);
		M=copydouble_longmatrix(ca);
	}
	write_map(files->co[ftca]+1, 1, par->format_out, M, UV);
	free_doublematrix(M);
	//calculate laplacian
	curv=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	nablaquadro(top->Z0,curv,UV->U,UV->V);
	top->pixel_type=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
	//copy drainage directions in top->pixel_type
	copy_shortmatrix(top->DD,top->pixel_type);
	//see geomorphologic library (overwrites pixel_type=10 for channels)
	//select_hillslopes_mod(ca,top->i_DD,curv,top->pixel_type,par->channel_thres,UV->U);
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			if(ca->co[r][c]>=par->channel_thres) top->pixel_type->co[r][c]=10;
		}
	}
	free_longmatrix(ca);
	free_doublematrix(curv);
}

// Creation of pixel-type matrix "top->pixel_type" on the basis channel network (already in top->pixel_type matrix):
M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
for(r=1;r<=top->Z0->nrh;r++){
	for(c=1;c<=top->Z0->nch;c++){
		if (top->pixel_type->co[r][c]!=10) top->pixel_type->co[r][c]=0;
		if(land->LC->co[r][c]==UV->V->co[2]) top->pixel_type->co[r][c]=9;
		M->co[r][c]=(double)top->pixel_type->co[r][c];
	}
}
write_map(files->co[fnet]+1, 1, par->format_out, M, UV);
free_doublematrix(M);

/****************************************************************************************************/

//Outlet distances
if(existing_file(files->co[fdist]+1)>0){
	top->pixel_distance=read_map(1, files->co[fdist]+1, land->LC, UV);
	//Check values
	for(r=1;r<=top->pixel_distance->nrh;r++){
		for(c=1;c<=top->pixel_distance->nch;c++){
			if(top->pixel_distance->co[r][c]==UV->V->co[2]){
				if(land->LC->co[r][c]!=UV->V->co[2]){
					f=fopen(error_file_name,"a");
					fprintf(f,"\nPixeldistance has NOVALUE where DTM does not in pixel r=%4ld c=%4ld, corrected, but inconsistently\n",r,c);
					fclose(f);
					top->pixel_distance->co[r][c]=0.0;
				}
			}
		}
	}
	//A null top->pixel_distance is changed in half pixel size
	for(r=1;r<=land->LC->nrh;r++){
		for(c=1;c<=land->LC->nch;c++){
			if(land->LC->co[r][c]!=UV->V->co[2]){
				if (top->pixel_distance->co[r][c]<=0.0) top->pixel_distance->co[r][c]=UV->U->co[2]/2.0;
			}
		}
	}
}else{
	top->pixel_distance=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	initialize_doublematrix(top->pixel_distance,UV->V->co[2]);
	//se geomorphologic libraries
	outletdistance(top->DD,top->pixel_distance,UV->U);
	for(r=1;r<=land->LC->nrh;r++){
		for(c=1;c<=land->LC->nch;c++){
			if(top->pixel_distance->co[r][c]==UV->V->co[2]){
				if(land->LC->co[r][c]!=UV->V->co[2]){
					f=fopen(error_file_name,"a");
					fprintf(f,"\nPixeldistance has NOVALUE where DTM does not in pixel r=%4ld c=%4ld, corrected, but inconsistently\n",r,c);
					fclose(f);
					top->pixel_distance->co[r][c]=0.0;
				}
			}
			if(land->LC->co[r][c]!=UV->V->co[2]){
				if (top->pixel_distance->co[r][c]<=0.0) top->pixel_distance->co[r][c]=UV->U->co[2]/2.0;
			}
		}
	}
}
write_map(files->co[fdist]+1, 0, par->format_out, top->pixel_distance, UV);



//Additional maps
M=copydoublematrix_const(1.0, land->LC, UV->V->co[2]);
land->LC2=copyshort_doublematrix(M);
par->nLC=1;
free_doublematrix(M);

if(files->index->nh!=nfiles) t_error("Error in file number");
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void read_parameterfile(char *name, PAR *par, INIT_TOOLS *itools){
/* reads the file __parameters.txt */
	FILE *f;
	long i,index;
	DOUBLEVECTOR *v;
	char *temp;
	temp=join_strings(name,textfile); /* added by Emanuele Cordano on 22 September 2008 */
	f=t_fopen(temp,"r");
	free(temp);
	index=read_index(f,PRINT);
	//1st block
	itools->land_classes=read_doublematrix(f,"a",PRINT);
	//2nd block
	v=read_doublearray(f,PRINT);
	if(v->co[1]<=0) t_error("ERROR: 0 is not admitted for nDtwater");
	par->nDt_water=(long)v->co[1];
	par->f_bound_Richards=v->co[2];	/*Parameter for the bottom boundary condition for the Richards' equation: =0 no flux, =1 free drainage*/
	par->imp=v->co[3];	/*Impedence factor for (partially) frozen soil*/
	par->psimin=v->co[4];
	par->psimin2=par->psimin;
	par->Esoil=v->co[5];
	par->MaxiterVWB=(long)v->co[6];
	par->TolVWb=v->co[7];
	par->Dpsi=v->co[8];
	par->dtmin=v->co[9];
	itools->u0=v->co[10]; /*THE MEAN VELOCITY IN CHANNELS*/
	itools->D=v->co[11];  /*THE HYDRODYNAMICAL DISPERSION IN CHANNELS*/
	par->gamma_m=v->co[12]; /*Exponent of the law of uniform motion on the surface*/
	free_doublevector(v);
	//3rd block
	v=read_doublearray(f,PRINT);
	par->latitude=v->co[1]*Pi/180.0;
	par->longitude=v->co[2]*Pi/180.0;
	par->Vis=v->co[3];
	par->Lozone=v->co[4];
	par->Vmin=v->co[5];		/*MINIMUM WIND VELOCITY VALUE (m/s)*/
	par->RHmin=v->co[6];		/*MINIMUM RELATIVE HUMIDITY VALUE (%)*/

	par->ifill=(int)v->co[7];
	par->iobsint=(int)v->co[8];
	par->dn=(float)v->co[9];
	par->curve_len_scale=(float)v->co[10];
	par->slopewt=(float)v->co[11];
	par->curvewt=(float)v->co[12];
	par->topoflag=(float)v->co[13];

	free_doublevector(v);
	/* 4th block Snow Parameters*/
	v=read_doublearray(f,PRINT);
	itools->Dsnow0=v->co[1];/* INITIAL SNOW DEPTH (initial snow depth) [mm] valid only if there is no snow map*/
	itools->rhosnow0=v->co[2]; /*SNOW INITIAL DENSITY [kg/mc]*/
	itools->agesnow0=v->co[3];/* INITIAL SNOW AGE (in days), valid only if there is no snow age map */
	par->T_rain=v->co[4];   /*TEMPERATURE ABOVE WICH ALL PRECIPITAION IS RAIN [C]*/
	par->T_snow=v->co[5];   /*TEMPERATURE BELOW WICH ALL PRECIPITAION IS SNOW [C]*/
	par->aep=v->co[6];      /*AßLBEDO EXTINCTION PARAMETER [m]*/
	par->avo=v->co[7];      /*NEW SNOW VISIBLE BAND REFLECTANCE*/
	par->airo=v->co[8];     /*NEW NEAR INFRARED BAND REFLECTANCE*/
	par->Sr=v->co[9];       /*IRREDUCIBLE WATER SATURATION [-]*/
	par->epsilon_snow=v->co[10];/* SNOW LONGWAVE EMISSIVITY [-]  */
	par->z0_snow=v->co[11]*0.001;/* roughness length over snow [mm]*/
	par->snowcorrfact=v->co[12];/* INCREASING FACTOR WHEN THE RAIN GAUGE IS SUPPOSED TO RECORD SNOW PRECIPITATION */
	par->raincorrfact=v->co[13];/* INCREASING FACTOR WHEN THE RAIN GAUGE IS SUPPOSED TO RECORD RAIN PRECIPITATION */
	par->snowlayer_max=(long)v->co[14]; /* MAXIMUM NUMBER OF SNOW LAYERS */
	par->snowlayer_inf=(long)v->co[15];
	par->snow_maxpor=v->co[16];/* MAXIMUM SNOW POROSITY ALLOWED [-]*/
	par->drysnowdef_rate=v->co[17];
	par->snow_density_cutoff=v->co[18];
	par->wetsnowdef_rate=v->co[19];
	par->snow_viscosity=v->co[20];
	par->snow_fetch=v->co[21];
	//if(par->snowlayer_max<3) t_error("Maximum number of snow layers too small, it must be not smaller than 3");
	//if( par->snowlayer_max!=1+2*floor(0.5*(par->snowlayer_max-1)) ) t_error("Maximum number of snow layers must be odd");
	free_doublevector(v);
	/* 5th-6th block MINIMUM and MAXIMUM SNOW LAYER THICKNESS*/
	par->Dmin=read_doublearray(f,PRINT);
	if(par->Dmin->nh!=par->snowlayer_max) t_error("Error in assigning max and min thickness to the snow layers");
	par->Dmax=read_doublearray(f,PRINT);
	if(par->Dmin->nh!=par->snowlayer_max) t_error("Error in assigning max and min thickness to the snow layers");
	/*if(par->Dmax->co[(long)(1+floor(0.5*par->snowlayer_max))]<1.0E10){
		printf("\nWARNING: Snow layer %ld must have infinite thickness!!! Assigned %f mm instead of %f mm\n",(long)(1+floor(0.5*par->snowlayer_max)),1.0E10,par->Dmax->co[(long)(1+floor(0.5*par->snowlayer_max))]);
		stop_execution();
		par->Dmax->co[(long)(1+floor(0.5*par->snowlayer_max))]=1.0E10;
	}*/
	//7th block
	v=read_doublearray(f,PRINT);
	itools->Dglac0=v->co[1];
	itools->rhoglac0=v->co[2];
	itools->Tglac0=v->co[3];
	par->Sr_glac=v->co[4];
	par->glaclayer_max=(long)v->co[5];
	if(par->glaclayer_max<0) t_error("Maximum number of glacier cannot be negative");
	if(par->glaclayer_max==0) printf("\nNo glacier is considered\n");
	free_doublevector(v);
	//8th block
	v=read_doublearray(f,PRINT);
	if(par->glaclayer_max>0){
		par->Dmin_glac=new_doublevector(par->glaclayer_max);
		par->Dmax_glac=new_doublevector(par->glaclayer_max);
		if(v->nh!=2*par->glaclayer_max-1) t_error("Error in assigning max and min thickness to the glacier layers");
		for(i=1;i<=par->glaclayer_max;i++){
			par->Dmin_glac->co[i]=v->co[i];
		}
		for(i=1;i<=par->glaclayer_max-1;i++){
			par->Dmax_glac->co[i]=v->co[i+par->glaclayer_max];
		}
		par->Dmax_glac->co[par->glaclayer_max]=1.0E10;
	}
	free_doublevector(v);
	//9th block
	v=read_doublearray(f,PRINT);
	par->state_turb=(short)v->co[1];
	par->state_lwrad=(short)v->co[2];
	par->monin_obukhov=(short)v->co[3];
	par->micromet=(short)v->co[4];
	par->blowing_snow=(short)v->co[5];
	if(v->nh>5) par->superfast=(short)v->co[6];// superfast version
	//printf("par->point_sim=%d, part->superfast=%ld",par->point_sim, par->superfast); stop_execution();
	if(par->point_sim!=1) par->superfast=0;
	//printf("par->point_sim=%d, part->superfast=%ld",par->point_sim, par->superfast); stop_execution();
	if(par->blowing_snow==1 && par->micromet==0){
		par->blowing_snow=0;
		printf("\nWarning: if you do not run Micromet, you can't run SnowTrans3D\n");
	}
	free_doublevector(v);
	//other par
	par->print=0;

	t_fclose(f);

	/*Parameters derived, already calculated or that have to be calculated because useful in the next subroutine:*/
	par->n_error=0;        /*Current number of error of the simulation*/
	par->max_error=10000; /*Maximum number of error for the simulation*/
	par->state_pixel=1;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
void read_optionsfile_distr(char *name, PAR *par, TIMES *times, DOUBLEMATRIX *Z0){

	FILE *f;
	long index, i;
	DOUBLEVECTOR *v;

	f=t_fopen(join_strings(name,textfile),"r");
	index=read_index(f,PRINT);

	//1. base par
	v=read_doublearray(f,PRINT);
	par->wat_balance=(short)v->co[1];
	par->en_balance=(short)v->co[2];
	par->state_px_coord=(short)v->co[3];
	par->distr_stat=(short)v->co[4];
	free_doublevector(v);

	//2. chkpt
	/* 2 block - COORDINATES of the points for which point values are plotted
	if state_px_coord==1 (E, N, layer) ---  if state_px_coord==0 (row,col,layer) (max 9999 points)*/
	par->chkpt=read_doublematrix(f,"a",PRINT);

	//3. saving points
	par->saving_points=read_doublearray(f,PRINT);

	//4. output par
	v=read_doublearray(f,PRINT);
	par->output_Txy=v->co[1];
	par->output_TETAxy=v->co[2];
	par->output_TETAICExy=v->co[3];
	par->output_PSIxy=v->co[4];
	par->output_snow=v->co[5];
	par->output_glac=v->co[6];
	par->output_h_sup=v->co[7];
	par->output_albedo=v->co[8];
	par->output_Rn=v->co[9];
	par->output_G=v->co[10];
	par->output_H=v->co[11];
	par->output_ET=v->co[12];
	par->output_Ts=v->co[13];
	par->output_P=v->co[14];
	par->output_Wr=v->co[15];
	par->output_balancesn=v->co[16];
	par->output_balancegl=v->co[17];
	par->output_Rswdown=v->co[18];
	par->output_meteo=v->co[19];
	free_doublevector(v);

	//5. special output
	/* 5 block - SPECIAL OUTPUT PARAMETERS
	It is possible to select some days for which energy balance and meteo data are plotted with a very short time step.
	This is useful if some careful analysis has to be done. If you do not want to use this possibility, just write 0 as first
	component of this vector */
	v=read_doublearray(f,PRINT);

	if(v->co[1]>0 && v->nh>1){
		if(fmod(v->nh-1,2)!=0) t_error("special output vector must have an odd number of components");
		times->n_plot=floor(v->co[1]*3600/par->Dt);

		par->JD_plots=new_longvector(v->nh-1);	//even number of components// Matteo on 6.10.09 instead of new_doublevector puts new_longvector
		for(i=1;i<=v->nh-1;i++){
			par->JD_plots->co[i]=v->co[i+1];
		}
	}else{
		par->JD_plots=new_longvector(1);// Matteo on 6.10.09 instead of new_doublevector puts new_longvector
		par->JD_plots->co[1]=-1.;
	}
	free_doublevector(v);

	t_fclose(f);

	if(par->chkpt->nch<2) t_error("Error in I_OPTIONS file: E and N coordinates missing");
	if(par->chkpt->nrh>9999) t_error("Error in I_OPTIONS file: no more than 9999 points allowed");

	if(par->state_pixel==1){
		par->rc=new_longmatrix(par->chkpt->nrh,2);
		for(i=1;i<=par->chkpt->nrh;i++){
			if(par->state_px_coord==1){
				par->rc->co[i][1]=row(par->chkpt->co[i][2], Z0->nrh, i, UV);
				par->rc->co[i][2]=col(par->chkpt->co[i][1], Z0->nch, i, UV);
			}else{
				par->rc->co[i][1]=(long)par->chkpt->co[i][1];
				par->rc->co[i][2]=(long)par->chkpt->co[i][2];
			}
			if(Z0->co[par->rc->co[i][1]][par->rc->co[i][2]]==UV->V->co[2]){
				printf("\nWARNING: Point #%4ld corresponds to NOVALUE pixel",i);
				stop_execution();
			}
		}
	}
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
void read_optionsfile_point(char *name, PAR *par, TOPO *top, LAND *land, SOIL *sl, INIT_TOOLS *IT, TIMES *times){

	FILE *f;
	long index, i, r, c;
	DOUBLEVECTOR *v;
	STRINGBIN *s;
	DOUBLEMATRIX *M, *Q;
	DOUBLEMATRIX *P=NULL; /* corrected by Emanuele Cordano on 24/9/9 */
	SHORTMATRIX *curv;
	short read_dem, read_lu, read_soil, read_sl, read_as, read_sk, read_dx, read_dy;
	char *temp;  /*added by Emanuele Cordano on 22 September 2009 */

	temp=join_strings(name,textfile);
	f=t_fopen(temp,"r");
	free(temp);
	index=read_index(f,PRINT);

	//1. base par
	v=read_doublearray(f,PRINT);
	par->wat_balance=(short)v->co[1];
	par->en_balance=(short)v->co[2];
	par->state_px_coord=(short)v->co[3];
	//ly=(short)v->co[4];


	//2. chkpt
	/** 2 block - COORDINATES of the points for which the simulations is run
	if state_px_coord==1 coordinate are (E, N) ---  if state_px_coord==0 coordinate are (row,col) (max 9999 points)
	#1	E(row)
	#2	N(col)
	#3	Z[m]
	#4	landuse
	#5	soiltype
	#6	slopes[rad]
	#7	aspect[rad]
	#8	sky[-]
	#9	dz/dE(net)
	#10	dz/dS(net)
	write -99 if you want to read the values from the distributed files */
	M=read_doublematrix(f,"a",PRINT);
	if(M->nch!=10) t_error("Error in block 2 option file");

	//3. saving points
	par->saving_points=read_doublearray(f,PRINT);
	if(par->saving_points->nh>1)  {
		printf("You have inserted %ld saving points when in 1D simulation they are not allowed",par->saving_points->nh);
		t_error("Saving points in block 3 of file option are not allowed! Set to {0} ");
	}else if (par->saving_points->co[1]>0.0) {
		printf("You have inserted one saving point equal to %f when in 1D simulation only 0 is allowed",par->saving_points->co[1]);
		t_error("Saving points in block 3 of file option are not allowed! Set to {0} ");
	}



	//4. horizon name
	s=read_stringarray(f,PRINT);

	t_fclose(f);

	//4. CALCULATE TOPOGRAPHIC PROPERTIES
	//a. read dem
	read_dem=0;
	for(i=1;i<=M->nrh;i++){ // for every given point
		if(M->co[i][3]==-99.0 || M->co[i][6]==-99.0 || M->co[i][7]==-99.0 || M->co[i][8]==-99.0 || M->co[i][9]==-99.0 || M->co[i][10]==-99.0) read_dem=1;
	}
	if(par->micromet==1 || par->recover==1) read_dem=1;
	if(read_dem==1){
		if(existing_file(files->co[fdem]+1)>0){
			Q=new_doublematrix(1,1);
			top->Z1=read_map(0, files->co[fdem]+1, Q, UV); //topography
			free_doublematrix(Q);
		}else{
			printf("Warning: Dem file not present - Cannot apply micromet\n");
			read_dem=0;
		}
	}
	if(read_dem==0 && par->state_px_coord==1){
		printf("\nWarning: if a dem file is not used or not present, you can't set state_px_coord at 1\n");
		printf("CORRECTED\n");
		par->state_px_coord=0;
	}
	if(read_dem==0){
		par->micromet=0;
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][3]==-99) M->co[i][3]=0.0;
		}
	}else{
		par->r_points=new_longvector(M->nrh);
		par->c_points=new_longvector(M->nrh);
		for(i=1;i<=M->nrh;i++){
			if(par->state_px_coord==0){
				par->r_points->co[i]=(long)M->co[i][1];
				par->c_points->co[i]=(long)M->co[i][2];
			}else{
				par->r_points->co[i]=row(M->co[i][2], top->Z1->nrh, i, UV);
				par->c_points->co[i]=col(M->co[i][1], top->Z1->nch, i, UV);
			}
			if(M->co[i][3]==-99) M->co[i][3]=top->Z1->co[par->r_points->co[i]][par->c_points->co[i]];
		}
	}


	//b. read land use
	read_lu=0;
	for(i=1;i<=M->nrh;i++){ if(M->co[i][4]==-99) read_lu=1; }
	if(par->micromet==1 || par->recover==1) read_lu=1;
	if(read_lu==1){
		if(existing_file(files->co[flu]+1)>0){
			if(read_dem==0){
				Q=new_doublematrix(1,1);
				IT->LU=read_map(0, files->co[flu]+1, Q, UV);
				free_doublematrix(Q);
			}else{
				IT->LU=read_map(1, files->co[flu]+1, top->Z1, UV);
			}
		}else{
			printf("Warning: Landuse file not present, uniform cover considered\n");
			if(read_dem==1){
				IT->LU=copydoublematrix_const(1.0, top->Z1, UV->V->co[2]);
			}else{
				read_lu=0;
			}
		}
	}
	printf("read land use 0=No, 1=Yes: %d\n",read_lu);
	if(read_lu==0){
		par->micromet=0;
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][4]==-99) M->co[i][4]=1.0;
		}
		printf("%f\n",M->co[1][4]);
	}else{
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][4]==-99){
				if(par->state_px_coord==0){
					r=(long)M->co[i][1];
					c=(long)M->co[i][2];
				}else{
					r=row(M->co[i][2], IT->LU->nrh, i, UV);
					c=col(M->co[i][1], IT->LU->nch, i, UV);
				}
				M->co[i][4]=IT->LU->co[r][c];
			}
		}
	}

	//c. read soil type
	read_soil=0;
	for(i=1;i<=M->nrh;i++){ if(M->co[i][5]==-99) read_soil=1; }
	if(read_soil==1){
		if(existing_file(files->co[fsoil]+1)>0){
			if(read_lu==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files->co[fsoil]+1, Q, UV);
				free_doublematrix(Q);
			}else{
				P=read_map(2, files->co[fsoil]+1, IT->LU, UV);
			}
		}else{
			printf("Warning: Soiltype file not present\n");
			read_soil=0;
		}
	}
	if(read_soil==0){
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][5]==-99) M->co[i][5]=1.0;
		}
	}else{
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][5]==-99){
				if(par->state_px_coord==0){
					r=(long)M->co[i][1];
					c=(long)M->co[i][2];
				}else{
					r=row(M->co[i][2], P->nrh, i, UV);
					c=col(M->co[i][1], P->nch, i, UV);
				}
				M->co[i][5]=P->co[r][c];
			}
		}
		free_doublematrix(P);
	}


	//d. read slopes
	read_sl=0;
	for(i=1;i<=M->nrh;i++){ if(M->co[i][6]==-99) read_sl=1; }
	if(read_sl==1){
		if(existing_file(files->co[fslp]+1)>0){
			if(read_lu==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files->co[fslp]+1, Q, UV);
				free_doublematrix(Q);
			}else{
				P=read_map(2, files->co[fslp]+1, IT->LU, UV);
			}
		}else{
			if(read_dem==0){
				printf("Warning: Slopes file not present\n");
				read_sl=0;
			}else{
				P=new_doublematrix(top->Z1->nrh,top->Z1->nch);
				slopes0875(top->Z1,UV->U,UV->V,P);
				for(r=1;r<=top->Z1->nrh;r++){
					for(c=1;c<=top->Z1->nch;c++){
						if(top->Z1->co[r][c]!=UV->V->co[2]){
							P->co[r][c]*=180.0/Pi;
						}
					}
				}
				write_map(files->co[fslp]+1, 0, par->format_out, P, UV);
			}
		}
	}

	if(read_sl==0){
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][6]==-99) M->co[i][6]=0.0;
		}
	}else{
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][6]==-99){
				if(par->state_px_coord==0){
					r=(long)M->co[i][1];
					c=(long)M->co[i][2];
				}else{
					r=row(M->co[i][2], P->nrh, i, UV);
					c=col(M->co[i][1], P->nch, i, UV);
				}
				M->co[i][6]=P->co[r][c];
			}
		}
		free_doublematrix(P);
	}

	//e. read aspect
	read_as=0;
	for(i=1;i<=M->nrh;i++){ if(M->co[i][7]==-99) read_as=1; }
	if(read_as==1){
		if(existing_file(files->co[fasp]+1)>0){
			if(read_lu==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files->co[fasp]+1, Q, UV);
				free_doublematrix(Q);
			}else{
				P=read_map(2, files->co[fasp]+1, IT->LU, UV);
			}
		}else{
			if(read_dem==0){
				printf("Warning: Aspect file not present\n");
				read_as=0;
			}else{
				P=new_doublematrix(top->Z1->nrh,top->Z1->nch);
				aspect0875(top->Z1,UV->U,UV->V,P);
				for(r=1;r<=top->Z1->nrh;r++){
					for(c=1;c<=top->Z1->nch;c++){
						if(top->Z1->co[r][c]!=UV->V->co[2]){
							P->co[r][c]*=180.0/Pi;
						}
					}
				}
				write_map(files->co[fasp]+1, 0, par->format_out, P, UV);
			}
		}
	}
	if(read_as==0){
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][7]==-99) M->co[i][7]=0.0;
		}
	}else{
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][7]==-99){
				if(par->state_px_coord==0){
					r=(long)M->co[i][1];
					c=(long)M->co[i][2];
				}else{
					r=row(M->co[i][2], P->nrh, i, UV);
					c=col(M->co[i][1], P->nch, i, UV);
				}
				M->co[i][7]=P->co[r][c];
			}
		}
		free_doublematrix(P);
	}

	//f. sky view factor file
	read_sk=0;
	for(i=1;i<=M->nrh;i++){ if(M->co[i][8]==-99) read_sk=1; }
	if(read_sk==1){
		if(existing_file(files->co[fsky]+1)>0){
			if(read_lu==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files->co[fsky]+1, Q, UV);
				free_doublematrix(Q);
			}else{
				P=read_map(2, files->co[fsky]+1, IT->LU, UV);
			}
		}else{
			if(read_dem==0){
				printf("Warning: Sky view factor file not present\n");
				read_sk=0;
			}else{
				P=new_doublematrix(top->Z1->nrh,top->Z1->nch);
				curv=new_shortmatrix(top->Z1->nrh,top->Z1->nch);
				nablaquadro_mask(top->Z1, curv, UV->U, UV->V);
				sky_view_factor(P, 36, UV, top->Z1, curv);
				free_shortmatrix(curv);
				write_map(files->co[fsky]+1, 0, par->format_out, P, UV);
			}
		}
	}
	if(read_sk==0){
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][8]==-99) M->co[i][8]=1.0;
		}
	}else{
		for(i=1;i<=M->nrh;i++){
			if(M->co[i][8]==-99){
				if(par->state_px_coord==0){
					r=(long)M->co[i][1];
					c=(long)M->co[i][2];
				}else{
					r=row(M->co[i][2], P->nrh, i, UV);
					c=col(M->co[i][1], P->nch, i, UV);
				}
				M->co[i][8]=P->co[r][c];
			}
		}
		free_doublematrix(P);
	}

	//g. dz/dx(east)
	read_dx=0;
	for(i=1;i<=M->nrh;i++){ if(M->co[i][9]==-99) read_dx=1; }
	if(read_dx==1){
		if(read_dem==0){
			printf("Impossible to calculate net E slope: set at 0");
			read_dx=0;
		}
	}
	for(i=1;i<=M->nrh;i++){
		if(par->state_px_coord==0){
			r=(long)M->co[i][1];
			c=(long)M->co[i][2];
		}else{
			r=row(M->co[i][2], top->Z1->nrh, i, UV);
			c=col(M->co[i][1], top->Z1->nch, i, UV);
		}
		if(read_dx==1){
			if(c>1 && c<top->Z1->nch){
				if(top->Z1->co[r][c+1]!=UV->V->co[2] && top->Z1->co[r][c-1]!=UV->V->co[2]){
					M->co[i][9]=( (top->Z1->co[r][c]-top->Z1->co[r][c+1]) - (top->Z1->co[r][c-1]-top->Z1->co[r][c])  )/UV->U->co[1];
				}
			}else{
				printf("Impossible to calculate net E slope for point %ld\n: set at 0",i);
				M->co[i][9]=0.0;
			}
		}else{
			if(M->co[i][9]==-99) M->co[i][9]=0.0;
		}
	}


	//h. dz/dy(north)
	read_dy=0;
	for(i=1;i<=M->nrh;i++){ if(M->co[i][10]==-99) read_dy=1; }
	if(read_dy==1){
		if(read_dem==0){
			printf("Impossible to calculate net S slope: set at 0");
			read_dy=0;
		}
	}
	for(i=1;i<=M->nrh;i++){
		if(par->state_px_coord==0){
			r=(long)M->co[i][1];
			c=(long)M->co[i][2];
		}else{
			r=row(M->co[i][2], top->Z1->nrh, i, UV);
			c=col(M->co[i][1], top->Z1->nch, i, UV);
		}
		if(read_dy==1){
			if(r>1 && r<top->Z1->nrh){
				if(top->Z1->co[r+1][c]!=UV->V->co[2] && top->Z1->co[r-1][c]!=UV->V->co[2]){
					M->co[i][10]=( (top->Z1->co[r][c]-top->Z1->co[r+1][c]) - (top->Z1->co[r-1][c]-top->Z1->co[r][c])  )/UV->U->co[2];
				}
			}else{
				printf("Impossible to calculate net S slope for point %ld\n: set at 0",i);
				M->co[i][10]=0.0;
			}
		}else{
			if(M->co[i][10]==-99) M->co[i][10]=0.0;
		}
	}


	//i.show results
	printf("\nPOINTS:\n");
	f=fopen(error_file_name,"a");
	fprintf(f,"\nPOINTS:\n");
	for(r=1;r<=M->nrh;r++){// # of points
		for(c=1;c<=10;c++){// # of properties
			printf("%f  ",M->co[r][c]);
			fprintf(f,"%f  ",M->co[r][c]);
		}
		printf("\n");
		fprintf(f,"\n");
	}
	fclose(f);

	//l. set UV
	if(read_dem==0 && read_lu==0 && read_soil==0 && read_sl==0 && read_as==0 && read_sk==0){
		UV=(T_INIT *)malloc(sizeof(T_INIT));
		if(!UV) t_error("UV was not allocated in input");
		UV->U=new_doublevector(4);
		UV->V=new_doublevector(2);
		UV->U->co[2]=1.0;
		UV->U->co[1]=1.0;
		UV->U->co[4]=0.0;
		UV->U->co[3]=0.0;
		UV->V->co[1]=-1;
		UV->V->co[2]=-9999.0;
	}

	//m. deallocation
	if(read_dem==1 && (par->micromet!=1 && par->recover!=1)){
		free_doublematrix(top->Z1);
		free_longvector(par->r_points);
		free_longvector(par->c_points);
	}

	//5. SET CHECKPOINT
	par->chkpt=new_doublematrix(M->nrh,2);
	par->rc=new_longmatrix(M->nrh,2);
	for(i=1;i<=M->nrh;i++){
		par->chkpt->co[i][1]=M->co[i][1];
		par->chkpt->co[i][2]=M->co[i][2];
		//par->chkpt->co[i][3]=(double)ly;
		par->rc->co[i][1]=1;
		par->rc->co[i][2]=i;
	}
	printf("%f\n",M->co[1][4]);

	//6. SET PROPERTIES
	top->Z0=new_doublematrix(1,M->nrh);// matrice 1x(#punti di misura)
	land->LC=new_doublematrix(1,M->nrh);
	sl->type=new_shortmatrix(1,M->nrh);
	top->slopes=new_doublematrix(1,M->nrh);
	top->aspect=new_doublematrix(1,M->nrh);
	top->sky=new_doublematrix(1,M->nrh);
	top->dz_dx=new_doublematrix(1,M->nrh);
	top->dz_dy=new_doublematrix(1,M->nrh);
	top->DD=new_shortmatrix(1,M->nrh);
	top->curv=new_shortmatrix(1,M->nrh);
	top->i_DD=new_doublematrix(1,M->nrh);
	/*top->area=new_doublematrix(1,M->nrh);// commented by Matteo on 3/11/09*/
	top->pixel_type=new_shortmatrix(1,M->nrh);
	top->pixel_distance=new_doublematrix(1,M->nrh);
	for(i=1;i<=M->nrh;i++){
		top->Z0->co[1][i]=M->co[i][3];
		land->LC->co[1][i]=M->co[i][4];
		sl->type->co[1][i]=(short)M->co[i][5];
		top->slopes->co[1][i]=M->co[i][6]*Pi/180.0;
		top->aspect->co[1][i]=M->co[i][7]*Pi/180.0;
		top->sky->co[1][i]=M->co[i][8];
		top->dz_dx->co[1][i]=M->co[i][9];
		top->dz_dy->co[1][i]=M->co[i][10];
		//these values are not used in 1d simulations
		top->DD->co[1][i]=10;	//outlet
		top->curv->co[1][i]=0;
		top->i_DD->co[1][i]=0.0;
		/*top->area->co[1][i]=UV->U->co[1]*UV->U->co[2];// commented by Matteo on 3/11/09 */
		top->pixel_type->co[1][i]=0;
		top->pixel_distance->co[1][i]=1.0;
	}

	//7. SET PAR
	par->output_Txy=0;
	par->output_TETAxy=0;
	par->output_TETAICExy=0;
	par->output_PSIxy=0;
	par->output_snow=0;
	par->output_glac=0;
	par->output_h_sup=0;
	par->output_albedo=0;
	par->output_Rn=0;
	par->output_G=0;
	par->output_H=0;
	par->output_ET=0;
	par->output_Ts=0;
	par->output_P=0;
	par->output_Wr=0;
	par->output_balancesn=0;
	par->output_balancegl=0;
	par->output_Rswdown=0;
	par->output_meteo=0;
	par->JD_plots=new_longvector(1); /* corrected in longvector by Emanuele Cordano on 24/9/9 */
	par->JD_plots->co[1]=-1.;
	par->blowing_snow=0;
	par->topoflag=0;

	//8. READ HORIZONS
	top->horizon_height=(double ****)malloc(top->Z0->nrh*sizeof(double***));
	for(r=1;r<=top->Z0->nrh;r++){// top->Z0->nrh=1 therefore there is only one line...I don't know why
		top->horizon_height[r-1]=(double ***)malloc(top->Z0->nch*sizeof(double**));// as many as the number of points
		for(c=1;c<=top->Z0->nch;c++){// for every point of the simulation 1D
			i=c;
			temp=join_strings(WORKING_DIRECTORY,s->co[1]+1);
			top->horizon_height[r-1][c-1]=read_horizon(temp,i);
			free(temp);
		}
	}
	free_stringbin(s); /* modified by Emanuele Cordano on 23 September 2009 again */

	free_doublematrix(M);

}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
double **read_horizon(char *name, long i){
	/* Author: Stefano Endrizzi, 2008
	* reads the horizon file
	* Input:
	* name: name of the horizon file
	* i: number of the point
	* Output:
	* hor: matrix of double containing the elevation angle of the obstacle at each azimut angle
	* Modified and commented by: Matteo Dall'Amico, September 2009 */

	FILE *f;
	DOUBLEMATRIX *M;
	long j, index;
	double **hor;
	char *temp;
	temp=namefile_i(name,i);// adds the suffix "0001" etc. to the file name
	f=fopen(temp,"r");
	free(temp);
	if(f!=NULL){
		index=read_index(f,PRINT);
		M=read_doublematrix(f,"a",PRINT);
		fclose(f);

		if(M->nch>2){
			printf("Warning: in file %s columns after the 2nd will be neglected\n",temp);
		}else if(M->nch<2){
			printf("Error: in file %s insufficient number of columns\n",temp);
			t_error("ERROR!!");
		}
		// allocates the matrix of horizon
		hor=alloc2(M->nrh,2);/* M->rh is the number of azimuth in which the horizon is divided */
		for(j=1;j<=M->nrh;j++){
			hor[j-1][0]=M->co[j][1];/* azimuth angle */
			hor[j-1][1]=M->co[j][2];/* horizon elevation angle */
		}
		free_doublematrix(M);
	}else{
		f=fopen(error_file_name,"a");
		fprintf(f,"\nNo horizon file of the Meteo Station #%ld is present. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
		printf("\nNo horizon file of the Meteo Station #%ld is present. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
		fprintf(f,"/** Horizon file for met station or point #%ld */\n",i);
		fprintf(f,"\n");
		fprintf(f,"1: double matrix horizon {4,2}\n");

		hor=alloc2(4,2);
		for(j=1;j<=4;j++){
			hor[j-1][0]=45.0+(j-1)*90.0;
			hor[j-1][1]=0.0;
			fprintf(f,"%f %f\n",hor[j-1][0],hor[j-1][1]);
		}
		fclose(f);
	}
	return(hor);
}




/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
void init_meteo_stations(DOUBLEMATRIX *INPUTmeteo, METEO_STATIONS *st){

	long i;

	if(INPUTmeteo->nch!=13) t_error("ERROR IN INPUT OF METEO STATIONS");

	st->E=new_doublevector(INPUTmeteo->nrh); /* East coordinate [m] of the meteo station */
	st->N=new_doublevector(INPUTmeteo->nrh);/* North coordinate [m] of the meteo station */
	st->lat=new_doublevector(INPUTmeteo->nrh);/* Latitude [rad] of the meteo station */
	st->lon=new_doublevector(INPUTmeteo->nrh);/* Longitude [rad] of the meteo station */
	st->Z=new_doublevector(INPUTmeteo->nrh);/* Elevation [m] of the meteo station */
	st->sky=new_doublevector(INPUTmeteo->nrh);/* Sky-view-factor [-] of the meteo station */
	st->ST=new_doublevector(INPUTmeteo->nrh);/* Standard time minus UTM [hours] of the meteo station */
	st->Vheight=new_doublevector(INPUTmeteo->nrh); /* Wind velocity measurement height [m] (a.g.l.)  */
	st->Theight=new_doublevector(INPUTmeteo->nrh); /* Air temperature measurement height [m] (a.g.l.)  */
	st->JD0=new_doublevector(INPUTmeteo->nrh);/* Decimal Julian Day of the first data */
	st->Y0=new_longvector(INPUTmeteo->nrh);/* Year of the first data */
	st->Dt=new_doublevector(INPUTmeteo->nrh);/* Dt of sampling of the data [sec]*/
	st->offset=new_longvector(INPUTmeteo->nrh);/* offset column*/

	for(i=1;i<=INPUTmeteo->nrh;i++){
		st->E->co[i]=INPUTmeteo->co[i][1];
		st->N->co[i]=INPUTmeteo->co[i][2];
		st->lat->co[i]=INPUTmeteo->co[i][3]*Pi/180.0;/* from deg to [rad] */
		st->lon->co[i]=INPUTmeteo->co[i][4]*Pi/180.0;/* from deg to [rad] */
		st->Z->co[i]=INPUTmeteo->co[i][5];
		st->sky->co[i]=INPUTmeteo->co[i][6];
		st->ST->co[i]=INPUTmeteo->co[i][7];
		st->Vheight->co[i]=INPUTmeteo->co[i][8];
		st->Theight->co[i]=INPUTmeteo->co[i][9];
		st->JD0->co[i]=INPUTmeteo->co[i][10];
		st->Y0->co[i]=(long)(INPUTmeteo->co[i][11]);
		st->Dt->co[i]=INPUTmeteo->co[i][12]*3600;		//from hours to seconds
		st->offset->co[i]=(long)(INPUTmeteo->co[i][13]);
	}

}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
void ReadMeteoHeader(FILE *f, STRINGBIN *ColDescr, long offset, long *ncols, long *MeteoCont){
	/* Author: Stefano Endrizzi, 2008
	 * reads the header of the meteo file and finds the number of variables (ncols)
	 * Input:
	 * f: file of the meteo station
	 * ColDescr: vector of strings read in the parameter file, block n.12
	 * offset: offset column as read in the parameter file, block 11
	 * Output:
	 * ncols: number of strings in the header of the meteo file
	 * MeteoCont: vector of long representing the vector of the columns of the met structure (met->column): is equal
	 * to the number of strings in parameters file block 12. The vector is initialized with -1. If the name of the meteo variable
	 * read in the meteo file corresponds to the allowed meteo variables given in ColDescr, then the corresponding column is
	 * set with the column of the meteo file
	 * Comment: Matteo Dall'Amico, april 2009 */
	char **a;
	long i,j,num_meteo_var;
	num_meteo_var=dim1l(MeteoCont);
	for(i=0;i<num_meteo_var;i++){
		MeteoCont[i]=-1;
	}

	a=readline_textarray(f, offset); /* vector of strings containing the names of meteo variables read in the header of the meteo file */

	*ncols=dim_vect_strings(a); /* number of elements (strings) in the header of the meteo file */

	for(j=1;j<=num_meteo_var;j++){/* for every allowed meteo string given in parameter file block 12 */
		for(i=1;i<=*ncols;i++){/* for every meteo strings in the meteo file */
			//printf("%s %s\n",ColDescr->co[j]+1, a[i-1]);stop_execution();
			if(compare_strings(ColDescr->co[j]+1, a[i-1])==1 && MeteoCont[j-1]==-1){
				/* checks that the description string in the meteo file corresponds with the string in parameters file block 12 which is the same in constant.h  */
				MeteoCont[j-1]=i-1;// it means that if the SW in the meteo file is the third string, then at the 7th position of the vector (it is the 7th string in constant.h) there will be written 2
				//printf("j:%ld dim:%ld ncol:%ld col:%ld\n",j,dim1l(MeteoCont),*ncols,MeteoCont[j-1]);
			}else if(compare_strings(ColDescr->co[j]+1, a[i-1])==1 && MeteoCont[j-1]!=-1){ /*it means that in the meteo file there is a string that has been repeated */
				printf("Column '%s' is present twice\n",ColDescr->co[j]+1);
				t_error("Meteo Column name DUPLICATED!");
			}

		}
	}

	/* added by Emanuele Cordano on 22 September 2009 (bad deallocation)
	i=0;
	while (a[i]!=NULL) {
		free(a[i]);
		i++;
	}*/
	free(a);
}




/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
double **read_datameteo(FILE *f, long offset, long ncols, double ndef, int numlines){
	/* Author: Stefano Endrizzi, 2008
	 * reads the data of the meteo file
	 * Input:
	 * f: 		file of the meteo station
	 * offset:	offset column as read in the parameter file, block 11
	 * ncols:	number of columns present in the meteo file
	 * ndef:	novalue
	 * numlines: number of lines of meteo data
	 * Output:
	 * a: matrix (n X ncols) with the meteo values with n=# of rows in the file of the meteo station
	 * Comment: Matteo Dall'Amico, april 2009
	 * Modified by Matteo Dall'Amico on september 2009 in Zurich*/
	double **a=NULL;
	long i;//number of lines of data in the meteo file
	//long k;// counter
	short end;// novalueend=1;

	//printf("num lines=%d, offset=%ld,ncols=%ld, ndef=%f",numlines,offset,ncols,ndef);stop_execution();
	a=alloc2(numlines,ncols);
	for (i=0;i<numlines;i++) {
		readline_array(f, a[i], offset, ncols, ndef, &end);
	}
	//for(k=0;k<=ncols;k++){printf("a[%ld][%ld]=%f\n",numlines-1,k,a[numlines-1][k]);}printf("\na[%ld][0]=%f",numlines,a[numlines][0]);stop_execution();
	return(a);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

/***********************************************************/
/***********************************************************/

void read_inpts_par(PAR *par, TIMES *times, char *filename, char *ext, char *pos){
// reads __ocontrol_parameters.txt
DOUBLEVECTOR *V;
double Dt_output;
FILE *fd;
short  index;
char *temp; /* added by Emanuele Cordano on 21 September 2009 */
printf("\nENTERING SEVERAL CONTROL PROGRAM PAR \n");
temp=join_strings(filename,ext);
fd=t_fopen(temp,"r");
free(temp);
index=read_index(fd,PRINT);
V=read_doublearray(fd,PRINT);
t_fclose(fd);

//	V=read_parameters("",program, ext, pos);

printf("ENTER THE INTEGRATION INTERVAL [s]: %f\n",V->co[1]);
par->Dt=(double)V->co[1];

printf("ENTER THE Decimal julian day at the beginning of simulation (0.0 - 365.99): %f\n",V->co[2]);
par->JD0=(double)V->co[2];

printf("ENTER THE YEAR at the beginning of simulation (0.0 - 365.99): %f\n",V->co[3]);
par->year0=(long)V->co[3];

printf("ENTER THE NUMBER OF DAYS OF SIMULATION: %f\n",V->co[4]);
times->TH=(double)V->co[4];
times->TH*=24; /*TH in hour*/

printf("ENTER THE Standard time to which all the output data are referred (difference respect UMT, in hour): %f\n",V->co[5]);
par->ST=(double)V->co[5];

printf("ENTER THE Delta TIME (in hour) with which THE OUTPUT FOR SPECIFIED PIXELS IS PRINTED: %f\n",V->co[6]);
Dt_output=(double)V->co[6];
if(Dt_output*3600<par->Dt) Dt_output=par->Dt/3600.0;
times->n_pixel=(long)(Dt_output*3600/(long)par->Dt);
times->i_pixel=0;/*counter for the output of a pixel*/
times->count_super_printed=0;

printf("ENTER THE NUMBER OF Dt AFTER WHICH THE OUTPUT FOR THE BASIN ARE PRINTED: %f\n",V->co[7]);
Dt_output=(double)V->co[7];
if(Dt_output*3600<par->Dt) Dt_output=par->Dt/3600.0;
times->n_basin=(long)(Dt_output*3600/(long)par->Dt);
times->i_basin=0;/*counter for the output of a pixel*/

printf(">=1 if you want to display RESULTS for altimetric stripes, it is how many altimetric stripes you want to consider (up to 99): %f\n",V->co[8]);
par->ES_num=(short)V->co[8];
if(par->ES_num>99 || par->ES_num<-99) t_error("Number of altimetric stripes not valid");
if(par->ES_num<0){
par->glac_thr=0.1;
par->ES_num*=(-1);
}else if(par->ES_num>0){
par->glac_thr=-0.1;
}

printf("Multiplying factor decreasing the dem resolution for the calculation of the sky view factor: %f\n",V->co[9]);
par->nsky=(short)V->co[9];
if(par->nsky<1) t_error("Multiplying factor for sky view factor calculation must be greater than or equal to 1: %f");

printf("Thershold for the definition of channel network (in pixel number draining): %f\n",V->co[10]);
par->channel_thres=(double)V->co[10];

printf("OUTPUT MAPS in fluidturtle format (=1), GRASS ASCII (=2), ESRI ASCII (=3): %f\n",V->co[11]);
par->format_out=(short)V->co[11];

printf("=1 for one point simulation, =0 for distributed simulation. The parameter file is different in these cases: %f\n",V->co[12]);
par->point_sim=(short)V->co[12];

printf("=1 if you want to recover a simulation, 0 otherwise: %f\n",V->co[13]);
par->recover=(short)V->co[13];

printf("\n");

	free_doublevector(V);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

void read_soil_parameters(char *name, long *nsoil, SOIL *sl){
/* function that reads the soil parameter file _soil.txt and create the doubletensor of soil parameters
* the doubletensor is composed by three indexes: [sy][sprop][slay]:
* sy is the number of soil types
* sprop is the number of soil properties considered (see constant.h)
* [slay] is the number of soil layer considered
* Input: name of the file to read
* Output: number of soil layers (*nsoil) and doubletensor sl->pa
* Author: Stefano Endrizzi
* Comments: Matteo Dall'Amico, April 2009 */
	FILE *f;
	long index, i, j, k;
	DOUBLEMATRIX *M;
	char *temp; /* added by Emanuele Cordano on 22 September 2009 */
	temp=join_strings(name,textfile);
	f=t_fopen(temp,"r");
	free(temp);
	index=read_index(f,PRINT);

	printf("%ld Different soil types given\n",index);

	M=read_doublematrix(f,"a",PRINT);
	Nl=M->nrh;

	sl->pa=new_doubletensor(index, nsoilprop, Nl);
	initialize_doubletensor(sl->pa,0.0);

	for(i=1;i<=index;i++){// for each soil type

		if(i>1){
			M=read_doublematrix(f,"a",PRINT);
			if(M->nrh!=Nl) t_error("The number of sl layers must be the same for each sl type");
		}

		if(M->nch!=nsoilprop)t_error("Wrong column number in the parameter file");

		for(j=1;j<=nsoilprop;j++){// for each soil property
			for(k=1;k<=M->nrh;k++){// for each layer
				sl->pa->co[i][j][k]=M->co[k][j];
				if(j==jdz && i!=1){
					if(sl->pa->co[i][j][k]!=sl->pa->co[i-1][j][k]) t_error("Soil layer thicknesses must be the same for each sl type");
				}
				if(sl->pa->co[i][jlatfl][k]<0 || sl->pa->co[i][jlatfl][k]>2) t_error("Value not admitted of jlatfl in soil parameters");
			}
		}
		free_doublematrix(M);

	}

	t_fclose(f);

	*nsoil=index;
}


//***************************************************************************
//Check DD

DOUBLEMATRIX *depitted(SHORTMATRIX *DD, DOUBLEMATRIX *Z){

	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; /* differential of number-pixel for rows and*/
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; /* columns, depending on Drainage Directions*/
	DOUBLEMATRIX *M;
	long r,c,cont=0;
	short a;

	M=new_doublematrix(Z->nrh,Z->nch);
	for(r=1;r<=Z->nrh;r++){
		for(c=1;c<=Z->nch;c++){
			M->co[r][c]=Z->co[r][c];
		}
	}

	do{
		a=0;
		cont++;
		for(r=1;r<=Z->nrh;r++){
			for(c=1;c<=Z->nch;c++){
				if(DD->co[r][c]>=1 && DD->co[r][c]<=8){
					if(M->co[r][c]<M->co[r+r_DD[DD->co[r][c]]][c+c_DD[DD->co[r][c]]]){
						M->co[r+r_DD[DD->co[r][c]]][c+c_DD[DD->co[r][c]]]=M->co[r][c];
						a=1;
					}
				}
			}
		}
		if(cont>50) printf("DEPITTING: iteration #%ld\n",cont);
	}while(a==1);

	return(M);

}

//***************************************************************************
void assign_recovered(char *name, double **assign, PAR *par, DOUBLEMATRIX *Z1, DOUBLEMATRIX *Z2){

	long r, c, i;
	DOUBLEMATRIX *M;

	if(par->point_sim==0){

		M=read_map(2, name, Z1, UV);

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				assign[r][c]=M->co[r][c];
			}
		}
	}else{
		printf("%ld %ld %f %f %f %f %f %f\n",Z2->nrh,Z2->nch,UV->U->co[1],UV->U->co[2], UV->U->co[3],UV->U->co[4],UV->V->co[1],UV->V->co[2]);
		M=read_map(2, name, Z2, UV);
		for(i=1;i<=par->r_points->nh;i++){
			r=par->r_points->co[i];
			c=par->c_points->co[i];
			assign[1][i]=M->co[r][c];
		}
	}
	free_doublematrix(M);
}


//***************************************************************************
void assign_recovered_long(char *name, long **assign, PAR *par, DOUBLEMATRIX *Z1, DOUBLEMATRIX *Z2){

	long r, c, i;
	DOUBLEMATRIX *M;

	if(par->point_sim==0){
		M=read_map(2, name, Z1, UV);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				assign[r][c]=(long)M->co[r][c];
			}
		}
	}else{
		M=read_map(2, name, Z2, UV);
		for(i=1;i<=par->r_points->nh;i++){
			r=par->r_points->co[i];
			c=par->c_points->co[i];
			assign[1][i]=(long)M->co[r][c];
		}
	}
	free_doublematrix(M);
}
//***************************************************************************

void i_lrc_cont(DOUBLEMATRIX *LC, long ***i, LONGMATRIX *lrc){

	long cont=0;
	long l, r, c;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(LC->co[r][c]!=NoV){
				for(l=1;l<=Nl;l++){
					cont++;
					i[l][r][c]=cont;
					lrc->co[cont][1]=l;
					lrc->co[cont][2]=r;
					lrc->co[cont][3]=c;
				}
			}
		}
	}
}
