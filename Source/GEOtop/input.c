
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
#include "input.h"
#include "geomorphology.0875.h"
#include "pedo.funct.h"
#include "geo_statistic.h"
#include "networks.h"
#include "constant.h"
#include "keywords_file.h"
#include "dtm_resolution.h"
#include "rw_maps.h"
#include "extensions.h"
#include "tabs.h"
#include "snow.h"
#include "micromet.h"
#include "vegetation.h"
#include "get_filenames.h"
#include "output.h"

#include <unistd.h>

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;
extern char *MISSING_FILE;


/*--------------------------------------------------------------------------------------------------*/
//! Subroutine which reads input data, performs  geomporphological analisys and allocates data
void get_all_input(int argc, char *argv[], TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet, 
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times)

{

 /*counters of layers(l),rows(r),columns(c) and other things(i) (internal variables):*/
 long l,r,c,i,j,index,ncols,iland,cont,rnext,cnext;
 /*auxiliary internal variables:*/
// static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0};/*differential of number-pixel for rows and*/
// static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0};/*columns, depending on Drainage Directions*/
 FILE *f;
 /*auxiliary internal vectors-variables:*/
 DOUBLEMATRIX *SNOW0, *GLACIER0, *h_channel, *bedrock_depth;
 DOUBLETENSOR *T, *pa_bed;
 //Checking variables
 long sy, synew, offset;
 short a;
 double z, zlim, th_oversat, Zaverage;
 INIT_TOOLS *IT;
 char *temp, *temp2, **temp3;
 char SSSS[ ]={"SSSS"};

 
 IT=(INIT_TOOLS *)malloc(sizeof(INIT_TOOLS));
 if(!IT) t_error("IT was not allocated");

/****************************************************************************************************/
/*! Reading of program control par from geo_top.inpts file											*/
/*  (if it is in WORKING_DIRECTORY), otherwise from prompt        									*/
/****************************************************************************************************/
	
printf("STATEMENT:\n\n");
	 
printf("GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE\n");
printf("GEOtop 1.0 Public - Version Montebello - Update 2 (29 April 2010)\n\n");
	 
printf("Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon\n\n");
	 
printf("This file is part of GEOtop 1.0 Public\n\n");
	 
printf("GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
printf("WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE\n\n");
	 
printf("GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.\n");
printf("If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/\n");
printf("Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.\n\n");
	 
printf("If you have satisfactorily used the code, please acknowledge the authors\n\n");
	 
i=usleep(10);	
	
	
if(!argv[1]){
	WORKING_DIRECTORY=get_workingdirectory();
}else{
	WORKING_DIRECTORY=argv[1];
}

printf("\nWORKING DIRECTORY: %s\n",WORKING_DIRECTORY);

//READ FILES
//reads the keywords in __geotop.init and reads the filenames from __geotop.inpts 
files=get_filenames_from_keys(WORKING_DIRECTORY, PROGRAM_NAME, PRINT); //READ INPUT PAR
MISSING_FILE = join_strings(WORKING_DIRECTORY, LOCAL_MISSING_FILE);
		
//reads the parameters in __control_parameters
read_inpts_par(par, times, files->co[fcontrol]+1, textfile); 

if(files->index->nh!=nfiles){
	for(i=1;i<=files->index->nh;i++){
		printf("i:%ld file:%s\n",i,files->co[i]+1);
	}
	printf("files given:%ld files required %d\n",files->index->nh,nfiles);
	t_error("Wrong number of files in .inpts file");
}else if(files->index->nh>nfiles){
	printf("Warning: Excess files given:%ld files required:%d\n",files->index->nh,nfiles);
}
	
files->co[ferr]=join_strings(files->co[ferr]+1,textfile)-1;
f=t_fopen(files->co[ferr]+1,"w");
fprintf(f,"Files of the errors in the simulation\n");
t_fclose(f);

/****************************************************************************************************/
/*! Reading of the Input files:                                                                      */
/****************************************************************************************************/
read_parameterfile(files->co[fpar]+1, par, IT);	
sl->pa = read_soil_parameters(files->co[fspar]+1);
	
if(par->point_sim!=1){ /* distributed simulation */	
	read_inputmaps(top, land, sl, par);
	read_optionsfile_distr(files->co[fopt]+1, par, times, land->LC);
}else{
	read_optionsfile_point(files->co[fopt]+1, par, top, land, sl, IT);
}
	

times->time=0.0;
Nr=top->Z0->nrh;
Nc=top->Z0->nch;
NoV=UV->V->co[2];
Nl = sl->pa->nch;

par->total_pixel=0;
for(r=1;r<=Nr;r++){
    for(c=1;c<=Nc;c++){
       if (land->LC->co[r][c]!=NoV) par->total_pixel++;
    }
}

//bedrock
if( strcmp(files->co[fbed]+1 , MISSING_FILE) != 0 ){
	
	if(existing_file(files->co[fbed]+1)==0){
		printf("File %s is missing. Please check if you have a bedrock topography map. If it is not available, remove the file name and keyword from .inpts file\n",files->co[fbed]+1);
		t_error("Check input files");
	}
	
	printf("A bedrock depth map has been assigned and read from %s\n\n",files->co[fbed]+1);

	f=fopen(files->co[ferr]+1, "a");
	fprintf(f,"A bedrock depth map has been assigned and read from %s\n\n",files->co[fbed]+1);
	fclose(f);
	
	par->bedrock = 1;
	bedrock_depth = read_map(1, files->co[fbed]+1, land->LC, UV);
	
	if(existing_file_text(files->co[fspar2]+1)==0){
		printf("File %s is missing. It is necessary if you give the bedrock topography.\n",files->co[fspar2]+1);
		t_error("Check input files");
	}
	
	pa_bed = read_soil_parameters(files->co[fspar2]+1);
	
	if(pa_bed->ndh != sl->pa->ndh) t_error("Soil parameters bedrock must have the same classes of Soil parameters");
	if(pa_bed->nch != sl->pa->nch) t_error("Soil parameters bedrock must have the same number of layers of Soil parameters");
	for(i=1;i<=pa_bed->ndh;i++){
		for(l=1;l<=Nl;l++){
			if(pa_bed->co[i][jdz][l] != sl->pa->co[i][jdz][l]){
				printf("Soil parameters bedrock must have the same layer thicknesses of Soil parameters,l:%ld soiltype:%ld in the soil d:%f in the bedrock d:%f\n",l,i,sl->pa->co[i][jdz][l],pa_bed->co[i][jdz][l]);
				t_error("Review soil files");
			}
		}
	}
	
	T=new_doubletensor(sl->pa->ndh, nsoilprop, Nl);
	for(i=1;i<=sl->pa->ndh;i++){
		for(j=1;j<=nsoilprop;j++){
			for(l=1;l<=Nl;l++){		
				T->co[i][j][l]=sl->pa->co[i][j][l];
			}
		}
	}
	free_doubletensor(sl->pa);

	sl->pa=new_doubletensor(par->total_pixel, nsoilprop, Nl);

	//assign jdz (is needed later)
	for(i=1;i<=sl->pa->ndh;i++){
		for(l=1;l<=Nl;l++){		
			sl->pa->co[i][jdz][l]=T->co[1][jdz][l];
		}
	}	
	
		
	cont=0;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				
				cont++;
			
				zlim = 1.E3*(top->Z0->co[r][c] - bedrock_depth->co[r][c]);
							
				sy = sl->type->co[r][c];
				synew = cont;
				sl->type->co[r][c] = synew;
		
				z=0.0;
				for(l=1;l<=Nl;l++){			
									
					z += 0.5*sl->pa->co[synew][jdz][l]*cos(top->slopes->co[r][c]);
																			
					if(z <= zlim){
					
						for(j=1;j<=nsoilprop;j++){
							sl->pa->co[synew][j][l] = T->co[sy][j][l];
						}
						
					}else{
					
						for(j=1;j<=nsoilprop;j++){
							sl->pa->co[synew][j][l] = pa_bed->co[sy][j][l] ;
						}				
						
					}
					
					z += 0.5*sl->pa->co[synew][jdz][l]*cos(top->slopes->co[r][c]);
														
				}
			}
		}
	}

	free_doubletensor(T);
	free_doublematrix(bedrock_depth);
	free_doubletensor(pa_bed);
		
}else{

	par->bedrock = 0;

}

/****************************************************************************************************/
//Reading of RAIN data file,	METEO data file, 	and CLOUD data file

//METEO DATA
//meteo parameters
f=t_fopen(join_strings(files->co[fmet]+1,textfile),"r");
index=read_index(f,PRINT);
IT->met=read_doublematrix(f,"a",PRINT);
if(IT->met->nch!=13) t_error("Error in meteo parameters");
IT->met_col_names=read_stringarray(f,PRINT);
t_fclose(f);

//meteo_stations
met->st=(METEO_STATIONS *)malloc(sizeof(METEO_STATIONS));	
if(!met->st) t_error("meteo_stations was not allocated"); 
init_meteo_stations(IT->met, met->st);

//meteo data
met->column=alloc_long2(IT->met->nrh);
met->data=(double ***)malloc(IT->met->nrh*sizeof(double**));
met->horizon=(double ***)malloc(IT->met->nrh*sizeof(double**));
met->var=alloc2(met->st->Z->nh,nmet);	//meteo variables for the current instant

for(i=1;i<=IT->met->nrh;i++){
	temp=namefile_i(files->co[fmet]+1, i);
	f=t_fopen(temp,"r");
	met->column[i-1]=alloc_long1(nmet);
	ReadMeteoHeader(f, IT->met_col_names, met->st->offset->co[i], &ncols, met->column[i-1]);
	met->data[i-1]=read_datameteo(f, met->st->offset->co[i], ncols, UV->V->co[2]);
	t_fclose(f);
	free(temp);
	met->horizon[i-1]=read_horizon(files->co[fhor]+1, i);
}


//read LAPSE RATES FILE  //0.Year 1.JD 2.Tair 3.Tdew 4.precipitation
met->LRp=new_longvector(3);	//column referring to the LRs - p stands for position
//year column 0 - JD column 1
met->LRp->co[1]=2;	//Column 2: Lapse rate for Tair (dT/dz in [C/m])
met->LRp->co[2]=3;	//Column 3: Lapse rate for Tdew (dT/dz in [C/m])
met->LRp->co[3]=4;	//Column 4: Lapse rate for prec (dP/dz in [mm/m])
if(strcmp(files->co[fLRs]+1 , MISSING_FILE) != 0){   //s stands for string
	if(existing_file_text(files->co[fLRs]+1)==0) printf("Lapse rate file unavailable. Check input files. If you do not have a lapse rate file, remove its name and keywords from .inpts file\n");
	f=t_fopen(join_strings(files->co[fLRs]+1,textfile),"r");
	temp3=readline_textarray(f, 0);
	met->LRs=read_datameteo(f, 0, 5, UV->V->co[2]);
	t_fclose(f);
	par->LRflag=1;
	for(i=1;i<=5;i++){
		free(temp3[i-1]);
	}
	free(temp3);
	printf("\nLapse rate file read\n");
}else{
	par->LRflag=0;
	printf("\nWARNING: LAPSE RATE FILE IS MISSING!!!\n");
	printf("BE AWARE OF THAT!!!! The default values will be used\n");
}
met->LRv=alloc1(5);
for(i=1;i<=3;i++){
	met->LRv[met->LRp->co[i]]=NoV;
}


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

par->n_landuses=0;
for(i=1;i<=IT->land_classes->nch;i++){
	if(par->n_landuses<(long)IT->land_classes->co[1][i]) par->n_landuses=(long)IT->land_classes->co[1][i];
}
printf("\nNumber of land use types: %ld\n",par->n_landuses);

//properties for each land use
land->ty=new_doublematrix(par->n_landuses,nlandprop);
initialize_doublematrix(land->ty, 0.0);

z = 0.;
l = 0;
do{	
	l++;
	z += sl->pa->co[1][jdz][l];
}while(l<Nl && z < z_transp);
land->root_fraction=new_doublematrix(par->n_landuses, l);
initialize_doublematrix(land->root_fraction, 0.0);

land->vegparp=new_longmatrix(par->n_landuses,jdvegprop);
land->vegpars=(double ***)malloc((par->n_landuses+1)*sizeof(double**));
land->vegparv=(double **)malloc((par->n_landuses+1)*sizeof(double*));
land->vegpar=new_doublevector(jdvegprop);

//check vegetation variable consistency
if(jHveg!=jdHveg+jHveg-1) t_error("Vegetation variables not consistent");
if(jz0thresveg!=jdz0thresveg+jHveg-1) t_error("Vegetation variables not consistent");
if(jz0thresveg2!=jdz0thresveg2+jHveg-1) t_error("Vegetation variables not consistent");
if(jLSAI!=jdLSAI+jHveg-1) t_error("Vegetation variables not consistent");
if(jcf!=jdcf+jHveg-1) t_error("Vegetation variables not consistent");
if(jdecay0!=jddecay0+jHveg-1) t_error("Vegetation variables not consistent");
if(jexpveg!=jdexpveg+jHveg-1) t_error("Vegetation variables not consistent");
if(jroot!=jdroot+jHveg-1) t_error("Vegetation variables not consistent");
if(jrs!=jdrs+jHveg-1) t_error("Vegetation variables not consistent");

par->vegflag=new_shortvector(par->n_landuses);	

offset=2;//Julian day and year
	
//time dependent vegetation parameters	
for(i=1;i<=par->n_landuses;i++){
	
	//vegparp[vegetationtype][vegetationvariable] represents the column of the veg file, so vegpars[vegetationtype][time][vegparp]
	for(j=1;j<=jdvegprop;j++){
		land->vegparp->co[i][j]=j+offset-1;
	}
	
	if(strcmp(files->co[fvegpar]+1 , MISSING_FILE) != 0){   //s stands for string
		
		write_suffix(SSSS, i, 0);
		temp=join_strings(files->co[fvegpar]+1, SSSS);	
		
		if(existing_file_text(temp)==1){ 
			printf("There is a specific vegetation parameter file for land cover type = %ld\n",i);
			temp2=join_strings(temp,textfile);
			f=t_fopen(temp2,"r");
			temp3=readline_textarray(f, 0);
			land->vegpars[i]=read_datameteo(f, 0, jdvegprop+2, UV->V->co[2]);
			t_fclose(f);
			for(j=1;j<=jdvegprop+2;j++){
				free(temp3[j-1]);
			}
			free(temp3);
			free(temp2);
			par->vegflag->co[i]=1;
		}else{
			printf("There is NOT a specific vegetation parameter file for land cover type = %ld\n",i);
			par->vegflag->co[i]=0;
		}
		
		free(temp);
	
	}else{
		
		par->vegflag->co[i]=0;
		
	}
	
	land->vegparv[i]=alloc1(jdvegprop+offset);
	
	for(c=1;c<=jdvegprop;c++){
		land->vegparv[i][land->vegparp->co[i][c]]=NoV;
	}

	for(j=1;j<=IT->land_classes->nch;j++){
		
		if(i==(long)IT->land_classes->co[1][j]){
		
			for(l=1;l<=nlandprop;l++){
				land->ty->co[i][l]=IT->land_classes->co[l+1][j];
			}

			//z0 (convert in m)
			land->ty->co[i][jz0]*=0.001;
						
			//find root fraction
			root(land->root_fraction->nch, land->ty->co[i][jroot], 0.0, sl->pa->co[1][jdz], land->root_fraction->co[i]);

			//check on vegetation height threshold
			/*if(land->ty->co[i][jz0thresveg2]>=land->ty->co[i][jz0thresveg]){
				printf("The thresholds on snow depth regarding vegetation snow burial are not set in the right way in land cover type %ld: thresveg:%f thresveg2:%f\n",
				   i,land->ty->co[i][jz0thresveg],land->ty->co[i][jz0thresveg2]);
				t_error("It has to be thresveg>thresveg2");
			}*/
			
			//error messages
			for(l=1;l<=met->st->Z->nh;l++){
				if(0.001*land->ty->co[i][jHveg]>met->st->Vheight->co[l] || 0.001*land->ty->co[i][jHveg]>met->st->Theight->co[l]){
					printf("hc:%f m, zmu:%f m, zmt:%f m - set hc lower than measurement height - land cover %ld, meteo station %ld\n",
						0.001*land->ty->co[i][jHveg],met->st->Vheight->co[l],met->st->Theight->co[l],i,l);
					t_error("ERROR 1");
				}
			}
		}
	}
}
/****************************************************************************************************/
// Completing of "top" (of the type TOPO):    


f=fopen(files->co[ferr]+1,"a");
fprintf(f,"Valid pixels: %ld\n",par->total_pixel);
fprintf(f,"Number of nodes: %ld\n",(Nl+1)*par->total_pixel);
fprintf(f,"Novalue pixels: %ld\n",(Nr*Nc-par->total_pixel));
fprintf(f,"Basin area: %f km2\n",(double)par->total_pixel*UV->U->co[1]*UV->U->co[2]/1.E6);
fclose(f);

Zaverage=0.0;
for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			Zaverage += top->Z0dp->co[r][c]/(double)par->total_pixel;
		}
	}
}		

top->Z=new_doubletensor0(Nl,Nr,Nc);
initialize_doubletensor(top->Z,NoV);

for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			
			top->Z0dp->co[r][c] -= Zaverage;
			
			sy=sl->type->co[r][c];
			z=1.E3*top->Z0dp->co[r][c];
			
			l=0;
			top->Z->co[l][r][c]=z;
			
			do{
				l++;
				z -= 0.5*sl->pa->co[sy][jdz][l]*cos(top->slopes->co[r][c]);
				top->Z->co[l][r][c]=z;
				z -= 0.5*sl->pa->co[sy][jdz][l]*cos(top->slopes->co[r][c]);
			}while(l<Nl);
		}
	}
}
	
//if(par->max_Courant_sup_sub<=1) searchDD(top->j_cont, top->rc_cont, top->DD, land->LC, top->DDup, top->DDdown);
	
/****************************************************************************************************/
/*! Filling up of the struct "channel" (of the type CHANNEL):                                        */

/*The number of channel-pixel are counted:*/
i=0;
for(r=1;r<=Nr;r++){
   for(c=1;c<=Nc;c++){
      if (top->pixel_type->co[r][c]>=10) i++;
   }
}
f=fopen(files->co[ferr]+1,"a");
fprintf(f,"Channel pixels: %ld\n",i);
fclose(f);
par->total_channel = i;

/* Creation of the vectors with the position of the channel-pixels ("cnet->r" and "cnet->c"):*/
if(i==0) i=1;

cnet->r=new_longvector(i);
initialize_longvector(cnet->r, 0);
cnet->c=new_longvector(i);
initialize_longvector(cnet->c, 0);
cnet->ch=new_longmatrix(Nr,Nc);
initialize_longmatrix(cnet->ch, 0);
	
cnet->length=new_doublevector(i);
initialize_doublevector(cnet->length, 0.);
	
if(i>1){
	i=0;
	do{
			
		find_max_constraint( top->Z0dp, land->LC, top->pixel_type, cnet->ch, &r, &c);
		if(r>0){

			i++;
			cnet->r->co[i]=r;
			cnet->c->co[i]=c;
			cnet->ch->co[r][c]=i;
			//printf("__________________________________________________________________________________\n");
			//printf("r:%ld c:%ld Z:%f i:%ld p:%d ",r,c,top->Z0dp->co[r][c],i,top->pixel_type->co[r][c]);

			do{
		
				next_down_channel_pixel( r, c, top->Z0dp, land->LC, top->pixel_type, cnet->ch, &rnext, &cnext);
				if(rnext>0){
					i++;
					if(fabs(rnext-r)==1 && fabs(cnext-c)==1){
						cnet->length->co[i-1]+=0.5*sqrt(2.);
						cnet->length->co[i]+=0.5*sqrt(2.);
					}else{
						cnet->length->co[i-1]+=0.5;
						cnet->length->co[i]+=0.5;				
					}
					//printf("length:%f\n",cnet->length->co[i-1]);
					cnet->r->co[i]=rnext;
					cnet->c->co[i]=cnext;
					cnet->ch->co[rnext][cnext]=i;
					r=rnext;
					c=cnext;
					//printf("r:%ld c:%ld Z:%f i:%ld p:%d ",r,c,top->Z0dp->co[r][c],i,top->pixel_type->co[r][c]);
				}
			
			}while(rnext>0);
			
			if(rnext<0){
				if(fabs(-rnext-r)==1 && fabs(-cnext-c)==1){
					cnet->length->co[i]+=0.5*sqrt(2.);
				}else{
					cnet->length->co[i]+=0.5;				
				}
			}else{
				cnet->length->co[i]+=0.5;				
			}
			//printf("length:%f\n",cnet->length->co[i]);

		
		}
				
	}while(r>0);
	
	for(i=1;i<=cnet->r->nh;i++){
		cnet->length->co[i] *= UV->U->co[1];
		//printf("%ld %f\n",i,cnet->length->co[i]);
	}
	

}else{

	cnet->r->co[1]=0;
	cnet->c->co[1]=0;

}

/* Initialization of the matrix with the sub/superficial flows which goes in a channel-cell ("cnet->q_sup"):*/
cnet->Qsup=new_doublevector(cnet->r->nh);
cnet->Qsub=new_doublevector(cnet->r->nh);

cnet->h_sup=new_doublevector(cnet->r->nh);
initialize_doublevector(cnet->h_sup, 0.0);
	
cnet->dh_sup=new_doublevector(cnet->r->nh);
	
if(par->output_h_sup>0){
	cnet->hsupav=new_doublematrix(Nr, Nc);	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(cnet->ch->co[r][c]>0){
				cnet->hsupav->co[r][c] = 0.;
			}else{
				cnet->hsupav->co[r][c] = NoV;
			}
		}
	}
}
	
if(par->recover==1){
	h_channel=new_doublematrix(Nr, Nc);		
	assign_recovered(files->co[rQch]+1, h_channel->co, par, land->LC, IT->LU);
	for(i=1;	i<=cnet->r->nh;i++){
		if(cnet->r->co[i]>0) cnet->h_sup->co[i] = h_channel->co[cnet->r->co[i]][cnet->c->co[i]];
	}
	free_doublematrix(h_channel);
}

	
/* Extraction of the vector of channel-pixels distances ("cnet->s0") from the matrix of the distances of
   each pixel from the outlet ("top->pixel_distance")*/
//cnet->s0=new_doublevector(cnet->r->nh);
/*for(i=1;i<=cnet->r->nh;i++){
	if(cnet->r->co[i]>0){
		cnet->s0->co[i]=top->pixel_distance->co[cnet->r->co[i]][cnet->c->co[i]];
	}else{
		cnet->s0->co[i]=UV->U->co[1];
	}
   
	if(cnet->s0->co[i]<=0.0){
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"Warning: negative distance from outlet of the channel pixel %ld, %ld\n",cnet->r->co[i],cnet->c->co[i]);
		fclose(f);
	}
}
if(par->print==1){
   f=fopen(files->co[ferr]+1,"a");
   for(i=1;i<=cnet->r->nh;i++){
      fprintf(f,"cnet->s0->co[i] = %f con i = %ld r = %ld c = %ld\n",cnet->s0->co[i],i,cnet->r->co[i],cnet->c->co[i]);
   }
   fclose(f);
}*/

/* Creation of the matrix with the coeficent to spread the channel-flow for each channel-pixel ("cnet->fraction_spread"):*/

//cnet->fraction_spread=De_Saint_Venant(cnet->s0, IT->u0, IT->D, par->Dt);

/* Initialization of the vector with channel-flow (derived from q_sup) for each virtual channel-pixel with the
   same distance from outlet ("cnet->Q_sup_s"); note: vector's dimension is the number of virtual stretches
   of channel (starting from outlet):*/
//cnet->Q_sup_s=new_doublevector(cnet->fraction_spread->nch);
//initialize_doublevector(cnet->Q_sup_s,0.0);

//cnet->Qsup_spread=new_doublevector(cnet->fraction_spread->nch);

/* Initialization of the vector with channel-flow (derived from q_sub) for each virtual channel-pixel with the
   same distance from outlet ("cnet->Q_sub_s"); note: vector's dimension is the number of virtual stretches
   of channel (starting from outlet):*/
//cnet->Q_sub_s=new_doublevector(cnet->fraction_spread->nch);
//initialize_doublevector(cnet->Q_sub_s,0.0);

//cnet->Qsub_spread=new_doublevector(cnet->fraction_spread->nch);

/*if(par->recover==1){
	if(existing_file_text(files->co[rQch]+1)==1){
		f=t_fopen(join_strings(files->co[rQch]+1,textfile),"r");
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
}*/


/*
distance_from_channel(M, top->DD, top->pixel_type);
distance_from_channel2(top->dist_channel, top->pixel_type, cnet->r, cnet->c);
for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			top->dist_channel->co[r][c]*=UV->U->co[1];
		}
	}
}
temp=join_strings(WORKING_DIRECTORY,"0dist_from_channel");
write_map(temp, 0, par->format_out, top->dist_channel, UV);	
free(temp);*/

//Cont for Richards 3D
top->i_cont=(long ***)malloc((Nl+1)*sizeof(long**));
for(l=0;l<=Nl;l++){
	top->i_cont[l]=(long **)malloc((Nr+1)*sizeof(long*));
	for(r=1;r<=Nr;r++){
		top->i_cont[l][r]=(long *)malloc((Nc+1)*sizeof(long));
	}
}
top->lrc_cont=new_longmatrix( (Nl+1)*par->total_pixel + cnet->r->nh , 3);
initialize_longmatrix(top->lrc_cont, 0);

i_lrc_cont(land->LC, top->i_cont, top->lrc_cont);

//Cont for overland flow
top->j_cont=(long **)malloc((Nr+1)*sizeof(long*));
for(r=1;r<=Nr;r++){
	top->j_cont[r]=(long *)malloc((Nc+1)*sizeof(long));
}

top->rc_cont=new_longmatrix(par->total_pixel,2);
j_rc_cont(land->LC, top->j_cont, top->rc_cont);

/****************************************************************************************************/
/*! Completing of the initialization of SOIL structure                               */
/****************************************************************************************************/

/****************************************************************************************************/
/*! Initialization of the sl temperature tensor "T" [C]: */

/****************************************************************************************************/
/*! Completing of "sl" (of the type SOIL):                                                    */

sl->th=new_doubletensor(Nl,Nr,Nc);
initialize_doubletensor(sl->th,0.0);

sl->thice=new_doubletensor(Nl,Nr,Nc);
initialize_doubletensor(sl->thice,0.0);
	
sl->ET=new_doubletensor(Nl,Nr,Nc);
initialize_doubletensor(sl->ET,0.0);

sl->P=new_doubletensor0(Nl,Nr,Nc);
initialize_doubletensor(sl->P,0.0);

sl->Ptot=new_doubletensor(Nl,Nr,Nc);
initialize_doubletensor(sl->Ptot,0.0);
	
sl->T=new_doubletensor(Nl,Nr,Nc);
initialize_doubletensor(sl->T,0.0);

sl->Tv=new_doublematrix(Nr,Nc);
initialize_doublematrix(sl->Tv,0.0);

if(par->state_pixel==1){
	sl->T_av=new_doublematrix(Nl, par->chkpt->nrh);
	sl->th_av=new_doublematrix(Nl, par->chkpt->nrh);
	sl->thice_av=new_doublematrix(Nl, par->chkpt->nrh);
}

for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			for(i=1;i<=sl->pa->ndh;i++){
				if(sl->type->co[r][c]==i){
					sy=sl->type->co[r][c];
										
					//sl->Tv->co[r][c]=sl->pa->co[sy][jT][1];
					
					for(l=1;l<=Nl;l++){
						
						sl->P->co[l][r][c]=sl->pa->co[sy][jpsi][l];
						sl->Ptot->co[l][r][c]=sl->pa->co[sy][jpsi][l];
						sl->T->co[l][r][c]=sl->pa->co[sy][jT][l];
						
						sl->th->co[l][r][c] = teta_psi(sl->P->co[l][r][c], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], 
														sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
						
						th_oversat = Fmax( sl->P->co[l][r][c] , 0.0 ) * sl->pa->co[sy][jss][l];
						sl->th->co[l][r][c] -= th_oversat;
														
						if(sl->T->co[l][r][c]<=Tfreezing){
							//Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
							sl->thice->co[l][r][c] = sl->th->co[l][r][c] - teta_psi(Psif(sl->T->co[l][r][c]), 0.0, sl->pa->co[sy][jsat][l], 
								sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);

							//if Theta(without freezing)<Theta_unfrozen(in equilibrium with T) Theta_ice is set at 0
							if(sl->thice->co[l][r][c]<0) sl->thice->co[l][r][c]=0.0;
							
							//Psi is updated taking into account the freezing
							sl->th->co[l][r][c] -= sl->thice->co[l][r][c];
							sl->P->co[l][r][c] = psi_teta(sl->th->co[l][r][c] + th_oversat, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l],
														  sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 
														  1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
						}
					}
				}
			}
			
		}else{
			
			sl->P->co[0][r][c]=NoV;
			
			for(l=1;l<=Nl;l++){
				sl->P->co[l][r][c]=NoV;
				sl->T->co[l][r][c]=NoV;
				sl->thice->co[l][r][c]=NoV;
				sl->th->co[l][r][c]=NoV;	
				sl->Ptot->co[l][r][c]=NoV;
			}
			
		}
	}	
}

	
if(par->state_pixel==1){
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				//average T theta and thetaice
				for(i=1;i<=par->chkpt->nrh;i++){		
				
					if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){
						for(l=1;l<=Nl;l++){
							sl->th_av->co[l][i] = sl->th->co[l][r][c];
							sl->thice_av->co[l][i] = sl->thice->co[l][r][c];
							sl->T_av->co[l][i] = sl->T->co[l][r][c];
						}
					}
				}
			}
		}
	}
}

if(par->recover==1){
	assign_recovered(files->co[rTv]+1, sl->Tv->co, par, land->LC, IT->LU);
	
	for(l=0;l<=Nl;l++){
		assign_recovered(namefile_i_we(files->co[rpsi]+1, l), sl->P->co[l], par, land->LC, IT->LU);
	}
	
	for(l=1;l<=Nl;l++){
		assign_recovered(namefile_i_we(files->co[riceg]+1, l), sl->thice->co[l], par, land->LC, IT->LU);
		assign_recovered(namefile_i_we(files->co[rTg]+1, l), sl->T->co[l], par, land->LC, IT->LU);

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
					sl->th->co[l][r][c] = theta_from_psi( Fmin ( sl->P->co[l][r][c], psisat_from(l, r, c, sl) ), l, r, c, sl, PsiMin );
				}
			}
		}

	}
}

//WRITE INITIAL CONDITION
write_output_headers(met->st->Z->nh, times, wat, par, top, land, sl, egy, snow, glac);

if(par->state_pixel==1){
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){//if the pixel is not a novalue				
				for(j=1;j<=par->rc->nrh;j++){
					if(r==par->rc->co[j][1] && c==par->rc->co[j][2]) write_soil_output(0, j, times->time, 0.0, par->year0, par->JD0, par->rc, sl, PsiMin);
				}
			}
		}
	}
}
	
//z boundary condition
for(l=1;l<=Nl;l++){
	par->Zboundary -= sl->pa->co[1][jdz][l];
}
if(par->Zboundary < sl->pa->co[1][jdz][Nl]/2.){
	printf("Z at which 0 annual temperature takes place is not lower than the soil column in land cover type %ld\n",i);
	t_error("Please correct");
}
par->Zboundary *= 1.E-3;	//convert in [m]
	
/****************************************************************************************************/
/*! Initialization of the struct "egy" (of the type ENERGY):*/

 if(par->output_Rn>0){
	egy->Rn_mean=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Rn_min=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->Rn_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->LW_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->LW_min=new_doublematrix(Nr,Nc);
	egy->LW_in=new_doublematrix(Nr,Nc);
	egy->LW=new_doublematrix(Nr,Nc);
	egy->SW=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->SW_max=new_doublematrix(Nr,Nc);
 }
 if(par->output_ET>0){
	egy->ET_mean=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->ET_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->ET_min=new_doublematrix(Nr,Nc);
 }
 if(par->output_H>0){
	egy->H_mean=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->H_max=new_doublematrix(Nr,Nc);
	if(par->distr_stat==1)egy->H_min=new_doublematrix(Nr,Nc); 
 }
 if(par->output_G>0){
	egy->SEB_mean=new_doublematrix(Nr,Nc);
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
				egy->LW->co[r][c]=0.0;
				egy->SW->co[r][c]=0.0;
				if(par->distr_stat==1)egy->SW_max->co[r][c]=-9.0E99;
			}
			if(par->output_ET>0){
				egy->ET_mean->co[r][c]=0.0;
				if(par->distr_stat==1)egy->ET_max->co[r][c]=-9.0E99;
				if(par->distr_stat==1)egy->ET_min->co[r][c]=+9.0E99;
			}
			if(par->output_H>0){
				egy->H_mean->co[r][c]=0.0;
				if(par->distr_stat==1)egy->H_max->co[r][c]=-9.0E99;
				if(par->distr_stat==1)egy->H_min->co[r][c]=+9.0E99;
			}
			if(par->output_G>0){
				egy->SEB_mean->co[r][c]=0.0;
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
				egy->LW->co[r][c]=UV->V->co[2];
				egy->SW->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->SW_max->co[r][c]=UV->V->co[2];			
			}
			if(par->output_ET>0){
				egy->ET_mean->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->ET_max->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->ET_min->co[r][c]=UV->V->co[2];
			}
			if(par->output_H>0){
				egy->H_mean->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->H_max->co[r][c]=UV->V->co[2];
				if(par->distr_stat==1)egy->H_min->co[r][c]=UV->V->co[2];
			}
			if(par->output_G>0){
				egy->SEB_mean->co[r][c]=UV->V->co[2];
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

//vectors used in energy_balance()
egy->Dlay = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->wliq = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->wice = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Temp = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->deltaw = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
	
egy->SWlayer = new_doublevector( par->snowlayer_max+1 );
	
egy->soil_transp_layer = new_doublevector(land->root_fraction->nch);
initialize_doublevector(egy->soil_transp_layer, 0.);
	
egy->dFenergy = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Kth0=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Kth1=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Fenergy=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Newton_dir=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->T0=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->T1=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Tstar=new_doublevector(Nl); //soil temperature at which freezing begins
egy->THETA=new_doublevector(Nl);	//water content (updated in the iterations)
	
//allocate vector	of soil layer contributions to evaporation (up to z_evap)
z = 0.;
l = 0;
do{	
	l++;
	z += sl->pa->co[1][jdz][l];
}while(l<Nl && z < z_evap);
egy->soil_evap_layer_bare = new_doublevector(l);
egy->soil_evap_layer_veg = new_doublevector(l);
initialize_doublevector(egy->soil_evap_layer_bare, 0.);
initialize_doublevector(egy->soil_evap_layer_veg, 0.);

printf("Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare->nh);
printf("Soil water transpires from the first %ld layers\n",egy->soil_transp_layer->nh);

/****************************************************************************************************/
/*! Completing of the struct "water" (of the type WATER) with the initializations of the remanent
    matrices (wat->rain, wat->Pnet, wat->wcan_rain):        */

/* Initialization of wat->Pnet (liquid precipitation that reaches the sl surface in mm):*/
wat->Pnet=new_doublematrix(Nr,Nc);
initialize_doublematrix(wat->Pnet,0.0);

/* Initialization of wat->wcan_rain: (liquid precipitation intercepted by vegetation in mm):*/
wat->wcan_rain=new_doublematrix(Nr,Nc);
wat->wcan_snow=new_doublematrix(Nr,Nc);


/* Initialization of wat->PrecTot (total precipitation (rain+snow) precipitation):*/
if(par->point_sim==1 && par->micromet==1){
	wat->PrecTot=new_doublematrix(top->Z1->nrh,top->Z1->nch);	
}else{
	wat->PrecTot=new_doublematrix(Nr,Nc);
}
initialize_doublematrix(wat->PrecTot,0.0);

/* Initialization of the matrices with the output of total precipitation and interception:*/
if(par->output_P>0){
	wat->PrTOT_mean=new_doublematrix(Nr,Nc);
	wat->PrSNW_mean=new_doublematrix(Nr,Nc);
}

if(par->output_h_sup>0) wat->hsupav=new_doublematrix(Nr,Nc);

wat->dh_sup=new_doublematrix(Nr, Nc);
initialize_doublematrix(wat->dh_sup, NoV);	

for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			wat->wcan_rain->co[r][c]=0.0;
			wat->wcan_snow->co[r][c]=0.0;
			if(par->output_P>0){
				wat->PrTOT_mean->co[r][c]=0.0;
				wat->PrSNW_mean->co[r][c]=0.0;
			}
			if(par->output_h_sup>0) wat->hsupav->co[r][c]=0.0;
		}else{
			wat->wcan_rain->co[r][c]=UV->V->co[2];
			wat->wcan_snow->co[r][c]=UV->V->co[2];
			if(par->output_P>0){
				wat->PrTOT_mean->co[r][c]=UV->V->co[2];
				wat->PrSNW_mean->co[r][c]=UV->V->co[2];
			}
			if(par->output_h_sup>0) wat->hsupav->co[r][c]=UV->V->co[2];
		}
	}
}

if(par->recover==1){
	assign_recovered(files->co[rwcrn]+1, wat->wcan_rain->co, par, land->LC, IT->LU);
	assign_recovered(files->co[rwcsn]+1, wat->wcan_rain->co, par, land->LC, IT->LU);
}



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




/****************************************************************************************************/
/*! Initialization of the struct "snow" (of the type SNOW):*/

/***************************************************************************************************/
/*! Optional reading of initial real snow thickness map in the whole basin ("SNOW0"):    */
if( strcmp(files->co[fsn0]+1 , MISSING_FILE) != 0 ){
	printf("Snow initial condition from file %s\n",files->co[fsn0]+1);
	SNOW0=read_map(2, files->co[fsn0]+1, land->LC, UV);
}else{
	SNOW0=copydoublematrix_const(IT->snow0*rho_w/IT->rhosnow0, land->LC, UV->V->co[2]);
}

/*! Optional reading of snow age in the whole basin     */
if( strcmp(files->co[fsnag0]+1 , MISSING_FILE) != 0 ){
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

if(par->output_balancesn>0){
	snow->MELTED=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->MELTED,NoV);
	snow->SUBL=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->SUBL,NoV);
	snow->t_snow=new_doublematrix(Nr,Nc);
	initialize_doublematrix(snow->t_snow,NoV);
	//snow->totav_snow=new_doublematrix(Nr,Nc);	
	//initialize_doublematrix(snow->totav_snow,NoV);
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
					snow->T->co[l][r][c]=IT->Tsnow0;
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
				//snow->totav_snow->co[r][c]=0.0;					
			}
			
			snow->nondimens_age->co[r][c]=snow->dimens_age->co[r][c];
			non_dimensionalize_snowage(&(snow->nondimens_age->co[r][c]), IT->Tsnow0);		
				
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
		assign_recovered(namefile_i_we(files->co[rDzs]+1, l), snow->Dzl->co[l], par, land->LC, IT->LU);
		assign_recovered(namefile_i_we(files->co[rwls]+1, l), snow->w_liq->co[l], par, land->LC, IT->LU);
		assign_recovered(namefile_i_we(files->co[rwis]+1, l), snow->w_ice->co[l], par, land->LC, IT->LU);
		assign_recovered(namefile_i_we(files->co[rTs]+1, l), snow->T->co[l], par, land->LC, IT->LU);
	}

}




/****************************************************************************************************/
/*! Initialization of the struct "glac" (of the type GLACIER):*/

/***************************************************************************************************/ 
/*! Optional reading of glacier depth in the whole basin ("GLACIER0"):    */
if( par->point_sim!=1 && strcmp(files->co[fgl0]+1 , MISSING_FILE) != 0 ){
	if(par->glaclayer_max==0){
		printf("Warning: Glacier map present, but glacier represented with 0 layers\n");
		stop_execution();
	}
}

if(par->glaclayer_max>0){
	if( par->point_sim!=1 && strcmp(files->co[fgl0]+1 , MISSING_FILE) != 0 ){
		printf("Glacier initial condition from file %s\n",files->co[fgl0]+1);
		GLACIER0=read_map(2, files->co[fgl0]+1, land->LC, UV);
	}else{
		GLACIER0=copydoublematrix_const(IT->Dglac0, land->LC, UV->V->co[2]);
	}

}else{
	//check
	if(IT->Dglac0>0){
		f=fopen(files->co[ferr]+1,"a");
		printf("\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
		fprintf(f,"\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
		fclose(f);	
	}	
}

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
	
	if(par->recover==1){
	
		assign_recovered_long(files->co[rni]+1, glac->lnum->co, par, land->LC, IT->LU);

		for(l=1;l<=par->glaclayer_max;l++){

			assign_recovered(namefile_i_we(files->co[rDzi]+1, l),  glac->Dzl->co[l], par, land->LC, IT->LU);
			assign_recovered(namefile_i_we(files->co[rwli]+1, l),  glac->w_liq->co[l], par, land->LC, IT->LU);
			assign_recovered(namefile_i_we(files->co[rwii]+1, l),  glac->w_ice->co[l], par, land->LC, IT->LU);
			assign_recovered(namefile_i_we(files->co[rTi]+1, l),  glac->T->co[l], par, land->LC, IT->LU);
	
		}
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
if(par->point_sim!=1 && strcmp(files->co[fSCA]+1 , MISSING_FILE) != 0){
	temp=join_strings(files->co[fSCA]+1,textfile);
	f=t_fopen(temp,"w");
	fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,snowDav,SWEav,Tav,Tsav,perc.SFA,perc.SCA\n");
	t_fclose(f);
	free(temp);
}


/****************************************************************************************************/

/*Free the struct allocated in this subroutine:*/
free_doublematrix(SNOW0);
if(par->glaclayer_max>0) free_doublematrix(GLACIER0);

free_doublematrix(IT->land_classes);
free_doublematrix(IT->met);
free(IT->met_col_names);
free(IT);

//cont_nonzero_values_matrix(wat, land->LC, top->i_cont, (Nl+1)*par->total_pixel, par->point_sim);

cont_nonzero_values_matrix2(&i, land->LC, top->lrc_cont, top->i_cont, (Nl+1)*par->total_pixel, par->point_sim);
top->Li = new_longvector(i);
top->Lp = new_longvector((Nl+1)*par->total_pixel);
wat->Lx = new_doublevector(i);	
cont_nonzero_values_matrix3(top->Lp, top->Li, land->LC, top->lrc_cont, top->i_cont, (Nl+1)*par->total_pixel, par->point_sim);
	
//top->Ui = new_longvector(i);
//top->Up = new_longvector((Nl+1)*par->total_pixel);
//wat->Ux = new_doublevector(i);
//cont_nonzero_values_matrix4(top->Lp, top->Li, top->Up, top->Ui, land->LC, top->lrc_cont, top->i_cont, (Nl+1)*par->total_pixel, par->point_sim);

wat->H0 = new_doublevector((Nl+1)*par->total_pixel);
wat->H1 = new_doublevector((Nl+1)*par->total_pixel);
wat->dH = new_doublevector((Nl+1)*par->total_pixel);
wat->B = new_doublevector((Nl+1)*par->total_pixel);
wat->f = new_doublevector((Nl+1)*par->total_pixel);
wat->df = new_doublevector((Nl+1)*par->total_pixel);
	
	
	
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

	f=fopen(files->co[ferr]+1,"a");
	
	printf("Dt:%f u0:%f\n",Dt,u0);

	/* Finding out the distance of the farthest virtual channel-stretch; it has to be enough to allow an
    adeguate diffusion (to do this it is fixed a minimum fraction_spread=0.00005):*/
	smax=(long)(s0max/(u0*Dt))+1;
	
	do{ 
		smax++; 
	}while(((s0max*Dt)/pow((4.0*Pi*D*pow(smax*Dt,3.0)),0.5)*exp(-u0*pow(s0max-smax*u0*Dt,2.0)/(4.0*D*smax*u0*Dt)))>0.00001);

	fprintf(f,"s0max(farthest channel distance)= %f",s0max);
	fprintf(f,"\nsmax(number channel-pixels)=%ld\n\n",smax);
 
	printf("s0max(farthest channel distance)= %f",s0max);
	printf("\nsmax(number channel-pixels)=%ld\n\n",smax);

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
			fraction_spread->co[ch][s]=(s0->co[ch]*Dt)/pow((4.0*Pi*D*pow(s*Dt,3.0)),0.5)*exp(-u0*pow(s0->co[ch]-s*u0*Dt,2.0)/(4.0*D*s*u0*Dt));
										

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
double temp,temp_max,temp_min;

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

long r, c, i, rtot, ctot, cont;					
DOUBLEMATRIX *M, *Q, *curv;
SHORTMATRIX *P;
LONGMATRIX *ca;
T_INIT *UV2;

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
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if (land->LC->co[r][c]!=NoV){					
				if(top->DD->co[r][c]<1 || top->DD->co[r][c]>8){
					printf("cell %ld %ld has DD=%d, LC:%f Z:%f\n",r,c,top->DD->co[r][c],land->LC->co[r][c],top->Z0->co[r][c]);
					t_error("check again the drainage directions");
				}
			}
		}
	}
	
	write_map(files->co[fdd]+1, 1, par->format_out, M, UV);
	free_doublematrix(M);
	top->Z0dp=depitted(top->DD, top->Z0);
	par->channel_network = 1;
	
	top->i_DD=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	initialize_doublematrix(top->i_DD,UV->V->co[2]);
	gradients(top->Z0dp,top->DD,top->i_DD,UV);	
	
}else{
	printf("YOU HAVE CHOSEN NOT TO GIVE DRAINAGE DIRECTIONS, AND THEN CHANNEL NETWORK\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	t_error("This possibility has not been implemented yet. You should include the map of drainage directions");
	top->DD=new_shortmatrix(top->Z0->nrh, top->Z0->nch);
	initialize_shortmatrix(top->DD, 0);
	par->channel_network = 0;
	
	top->Z0dp=new_doublematrix(top->Z0->nrh, top->Z0->nch);
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			top->Z0dp->co[r][c] = top->Z0->co[r][c];
		}
	}	
}

/****************************************************************************************************/
//reading SOIL MAP
if(existing_file(files->co[fsoil]+1)>0){
	M=read_map(2, files->co[fsoil]+1, land->LC, UV);
	sl->type=copylong_doublematrix(M);

}else{//default value (99)
	M=copydoublematrix_const(1.0, land->LC, UV->V->co[2]);
	sl->type=copylong_doublematrix(M);
}
write_map(files->co[fsoil]+1, 1, par->format_out, M, UV);
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
//Channel network (in top->pixel_type)
if(existing_file(files->co[fnet]+1)>0){	
	M=read_map(2, files->co[fnet]+1, land->LC, UV);
	top->pixel_type=copyshort_doublematrix(M);
	free_doublematrix(M);
	cont = 0;
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			if(land->LC->co[r][c]!=UV->V->co[2]){
				if(par->channel_network == 1){
					if(top->pixel_type->co[r][c]!=0 && top->pixel_type->co[r][c]!=10){
						t_error("if you give drainage directions (par->channel_network=1) only values of 0 and 10 in the network map are admitted");
					}
					if(top->pixel_type->co[r][c]==10) cont++;
				}else if(par->channel_network == 0){
					if(top->pixel_type->co[r][c]!=0 && top->pixel_type->co[r][c]!=1){
						printf("r:%ld c:%ld %d\n",r,c,top->pixel_type->co[r][c]);
						t_error("if you don't give drainage directions (par->channel_network=0) only values of 0 and 1 in the network map are admitted");
					}
					if(top->pixel_type->co[r][c]==1) cont++;
				}
			}
		}
	}
	if(par->channel_network == 1){
		printf("Channel networks has %ld pixels set to channel\n",cont);
		printf("channel_network:%d\n",par->channel_network);
		printf("_________________________________________________________________________\n");
		printf("_________________________________________________________________________\n");
		printf("_________________________________________________________________________\n");
		printf("_________________________________________________________________________\n");
		
		//remove subbasins draining out of the network
		if(cont>0) remove_outdraining_subbasins(top->DD, top->pixel_type, land->LC, UV->V->co[2]);
		
	}else if(par->channel_network == 0){
		printf("No channel network and %ld pixels set to outlet\n",cont);
		printf("channel_network:%d\n",par->channel_network);
		printf("_________________________________________________________________________\n");
		printf("_________________________________________________________________________\n");
		printf("_________________________________________________________________________\n");
		printf("_________________________________________________________________________\n");
	}
	
}else if(par->channel_network == 1){
	//total contributing areas
	ca=new_longmatrix(top->Z0->nrh,top->Z0->nch);
	initialize_longmatrix(ca,0.0);
	tca(top->DD,ca);
	
	//calculate laplacian
	curv=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	nablaquadro(top->Z0,curv,UV->U,UV->V);
	
	//allocate and initialize
	top->pixel_type=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
	initialize_shortmatrix(top->pixel_type, 0);	
	
	//see geomorphologic library (overwrites pixel_type=10 for channels)
	//select_hillslopes_mod(ca,top->i_DD,curv,top->pixel_type,par->channel_thres,UV->U);
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			if(ca->co[r][c]>=par->channel_thres) top->pixel_type->co[r][c]=10; 
		}
	}
	free_longmatrix(ca);
	free_doublematrix(curv);

	// Creation of pixel-type matrix "top->pixel_type" on the basis channel network (already in top->pixel_type matrix): 
	cont=0;
	M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			M->co[r][c] = (double)top->pixel_type->co[r][c];
			if(land->LC->co[r][c]==UV->V->co[2]) M->co[r][c]=UV->V->co[2];
			if(top->pixel_type->co[r][c]==10) cont++;
		}
	}	
	write_map(files->co[fnet]+1, 1, par->format_out, M, UV);
	free_doublematrix(M);
	
	printf("Channel networks has been created with %ld pixels set to channel\n",cont);
	printf("channel_network:%d\n",par->channel_network);
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	printf("_________________________________________________________________________\n");
	
	//remove subbasins draining out of the network
	if(cont>0) remove_outdraining_subbasins(top->DD, top->pixel_type, land->LC, UV->V->co[2]);
	

}else if(par->channel_network == 0){
	top->pixel_type=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
	initialize_shortmatrix(top->pixel_type, 9);
	set_boundary_condition(top->Z0, land->LC, 1, top->pixel_type, UV->V->co[2]);
	
	M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	initialize_doublematrix(M, UV->V->co[2]);
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			if(land->LC->co[r][c]!=UV->V->co[2]){
				M->co[r][c] = (double)top->pixel_type->co[r][c];
			}
		}
	}
		
	write_map(files->co[fnet]+1, 1, par->format_out, M, UV);
	free_doublematrix(M);	
}

/****************************************************************************************************/
//Others

//distance_from_channel(M, top->DD, top->pixel_type);
//distance_from_channel2(top->dist_channel, top->pixel_type, cnet->r, cnet->c);
/*for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){
		if(land->LC->co[r][c]!=NoV){
			top->dist_channel->co[r][c]*=UV->U->co[1];
		}
	}
}*/

}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void read_parameterfile(char *name, PAR *par, INIT_TOOLS *itools){

	FILE *f;
	long i,index;
	DOUBLEVECTOR *v;
	char *temp;

	temp=join_strings(name,textfile);
	f=t_fopen(temp,"r");
	index=read_index(f,PRINT);
	
	//1st block
	itools->land_classes=read_doublematrix(f,"a",PRINT);
	if(itools->land_classes->nrh!=nlandprop+1) t_error("Error in block 1 parameter files, check properties");
	
	//2nd block
	v=read_doublearray(f,PRINT);
	if(v->nh!=11) t_error("Block 2 of parameter file requires 12 parameters - Please review them");
	par->imp=v->co[1];	/*Impedence factor for (partially) frozen soil*/
	par->TolVWb=v->co[2];
	par->RelTolVWb=RelativeErrorRichards;
	par->MaxErrWb=1.E99;
	par->MaxiterTol=(long)v->co[3];
	par->TolCG=v->co[4];
	par->UpdateK=0;
	par->DtminWb=par->Dt/pow(2.,max_time_reduction_time);
	par->nredDtWb=2;
	par->harm_or_arit_mean=0;	
	par->gamma_m=v->co[5]; /*Exponent of the law of uniform motion on the surface*/
	par->thres_hsup=v->co[6];
	par->Ks_channel=v->co[7];
	par->thres_hchannel=v->co[8];
	par->Kch_b=v->co[9];
	par->w_dx=v->co[10];
	if(par->w_dx<=0) t_error("the ratio channel width/pixel size w_dx must be > 0");
	par->depr_channel=v->co[11];
	//itools->u0=v->co[11]; /*THE MEAN VELOCITY IN CHANNELS*/
	//itools->D=v->co[12];  /*THE HYDRODYNAMICAL DISPERSION IN CHANNELS*/
	free_doublevector(v);
	
	//3rd block
	v=read_doublearray(f,PRINT);
	if(v->nh!=21) t_error("Block 3 of parameter file requires 21 parameters - Please review them");
	par->latitude=v->co[1]*Pi/180.0;
	par->longitude=v->co[2]*Pi/180.0;
	par->Vmin=v->co[3];		/*MINIMUM WIND VELOCITY VALUE (m/s)*/
	par->RHmin=v->co[4];		/*MINIMUM RELATIVE HUMIDITY VALUE (%)*/
	par->alpha_snow=v->co[5];
	par->tol_energy=v->co[6];
	par->maxiter_energy=(long)v->co[7];
	par->maxiter_canopy=(long)v->co[8];
	par->maxiter_Businger=(long)v->co[9];
	par->maxiter_Ts=(long)v->co[10];
	par->maxiter_Loc=(long)v->co[11];
	par->stabcorr_incanopy=(short)v->co[12];
	par->ifill=(int)v->co[13];
	par->iobsint=(int)v->co[14];
	par->dn=(float)v->co[15];
	par->curve_len_scale=(float)v->co[16];
	par->slopewt=(float)v->co[17];
	par->curvewt=(float)v->co[18];	
	par->topoflag=(float)v->co[19];
	par->Zboundary=v->co[20];
	par->Tboundary=v->co[21];	
	free_doublevector(v);

	/* 4th block Snow Parameters*/
	v=read_doublearray(f,PRINT);
	if(v->nh!=26) t_error("Block 4 of parameter file requires 26 parameters - Please review them");
	itools->snow0=v->co[1];/* INITIAL SNOW DEPTH (initial snow depth) [mm] valid only if there is no snow map*/
	itools->rhosnow0=v->co[2]; /*SNOW INITIAL DENSITY [kg/mc]*/
	itools->Tsnow0=v->co[3];
	itools->agesnow0=v->co[4];/* INITIAL SNOW AGE (in days), valid only if there is no snow age map */
	par->T_rain=v->co[5];   /*TEMPERATURE ABOVE WICH ALL PRECIPITAION IS RAIN [C]*/
	par->T_snow=v->co[6];   /*TEMPERATURE BELOW WICH ALL PRECIPITAION IS SNOW [C]*/
	par->aep=v->co[7];      /*AßLBEDO EXTINCTION PARAMETER [m]*/
	par->avo=v->co[8];      /*NEW SNOW VISIBLE BAND REFLECTANCE*/
	par->airo=v->co[9];     /*NEW NEAR INFRARED BAND REFLECTANCE*/
	par->Sr=v->co[10];       /*IRREDUCIBLE WATER SATURATION [-]*/
	par->epsilon_snow=v->co[11];/* SNOW LONGWAVE EMISSIVITY [-]  */
	par->z0_snow=v->co[12]*0.001;
	par->snowcorrfact=v->co[13];/* INCREASING FACTOR WHEN THE RAIN GAUGE IS SUPPOSED TO RECORD SNOW PRECIPITATION */
	par->raincorrfact=v->co[14];/* INCREASING FACTOR WHEN THE RAIN GAUGE IS SUPPOSED TO RECORD RAIN PRECIPITATION */
	par->snowlayer_max=(long)v->co[15]; /* MAXIMUM NUMBER OF SNOW LAYERS */
	par->snowlayer_inf=(long)v->co[16];
	par->snow_maxpor=v->co[17];
	par->drysnowdef_rate=v->co[18];
	par->snow_density_cutoff=v->co[19];
	par->wetsnowdef_rate=v->co[20];
	par->snow_viscosity=v->co[21];
	par->snow_fetch=v->co[22];
	par->incr_V_factor=v->co[23];
	par->snow_smin=v->co[24];	//in degrees
	par->snow_smax=v->co[25];
	par->snow_curv=v->co[26];
	free_doublevector(v);
	
	/* 5th-6th block MINIMUM and MAXIMUM SNOW LAYER THICKNESS*/
	par->Dmin=read_doublearray(f,PRINT);
	if(par->Dmin->nh!=par->snowlayer_max) t_error("Error in assigning max and min thickness to the snow layers");
	par->Dmax=read_doublearray(f,PRINT);
	if(par->Dmax->nh!=par->snowlayer_max) t_error("Error in assigning max and min thickness to the snow layers");
	par->Dmax->co[par->snowlayer_max]=1.E10;
	
	//7th block
	v=read_doublearray(f,PRINT);
	if(v->nh!=5) t_error("Block 7 of parameter file requires 5 parameters - Please review them");
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
	if(v->nh!=4) t_error("Block 9 of parameter file requires 4 parameters (lwrad, monin_obukhov, micromet and blowing snow) - Please review them");
	par->state_turb=1;
	par->state_lwrad=(short)v->co[1];
	par->monin_obukhov=(short)v->co[2];
	par->micromet=(short)v->co[3];
	par->blowing_snow=(short)v->co[4];
	if(par->blowing_snow==1 && par->micromet==0){
		par->blowing_snow=0;
		printf("\nWarning: if you do not run Micromet, you can't run SnowTrans3D\n");
	}
	
	free_doublevector(v);	
	
	//other par
	par->print=0;
	
	t_fclose(f);

	
	free(temp);
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
	//par->distr_stat=(short)v->co[4];	
	par->distr_stat=0;
	free_doublevector(v);
	
	if(par->wat_balance!=3 && par->channel_network==0){
		t_error("You have give to give the drainage directions as an input maps unless you solve Richard 3D");
	}

	//2. chkpt
	par->chkpt=read_doublematrix(f,"a",PRINT);

	//3. saving points
	par->saving_points=read_doublearray(f,PRINT);
 
	//4. output par
	v=read_doublearray(f,PRINT);
	i=0;
	if(v->nh==19) i++;
	par->output_Txy=v->co[1];
	par->output_TETAxy=v->co[2];
	par->output_TETAICExy=v->co[3];	
	par->output_PSIxy=v->co[4];
	par->output_snow=v->co[5];
	par->output_glac=v->co[6];
	par->output_h_sup=v->co[7];
	par->output_Rn=v->co[8+i];
	par->output_G=v->co[9+i];
	par->output_H=v->co[10+i];
	par->output_ET=v->co[11+i];
	par->output_Ts=v->co[12+i];
	par->output_P=v->co[13+i];
	par->output_Wr=v->co[14+i];
	par->output_balancesn=v->co[15+i];
	par->output_balancegl=v->co[16+i];
	par->output_Rswdown=v->co[17+i];
	par->output_meteo=v->co[18+i];
	free_doublevector(v);
	
	//5. special output
	v=read_doublearray(f,PRINT);
	
	if(v->co[1]>0 && v->nh>1){
		if(fmod(v->nh-1,2)!=0) t_error("special output vector must have an odd number of components");
		times->n_plot=floor(v->co[1]*3600/par->Dt);
		
		par->JD_plots=new_doublevector(v->nh-1);	//even number of components
		for(i=1;i<=v->nh-1;i++){
			par->JD_plots->co[i]=v->co[i+1];
		}
	}else{
		par->JD_plots=new_doublevector(1);
		par->JD_plots->co[1]=-1.;
	}
	free_doublevector(v);
	
	t_fclose(f);
  
	if(par->chkpt->nch<2) t_error("Error in block 16 in parameter files: E and N coordinates missing");
	if(par->chkpt->nrh>9999) t_error("Error in block 16 in parameter files: no more than 9999 points allowed");

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
void read_optionsfile_point(char *name, PAR *par, TOPO *top, LAND *land, SOIL *sl, INIT_TOOLS *IT){

	FILE *f;
	long index, i, r, c;
	DOUBLEVECTOR *v;
	STRINGBIN *s;
	DOUBLEMATRIX *M, *Q, *P;
	SHORTMATRIX *curv;
	short read_dem, read_lu, read_soil, read_sl, read_as, read_sk;
	//short read_dx, read_dy;
	
	f=t_fopen(join_strings(name,textfile),"r");
	index=read_index(f,PRINT);

	//1. base par
	v=read_doublearray(f,PRINT);
	par->wat_balance=(short)v->co[1];
	//if(par->wat_balance == 3) t_error("Water balance 3d not admitted in 1D simulation");
	par->en_balance=(short)v->co[2];
	par->state_px_coord=(short)v->co[3];	
	free_doublevector(v);
	
	//2. chkpt
	M=read_doublematrix(f,"a",PRINT);
	if(M->nch<8) t_error("Error in block 2 option file");


	//3. saving points
	par->saving_points=read_doublearray(f,PRINT);

	
	//4. horizon name
	s=read_stringarray(f,PRINT);

	t_fclose(f);
	

	//4. CALCULATE TOPOGRAPHIC PROPERTIES
	//a. read dem	
	read_dem=0;
	for(i=1;i<=M->nrh;i++){ 
		if(M->co[i][3]==-99.0 || M->co[i][6]==-99.0 || M->co[i][7]==-99.0 || M->co[i][8]==-99.0 || M->co[i][9]==-99.0 || M->co[i][10]==-99.0) read_dem=1; 
	}	
	if(par->micromet==1 || par->recover==1) read_dem=1;

	if(read_dem==1){
		if(existing_file(files->co[fdem]+1)>0){
			
			printf("Warning: Dem file %s present\n",files->co[fdem]+1);
			Q=new_doublematrix(1,1);
			top->Z1=read_map(0, files->co[fdem]+1, Q, UV); //topography
			free_doublematrix(Q);

		}else{
			printf("Warning: Dem file not present\n");
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
			if(M->co[i][3]==UV->V->co[2]) printf("Point corresponding to NoValue");
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
	if(read_lu==0){
		par->micromet=0;
		for(i=1;i<=M->nrh;i++){ 
			if(M->co[i][4]==-99) M->co[i][4]=1.0; 
		}
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
				if(M->co[i][4]==UV->V->co[2]) printf("Point corresponding to NoValue");
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
				if(M->co[i][5]==UV->V->co[2]) printf("Point corresponding to NoValue");

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
				if(M->co[i][6]==UV->V->co[2]) printf("Point corresponding to NoValue");

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
				if(M->co[i][7]==UV->V->co[2]) printf("Point corresponding to NoValue");

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
				if(M->co[i][8]==UV->V->co[2]) printf("Point corresponding to NoValue");
			}
		}
		free_doublematrix(P);
	}

	//g. dz/dx(east)
	/*read_dx=0;
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
	}*/
	

	//i.show results
	printf("\nPOINTS:\n");
	f=fopen(files->co[ferr]+1,"a");
	fprintf(f,"\nPOINTS:\n");
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=8;c++){
			printf("%f  ",M->co[r][c]);
			fprintf(f,"%f  ",M->co[r][c]);
		}
		printf("\n");
		fprintf(f,"\n\n");
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
	
	//6. SET PROPERTIES
	top->Z0=new_doublematrix(1,M->nrh);
	land->LC=new_doublematrix(1,M->nrh);	
	sl->type=new_longmatrix(1,M->nrh);	
	top->slopes=new_doublematrix(1,M->nrh);
	top->aspect=new_doublematrix(1,M->nrh);	
	top->sky=new_doublematrix(1,M->nrh);
	//top->dz_dx=new_doublematrix(1,M->nrh);
	//top->dz_dy=new_doublematrix(1,M->nrh);
	top->DD=new_shortmatrix(1,M->nrh);
	top->i_DD=new_doublematrix(1,M->nrh);
	top->area=new_doublematrix(1,M->nrh);
	top->pixel_type=new_shortmatrix(1,M->nrh);
	top->Z0dp=new_doublematrix(1,M->nrh);
	for(i=1;i<=M->nrh;i++){
		top->Z0->co[1][i]=M->co[i][3];
		land->LC->co[1][i]=M->co[i][4];
		sl->type->co[1][i]=(long)M->co[i][5];
		top->slopes->co[1][i]=M->co[i][6]*Pi/180.0;
		top->aspect->co[1][i]=M->co[i][7]*Pi/180.0;
		top->sky->co[1][i]=M->co[i][8];
		//top->dz_dx->co[1][i]=0.;
		//top->dz_dy->co[1][i]=0.;
		//these values are not used in 1d simulations
		top->DD->co[1][i]=10;	//outlet
		top->i_DD->co[1][i]=0.0;
		top->area->co[1][i]=UV->U->co[1]*UV->U->co[2];
		top->pixel_type->co[1][i]=0;
		top->Z0dp->co[1][i]=top->Z0->co[1][i];
	}
	
	//7. SET PAR
	par->output_Txy=0;
	par->output_TETAxy=0;
	par->output_TETAICExy=0;	
	par->output_PSIxy=0;
	par->output_snow=0;
	par->output_glac=0;
	par->output_h_sup=0;
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
	par->JD_plots=new_doublevector(1);
	par->JD_plots->co[1]=-1.;	
	par->blowing_snow=0;
	par->topoflag=0;

	//8. READ HORIZONS
	top->horizon_height=(double ****)malloc(top->Z0->nrh*sizeof(double***));	
	for(r=1;r<=top->Z0->nrh;r++){
		top->horizon_height[r-1]=(double ***)malloc(top->Z0->nch*sizeof(double**));	
		for(c=1;c<=top->Z0->nch;c++){
			i=c;
			top->horizon_height[r-1][c-1]=read_horizon(join_strings(WORKING_DIRECTORY,s->co[1]+1),i);
		}
	}
	free(s);
	
	free_doublematrix(M);

}
	

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
double **read_horizon(char *name, long i){
	
	FILE *f;
	DOUBLEMATRIX *M;
	long j, index;
	double **hor;
	char *temp;
	short fileyes;
	
	if(strcmp(name , MISSING_FILE) != 0){
		fileyes=1;
	}else{
		fileyes=0;
	}
	
	if(fileyes==1){
		temp=namefile_i(name,i);
		f=fopen(temp,"r");
		if(f==NULL){
			fileyes=0;
		}else{
			fclose(f);
		}
		free(temp);
	}
	
	if(fileyes==1){
		
		temp=namefile_i(name,i);
		f=t_fopen(temp,"r");
		index=read_index(f,PRINT);
		M=read_doublematrix(f,"a",PRINT);
		t_fclose(f);
		
		if(M->nch>2){
			printf("Warning: in file %s columns after the 2nd will be neglected\n",temp);
		}else if(M->nch<2){
			printf("Error: in file %s insufficient number of columns\n",temp);
			t_error("ERROR!!");
		}
		hor=alloc2(M->nrh,2);
		for(j=1;j<=M->nrh;j++){
			hor[j-1][0]=M->co[j][1];
			hor[j-1][1]=M->co[j][2];
		}
		
		free_doublematrix(M);
		free(temp);

	}else{
	
		if(strcmp(name , MISSING_FILE) != 0){
			f=fopen(files->co[ferr]+1,"a");
			fprintf(f,"\nNo horizon file of the Meteo Station #%ld is present. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
			fclose(f);
			printf("\nNo horizon file of the Meteo Station #%ld is present. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
			
			temp=namefile_i(name,i);
			f=t_fopen(temp,"w");
			fprintf(f,"/** Horizon file for met station or point #%ld */\n",i);
			fprintf(f,"\n");
			fprintf(f,"1: double matrix horizon {4,2}\n");
		
			hor=alloc2(4,2);
			for(j=1;j<=4;j++){
				hor[j-1][0]=45.0+(j-1)*90.0;
				hor[j-1][1]=0.0;
				fprintf(f,"%f %f\n",hor[j-1][0],hor[j-1][1]);
			}
			t_fclose(f);
			free(temp);
		}else{
			hor=alloc2(4,2);
			for(j=1;j<=4;j++){
				hor[j-1][0]=45.0+(j-1)*90.0;
				hor[j-1][1]=0.0;
			}
		}			
			
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
	
	st->E=new_doublevector(INPUTmeteo->nrh);
	st->N=new_doublevector(INPUTmeteo->nrh);
	st->lat=new_doublevector(INPUTmeteo->nrh);
	st->lon=new_doublevector(INPUTmeteo->nrh);
	st->Z=new_doublevector(INPUTmeteo->nrh);
	st->sky=new_doublevector(INPUTmeteo->nrh);
	st->ST=new_doublevector(INPUTmeteo->nrh);
	st->Vheight=new_doublevector(INPUTmeteo->nrh);
	st->Theight=new_doublevector(INPUTmeteo->nrh);
	st->JD0=new_doublevector(INPUTmeteo->nrh);
	st->Y0=new_longvector(INPUTmeteo->nrh);
	st->Dt=new_doublevector(INPUTmeteo->nrh);
	st->offset=new_longvector(INPUTmeteo->nrh);	
	
	for(i=1;i<=INPUTmeteo->nrh;i++){
		st->E->co[i]=INPUTmeteo->co[i][1];
		st->N->co[i]=INPUTmeteo->co[i][2];
		st->lat->co[i]=INPUTmeteo->co[i][3]*Pi/180.0;
		st->lon->co[i]=INPUTmeteo->co[i][4]*Pi/180.0;
		st->Z->co[i]=INPUTmeteo->co[i][5];
		st->sky->co[i]=INPUTmeteo->co[i][6];
		st->ST->co[i]=INPUTmeteo->co[i][7];
		st->Vheight->co[i]=INPUTmeteo->co[i][8];
		st->Theight->co[i]=INPUTmeteo->co[i][9];
		st->JD0->co[i]=INPUTmeteo->co[i][10];
		st->Y0->co[i]=(long)(INPUTmeteo->co[i][11]);
		st->Dt->co[i]=INPUTmeteo->co[i][12]*3600;		//seconds
		st->offset->co[i]=(long)(INPUTmeteo->co[i][13]);
	}
	
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		
void ReadMeteoHeader(FILE *f, STRINGBIN *ColDescr, long offset, long *ncols, long *MeteoCont){
	
	char **a;	
	long i,j;
	
	for(i=0;i<dim1l(MeteoCont);i++){
		MeteoCont[i]=-1;
	}
			
	a=readline_textarray(f, offset);
	
	*ncols=dim_vect_strings(a);
	
	for(j=1;j<=dim1l(MeteoCont);j++){
		for(i=1;i<=*ncols;i++){
			if(compare_strings(ColDescr->co[j]+1, a[i-1])==1 && MeteoCont[j-1]==-1){
				MeteoCont[j-1]=i-1;
			}else if(compare_strings(ColDescr->co[j]+1, a[i-1])==1 && MeteoCont[j-1]!=-1){
				printf("Column '%s' is present twice\n",ColDescr->co[j]+1);
				t_error("Meteo Column name DUPLICATED!");
			}
		}
	}
		
}
	
	
	
		
/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		
double **read_datameteo(FILE *f, long offset, long ncols, double ndef){

	double **a;
	long i, j;
	short end, novalueend=1;
		
	i=0;
	do{
		if(i==0){
			a=(double **)malloc(sizeof(double*));
		}else{
			a=(double **)realloc(a,(i+1)*sizeof(double*));
		}
		a[i]=(double *)malloc((ncols+1)*sizeof(double));
		a[i][ncols]=end_vector;
		readline_array(f, a[i], offset, ncols, ndef, &end);
		if(end==0)i++;
	}while(end==0);


	for(j=0;j<ncols;j++){
		if(a[i][j]!=ndef) novalueend=0;
	}
	
	if(novalueend==0){
		i++;
		a=(double **)realloc(a,(i+1)*sizeof(double*));
		a[i]=(double *)malloc(sizeof(double));
	}
	a[i][0]=end_vector;
		
	return(a);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

void read_inpts_par(PAR *par, TIMES *times, char *program, char *ext){
	
	DOUBLEVECTOR *V;
	double Dt_output;
	long index;
	FILE *f;
	char *temp;
	
	temp = join_strings(program, ext);
	f = fopen(temp, "r");
	
	index=read_index(f,PRINT);
		
	printf("\nENTERING SEVERAL CONTROL PROGRAM PAR %s %s\n",program,ext);

	V=read_doublearray(f, PRINT);
	if(V->nh != 13) t_error("GEOtop requires 13 parameters in the base parameter file");

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

	printf("ENTER THE NUMBER OF Dt AFTER WHICH THE DISCHARGE IS PRINTED: %f\n",V->co[6]);
	Dt_output=(double)V->co[6];
	if(Dt_output>0 && Dt_output*3600<par->Dt) Dt_output=par->Dt/3600.0;
	times->n_discharge=(long)(Dt_output*3600/(long)par->Dt);
	times->i_discharge=0;/*counter for the output of a pixel*/
	if(times->n_discharge>0){
		par->state_discharge=1;
	}else{
		par->state_discharge=0;
	}

	
	printf("ENTER THE NUMBER OF Dt AFTER WHICH THE OUTPUT FOR A SPECIFIED PIXEL ARE PRINTED: %f\n",V->co[7]);
	Dt_output=(double)V->co[7];
	if(Dt_output>0 && Dt_output*3600<par->Dt) Dt_output=par->Dt/3600.0;
	times->n_pixel=(long)(Dt_output*3600/(long)par->Dt);
	times->i_pixel=0;/*counter for the output of a pixel*/
	if(times->n_pixel>0){
		par->state_pixel=1;
	}else{
		par->state_pixel=0;
	}	

	printf("ENTER THE NUMBER OF Dt AFTER WHICH THE OUTPUT FOR THE BASIN ARE PRINTED: %f\n",V->co[8]);
	Dt_output=(double)V->co[8];
	if(Dt_output>0 && Dt_output*3600<par->Dt) Dt_output=par->Dt/3600.0;
	times->n_basin=(long)(Dt_output*3600/(long)par->Dt);
	times->i_basin=0;/*counter for the output of a pixel*/
	if(times->n_basin>0){
		par->state_basin=1;
	}else{
		par->state_basin=0;
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
	
	fclose(f);
	free(temp);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
	
DOUBLETENSOR *read_soil_parameters(char *name){

	FILE *f;
	long index, i, j, k, n;
	DOUBLEMATRIX *M;
	DOUBLETENSOR *Soil_parameters;
	char *temp;
	
	temp = join_strings(name,textfile);
	f = t_fopen(temp,"r");
	index = read_index(f,PRINT);	
	
	printf("%ld Different soil types given in %s\n",index,temp);
	
	M = read_doublematrix(f,"a",PRINT);
	n = M->nrh;
	
	if(M->nch != nsoilprop){
		printf("GEOtop requires %d soil parameters, you have given %ld\n",nsoilprop,M->nch);
		t_error("Please check that");
	}
	
	Soil_parameters = new_doubletensor(index, nsoilprop, n);
	initialize_doubletensor(Soil_parameters, 0.0);
	
	for(i=1;i<=index;i++){
		
		if(i>1){
			M = read_doublematrix(f,"a",PRINT);
			if(M->nrh != n) t_error("The number of soil layers must be the same for each sl type");
		}
		
		for(j=1;j<=nsoilprop;j++){
			for(k=1;k<=n;k++){
				Soil_parameters->co[i][j][k] = M->co[k][j];
				if(j==jdz && i!=1){
					if(Soil_parameters->co[i][j][k] != Soil_parameters->co[i-1][j][k]) t_error("Soil layer thicknesses must be the same for each soil type");
				}
			}
		}
		free_doublematrix(M);
		
	}
	
	t_fclose(f);
	
	free(temp);
		
	return(Soil_parameters);
		
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
	
long n_entries_matrix(TOPO *top, DOUBLEMATRIX *LC, long n){	

	int i, j=0, l, r, c;
	
	for(i=1;i<=n;i++){
		l=top->lrc_cont->co[i][1];
		r=top->lrc_cont->co[i][2];
		c=top->lrc_cont->co[i][3];
		j++;
		if(l>1) j++;
		if(l<Nl) j++;
		if(r>1){
			if(LC->co[r-1][c]!=NoV) j++;
		}
		if(r<Nr){
			if(LC->co[r+1][c]!=NoV) j++;
		}
		if(c>1){
			if(LC->co[r][c-1]!=NoV) j++;
		}
		if(c<Nc){
			if(LC->co[r][c+1]!=NoV) j++;
		}
	}
	return(j);
}

//***************************************************************************

/*void complete_Ap_Ai(TOPO *top, DOUBLEMATRIX *LC, long n){	

	long i, j, cont=0, l, r, c, L, R, C;

	for(i=1;i<=n;i++){
		l=top->lrc_cont->co[i][1];
		r=top->lrc_cont->co[i][2];
		c=top->lrc_cont->co[i][3];
		for(j=1;j<=n;j++){
			L=top->lrc_cont->co[j][1];
			R=top->lrc_cont->co[j][2];
			C=top->lrc_cont->co[j][3];
			if( (R==r && C==c && fabs(L-l)<=1) || (fabs(R-r)==1 && C==c && L==l && l>0) || (fabs(C-c)==1 && R==r && L==l && l>0) ) cont++;
		}
	}
	
	top->Ai=new_longvector(cont);
	top->Ap=new_longvector(n);
		
	cont=0;
	
	for(i=1;i<=n;i++){	//column
		l=top->lrc_cont->co[i][1];
		r=top->lrc_cont->co[i][2];
		c=top->lrc_cont->co[i][3];
		
		for(j=1;j<=n;j++){	//row
			L=top->lrc_cont->co[j][1];
			R=top->lrc_cont->co[j][2];
			C=top->lrc_cont->co[j][3];
						
			if( (R==r && C==c && fabs(L-l)<=1) || (fabs(R-r)==1 && C==c && L==l && l>0) || (fabs(C-c)==1 && R==r && L==l && l>0) ){
				cont++;
				top->Ai->co[cont] = j;
			}
			
		}

		top->Ap->co[i]=cont;
	}
		
}*/
	
//***************************************************************************
	
void find_max_constraint( DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long *R, long *C){
	
	long r, c;
	double z = -1.E99;
	
	*R=0;
	*C=0;
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(LC->co[r][c]!=NoV){
				if(pixel_type->co[r][c]>=10 && CH->co[r][c]==0){
					if(Z->co[r][c]>z){
						z=Z->co[r][c];
						*R=r;
						*C=c;
					}
				}
			}
		}
	}
}

//***************************************************************************

void next_down_channel_pixel( long r, long c, DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long *R, long *C){

	*R=0;
	*C=0;

	if(neighboring_down_channel_pixel(r, c, 1, 1, Z, LC, pixel_type, CH) == -1){
		*R=r+1;
		*C=c+1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, -1, Z, LC, pixel_type, CH) == -1){
		*R=r+1;
		*C=c-1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, -1, 1, Z, LC, pixel_type, CH) == -1){
		*R=r-1;
		*C=c+1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}																																																														
	
	if(neighboring_down_channel_pixel(r, c, -1, -1, Z, LC, pixel_type, CH) == -1){
		*R=r-1;
		*C=c-1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, 0, Z, LC, pixel_type, CH) == -1){
		*R=r+1;
		*C=c;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 0, 1, Z, LC, pixel_type, CH) == -1){
		*R=r;
		*C=c+1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, -1, 0, Z, LC, pixel_type, CH) == -1){
		*R=r-1;
		*C=c;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 0, -1, Z, LC, pixel_type, CH) == -1){
		*R=r;
		*C=c-1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, 0, Z, LC, pixel_type, CH) == 1){
		*R=r+1;
		*C=c;
	}
	
	if(neighboring_down_channel_pixel(r, c, 0, 1, Z, LC, pixel_type, CH) == 1){
		*R=r;
		*C=c+1;
	}
	
	if(neighboring_down_channel_pixel(r, c, -1, 0, Z, LC, pixel_type, CH) == 1){
		*R=r-1;
		*C=c;
	}
	
	if(neighboring_down_channel_pixel(r, c, 0, -1, Z, LC, pixel_type, CH) == 1){
		*R=r;
		*C=c-1;
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, 1, Z, LC, pixel_type, CH) == 1){
		*R=r+1;
		*C=c+1;
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, -1, Z, LC, pixel_type, CH) == 1){
		*R=r+1;
		*C=c-1;
	}
	
	if(neighboring_down_channel_pixel(r, c, -1, 1, Z, LC, pixel_type, CH) == 1){
		*R=r-1;
		*C=c+1;
	}																																																														

	if(neighboring_down_channel_pixel(r, c, -1, -1, Z, LC, pixel_type, CH) == 1){
		*R=r-1;
		*C=c-1;
	}
}

//***************************************************************************	

short neighboring_down_channel_pixel( long r, long c, long ir, long ic, DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH){
	
	short yes=0;
	long R=r+ir, C=c+ic;
	
	if(R>=1 && R<=Nr && C>=1 && C<=Nc){
		if(LC->co[R][C]!=NoV){
			if(Z->co[R][C]<=Z->co[r][c] && pixel_type->co[R][C]>=10) yes=-1;
			if(Z->co[R][C]<=Z->co[r][c] && pixel_type->co[R][C]>=10 && CH->co[R][C]==0) yes=1;
		}
	}
	
	return(yes);
}


//***************************************************************************

void i_lrc_cont(DOUBLEMATRIX *LC, long ***i, LONGMATRIX *lrc){
	
	long cont=0;
	long l, r, c;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(LC->co[r][c]!=NoV){
				for(l=0;l<=Nl;l++){
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

//***************************************************************************

void j_rc_cont(DOUBLEMATRIX *LC, long **j, LONGMATRIX *rc){
	
	long cont=0;
	long r, c;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(LC->co[r][c]!=NoV){
				cont++;
				j[r][c]=cont;
				rc->co[cont][1]=r;
				rc->co[cont][2]=c;
			}
		}
	}
}

//***************************************************************************

/*void searchDD(long ***I, LONGMATRIX *LRC, SHORTMATRIX *DD, DOUBLEMATRIX *LC, LONGMATRIX *up, LONGVECTOR *down){
	
	short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	
	long r, c, j, R, C;
	long n = RC->nrh;
	long cont;
	
	up = new_longmatrix(n, 8);
	initialize_longmatrix(up, 0);
	
	down = new_longvector(n);
	initialize_longvector(down, 0);
	
	for(j=1;j<=n;j++){
		
		r = RC->co[j][1];
		c = RC->co[j][2];
		
		cont=0;
		
		down->co[j] = J[r+r_DD[DD->co[r][c]]][c+c_DD[DD->co[r][c]]];
		
		R=r-1;
		C=c;
		if(LC->co[R][C]!=NoV){
			if(R+r_DD[DD->co[R][C]]==r && C+c_DD[DD->co[R][C]]==c){
				cont++;
				up->co[j][cont] = J[R][C];
			}
		}
		
		R=r+1;
		C=c;
		if(LC->co[R][C]!=NoV){
			if(R+r_DD[DD->co[R][C]]==r && C+c_DD[DD->co[R][C]]==c){
				cont++;
				up->co[j][cont] = J[R][C];
			}
		}
		
		R=r;
		C=c+1;
		if(LC->co[R][C]!=NoV){
			if(R+r_DD[DD->co[R][C]]==r && C+c_DD[DD->co[R][C]]==c){
				cont++;
				up->co[j][cont] = J[R][C];
			}
		}
		
		R=r;
		C=c-1;
		if(LC->co[R][C]!=NoV){
			if(R+r_DD[DD->co[R][C]]==r && C+c_DD[DD->co[R][C]]==c){
				cont++;
				up->co[j][cont] = J[R][C];
			}
		}
		
		R=r-1;
		C=c-1;
		if(LC->co[R][C]!=NoV){
			if(R+r_DD[DD->co[R][C]]==r && C+c_DD[DD->co[R][C]]==c){
				cont++;
				up->co[j][cont] = J[R][C];
			}
		}
		
		R=r+1;
		C=c+1;
		if(LC->co[R][C]!=NoV){
			if(R+r_DD[DD->co[R][C]]==r && C+c_DD[DD->co[R][C]]==c){
				cont++;
				up->co[j][cont] = J[R][C];
			}
		}
		
		R=r-1;
		C=c+1;
		if(LC->co[R][C]!=NoV){
			if(R+r_DD[DD->co[R][C]]==r && C+c_DD[DD->co[R][C]]==c){
				cont++;
				up->co[j][cont] = J[R][C];
			}
		}
		
		R=r+1;
		C=c-1;
		if(LC->co[R][C]!=NoV){
			if(R+r_DD[DD->co[R][C]]==r && C+c_DD[DD->co[R][C]]==c){
				cont++;
				up->co[j][cont] = J[R][C];
			}
		}
		
	}
}
*/
//***************************************************************************

/*void cont_nonzero_values_matrix(WATER *wat, DOUBLEMATRIX *LC, long ***i, long n, short point){
			
	long l,r,c;
	long cnt = 0;

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(LC->co[r][c]!=NoV){
				for(l=0;l<=Nl;l++){
					cnt++;//the cell itself
					if(l<Nl) cnt++; //the cell below
					if(l>0) cnt++; //the cell above
					if(l>0){
						if(LC->co[r-1][c]!=NoV && point!=1) cnt++;
						if(LC->co[r+1][c]!=NoV && point!=1) cnt++;
						if(LC->co[r][c-1]!=NoV && point!=1) cnt++;
						if(LC->co[r][c+1]!=NoV && point!=1) cnt++;
					}
				}
			}
		}
	}
	
	wat->Jtriplet = new_UMFPACK_REAL_TRIPLET(cnt);
	wat->Jmatrix = new_UMFPACK_REAL_MATRIX(n, n, cnt);
	
	cnt = 0;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(LC->co[r][c]!=NoV){
				for(l=0;l<=Nl;l++){
					
					cnt++;//the cell itself
					wat->Jtriplet->Ti[cnt]=(UF_long)(i[l][r][c]-1);
					wat->Jtriplet->Tj[cnt]=(UF_long)(i[l][r][c]-1);

					if(l<Nl){
						cnt++; //the cell below
						wat->Jtriplet->Ti[cnt]=(UF_long)(i[l][r][c]-1);
						wat->Jtriplet->Tj[cnt]=(UF_long)(i[l+1][r][c]-1);
					}
					
					if(l>0){
						cnt++; //the cell above
						wat->Jtriplet->Ti[cnt]=(UF_long)(i[l][r][c]-1);
						wat->Jtriplet->Tj[cnt]=(UF_long)(i[l-1][r][c]-1);
					}
					
					if(l>0){
						if(LC->co[r-1][c]!=NoV && point!=1){
							cnt++;
							wat->Jtriplet->Ti[cnt]=(UF_long)(i[l][r][c]-1);
							wat->Jtriplet->Tj[cnt]=(UF_long)(i[l][r-1][c]-1);
						}
						
						if(LC->co[r+1][c]!=NoV && point!=1){
							cnt++;
							wat->Jtriplet->Ti[cnt]=(UF_long)(i[l][r][c]-1);
							wat->Jtriplet->Tj[cnt]=(UF_long)(i[l][r+1][c]-1);
						}
						
						if(LC->co[r][c-1]!=NoV && point!=1){
							cnt++;
							wat->Jtriplet->Ti[cnt]=(UF_long)(i[l][r][c]-1);
							wat->Jtriplet->Tj[cnt]=(UF_long)(i[l][r][c-1]-1);
						}
						
						if(LC->co[r][c+1]!=NoV && point!=1){
							cnt++;
							wat->Jtriplet->Ti[cnt]=(UF_long)(i[l][r][c]-1);
							wat->Jtriplet->Tj[cnt]=(UF_long)(i[l][r][c+1]-1);
						}
					}
				}
			}
		}
	}
}*/

//***************************************************************************


void cont_nonzero_values_matrix2(long *tot, DOUBLEMATRIX *LC, LONGMATRIX *lrc, long ***i, long n, short point){

	long j, l, r, c, cnt=0;
	
	for(j=1;j<=n;j++){
		
		l=lrc->co[j][1];
		r=lrc->co[j][2];
		c=lrc->co[j][3];
		
		//the cell itself
		//cnt ++;		
		
		//the cell below
		if(l<Nl) cnt ++;
		
		if(l>0 && point!=1 && LC->co[r-1][c]!=NoV){
			if(i[l][r-1][c]>j) cnt ++;
		}
		
		if(l>0 && point!=1 && LC->co[r+1][c]!=NoV){
			if(i[l][r+1][c]>j) cnt ++;
		}
		
		if(l>0 && point!=1 && LC->co[r][c-1]!=NoV){
			if(i[l][r][c-1]>j) cnt ++;
		}
		
		if(l>0 && point!=1 && LC->co[r][c+1]!=NoV){
			if(i[l][r][c+1]>j) cnt ++;
		}
	}
	
	*tot = cnt;
			
}

//***************************************************************************


void cont_nonzero_values_matrix3(LONGVECTOR *Lp, LONGVECTOR *Li, DOUBLEMATRIX *LC, LONGMATRIX *lrc, long ***i, long n, short point){
	
	//Ai = line index
	//Ap = number of values for each row
	long j,l,r,c;
	long cnt = 0;
		
	for(j=1;j<=n;j++){
		
		l=lrc->co[j][1];
		r=lrc->co[j][2];
		c=lrc->co[j][3];
		
		//the cell itself
		//cnt++;
		//Li->co[cnt] = j;
		
		//the cell below
		if(l<Nl){
			cnt++;
			Li->co[cnt] = j+1;
		}
		
		if(l>0 && point!=1 && LC->co[r-1][c]!=NoV){
			if(i[l][r-1][c]>j){
				cnt++;
				Li->co[cnt] = i[l][r-1][c];
			}
		}
		
		if(l>0 && point!=1 && LC->co[r+1][c]!=NoV){
			if(i[l][r+1][c]>j){
				cnt++;
				Li->co[cnt] = i[l][r+1][c];
			}
		}
		
		if(l>0 && point!=1 && LC->co[r][c-1]!=NoV){
			if(i[l][r][c-1]>j){
				cnt++;
				Li->co[cnt] = i[l][r][c-1];
			}
		}
		
		if(l>0 && point!=1 && LC->co[r][c+1]!=NoV){
			if(i[l][r][c+1]>j){
				cnt++;
				Li->co[cnt] = i[l][r][c+1];
			}
		}	
		
		Lp->co[j] = cnt;
	}

}

//***************************************************************************

void cont_nonzero_values_matrix4(LONGVECTOR *Lp, LONGVECTOR *Li, LONGVECTOR *Up, LONGVECTOR *Ui, DOUBLEMATRIX *LC, 
								 LONGMATRIX *lrc, long ***i, long n, short point){
	
	//Li = line index
	//Lp = number of values for each row
	//Ui = line index transposed
	//Up = number of values for each row trasnposed
	//Axt such that Ax[Axt[i]] is the transposed
	
	long j,l,r,c;
	long cnt = 0, cntt = 0;
		
	for(j=1;j<=n;j++){
		
		l=lrc->co[j][1];
		r=lrc->co[j][2];
		c=lrc->co[j][3];
		
		//the cell itself
		cnt++;
		Li->co[cnt] = j;
		
		//the cell below
		if(l<Nl){
			cnt++;
			Li->co[cnt] = j+1;
		}
		
		if(l>0 && point!=1 && LC->co[r-1][c]!=NoV){
			if(i[l][r-1][c]>j){
				cnt++;
				Li->co[cnt] = i[l][r-1][c];
			}else if(i[l][r-1][c]<j){
				cntt++;
				Ui->co[cntt] = i[l][r-1][c];
			}
		}
		
		if(l>0 && point!=1 && LC->co[r+1][c]!=NoV){
			if(i[l][r+1][c]>j){
				cnt++;
				Li->co[cnt] = i[l][r+1][c];
			}else if(i[l][r+1][c]<j){
				cntt++;
				Ui->co[cntt] = i[l][r+1][c];
			}		
		}
		
		if(l>0 && point!=1 && LC->co[r][c-1]!=NoV){
			if(i[l][r][c-1]>j){
				cnt++;
				Li->co[cnt] = i[l][r][c-1];
			}else if(i[l][r][c-1]<j){
				cntt++;
				Ui->co[cntt] = i[l][r][c-1];
			}		
		}
		
		if(l>0 && point!=1 && LC->co[r][c+1]!=NoV){
			if(i[l][r][c+1]>j){
				cnt++;
				Li->co[cnt] = i[l][r][c+1];
			}else if(i[l][r][c+1]<j){
				cntt++;
				Ui->co[cntt] = i[l][r][c+1];
			}		
		}
		
		//the cell above
		if(l>0){
			cntt++;
			Ui->co[cntt] = j-1;
		}
		
		//the cell itself
		cntt++;
		Ui->co[cntt] = j;
		
		Lp->co[j] = cnt;
		Up->co[j] = cntt;
	}
	
}

//***************************************************************************

void remove_outdraining_subbasins(SHORTMATRIX *DD, SHORTMATRIX *net, DOUBLEMATRIX *LC, double novalue){
	
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	long r, c, R, C, i=0, max_i=50, j=0;
	short ok;
	
	do{
		
		i++;
		ok=1;
		
		for(r=1;r<=LC->nrh;r++){
			for(c=1;c<=LC->nch;c++){
				if(LC->co[r][c]!=novalue){
					
					R=r+r_DD[DD->co[r][c]];
					C=c+c_DD[DD->co[r][c]];
									
					if(LC->co[R][C]==novalue && net->co[r][c]==0){
																		
						LC->co[r][c]=novalue;
						ok=0;
						j++;
						
					}
				}
			}
		}
				
	}while(ok==0 && i<max_i);
	
	if(i==max_i) t_error("There are part of the basin not draining on the channel network that it is not possible to remove. Please review the drainage direction map.");

	printf("Removed %ld pixels draining out of the channel network\n",j);
}

//***************************************************************************


double peat_thickness(double dist_from_channel){
	
	double D;
	
	if(dist_from_channel<45.23){
		D = 10.*(47.383 - 0.928*dist_from_channel + 0.010*pow(dist_from_channel,2.));
	}else{
		D = 10.*26.406;
	}
	
	return(D);
}

//***************************************************************************
