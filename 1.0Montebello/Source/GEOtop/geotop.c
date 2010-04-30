
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
    
    
    /*--------  1.  Include File, Prototype of the subroutine "time_loop", global variables  -------*/

#include <sys/stat.h>
#include "struct.geotop.h"
#include "input.h"
#include "output.h"
#include "times.h"
#include "constant.h"
#include "keywords_file.h"
#include "energy.balance.h"
#include "meteo.h"
#include "water.balance.h"

void time_loop(ALLDATA *all);
               

/*----------   1. Global variables  ------------*/
T_INIT *UV;
STRINGBIN *files;
long Nl,Nr,Nc;
double NoV;
DOUBLEMATRIX *outdata_point;
DOUBLEVECTOR *outdata_basin;
char *MISSING_FILE;

//double Tgmean;

/*----------   2.  Begin of main and declaration of its variables (several structs)   ----------*/
int main(int argc,char *argv[]){

	
	/*structs' declarations:*/
	ALLDATA *adt;
   
	/*dinamic allocations:*/
	UV=(T_INIT *)malloc(sizeof(T_INIT));
	if(!UV) t_error("UV was not allocated");
 
	adt=(ALLDATA *)malloc(sizeof(ALLDATA));
	if(!adt) t_error("adt was not allocated");

//---------------------------------------------------
//---------------------------------------------------

	adt->T=(TOPO *)malloc(sizeof(TOPO));
	if(!(adt->T)) t_error("top was not allocated");
 
	adt->S=(SOIL *)malloc(sizeof(SOIL));
	if(!(adt->S)) t_error("sl was not allocated");
 
	adt->L=(LAND *)malloc(sizeof(LAND));
	if(!(adt->L)) t_error("land was not allocated");
  
	adt->W=(WATER *)malloc(sizeof(WATER));
	if(!(adt->W)) t_error("water was not allocated");

	adt->P=(PAR *)malloc(sizeof(PAR));
	if(!(adt->P)) t_error("par was not allocated");

	adt->C=(CHANNEL *)malloc(sizeof(CHANNEL));
	if(!(adt->C)) t_error("channel was not allocated"); 
 
	adt->E=(ENERGY *)malloc(sizeof(ENERGY));
	if(!(adt->E)) t_error("egy was not allocated");

	adt->N=(SNOW *)malloc(sizeof(SNOW));	
	if(!(adt->N)) t_error("snow was not allocated");	

	adt->G=(GLACIER *)malloc(sizeof(GLACIER));	
	if(!(adt->G)) t_error("glac was not allocated"); 

	adt->M=(METEO *)malloc(sizeof(METEO));	
	if(!(adt->M)) t_error("met was not allocated"); 
 
	adt->I=(TIMES *)malloc(sizeof(TIMES));
	if(!(adt->I)) t_error("times was not allocated");	
	

/*------------------    3.  Acquisition of input data and initialisation    --------------------*/
get_all_input(argc, argv, adt->T, adt->S, adt->L, adt->M, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->I);
		
/*-----------------   4. Time-loop for the balances of water-mass and egy   -----------------*/
time_loop(adt);

/*--------------------   5.Completion of the output files and deallocaions  --------------------*/
dealloc_all(adt->T, adt->S, adt->L, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->M);
free(adt);
            
printf("End of simulation!\n");

return 0;
}


/*----------------   6. The most important subroutine of the main: "time_loop"   ---------------*/
void time_loop(ALLDATA *all)
{ 

 	
 do{
	
    updates_times(all->I, all->P);
		
	meteo_distr(all->M, all->E, all->W, all->T, all->N, all->I->time, all->P);
	 	
    if(all->P->en_balance==1) energy_balance(all->I, all->P, all->L, all->T, all->S, all->M, all->W, all->E, all->N, all->G);
    			
	if(all->P->wat_balance==1) water_balance(all);

	write_output(all->I, all->W, all->C, all->P, all->T, all->L, all->S, all->E, all->N, all->G, all->M);
				
	all->I->time+=all->P->Dt;//Increase TIME
	
 }while(all->I->time<(all->I->TH*3600.0));//end of time-cycle

 free(all->I);
 
}

/*--------------------------------------------------------------------------------------------*/

