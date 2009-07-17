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


    /*--------  1.  Include File, Prototype of the subroutine "time_loop", global variables  -------*/
#include <sys/stat.h>
#include "struct.geotop.09375.h"
#include "liston.h"
#include "input.09375.h"
#include "output.09375.h"
#include "times.h"
#include "constant.h"
#include "energy.balance.h"
#include "water.balance.h"
#include "meteo.09375.h"

void time_loop(TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet, PAR *par, ENERGY *egy,
			SNOW *snow, GLACIER *glac, TIMES *times, LISTON *liston);

void checkErrorSize(char *errfilepath);


/*----------   1. Global variables  ------------*/
T_INIT *UV;
STRINGBIN *files;
long Nl; // total number of soil layers (constant in the whole basin)
long Nr; // total number of rows (of the map)
long Nc;// total number of columns (of the map)
double NoV;

clock_t init, end, start, stop, start_loop, stop_loop;


/*----------   2.  Begin of main and declaration of its variables (several structs)   ----------*/
int main(int argc,char *argv[]){

 /*structs' declarations:*/
 TOPO *top; /* topographical characteristics */
 SOIL *sl; /* soil characteristics */
 LAND *land; /* land characteristics */
 WATER *wat;  /* water infiltrating */
 CHANNEL *cnet; /* channel routing characteristics */
 PAR *par; /* various parameters */
 ENERGY *egy; /* energy radiation characteristics */
 SNOW *snow; /* snow characteristics */
 GLACIER *glac; /* glacier characteristics */
 METEO *met; /* meteo data characteristics */
 TIMES *times; /* time variables */
 LISTON *liston; /* structure for Micromet */

 init=clock();

 /*dinamic allocations:*/
 UV=(T_INIT *)malloc(sizeof(T_INIT));
 if(!UV) t_error("UV was not allocated");

 top=(TOPO *)malloc(sizeof(TOPO));
 if(!top) t_error("top was not allocated");

 sl=(SOIL *)malloc(sizeof(SOIL));
 if(!sl) t_error("sl was not allocated");

 land=(LAND *)malloc(sizeof(LAND));
 if(!land) t_error("land was not allocated");

 wat=(WATER *)malloc(sizeof(WATER));
 if(!wat) t_error("water was not allocated");

 cnet=(CHANNEL *)malloc(sizeof(CHANNEL));
 if(!cnet) t_error("channel was not allocated");

 par=(PAR *)malloc(sizeof(PAR));
 if(!par) t_error("par was not allocated");

 egy=(ENERGY *)malloc(sizeof(ENERGY));
 if(!egy) t_error("egy was not allocated");

 snow=(SNOW *)malloc(sizeof(SNOW));
 if(!snow) t_error("snow was not allocated");

 glac=(GLACIER *)malloc(sizeof(GLACIER));
 if(!glac) t_error("glac was not allocated");

 met=(METEO *)malloc(sizeof(METEO));
 if(!met) t_error("met was not allocated");

 times=(TIMES *)malloc(sizeof(TIMES));
 if(!times) t_error("times was not allocated");

 liston=(LISTON *)malloc(sizeof(LISTON));
 if(!liston) t_error("liston was not allocated");


/*------------------    3.  Acquisition of input data and initialisation    --------------------*/
get_all_input(argc, argv, top, sl, land, met, wat, cnet, par, egy, snow, glac, times, liston);

write_init_condit(met->st->Z->nh, times, wat, par, top, land, sl, egy, snow, glac);


/*-----------------   4. Time-loop for the balances of water-mass and egy   -----------------*/
time_loop(top,sl,land,met,wat,cnet,par,egy,snow,glac,times,liston);

/*--------------------   5.Completion of the output files and deallocaions  --------------------*/
dealloc_all(top,sl,land,wat,cnet,par,egy,snow,glac,met,liston);

end=clock();
/* Print on the screen the time of water-lateral-distribution subroutine: */
printf("%10.2f second time of simulation \n",((double)end-(double)init)/CLOCKS_PER_SEC);
printf("End of simulation!\n");
//stop_execution();
return 0;
}


/*----------------   6. The most important subroutine of the main: "time_loop"   ---------------*/
void time_loop( TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet, PAR *par, ENERGY *egy,
			SNOW *snow, GLACIER *glac, TIMES *times, LISTON *liston)
{


 /*begin of time-cycle:*/
 do{
    /*update of time vector:*/
	/*--------------------------------------------------------------------------------------------*/
	/*! Call updates_time(), which uptates time counters every time step*/

    updates_times(times,par);
    //printf("time=%f",times->time); stop_execution();
	meteo_distr(met, liston, egy, wat, land, top, snow, times->time, par);

	//printf("\n ENERGY START %10.2f\n",times->time);
	//stop_execution();

    if(par->en_balance==1) energy_balance(times, par, land, top, sl, met, wat, egy, snow, glac, liston);

	//printf("\n ENERGY END %10.2f\n",times->time);
	//stop_execution();

	//printf("\n MASS START %10.2f\n",times->time);
	//stop_execution();

	if(par->wat_balance==1) water_balance(top, sl, land, wat, cnet, par, times->time);

	//printf("\n MASS END %10.2f\n",times->time);
	//stop_execution();

	//printf("\n WRITE OUTPUT b %10.2f\n",times->time);
	//stop_execution();

	/*! Call write_output(), which writes output data */
    write_output(times,wat,cnet,par,top,land,sl,egy,snow,glac,met);

	//printf("\n WRITE OUTPUT e %10.2f\n",times->time);
	//stop_execution();

	//Increase TIME!
	times->time+=(long)par->Dt;

 }while(times->time<(times->TH*3600.0));/*end of time-cycle*/

 free(times);

}

/*--------------------------------------------------------------------------------------------*/

void checkErrorSize(char *errfilepath)
{
 struct stat testfile;
 double errordimension;
 long errsize;

	if(!stat(errfilepath, &testfile))
	{
		errsize = testfile.st_size;
		errordimension = errsize/1024000.0;
/* 		printf("Error file '%s' has grown to dimension = %f Mbytes\n", errfilepath, errordimension); */
		if(errordimension > 100.0)
		{
			printf("\n*****************************************************************************************************");
			printf("\n*****************************************************************************************************");
			printf("\n**                                                                         ");
			printf("\n**                                                                         ");
			printf("\n**  Error file '%s' has grown TOO MUCH!!!", errfilepath);
			printf("\n**  Check the parameter choice!");
			printf("\n**                                                                         ");
			printf("\n**                                                                         ");
			printf("\n*****************************************************************************************************");
			printf("\n*****************************************************************************************************\n");
			exit(1);
		}
	}

}


