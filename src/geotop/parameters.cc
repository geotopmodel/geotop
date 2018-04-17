/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */


#include "struct.geotop.h"
#include "input.h"
#include "parameters.h"
#include "constants.h"
#include "tabs.h"
#include "extensions.h"
#include "rw_maps.h"
#include "times.h"
#include "pedo.funct.h"

extern long number_novalue, number_absent;
extern char *string_novalue;

extern char *WORKING_DIRECTORY;

extern char **files;

extern long *opnt, nopnt, *obsn, nobsn, *osnw, nosnw, *oglc, noglc, *osl,
       nosl;
extern short *ipnt, *ibsn;
extern char **hpnt, * *hbsn, * *hsnw, * *hglc, * *hsl;
extern char *keywords_num[num_par_number], *keywords_char[num_par_char];

extern char *SuccessfulRunFile, *FailedRunFile;

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

short read_inpts_par(PAR *par, LAND *land, TIMES *times, SOIL *sl, METEO *met,
                     INIT_TOOLS *itools, char *filename, FILE *flog)
{

  //variables
  FILE *f;

  short res;
  long i, j, inum, istr, n;

  long keylength, stringlength, numberlength;
  long *key, *string;
  double *number;

  char **keywords_num_read, **keywords_str_read;
  long **string_read;
  long *string_length_read;
  double **number_read;
  long *number_comp_read;

  long *num_param_components;
  double **num_param;

  char **string_param;

  short endoffile, ok;

  char *temp;
  char *path_rec_files;

  char **keywords_num_lower_case, **keywords_char_lower_case;

  long beg=0, end=0;

  //convert keyword listed on top of the file in lower case
  n = (long)num_par_number;
  keywords_num_lower_case = (char **)malloc(n*sizeof(char *));
  for (i=0; i<n; i++)
    {
      //printf("%ld,%s\n",i,keywords_num[i]);
      keywords_num_lower_case[i] = assign_string(keywords_num[i]);
      convert_string_in_lower_case(keywords_num_lower_case[i]);
    }

  n = (long)num_par_char;
  keywords_char_lower_case = (char **)malloc(n*sizeof(char *));
  for (i=0; i<n; i++)
    {
      //printf("%ld,%s\n",i,keywords_char[i]);
      keywords_char_lower_case[i] = assign_string(keywords_char[i]);
      convert_string_in_lower_case(keywords_char_lower_case[i]);
    }

  //Allocation
  n = (long)max_charstring;
  key = (long *)malloc(n*sizeof(long));
  string = (long *)malloc(n*sizeof(long));

  n = (long)max_numvect;
  number = (double *)malloc(n*sizeof(double));

  //read how many lines there are
  inum=0;
  istr=0;
  f = t_fopen(filename, "r");
  do
    {
      res=readline_par(f, 33, 61, 44, max_charstring, max_numvect, key, &keylength,
                       string, &stringlength, number, &numberlength, &endoffile);
      if (res==1) inum++;
      if (res==2) istr++;
    }
  while (endoffile==0);
  t_fclose(f);

  string_length_read=(long *)malloc(istr*sizeof(long));
  number_comp_read=(long *)malloc(inum*sizeof(long));

  keywords_str_read=(char **)malloc(istr*sizeof(char *));
  keywords_num_read=(char **)malloc(inum*sizeof(char *));

  string_read=(long **)malloc(istr*sizeof(long *));
  number_read=(double **)malloc(inum*sizeof(double *));

  //read single lines
  f = fopen(filename, "r");
  inum=0;
  istr=0;
  do
    {
      res=readline_par(f, 33, 61, 44, max_charstring, max_numvect, key, &keylength,
                       string, &stringlength, number, &numberlength, &endoffile);
      if (res==1)
        {
          inum++;
          keywords_num_read[inum-1]=find_string(key, keylength);
          convert_string_in_lower_case(keywords_num_read[inum-1]);
          number_read[inum-1]=find_number_vector(number, numberlength);
          number_comp_read[inum-1]=numberlength;
        }
      else if (res==2)
        {
          istr++;
          keywords_str_read[istr-1]=find_string(key, keylength);
          convert_string_in_lower_case(keywords_str_read[istr-1]);
          string_read[istr-1]=find_string_int(string, stringlength);
          string_length_read[istr-1]=stringlength;
        }
    }
  while (endoffile==0);
  fclose(f);
  free(key);
  free(string);
  free(number);

  //compare keywords number
  n = (long)num_par_number;
  num_param_components = (long *)malloc(n*sizeof(long));
  num_param = (double **)malloc(n*sizeof(double *));
  for (i=0; i<num_par_number; i++)
    {
      ok=0;
      for (j=0; j<inum; j++)
        {
          if ( strcmp (keywords_num_lower_case[i], keywords_num_read[j]) == 0)
            {
              ok=1;
              num_param_components[i] = number_comp_read[j];
              num_param[i] = find_number_vector(number_read[j], number_comp_read[j]);
            }
        }
      if (ok==0) //parameter not read
        {
          num_param_components[i] = 1;
          num_param[i] = (double *)malloc(sizeof(double));
          num_param[i][0] = (double)number_novalue;
        }
    }

  //deallocate read arrays
  for (j=0; j<inum; j++)
    {
      free(number_read[j]);
      free(keywords_num_read[j]);
    }
  free(number_read);
  free(keywords_num_read);
  free(number_comp_read);

  //assign parameter
  assign_numeric_parameters(par, land, times, sl, met, itools, num_param,
                            num_param_components, keywords_num, flog);

  //deallocate keyword arrays
  for (i=0; i<num_par_number; i++)
    {
      free(num_param[i]);
    }
  free(num_param);
  free(num_param_components);

  //compare keywords string
  n = (long)num_par_char;
  string_param = (char **)malloc(n*sizeof(char *));
  for (i=0; i<num_par_char; i++)
    {
      ok=0;
      for (j=0; j<istr; j++)
        {
          if ( strcmp (keywords_char_lower_case[i], keywords_str_read[j] ) == 0)
            {
              ok=1;
              string_param[i] = find_string(string_read[j], string_length_read[j]);
            }
        }
      if (ok==0)
        {
          n = strlen(string_novalue)+1;
          string_param[i] = (char *)malloc(n*sizeof(char));
          string_param[i] = strcpy(string_param[i], string_novalue);
        }
    }

  //deallocate read arrays
  for (j=0; j<istr; j++)
    {
      free(string_read[j]);
      free(keywords_str_read[j]);
    }
  free(string_read);
  free(keywords_str_read);
  free(string_length_read);

  //assign parameter
  end += nmet;
  itools->met_col_names = assign_string_parameter(flog, beg, end, string_param,
                                                  keywords_char);
  beg = end;
  end += nsoilprop;
  itools->soil_col_names = assign_string_parameter(flog, beg, end, string_param,
                                                   keywords_char);
  beg = end;
  end += 2;
  itools->horizon_col_names = assign_string_parameter(flog, beg, end,
                                                      string_param, keywords_char);
  beg = end;
  end += nfiles;
  files = assign_string_parameter(flog, beg, end, string_param, keywords_char);
  beg = end;
  end += otot;
  hpnt = assign_string_parameter(flog, beg, end, string_param, keywords_char);
  beg = end;
  end += ootot;
  hbsn = assign_string_parameter(flog, beg, end, string_param, keywords_char);
  beg = end;
  end += 10;
  hsnw = assign_string_parameter(flog, beg, end, string_param, keywords_char);
  beg = end;
  end += 10;
  hglc = assign_string_parameter(flog, beg, end, string_param, keywords_char);
  beg = end;
  end += 6;
  hsl = assign_string_parameter(flog, beg, end, string_param, keywords_char);
  beg = end;
  end += ptTOT;
  itools->point_col_names = assign_string_parameter(flog, beg, end,
                                                    string_param, keywords_char);
  beg = end;
  end += nlstot;
  itools->lapserates_col_names = assign_string_parameter(flog, beg, end,
                                                         string_param, keywords_char);
  beg = end;
  end += 8;
  itools->meteostations_col_names = assign_string_parameter(flog, beg, end,
                                                            string_param, keywords_char);

  beg = end;
  end += 1;
  temp = assignation_string(flog, beg, keywords_char, string_param);
  if (strcmp(string_novalue,temp) == 0)
    {
      free(temp);
      temp = assign_string("_SUCCESSFUL_RUN");
    }
  SuccessfulRunFile = join_strings(WORKING_DIRECTORY, temp);
  free(temp);

  beg = end;
  end += 1;
  temp = assignation_string(flog, beg, keywords_char, string_param);
  if (strcmp(string_novalue,temp) == 0)
    {
      free(temp);
      temp = assign_string("_FAILED_RUN");
    }
  FailedRunFile = join_strings(WORKING_DIRECTORY, temp);
  free(temp);

  beg = end;
  end += 1;
  path_rec_files = assignation_string(flog, beg, keywords_char,
                                      string_param);//path of recovery files

  //deallocate keyword arrays
  for (i=0; i<num_par_char; i++)
    {
      free(string_param[i]);
    }
  free(string_param);

  n = (long)num_par_number;
  for (i=0; i<n; i++)
    {
      free(keywords_num_lower_case[i]);
    }
  free(keywords_num_lower_case);

  n = (long)num_par_char;
  for (i=0; i<n; i++)
    {
      free(keywords_char_lower_case[i]);
    }
  free(keywords_char_lower_case);

  //replace none value with some default values

  //Horizon
  for (j=0; j<2; j++)
    {
      if (strcmp(itools->horizon_col_names[j], string_novalue) == 0)
        {
          free(itools->horizon_col_names[j]);
          if (j==0)
            {
              itools->horizon_col_names[j] = assign_string("AngleFromNorthClockwise");
            }
          else if (j==1)
            {
              itools->horizon_col_names[j] = assign_string("HorizonHeight");
            }
        }
    }

  //HeaderPointFile
  for (j=0; j<otot; j++)
    {
      if (strcmp(hpnt[j], string_novalue) == 0)
        {
          free(hpnt[j]);
          if (j==odate12)
            {
              hpnt[j] = assign_string("Date12[DDMMYYYYhhmm]");
            }
          else if (j==oJDfrom0)
            {
              hpnt[j] = assign_string("JulianDayFromYear0[days]");
            }
          else if (j==odaysfromstart)
            {
              hpnt[j] = assign_string("TimeFromStart[days]");
            }
          else if (j==operiod)
            {
              hpnt[j] = assign_string("Simulation_Period");
            }
          else if (j==orun)
            {
              hpnt[j] = assign_string("Run");
            }
          else if (j==opoint)
            {
              hpnt[j] = assign_string("IDpoint");
            }
          else if (j==osnowover)
            {
              hpnt[j] = assign_string("Psnow_over_canopy[mm]");
            }
          else if (j==orainover)
            {
              hpnt[j] = assign_string("Prain_over_canopy[mm]");
            }
          else if (j==oprecsnow)
            {
              hpnt[j] = assign_string("Psnow_under_canopy[mm]");
            }
          else if (j==oprecrain)
            {
              hpnt[j] = assign_string("Prain_under_canopy[mm]");
            }
          else if (j==orainonsnow)
            {
              hpnt[j] = assign_string("Prain_rain_on_snow[mm]");
            }
          else if (j==oV)
            {
              hpnt[j] = assign_string("Wind_speed[m/s]");
            }
          else if (j==oVdir)
            {
              hpnt[j] = assign_string("Wind_direction[deg]");
            }
          else if (j==oRH)
            {
              hpnt[j] = assign_string("Relative_Humidity[-]");
            }
          else if (j==oPa)
            {
              hpnt[j] = assign_string("Pressure[mbar]");
            }
          else if (j==oTa)
            {
              hpnt[j] = assign_string("Tair[C]");
            }
          else if (j==oTdew)
            {
              hpnt[j] = assign_string("Tdew[C]");
            }
          else if (j==oTg)
            {
              hpnt[j] = assign_string("Tsurface[C]");
            }
          else if (j==oTv)
            {
              hpnt[j] = assign_string("Tvegetation[C]");
            }
          else if (j==oTs)
            {
              hpnt[j] = assign_string("Tcanopyair[C]");
            }
          else if (j==oEB)
            {
              hpnt[j] = assign_string("Surface_Energy_balance[W/m2]");
            }
          else if (j==oG)
            {
              hpnt[j] = assign_string("Soil_heat_flux[W/m2]");
            }
          else if (j==oSWin)
            {
              hpnt[j] = assign_string("SWin[W/m2]");
            }
          else if (j==oSWb)
            {
              hpnt[j] = assign_string("SWbeam[W/m2]");
            }
          else if (j==oSWd)
            {
              hpnt[j] = assign_string("SWdiff[W/m2]");
            }
          else if (j==oLWin)
            {
              hpnt[j] = assign_string("LWin[W/m2]");
            }
          else if (j==ominLWin)
            {
              hpnt[j] = assign_string("LWin_min[W/m2]");
            }
          else if (j==omaxLWin)
            {
              hpnt[j] = assign_string("LWin_max[W/m2]");
            }
          else if (j==oSW)
            {
              hpnt[j] = assign_string("SWnet[W/m2]");
            }
          else if (j==oLW)
            {
              hpnt[j] = assign_string("LWnet[W/m2]");
            }
          else if (j==oH)
            {
              hpnt[j] = assign_string("H[W/m2]");
            }
          else if (j==oLE)
            {
              hpnt[j] = assign_string("LE[W/m2]");
            }
          else if (j==ofc)
            {
              hpnt[j] = assign_string("Canopy_fraction[-]");
            }
          else if (j==oLSAI)
            {
              hpnt[j] = assign_string("LSAI[m2/m2]");
            }
          else if (j==oz0v)
            {
              hpnt[j] = assign_string("z0veg[m]");
            }
          else if (j==od0v)
            {
              hpnt[j] = assign_string("d0veg[m]");
            }
          else if (j==oEcan)
            {
              hpnt[j] = assign_string("Estored_canopy[W/m2]");
            }
          else if (j==oSWv)
            {
              hpnt[j] = assign_string("SWv[W/m2]");
            }
          else if (j==oLWv)
            {
              hpnt[j] = assign_string("LWv[W/m2]");
            }
          else if (j==oHv)
            {
              hpnt[j] = assign_string("Hv[W/m2]");
            }
          else if (j==oLEv)
            {
              hpnt[j] = assign_string("LEv[W/m2]");
            }
          else if (j==oHg0)
            {
              hpnt[j] = assign_string("Hg_unveg[W/m2]");
            }
          else if (j==oLEg0)
            {
              hpnt[j] = assign_string("LEg_unveg[W/m2]");
            }
          else if (j==oHg1)
            {
              hpnt[j] = assign_string("Hg_veg[W/m2]");
            }
          else if (j==oLEg1)
            {
              hpnt[j] = assign_string("LEg_veg[W/m2]");
            }
          else if (j==oevapsur)
            {
              hpnt[j] = assign_string("Evap_surface[mm]");
            }
          else if (j==otrasp)
            {
              hpnt[j] = assign_string("Trasp_canopy[mm]");
            }
          else if (j==owcan_rain)
            {
              hpnt[j] = assign_string("Water_on_canopy[mm]");
            }
          else if (j==owcan_snow)
            {
              hpnt[j] = assign_string("Snow_on_canopy[mm]");
            }
          else if (j==oQv)
            {
              hpnt[j] = assign_string("Qvegetation[-]");
            }
          else if (j==oQg)
            {
              hpnt[j] = assign_string("Qsurface[-]");
            }
          else if (j==oQa)
            {
              hpnt[j] = assign_string("Qair[-]");
            }
          else if (j==oQs)
            {
              hpnt[j] = assign_string("Qcanopyair[-]");
            }
          else if (j==oLobuk)
            {
              hpnt[j] = assign_string("LObukhov[m]");
            }
          else if (j==oLobukcan)
            {
              hpnt[j] = assign_string("LObukhovcanopy[m]");
            }
          else if (j==outop)
            {
              hpnt[j] = assign_string("Wind_speed_top_canopy[m/s]");
            }
          else if (j==odecay)
            {
              hpnt[j] = assign_string("Decay_of_K_in_canopy[-]");
            }
          else if (j==oSWup)
            {
              hpnt[j] = assign_string("SWup[W/m2]");
            }
          else if (j==oLWup)
            {
              hpnt[j] = assign_string("LWup[W/m2]");
            }
          else if (j==oHup)
            {
              hpnt[j] = assign_string("Hup[W/m2]");
            }
          else if (j==oLEup)
            {
              hpnt[j] = assign_string("LEup[W/m2]");
            }
          else if (j==osnowdepth)
            {
              hpnt[j] = assign_string("snow_depth[mm]");
            }
          else if (j==oSWE)
            {
              hpnt[j] = assign_string("snow_water_equivalent[mm]");
            }
          else if (j==osnowdens)
            {
              hpnt[j] = assign_string("snow_density[kg/m3]");
            }
          else if (j==osnowT)
            {
              hpnt[j] = assign_string("snow_temperature[C]");
            }
          else if (j==omrsnow)
            {
              hpnt[j] = assign_string("snow_melted[mm]");
            }
          else if (j==osrsnow)
            {
              hpnt[j] = assign_string("snow_subl[mm]");
            }
          else if (j==oblowingsnowtrans)
            {
              hpnt[j] = assign_string("snow_blown_away[mm]");
            }
          else if (j==oblowingsnowsubl)
            {
              hpnt[j] = assign_string("snow_subl_while_blown[mm]");
            }
          else if (j==oglacdepth)
            {
              hpnt[j] = assign_string("glac_depth[mm]");
            }
          else if (j==oGWE)
            {
              hpnt[j] = assign_string("glac_water_equivalent[mm]");
            }
          else if (j==oglacdens)
            {
              hpnt[j] = assign_string("glac_density[kg/m3]");
            }
          else if (j==oglacT)
            {
              hpnt[j] = assign_string("glac_temperature[C]");
            }
          else if (j==omrglac)
            {
              hpnt[j] = assign_string("glac_melted[mm]");
            }
          else if (j==osrglac)
            {
              hpnt[j] = assign_string("glac_subl[mm]");
            }
          else if (j==othawedup)
            {
              hpnt[j] = assign_string("lowest_thawed_soil_depth[mm]");
            }
          else if (j==othaweddw)
            {
              hpnt[j] = assign_string("highest_thawed_soil_depth[mm]");
            }
          else if (j==owtableup)
            {
              hpnt[j] = assign_string("lowest_water_table_depth[mm]");
            }
          else if (j==owtabledw)
            {
              hpnt[j] = assign_string("highest_water_table_depth[mm]");
            }
        }
    }

  //HeaderBasinFile
  for (j=0; j<ootot; j++)
    {
      if (strcmp(hbsn[j], string_novalue) == 0)
        {
          free(hbsn[j]);
          if (j==oodate12)
            {
              hbsn[j] = assign_string("Date12[DDMMYYYYhhmm]");
            }
          else if (j==ooJDfrom0)
            {
              hbsn[j] = assign_string("JulianDayFromYear0[days]");
            }
          else if (j==oodaysfromstart)
            {
              hbsn[j] = assign_string("TimeFromStart[days]");
            }
          else if (j==ooperiod)
            {
              hbsn[j] = assign_string("Simulation_Period");
            }
          else if (j==oorun)
            {
              hbsn[j] = assign_string("Run");
            }
          else if (j==ooprecrain)
            {
              hbsn[j] = assign_string("Prain_below_canopy[mm]");
            }
          else if (j==ooprecsnow)
            {
              hbsn[j] = assign_string("Psnow_below_canopy[mm]");
            }
          else if (j==oorainover)
            {
              hbsn[j] = assign_string("Prain_above_canopy[mm]");
            }
          else if (j==oosnowover)
            {
              hbsn[j] = assign_string("Prain_above_canopy[mm]");
            }
          else if (j==oopnet)
            {
              hbsn[j] = assign_string("Pnet[mm]");
            }
          else if (j==ooTa)
            {
              hbsn[j] = assign_string("Tair[C]");
            }
          else if (j==ooTg)
            {
              hbsn[j] = assign_string("Tsurface[C]");
            }
          else if (j==ooTv)
            {
              hbsn[j] = assign_string("Tvegetation[C]");
            }
          else if (j==ooevapsur)
            {
              hbsn[j] = assign_string("Evap_surface[mm]");
            }
          else if (j==ootrasp)
            {
              hbsn[j] = assign_string("Transpiration_canopy[mm]");
            }
          else if (j==ooLE)
            {
              hbsn[j] = assign_string("LE[W/m2]");
            }
          else if (j==ooH)
            {
              hbsn[j] = assign_string("H[W/m2]");
            }
          else if (j==ooSW)
            {
              hbsn[j] = assign_string("SW[W/m2]");
            }
          else if (j==ooLW)
            {
              hbsn[j] = assign_string("LW[W/m2]");
            }
          else if (j==ooLEv)
            {
              hbsn[j] = assign_string("LEv[W/m2]");
            }
          else if (j==ooHv)
            {
              hbsn[j] = assign_string("Hv[W/m2]");
            }
          else if (j==ooSWv)
            {
              hbsn[j] = assign_string("SWv[W/m2]");
            }
          else if (j==ooLWv)
            {
              hbsn[j] = assign_string("LWv[W/m2]");
            }
          else if (j==ooSWin)
            {
              hbsn[j] = assign_string("SWin[W/m2]");
            }
          else if (j==ooLWin)
            {
              hbsn[j] = assign_string("LWin[W/m2]");
            }
          else if (j==oomasserror)
            {
              hbsn[j] = assign_string("Mass_balance_error[mm]");
            }
          else if (j==ootimestep)
            {
              hbsn[j] = assign_string("Mean_Time_Step[s]");
            }
        }
    }

  //HeaderSnowFile
  for (j=0; j<10; j++)
    {
      if (strcmp(hsnw[j], string_novalue) == 0)
        {
          free(hsnw[j]);
          if (j==0)
            {
              hsnw[j] = assign_string("Date12[DDMMYYYYhhmm]");
            }
          else if (j==1)
            {
              hsnw[j] = assign_string("JulianDayFromYear0[days]");
            }
          else if (j==2)
            {
              hsnw[j] = assign_string("TimeFromStart[days]");
            }
          else if (j==3)
            {
              hsnw[j] = assign_string("Simulation_Period");
            }
          else if (j==4)
            {
              hsnw[j] = assign_string("Run");
            }
          else if (j==5)
            {
              hsnw[j] = assign_string("IDpoint");
            }
        }
    }

  //HeaderGlacierFile
  for (j=0; j<10; j++)
    {
      if (strcmp(hglc[j], string_novalue) == 0)
        {

          free(hglc[j]);
          if (j==0)
            {
              hglc[j] = assign_string("Date12[DDMMYYYYhhmm]");
            }
          else if (j==1)
            {
              hglc[j] = assign_string("JulianDayFromYear0[days]");
            }
          else if (j==2)
            {
              hglc[j] = assign_string("TimeFromStart[days]");
            }
          else if (j==3)
            {
              hglc[j] = assign_string("Simulation_Period");
            }
          else if (j==4)
            {
              hglc[j] = assign_string("Run");
            }
          else if (j==5)
            {
              hglc[j] = assign_string("IDpoint");
            }
          else if (j==6)
            {
              hglc[j] = assign_string("Temperature[C]");
            }
          else if (j==7)
            {
              hglc[j] = assign_string("wice[kg/m2]");
            }
          else if (j==8)
            {
              hglc[j] = assign_string("wliq[kg/m2]");
            }
          else if (j==9)
            {
              hglc[j] = assign_string("Dz[mm]");
            }
        }
    }


  //HeaderGlacierFile
  for (j=0; j<6; j++)
    {
      if (strcmp(hsl[j], string_novalue) == 0)
        {
          free(hsl[j]);
          if (j==0)
            {
              hsl[j] = assign_string("Date12[DDMMYYYYhhmm]");
            }
          else if (j==1)
            {
              hsl[j] = assign_string("JulianDayFromYear0[days]");
            }
          else if (j==2)
            {
              hsl[j] = assign_string("TimeFromStart[days]");
            }
          else if (j==3)
            {
              hsl[j] = assign_string("Simulation_Period");
            }
          else if (j==4)
            {
              hsl[j] = assign_string("Run");
            }
          else if (j==5)
            {
              hsl[j] = assign_string("IDpoint");
            }
        }
    }

  //Recovery Files
  for (j=rpsi; j<=rsux; j++)
    {
      if (strcmp(files[j], string_novalue) == 0)
        {
          free(files[j]);
          if (j==rpsi)
            {
              files[j] = assign_string("SoilPressure");
            }
          else if (j==riceg)
            {
              files[j] = assign_string("SoilIceContent");
            }
          else if (j==rTg)
            {
              files[j] = assign_string("SoilTemperature");
            }
          else if (j==rDzs)
            {
              files[j] = assign_string("SnowThickness");
            }
          else if (j==rwls)
            {
              files[j] = assign_string("SnowLiqWaterContent");
            }
          else if (j==rwis)
            {
              files[j] = assign_string("SnowIceContent");
            }
          else if (j==rTs)
            {
              files[j] = assign_string("SnowTemperature");
            }
          else if (j==rDzi)
            {
              files[j] = assign_string("GlacThickness");
            }
          else if (j==rwli)
            {
              files[j] = assign_string("GlacLiqWaterContent");
            }
          else if (j==rwii)
            {
              files[j] = assign_string("GlacIceContent");
            }
          else if (j==rTi)
            {
              files[j] = assign_string("GlacTemperature");
            }
          else if (j==rns)
            {
              files[j] = assign_string("SnowLayersNumber");
            }
          else if (j==rni)
            {
              files[j] = assign_string("GlacLayersNumber");
            }
          else if (j==rsnag)
            {
              files[j] = assign_string("SnowAge");
            }
          else if (j==rwcrn)
            {
              files[j] = assign_string("RainOnCanopy");
            }
          else if (j==rwcsn)
            {
              files[j] = assign_string("SnowOnCanopy");
            }
          else if (j==rTv)
            {
              files[j] = assign_string("VegTemperature");
            }
          else if (j==rpsich)
            {
              files[j] = assign_string("SoilChannelPressure");
            }
          else if (j==ricegch)
            {
              files[j] = assign_string("SoilChannelIceContent");
            }
          else if (j==rTgch)
            {
              files[j] = assign_string("SoilChannelTemperature");
            }
          else if (j==rTrun)
            {
              files[j] = assign_string("RunMeanSoilTemperature");
            }
          else if (j==rwrun)
            {
              files[j] = assign_string("RunMeanSoilTotWater");
            }
          else if (j==rdUrun)
            {
              files[j] = assign_string("RunSoilInternalEnergy");
            }
          else if (j==rSWErun)
            {
              files[j] = assign_string("RunMeanSWE");
            }
          else if (j==rTmaxrun)
            {
              files[j] = assign_string("RunMaxSoilTemperature");
            }
          else if (j==rTminrun)
            {
              files[j] = assign_string("RunMinSoilTemperature");
            }
          else if (j==rwmaxrun)
            {
              files[j] = assign_string("RunMaxSoilTotWater");
            }
          else if (j==rwminrun)
            {
              files[j] = assign_string("RunMinSoilTotWater");
            }
          else if (j==rtime)
            {
              files[j] = assign_string("RecoveryTime");
            }
          else if (j==rsux)
            {
              files[j] = assign_string("SuccessfulRecovery");
            }
        }
    }

  //add path to recovery files
  if (strcmp(path_rec_files, string_novalue) != 0)
    {
      temp = assign_string(path_rec_files);
      free(path_rec_files);
      path_rec_files = join_strings(temp, "/");
      free(temp);
      for (i=rpsi; i<=rsux; i++)
        {
          temp = assign_string(files[i]);
          free(files[i]);
          files[i] = join_strings(path_rec_files, temp);
          free(temp);
        }
    }
  free(path_rec_files);

  //add working path to the file name
  for (i=0; i<nfiles; i++)
    {
      if (strcmp(files[i], string_novalue) != 0)
        {
          temp = assign_string(files[i]);
          free(files[i]);
          files[i] = join_strings(WORKING_DIRECTORY, temp);
          free(temp);
        }
    }


  return 1;

}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

void assign_numeric_parameters(PAR *par, LAND *land, TIMES *times, SOIL *sl,
                               METEO *met, INIT_TOOLS *itools, double **num_param,
                               long *num_param_components, char **keyword, FILE *flog)
{

  short occurring;
  long cod, codn, i, j, k, n, m, nsoillayers, nmeteo_stations, npoints;
  double a, minDt=1.E99;

  fprintf(flog,"\n");

  par->print=0;

  //find components of times->Dt_vector
  cod = 0;
  n = (long)max_cols_time_steps_file + 1;
  times->Dt_vector=(double *)malloc(n*sizeof(double));
  times->Dt_vector[0] =
    0.;//it is the space for the date in case of time variable time step
  times->Dt_vector[1] = assignation_number(flog, cod, 0, keyword, num_param,
                                           num_param_components, 0., 1);
  for (i=2; i<n; i++)
    {
      if (i <= num_param_components[cod])
        {
          times->Dt_vector[i] = assignation_number(flog, cod, i-1, keyword, num_param,
                                                   num_param_components, (double)number_novalue, 0);
        }
      else
        {
          times->Dt_vector[i] = (double)number_novalue;
        }
    }
  for (i=1; i<n; i++)
    {
      if ((long)times->Dt_vector[i] != number_novalue)
        {
          if (times->Dt_vector[i] < minDt) minDt = times->Dt_vector[i];
        }
    }

  //init date
  cod = 1;
  par->init_date->reinit(num_param_components[cod]);
  for (i=1; i<=par->init_date->nh; i++)
    {
      par->init_date->co[i] = assignation_number(flog, cod, i-1, keyword, num_param,
                                                 num_param_components, 010119000000., 0);
      par->init_date->co[i] = convert_dateeur12_JDfrom0(par->init_date->co[i]);
    }

  //simulation time
  cod = 398;
  par->simulation_hours = assignation_number(flog, cod, 0, keyword, num_param,
                                             num_param_components, 1., 0);

  //end date
  cod = 2;
  par->end_date->reinit(num_param_components[cod]);

  if (par->end_date->nh != par->init_date->nh)
    {
      fprintf(flog,
              "Error:: End date has a number of components different from Init Date");
      printf("Error:: End date has a number of components different from Init Date");
      t_error("Fatal Error! Geotop is closed. See failing report.");
    }

  for (i=1; i<=par->end_date->nh; i++)
    {
      par->end_date->co[i] = assignation_number(flog, cod, i-1, keyword, num_param,
                                                num_param_components, (double)number_novalue, 0);
      if ((long)par->end_date->co[i] == number_novalue)
        {
          par->end_date->co[i] = par->init_date->co[i] + par->simulation_hours/24.;
        }
      else
        {
          par->end_date->co[i] = convert_dateeur12_JDfrom0(par->end_date->co[i]);
        }
    }

  //run times
  cod = 3;
  par->run_times = new_longvector(par->init_date->nh);
  par->run_times->co[1] = (long)assignation_number(flog, cod, 0, keyword,
                                                   num_param, num_param_components, 1., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->run_times->co[i] = (long)assignation_number(flog, cod, i-1, keyword,
                                                       num_param, num_param_components, (double)par->run_times->co[i-1], 0);
    }

  par->ST = assignation_number(flog, 4, 0, keyword, num_param,
                               num_param_components, 0., 0);

  cod = 5;
  par->Dtplot_discharge.reset(new Vector<double> {par->init_date->nh});
  par->Dtplot_discharge->co[1] = assignation_number(flog, cod, 0, keyword,
                                                    num_param, num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->Dtplot_discharge->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                        num_param, num_param_components, par->Dtplot_discharge->co[i-1], 0);
    }
  par->plot_discharge_with_Dt_integration = new_shortvector(par->init_date->nh);
  par->state_discharge = 0;
  for (i=1; i<=par->init_date->nh; i++)
    {
      par->Dtplot_discharge->co[i] *= 3600.;
      if (par->Dtplot_discharge->co[i] > 1.E-5
          && par->Dtplot_discharge->co[i] <= minDt)
        {
          par->plot_discharge_with_Dt_integration->co[i]=1;
        }
      else
        {
          par->plot_discharge_with_Dt_integration->co[i]=0;
        }
      if (par->Dtplot_discharge->co[i] > 1.E-5) par->state_discharge = 1;
    }

  cod = 6;
  par->Dtplot_point.reset(new Vector<double>{par->init_date->nh});
  par->Dtplot_point->co[1] = assignation_number(flog, cod, 0, keyword,
                                                num_param, num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->Dtplot_point->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                    num_param, num_param_components, par->Dtplot_point->co[i-1], 0);
    }
  par->plot_point_with_Dt_integration = new_shortvector(par->init_date->nh);
  par->state_pixel = 0;
  for (i=1; i<=par->init_date->nh; i++)
    {
      par->Dtplot_point->co[i] *= 3600.;
      if (par->Dtplot_point->co[i] > 1.E-5 && par->Dtplot_point->co[i] <= minDt)
        {
          par->plot_point_with_Dt_integration->co[i]=1;
        }
      else
        {
          par->plot_point_with_Dt_integration->co[i]=0;
        }
      if (par->Dtplot_point->co[i] > 1.E-5) par->state_pixel = 1;
    }

  cod = 7;
  par->Dtplot_basin.reset(new Vector<double>{par->init_date->nh});
  par->Dtplot_basin->co[1] = assignation_number(flog, cod, 0, keyword,
                                                num_param, num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->Dtplot_basin->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                    num_param, num_param_components, par->Dtplot_basin->co[i-1], 0);
    }
  par->plot_basin_with_Dt_integration = new_shortvector(par->init_date->nh);
  par->state_basin = 0;
  for (i=1; i<=par->init_date->nh; i++)
    {
      par->Dtplot_basin->co[i] *= 3600.;
      if (par->Dtplot_basin->co[i] > 1.E-5 && par->Dtplot_basin->co[i] <= minDt)
        {
          par->plot_basin_with_Dt_integration->co[i]=1;
        }
      else
        {
          par->plot_basin_with_Dt_integration->co[i]=0;
        }
      if (par->Dtplot_basin->co[i] > 1.E-5) par->state_basin = 1;
    }


  par->lowpass = (long)assignation_number(flog, 8, 0, keyword, num_param,
                                          num_param_components, 0., 0);
  par->lowpass_curvatures = (long)assignation_number(flog, 9, 0, keyword,
                                                     num_param, num_param_components, 0., 0);
  par->sky = (short)assignation_number(flog, 10, 0, keyword, num_param,
                                       num_param_components, 0., 0);
  par->format_out = (short)assignation_number(flog, 11, 0, keyword, num_param,
                                              num_param_components, 3., 0);
  par->point_sim = (short)assignation_number(flog, 12, 0, keyword, num_param,
                                             num_param_components, 0., 0);
  par->recover = (short)assignation_number(flog, 13, 0, keyword, num_param,
                                           num_param_components, 0., 0);


  //land cover types
  par->n_landuses = (long)assignation_number(flog, 14, 0, keyword, num_param,
                                             num_param_components, 1., 0);

  land->ty = new_doublematrix(par->n_landuses, nlandprop);

  land->ty->co[1][jz0] = assignation_number(flog, 15, 0, keyword, num_param,
                                            num_param_components, 10., 0);
  land->ty->co[1][jz0thressoil] = assignation_number(flog, 16, 0, keyword,
                                                     num_param, num_param_components, land->ty->co[1][jz0], 0);
  land->ty->co[1][jHveg] = assignation_number(flog, 17, 0, keyword, num_param,
                                              num_param_components, 1000., 0);
  land->ty->co[1][jz0thresveg] = assignation_number(flog, 18, 0, keyword,
                                                    num_param, num_param_components, land->ty->co[1][jHveg], 0);
  land->ty->co[1][jz0thresveg2] = assignation_number(flog, 19, 0, keyword,
                                                     num_param, num_param_components, land->ty->co[1][jz0thresveg], 0);
  land->ty->co[1][jLSAI] = assignation_number(flog, 20, 0, keyword, num_param,
                                              num_param_components, 1., 0);
  land->ty->co[1][jcf] = assignation_number(flog, 21, 0, keyword, num_param,
                                            num_param_components, 0., 0);
  land->ty->co[1][jdecay0] = assignation_number(flog, 22, 0, keyword, num_param,
                                                num_param_components, 2.5, 0);
  land->ty->co[1][jexpveg] = assignation_number(flog, 23, 0, keyword, num_param,
                                                num_param_components, 1., 0);
  land->ty->co[1][jroot] = assignation_number(flog, 24, 0, keyword, num_param,
                                              num_param_components, 300., 0);
  land->ty->co[1][jrs] = assignation_number(flog, 25, 0, keyword, num_param,
                                            num_param_components, 60., 0);
  land->ty->co[1][jvR_vis] = assignation_number(flog, 26, 0, keyword, num_param,
                                                num_param_components, 0.2, 0);
  land->ty->co[1][jvR_nir] = assignation_number(flog, 27, 0, keyword, num_param,
                                                num_param_components, 0.2, 0);
  land->ty->co[1][jvT_vis] = assignation_number(flog, 28, 0, keyword, num_param,
                                                num_param_components, 0.2, 0);
  land->ty->co[1][jvT_nir] = assignation_number(flog, 29, 0, keyword, num_param,
                                                num_param_components, 0.2, 0);
  land->ty->co[1][jvCh] = assignation_number(flog, 30, 0, keyword, num_param,
                                             num_param_components, 0., 0);
  land->ty->co[1][jcd] = assignation_number(flog, 31, 0, keyword, num_param,
                                            num_param_components, 2., 0);
  land->ty->co[1][ja_vis_dry] = assignation_number(flog, 32, 0, keyword,
                                                   num_param, num_param_components, 0.2, 0);
  land->ty->co[1][ja_nir_dry] = assignation_number(flog, 33, 0, keyword,
                                                   num_param, num_param_components, land->ty->co[1][ja_vis_dry], 0);
  land->ty->co[1][ja_vis_sat] = assignation_number(flog, 34, 0, keyword,
                                                   num_param, num_param_components, land->ty->co[1][ja_vis_dry], 0);
  land->ty->co[1][ja_nir_sat] = assignation_number(flog, 35, 0, keyword,
                                                   num_param, num_param_components, land->ty->co[1][ja_nir_dry], 0);
  land->ty->co[1][jemg] = assignation_number(flog, 36, 0, keyword, num_param,
                                             num_param_components, 0.96, 0);
  land->ty->co[1][jcm] = assignation_number(flog, 37, 0, keyword, num_param,
                                            num_param_components, 0.5, 0);
  land->ty->co[1][jN] = assignation_number(flog, 38, 0, keyword, num_param,
                                           num_param_components, 0., 0);
  land->ty->co[1][jdv] = assignation_number(flog, 39, 0, keyword, num_param,
                                            num_param_components, 50., 0);

  for (i=2; i<=par->n_landuses; i++)
    {
      for (j=1; j<=nlandprop; j++)
        {
          land->ty->co[i][j] = assignation_number(flog, 15+j-1, i-1, keyword, num_param,
                                                  num_param_components, land->ty->co[i-1][j], 0);
        }
    }

  //former block 2
  par->imp = assignation_number(flog, 40, 0, keyword, num_param,
                                num_param_components, 7., 0);
  par->free_drainage_bottom = assignation_number(flog, 41, 0, keyword,
                                                 num_param, num_param_components, 0., 0);
  par->free_drainage_lateral = assignation_number(flog, 42, 0, keyword,
                                                  num_param, num_param_components, 1., 0);
  par->TolVWb = assignation_number(flog, 43, 0, keyword, num_param,
                                   num_param_components, 1.E-6, 0);
  par->RelTolVWb = RelativeErrorRichards;
  par->MaxErrWb = 1.E99;
  par->MaxiterTol = (long)assignation_number(flog, 44, 0, keyword, num_param,
                                             num_param_components, 100., 0);
  par->TolCG = assignation_number(flog, 45, 0, keyword, num_param,
                                  num_param_components, 0.01, 0);
  par->min_lambda_wat = assignation_number(flog, 46, 0, keyword, num_param,
                                           num_param_components, 1.E-7, 0);
  par->max_times_min_lambda_wat = (long)assignation_number(flog, 47, 0, keyword,
                                                           num_param, num_param_components, 0.0, 0);
  par->exit_lambda_min_wat = (short)assignation_number(flog, 48, 0, keyword,
                                                       num_param, num_param_components, 1., 0);
  par->min_Dt = assignation_number(flog, 49, 0, keyword, num_param,
                                   num_param_components, 10., 0);
  par->gamma_m = assignation_number(flog, 50, 0, keyword, num_param,
                                    num_param_components, 2./3., 0);
  par->thres_hsup_1 = assignation_number(flog, 51, 0, keyword, num_param,
                                         num_param_components, 0., 0);
  par->thres_hsup_2 = assignation_number(flog, 52, 0, keyword, num_param,
                                         num_param_components, 0., 0);
  par->Ks_channel = assignation_number(flog, 53, 0, keyword, num_param,
                                       num_param_components, 20., 0);
  par->thres_hchannel = assignation_number(flog, 54, 0, keyword, num_param,
                                           num_param_components, 50., 0);
  par->w_dx = assignation_number(flog, 55, 0, keyword, num_param,
                                 num_param_components, 0.1, 0);
  par->depr_channel = assignation_number(flog, 56, 0, keyword, num_param,
                                         num_param_components, 500., 0);
  par->max_courant_land = assignation_number(flog, 57, 0, keyword, num_param,
                                             num_param_components, 0.1, 0);
  par->max_courant_channel = assignation_number(flog, 58, 0, keyword, num_param,
                                                num_param_components, 0.1, 0);
  par->min_hsup_land = assignation_number(flog, 59, 0, keyword, num_param,
                                          num_param_components, 1., 0);
  par->min_hsup_channel = assignation_number(flog, 60, 0, keyword, num_param,
                                             num_param_components, 1., 0);
  par->min_dhsup_land_channel_in = assignation_number(flog, 61, 0, keyword,
                                                      num_param, num_param_components, 1., 0);
  par->dtmin_sup = assignation_number(flog, 62, 0, keyword, num_param,
                                      num_param_components, 0.01, 0);
  //former block 3
  par->latitude = assignation_number(flog, 63, 0, keyword, num_param,
                                     num_param_components, 45., 0);
  par->longitude = assignation_number(flog, 64, 0, keyword, num_param,
                                      num_param_components, 0., 0);
  par->Vmin = assignation_number(flog, 65, 0, keyword, num_param,
                                 num_param_components, 0.5, 0);
  par->RHmin = assignation_number(flog, 66, 0, keyword, num_param,
                                  num_param_components, 10., 0);
  par->alpha_snow = assignation_number(flog, 67, 0, keyword, num_param,
                                       num_param_components, 1.E5, 0);
  par->nsurface = (long)assignation_number(flog, 68, 0, keyword, num_param,
                                           num_param_components, 0., 0);
  par->tol_energy = assignation_number(flog, 69, 0, keyword, num_param,
                                       num_param_components, 1.E-4, 0);
  par->maxiter_energy = (long)assignation_number(flog, 70, 0, keyword,
                                                 num_param, num_param_components, 500., 0);
  par->min_lambda_en = assignation_number(flog, 71, 0, keyword, num_param,
                                          num_param_components, 1.E-5, 0);
  par->max_times_min_lambda_en = (long)assignation_number(flog, 72, 0, keyword,
                                                          num_param, num_param_components, 0.0, 0);
  par->exit_lambda_min_en = (short)assignation_number(flog, 73, 0, keyword,
                                                      num_param, num_param_components, 0., 0);
  par->dem_rotation = assignation_number(flog, 74, 0, keyword, num_param,
                                         num_param_components, 0., 0);
  par->maxiter_canopy = (long)assignation_number(flog, 75, 0, keyword,
                                                 num_param, num_param_components, 3., 0);
  par->maxiter_Businger = (long)assignation_number(flog, 76, 0, keyword,
                                                   num_param, num_param_components, 5., 0);
  par->maxiter_Ts = (long)assignation_number(flog, 77, 0, keyword, num_param,
                                             num_param_components, 2., 0);
  par->maxiter_Loc = (long)assignation_number(flog, 78, 0, keyword, num_param,
                                              num_param_components, 3., 0);
  par->stabcorr_incanopy = (short)assignation_number(flog, 79, 0, keyword,
                                                     num_param, num_param_components, 1., 0);
  par->iobsint = (short)assignation_number(flog, 80, 0, keyword, num_param,
                                           num_param_components, 1., 0);
  par->dn = assignation_number(flog, 81, 0, keyword, num_param,
                               num_param_components, 1., 0);
  par->slopewt = assignation_number(flog, 82, 0, keyword, num_param,
                                    num_param_components, 0., 0);
  par->curvewt = assignation_number(flog, 83, 0, keyword, num_param,
                                    num_param_components, 0., 0);
  par->slopewtD = assignation_number(flog, 84, 0, keyword, num_param,
                                     num_param_components, par->slopewt, 0);
  par->curvewtD = assignation_number(flog, 85, 0, keyword, num_param,
                                     num_param_components, par->curvewt, 0);
  par->slopewtI = assignation_number(flog, 86, 0, keyword, num_param,
                                     num_param_components, par->slopewt, 0);
  par->curvewtI = assignation_number(flog, 87, 0, keyword, num_param,
                                     num_param_components, par->curvewt, 0);
  par->Zboundary = assignation_number(flog, 88, 0, keyword, num_param,
                                      num_param_components, 1.E20, 0);
  par->Tboundary = assignation_number(flog, 89, 0, keyword, num_param,
                                      num_param_components, 20., 0);
  par->Fboundary = assignation_number(flog, 90, 0, keyword, num_param,
                                      num_param_components, 0., 0);
  //former block 4
  itools->swe0 = assignation_number(flog, 91, 0, keyword, num_param,
                                    num_param_components, 0., 0);
  itools->rhosnow0 = assignation_number(flog, 92, 0, keyword, num_param,
                                        num_param_components, 200., 0);
  itools->Tsnow0 = assignation_number(flog, 93, 0, keyword, num_param,
                                      num_param_components, -3., 0);
  itools->agesnow0 = assignation_number(flog, 94, 0, keyword, num_param,
                                        num_param_components, 0., 0);
  par->T_rain = assignation_number(flog, 95, 0, keyword, num_param,
                                   num_param_components, 3., 0);
  par->T_snow = assignation_number(flog, 96, 0, keyword, num_param,
                                   num_param_components, -1., 0);
  par->dew = (short)assignation_number(flog, 97, 0, keyword, num_param,
                                       num_param_components, 0., 0);
  par->aep = assignation_number(flog, 98, 0, keyword, num_param,
                                num_param_components, 10., 0);
  par->avo = assignation_number(flog, 99, 0, keyword, num_param,
                                num_param_components, 0.9, 0);
  par->airo = assignation_number(flog, 100, 0, keyword, num_param,
                                 num_param_components, 0.65, 0);
  par->Sr = assignation_number(flog, 101, 0, keyword, num_param,
                               num_param_components, 0.02, 0);
  par->epsilon_snow = assignation_number(flog, 102, 0, keyword, num_param,
                                         num_param_components, 0.98, 0);
  par->z0_snow = 0.001*assignation_number(flog, 103, 0, keyword, num_param,
                                          num_param_components, 0.1, 0);
  par->snowcorrfact = assignation_number(flog, 104, 0, keyword, num_param,
                                         num_param_components, 1., 0);
  par->raincorrfact = assignation_number(flog, 105, 0, keyword, num_param,
                                         num_param_components, 1., 0);
  par->snow_maxpor = assignation_number(flog, 106, 0, keyword, num_param,
                                        num_param_components, 0.7, 0);
  par->drysnowdef_rate = assignation_number(flog, 107, 0, keyword, num_param,
                                            num_param_components, 1., 0);
  par->snow_density_cutoff = assignation_number(flog, 108, 0, keyword,
                                                num_param, num_param_components, 100., 0);
  par->wetsnowdef_rate = assignation_number(flog, 109, 0, keyword, num_param,
                                            num_param_components, 1.5, 0);
  par->snow_viscosity = assignation_number(flog, 110, 0, keyword, num_param,
                                           num_param_components, 1.E6, 0);
  par->fetch_up = assignation_number(flog, 111, 0, keyword, num_param,
                                     num_param_components, 1000., 0);
  par->fetch_down = assignation_number(flog, 112, 0, keyword, num_param,
                                       num_param_components, 100., 0);
  par->Wice_PBSM = assignation_number(flog, 113, 0, keyword, num_param,
                                      num_param_components, 0., 0);
  par->Dt_PBSM = assignation_number(flog, 114, 0, keyword, num_param,
                                    num_param_components, times->Dt_vector[1], 0);
  par->snow_smin = assignation_number(flog, 115, 0, keyword, num_param,
                                      num_param_components, 30., 0);
  par->snow_smax = assignation_number(flog, 116, 0, keyword, num_param,
                                      num_param_components, 80., 0);
  par->snow_curv = assignation_number(flog, 117, 0, keyword, num_param,
                                      num_param_components, -200., 0);

  //former blocks 5/6
  par->max_weq_snow = assignation_number(flog, 118, 0, keyword, num_param,
                                         num_param_components, 5., 0);
  par->max_snow_layers = (long)assignation_number(flog, 119, 0, keyword,
                                                  num_param, num_param_components, 10., 0);

  cod = 120;
  par->max_weq_snow = assignation_number(flog, 118, 0, keyword, num_param,
                                         num_param_components, 5., 0);
  n = (long)assignation_number(flog, 119, 0, keyword, num_param,
                               num_param_components, 2., 0);
  if (n < 1)
    {
      fprintf(flog,"Error:: %s must be 1 or larger\n",keyword[119]);
      printf("Error:: %s must be 1 or larger\n",keyword[119]);
      t_error("Fatal Error! Geotop is closed.");
    }
  par->SWE_bottom = assignation_number(flog, 120, 0, keyword, num_param,
                                       num_param_components, 20., 0);
  par->SWE_top = assignation_number(flog, 121, 0, keyword, num_param,
                                    num_param_components, 20., 0);
  par->max_snow_layers = (long)floor(par->SWE_bottom/par->max_weq_snow) +
                         (long)floor(par->SWE_top/par->max_weq_snow) + n;
  par->inf_snow_layers = new_longvector(n);
  fprintf(flog,
          "Max snow layer number: %ld, of which %.0f at the bottom, %ld in the middle, and %.0f at the top.\n",
          par->max_snow_layers,floor(par->SWE_bottom/par->max_weq_snow),n,
          floor(par->SWE_top/par->max_weq_snow));
  fprintf(flog,"Infinite Snow layer numbers are numbers: ");
  for (i=1; i<=n; i++)
    {
      par->inf_snow_layers->co[i] = (long)floor(par->SWE_bottom/par->max_weq_snow) +
                                    i;
      fprintf(flog, "%ld ",par->inf_snow_layers->co[i]);
    }
  fprintf(flog,"\n");

  //former block 7
  itools->Dglac0 = assignation_number(flog, 122, 0, keyword, num_param,
                                      num_param_components, 0., 0);
  itools->rhoglac0 = assignation_number(flog, 123, 0, keyword, num_param,
                                        num_param_components, 800., 0);
  itools->Tglac0 = assignation_number(flog, 124, 0, keyword, num_param,
                                      num_param_components, -3., 0);
  par->Sr_glac = assignation_number(flog, 125, 0, keyword, num_param,
                                    num_param_components, 0.02, 0);

  //former block 8
  par->max_weq_glac = assignation_number(flog, 126, 0, keyword, num_param,
                                         num_param_components, 5., 0);
  n = (long)assignation_number(flog, 127, 0, keyword, num_param,
                               num_param_components, 0., 0);

  par->GWE_bottom = assignation_number(flog, 128, 0, keyword, num_param,
                                       num_param_components, 0., 0);
  par->GWE_top = assignation_number(flog, 129, 0, keyword, num_param,
                                    num_param_components, 0., 0);

  if (n < 1 && (par->GWE_bottom > 0 || par->GWE_top > 0))
    {
      fprintf(flog,"Error:: %s must be 1 or larger\n",keyword[127]);
      printf("Error:: %s must be 1 or larger\n",keyword[127]);
      t_error("Fatal Error! Geotop is closed.");
    }

  par->max_glac_layers = (long)floor(par->GWE_bottom/par->max_weq_glac) +
                         (long)floor(par->GWE_top/par->max_weq_glac) + n;
  par->inf_glac_layers = new_longvector(n);
  fprintf(flog,
          "Max glac layer number: %ld, of which %.0f at the bottom, %ld in the middle, and %.0f at the top.\n",
          par->max_glac_layers,floor(par->GWE_bottom/par->max_weq_glac),n,
          floor(par->GWE_top/par->max_weq_glac));
  fprintf(flog,"Infinite Glac layer numbers are numbers: ");
  for (i=1; i<=n; i++)
    {
      par->inf_glac_layers->co[i] = (long)floor(par->GWE_bottom/par->max_weq_glac) +
                                    i;
      fprintf(flog, "%ld ",par->inf_glac_layers->co[i]);
    }
  fprintf(flog,"\n");

  //former block 9
  par->state_turb = 1;
  par->state_lwrad = (short)assignation_number(flog, 130, 0, keyword, num_param,
                                               num_param_components, 9., 0);
  par->monin_obukhov = (short)assignation_number(flog, 131, 0, keyword,
                                                 num_param, num_param_components, 1., 0);
  par->surroundings = (short)assignation_number(flog, 132, 0, keyword,
                                                num_param, num_param_components, 0., 0);

  //distributed option file
  par->wat_balance = (short)assignation_number(flog, 133, 0, keyword, num_param,
                                               num_param_components, 0., 0);
  par->en_balance = (short)assignation_number(flog, 134, 0, keyword, num_param,
                                              num_param_components, 0., 0);
  par->blowing_snow = (short)assignation_number(flog, 135, 0, keyword,
                                                num_param, num_param_components, 0., 0);

  par->Wmin_BS = assignation_number(flog, 136, 0, keyword, num_param,
                                    num_param_components, 8., 0);

  cod = 137;
  npoints = 0;
  for (j=1; j<=16; j++)
    {
      if (npoints < num_param_components[cod + j-1]) npoints =
          num_param_components[cod + j-1];
    }

  if (par->point_sim == 1)
    {
      par->chkpt = new_doublematrix(npoints, ptTOT);
    }
  else
    {
      par->chkpt = new_doublematrix(npoints, 3);
    }

  for (i=1; i<=par->chkpt->nrh; i++)
    {
      for (j=1; j<=par->chkpt->nch; j++)
        {
          par->chkpt->co[i][j] = assignation_number(flog, cod + j-1, i-1, keyword,
                                                    num_param, num_param_components, (double)number_novalue, 0);
        }
    }

  cod = 156;
  par->saving_points.reset(new Vector<double>{num_param_components[cod]});
  for (i=1; i<=par->saving_points->nh; i++)
    {
      par->saving_points->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                     num_param, num_param_components, 0., 0);
    }

  par->output_soil.reset(new Vector<double>{par->init_date->nh});
  par->output_snow.reset(new Vector<double>{par->init_date->nh});
  par->output_glac.reset(new Vector<double>{par->init_date->nh});
  par->output_surfenergy.reset(new Vector<double>{par->init_date->nh});
  par->output_vegetation.reset(new Vector<double>{par->init_date->nh});
  par->output_meteo.reset(new Vector<double>{par->init_date->nh});

  cod = 157;
  par->output_soil->co[1] = assignation_number(flog, cod, 0, keyword, num_param,
                                               num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->output_soil->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                   num_param, num_param_components, par->output_soil->co[i-1], 0);
    }

  cod = 158;
  par->output_snow->co[1] = assignation_number(flog, cod, 0, keyword, num_param,
                                               num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->output_snow->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                   num_param, num_param_components, par->output_snow->co[i-1], 0);
    }

  cod = 159;
  par->output_glac->co[1] = assignation_number(flog, cod, 0, keyword, num_param,
                                               num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->output_glac->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                   num_param, num_param_components, par->output_glac->co[i-1], 0);
    }

  cod = 160;
  par->output_surfenergy->co[1] = assignation_number(flog, cod, 0, keyword,
                                                     num_param, num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->output_surfenergy->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                         num_param, num_param_components, par->output_surfenergy->co[i-1], 0);
    }

  cod = 161;
  par->output_vegetation->co[1] = assignation_number(flog, cod, 0, keyword,
                                                     num_param, num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->output_vegetation->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                         num_param, num_param_components, par->output_vegetation->co[i-1], 0);
    }

  cod = 162;
  par->output_meteo->co[1] = assignation_number(flog, cod, 0, keyword,
                                                num_param, num_param_components, 0., 0);
  for (i=2; i<=par->init_date->nh; i++)
    {
      par->output_meteo->co[i] = assignation_number(flog, cod, i-1, keyword,
                                                    num_param, num_param_components, par->output_meteo->co[i-1], 0);
    }

  par->output_soil_bin = 0;
  par->output_snow_bin = 0;
  par->output_glac_bin = 0;
  par->output_surfenergy_bin = 0;
  par->output_meteo_bin = 0;

  for (i=1; i<=par->init_date->nh; i++)
    {
      if (par->output_soil->co[i] > 0) par->output_soil_bin = 1;
      if (par->output_snow->co[i] > 0) par->output_snow_bin = 1;
      if (par->output_glac->co[i] > 0) par->output_glac_bin = 1;
      if (par->output_surfenergy->co[i] > 0) par->output_surfenergy_bin = 1;
      if (par->output_meteo->co[i] > 0) par->output_meteo_bin = 1;
    }

  cod = 163;
  codn = 164;

  if (num_param_components[cod] != num_param_components[codn])
    {
      fprintf(flog,
              "Error:: Number of components of parameters %s and %s must be equal\n",
              keyword[cod],keyword[codn]);
      printf("Error:: Number of components of parameters %s and %s must be equal\n",
             keyword[cod],keyword[codn]);
      t_error("Fatal Error! Geotop is closed. See failing report.");
    }

  times->JD_plots.reset(new Vector<double> {num_param_components[cod] +
                                     num_param_components[codn]});
  for (i=1; i<=(long)(times->JD_plots->nh/2.); i++)
    {
      times->JD_plots->co[2*i-1] = assignation_number(flog, cod, i-1, keyword,
                                                      num_param, num_param_components, 0., 0);
      times->JD_plots->co[2*i  ] = assignation_number(flog, codn, i-1, keyword,
                                                      num_param, num_param_components, 0., 0);
    }
  if (times->JD_plots->nh == 2 && times->JD_plots->co[1] < 1.E-5
      && times->JD_plots->co[2] < 1.E-5)
    {
      times->JD_plots.reset(new Vector<double>{1});
    }
  if (times->JD_plots->nh > 1)
    {
      for (i=1; i<=times->JD_plots->nh; i++)
        {
          times->JD_plots->co[i] = convert_dateeur12_JDfrom0(times->JD_plots->co[i]);
        }
    }

  //initial condition on the water pressure
  par->nsoiltypes = (long)assignation_number(flog, 165, 0, keyword, num_param,
                                             num_param_components, 1., 0);
  if (par->nsoiltypes < 1) par->nsoiltypes = 1;

  cod = 166;
  itools->init_water_table_depth = new_doublevector(par->nsoiltypes);
  itools->init_water_table_depth->co[1] = assignation_number(flog, cod, 0,
                                                             keyword, num_param, num_param_components, 5000., 0);
  for (k=2; k<=par->nsoiltypes; k++)
    {
      itools->init_water_table_depth->co[k] = assignation_number(flog, cod, k-1,
                                                                 keyword, num_param, num_param_components,
                                                                 itools->init_water_table_depth->co[k-1], 0);
    }

  //soil properties and discretization
  par->soil_type_land_default = (long)assignation_number(flog, 167, 0, keyword,
                                                         num_param, num_param_components, 1., 0);

  if (par->soil_type_land_default<1
      || par->soil_type_land_default>par->nsoiltypes)
    {
      fprintf(flog,
              "Error:: Soil_type_land_default lower than 0 or higher than soil types numbers");
      printf("Error:: Soil_type_land_default lower than 0 or higher than soil types numbers");
      t_error("Fatal Error! Geotop is closed. See failing report.");
    }

  par->soil_type_chan_default = (long)assignation_number(flog, 168, 0, keyword,
                                                         num_param, num_param_components, 1., 0);

  if (par->soil_type_chan_default<1
      || par->soil_type_chan_default>par->nsoiltypes)
    {
      fprintf(flog,
              "Error:  Soil_type_chan_default lower than 0 or higher than soil types numbers");
      printf("Error:  Soil_type_chan_default lower than 0 or higher than soil types numbers");
      t_error("Fatal Error! Geotop is closed. See failing report.");
    }

  cod = 169;
  a = assignation_number(flog, cod, 0, keyword, num_param, num_param_components,
                         (double)number_novalue, 0);

  //there is a specific layer discretization
  if ( (long)a != number_novalue && num_param_components[cod] > 1)
    {

      nsoillayers = num_param_components[cod];
      sl->pa = new_doubletensor(par->nsoiltypes, nsoilprop, nsoillayers);

      sl->pa->co[1][jdz][1] = a;
      for (i=2; i<=sl->pa->nch; i++)
        {
          sl->pa->co[1][jdz][i] = assignation_number(flog, cod, i-1, keyword, num_param,
                                                     num_param_components, sl->pa->co[1][jdz][i-1], 0);
        }

    }
  else
    {

      if ((long)a==number_novalue) a=100.;
      nsoillayers = (long)assignation_number(flog, 170, 0, keyword, num_param,
                                             num_param_components, 5., 0);

      sl->pa = new_doubletensor(par->nsoiltypes, nsoilprop, nsoillayers);

      for (i=1; i<=sl->pa->nch; i++)
        {
          sl->pa->co[1][jdz][i] = a;
        }

    }

  //first layer
  i = 1;
  sl->pa->co[1][jpsi][i] = assignation_number(flog, 171, i-1, keyword,
                                              num_param, num_param_components, (double)number_novalue, 0);
  sl->pa->co[1][jT][i] = assignation_number(flog, 172, i-1, keyword, num_param,
                                            num_param_components, 5., 0);
  sl->pa->co[1][jKn][i] = assignation_number(flog, 173, i-1, keyword, num_param,
                                             num_param_components, 1.E-4, 0);
  sl->pa->co[1][jKl][i] = assignation_number(flog, 174, i-1, keyword, num_param,
                                             num_param_components, 1.E-4, 0);
  sl->pa->co[1][jres][i] = assignation_number(flog, 175, i-1, keyword,
                                              num_param, num_param_components, 0.05, 0);
  sl->pa->co[1][jwp][i] = assignation_number(flog, 176, i-1, keyword, num_param,
                                             num_param_components, (double)number_novalue, 0);
  sl->pa->co[1][jfc][i] = assignation_number(flog, 177, i-1, keyword, num_param,
                                             num_param_components, (double)number_novalue, 0);
  sl->pa->co[1][jsat][i] = assignation_number(flog, 178, i-1, keyword,
                                              num_param, num_param_components, 0.5, 0);
  sl->pa->co[1][ja][i] = assignation_number(flog, 179, i-1, keyword, num_param,
                                            num_param_components, 0.004, 0);
  sl->pa->co[1][jns][i] = assignation_number(flog, 180, i-1, keyword, num_param,
                                             num_param_components, 1.3, 0);
  sl->pa->co[1][jv][i] = assignation_number(flog, 181, i-1, keyword, num_param,
                                            num_param_components, 0.5, 0);
  sl->pa->co[1][jkt][i] = assignation_number(flog, 182, i-1, keyword, num_param,
                                             num_param_components, 2.5, 0);
  sl->pa->co[1][jct][i] = assignation_number(flog, 183, i-1, keyword, num_param,
                                             num_param_components, 1.E6, 0);
  sl->pa->co[1][jss][i] = assignation_number(flog, 184, i-1, keyword, num_param,
                                             num_param_components, 1.E-7, 0);

  //other layers
  for (i=2; i<=sl->pa->nch; i++)
    {
      for (j=2; j<=sl->pa->nrh; j++)
        {
          sl->pa->co[1][j][i] = assignation_number(flog, 171 + j-2, i-1, keyword,
                                                   num_param, num_param_components, sl->pa->co[1][j][i-1], 0);
        }
    }

  //field capacity (-0.333 bar) and wilting point (-15 bar)
  for (i=1; i<=sl->pa->nch; i++)
    {
      if ( (long)sl->pa->co[1][jfc][i] == number_novalue)
        {
          sl->pa->co[1][jfc][i] = teta_psi( (-1./3.)*1.E5/g, 0., sl->pa->co[1][jsat][i],
                                            sl->pa->co[1][jres][i], sl->pa->co[1][ja][i],
                                            sl->pa->co[1][jns][i], 1.-1./sl->pa->co[1][jns][i], PsiMin,
                                            sl->pa->co[1][jss][i]);
        }

      if ( (long)sl->pa->co[1][jwp][i] == number_novalue)
        {
          sl->pa->co[1][jwp][i] = teta_psi( -15.*1.E5/g, 0., sl->pa->co[1][jsat][i],
                                            sl->pa->co[1][jres][i], sl->pa->co[1][ja][i],
                                            sl->pa->co[1][jns][i], 1.-1./sl->pa->co[1][jns][i], PsiMin,
                                            sl->pa->co[1][jss][i]);
        }
    }

  //other soil types
  for (k=2; k<=par->nsoiltypes; k++)
    {
      for (i=1; i<=sl->pa->nch; i++)
        {
          for (j=1; j<=sl->pa->nrh; j++)
            {
              sl->pa->co[k][j][i] = sl->pa->co[1][j][i];
            }
        }
    }

  //use water table for water pressure
  for (k=1; k<=par->nsoiltypes; k++)
    {
      occurring = 0;//check if psi initial has at least one novalue
      for (i=1; i<=sl->pa->nch; i++)
        {
          if ( (long)sl->pa->co[k][jpsi][i] == number_novalue) occurring = 1;
        }
      if (occurring == 0) itools->init_water_table_depth->co[k] =
          (double)number_novalue;
    }

  cod = 185;
  itools->pa_bed = new_doubletensor(1, nsoilprop, nsoillayers);
  for (i=1; i<=nsoillayers; i++)
    {
      itools->pa_bed->co[1][jdz][i] = sl->pa->co[1][jdz][i];
    }
  for (j=1; j<=nsoilprop; j++)
    {
      if (j != jdz)
        {
          itools->pa_bed->co[1][j][1] = assignation_number(flog, cod+j-2, 0, keyword,
                                                           num_param, num_param_components, (double)number_novalue, 0);
        }
    }
  for (i=2; i<=nsoillayers; i++)
    {
      for (j=1; j<=nsoilprop; j++)
        {
          if (j != jdz) itools->pa_bed->co[1][j][i] = assignation_number(flog, cod+j-2,
                                                        i-1, keyword, num_param, num_param_components, itools->pa_bed->co[1][j][i-1],
                                                        0);
        }
    }
  //field capacity (-0.333 bar) and wilting point (-15 bar)
  for (i=1; i<=sl->pa->nch; i++)
    {
      if ( (long)itools->pa_bed->co[1][jsat][i] != number_novalue
           && (long)itools->pa_bed->co[1][jres][i] != number_novalue &&
           (long)itools->pa_bed->co[1][ja][i] != number_novalue
           && (long)itools->pa_bed->co[1][jns][i] != number_novalue &&
           (long)itools->pa_bed->co[1][jss][i] )
        {
          if ( (long)itools->pa_bed->co[1][jfc][i] == number_novalue)
            {
              itools->pa_bed->co[1][jfc][i] = teta_psi( (-1./3.)*1.E5/g, 0.,
                                                        itools->pa_bed->co[1][jsat][i], itools->pa_bed->co[1][jres][i],
                                                        itools->pa_bed->co[1][ja][i],
                                                        itools->pa_bed->co[1][jns][i], 1.-1./itools->pa_bed->co[1][jns][i], PsiMin,
                                                        itools->pa_bed->co[1][jss][i]);
            }
          if ( (long)itools->pa_bed->co[1][jwp][i] == number_novalue)
            {
              itools->pa_bed->co[1][jwp][i] = teta_psi( -15.*1.E5/g, 0.,
                                                        itools->pa_bed->co[1][jsat][i], itools->pa_bed->co[1][jres][i],
                                                        itools->pa_bed->co[1][ja][i],
                                                        itools->pa_bed->co[1][jns][i], 1.-1./itools->pa_bed->co[1][jns][i], PsiMin,
                                                        itools->pa_bed->co[1][jss][i]);
            }
        }
    }

  //meteo stations
  cod = 199;
  met->imeteo_stations = new_longvector(num_param_components[cod]);
  met->imeteo_stations->co[1] = (long)assignation_number(flog, cod, 0, keyword,
                                                         num_param, num_param_components, (double)number_novalue, 0);
  if ( met->imeteo_stations->co[1] != number_novalue )
    {
      for (i=2; i<=num_param_components[cod]; i++)
        {
          met->imeteo_stations->co[i] = (long)assignation_number(flog, cod, i-1,
                                                                 keyword, num_param, num_param_components, 0., 1);
        }
      nmeteo_stations = num_param_components[cod];
    }
  else
    {
      nmeteo_stations = (long)assignation_number(flog, 200, 0, keyword, num_param,
                                                 num_param_components, 1., 0);
    }

  met->st=(METEO_STATIONS *)malloc(sizeof(METEO_STATIONS));
  if (!met->st) t_error("meteo_stations was not allocated");
  met->st->E.reset(new Vector<double>{nmeteo_stations});
  met->st->N.reset(new Vector<double>{nmeteo_stations});
  met->st->lat.reset(new Vector<double>{nmeteo_stations});
  met->st->lon.reset(new Vector<double>{nmeteo_stations});
  met->st->Z.reset(new Vector<double>{nmeteo_stations});
  met->st->sky.reset(new Vector<double>{nmeteo_stations});
  met->st->ST.reset(new Vector<double>{nmeteo_stations});
  met->st->Vheight.reset(new Vector<double>{nmeteo_stations});
  met->st->Theight.reset(new Vector<double>{nmeteo_stations});

  i=1;
  met->st->E->co[i] = assignation_number(flog, 201, i-1, keyword, num_param,
                                         num_param_components, (double)number_novalue, 0);
  met->st->N->co[i] = assignation_number(flog, 202, i-1, keyword, num_param,
                                         num_param_components, (double)number_novalue, 0);
  met->st->lat->co[i] = assignation_number(flog, 203, i-1, keyword, num_param,
                                           num_param_components, par->latitude, 0);
  met->st->lon->co[i] = assignation_number(flog, 204, i-1, keyword, num_param,
                                           num_param_components, par->longitude, 0);
  met->st->Z->co[i] = assignation_number(flog, 205, i-1, keyword, num_param,
                                         num_param_components, 0., 0);
  met->st->sky->co[i] = assignation_number(flog, 206, i-1, keyword, num_param,
                                           num_param_components, 1., 0);
  met->st->ST->co[i] = assignation_number(flog, 207, i-1, keyword, num_param,
                                          num_param_components, par->ST, 0);

  for (i=2; i<=nmeteo_stations; i++)
    {
      met->st->E->co[i] = assignation_number(flog, 201, i-1, keyword, num_param,
                                             num_param_components, met->st->E->co[i-1], 0);
      met->st->N->co[i] = assignation_number(flog, 202, i-1, keyword, num_param,
                                             num_param_components, met->st->N->co[i-1], 0);
      met->st->lat->co[i] = assignation_number(flog, 203, i-1, keyword, num_param,
                                               num_param_components, met->st->lat->co[i-1], 0);
      met->st->lon->co[i] = assignation_number(flog, 204, i-1, keyword, num_param,
                                               num_param_components, met->st->lon->co[i-1], 0);
      met->st->Z->co[i] = assignation_number(flog, 205, i-1, keyword, num_param,
                                             num_param_components, met->st->Z->co[i-1], 0);
      met->st->sky->co[i] = assignation_number(flog, 206, i-1, keyword, num_param,
                                               num_param_components, met->st->sky->co[i-1], 0);
      met->st->ST->co[i] = assignation_number(flog, 207, i-1, keyword, num_param,
                                              num_param_components, met->st->ST->co[i-1], 0);
    }

  a = assignation_number(flog, 208, 0, keyword, num_param, num_param_components,
                         10., 0);
  for (i=1; i<=nmeteo_stations; i++)
    {
      met->st->Vheight->co[i] = a;
    }

  a = assignation_number(flog, 209, 0, keyword, num_param, num_param_components,
                         2., 0);
  for (i=1; i<=nmeteo_stations; i++)
    {
      met->st->Theight->co[i] = a;
    }

  //lapse rates (cyclic)
  n = (long)nlstot;
  met->LRcnc = (long *)malloc(n*sizeof(long));
  met->LRcnc[ilsDate12] = 1;
  met->LRcnc[ilsTa] = num_param_components[210];
  met->LRcnc[ilsTdew] = num_param_components[211];
  met->LRcnc[ilsPrec] = num_param_components[212];
  met->LRc = (double **)malloc(n*sizeof(double *));
  for (i=0; i<nlstot; i++)
    {
      met->LRc[i] = (double *)malloc(met->LRcnc[i]*sizeof(double));
      for (j=0; j<met->LRcnc[i]; j++)
        {
          if (i==ilsDate12) met->LRc[i][j] = 0.;
          if (i==ilsTa) met->LRc[i][j] = assignation_number(flog, 210, j, keyword,
                                                              num_param, num_param_components, (double)number_novalue, 0);
          if (i==ilsTdew) met->LRc[i][j] = assignation_number(flog, 211, j, keyword,
                                                                num_param, num_param_components, (double)number_novalue, 0);
          if (i==ilsPrec) met->LRc[i][j] = assignation_number(flog, 212, j, keyword,
                                                                num_param, num_param_components, (double)number_novalue, 0);
        }
    }

  par->MinIncrFactWithElev = assignation_number(flog, 213, j, keyword,
                                                num_param, num_param_components, 0.1, 0);
  par->MaxIncrFactWithElev = assignation_number(flog, 214, j, keyword,
                                                num_param, num_param_components, 4.4, 0);


  //output point column
  n = (long)otot;
  opnt = (long *)malloc(n*sizeof(long));
  ipnt = (short *)malloc(n*sizeof(short));

  par->all_point = (short)assignation_number(flog, 294, 0, keyword, num_param,
                                             num_param_components, 0., 0);

  if (par->all_point == 1)
    {

      for (i=0; i<n; i++)
        {
          ipnt[i] = 1;
          opnt[i] = i;
        }

      nopnt = n;

    }
  else
    {

      for (i=0; i<n; i++)
        {
          ipnt[i] = 0;
          opnt[i] = -1;
        }
      for (i=0; i<n; i++)
        {
          j = (long)assignation_number(flog, 215+i, 0, keyword, num_param,
                                       num_param_components, -1., 0);
          if (j>=1 && j<=n)
            {
              opnt[j-1] = i;
              ipnt[i] = 1;
            }
        }
      nopnt = 0;
      for (i=0; i<n; i++)
        {
          if (opnt[i] > 0) nopnt = i+1;
        }
    }

  //output basin column
  n = (long)ootot;
  obsn = (long *)malloc(n*sizeof(long));
  ibsn = (short *)malloc(n*sizeof(short));

  par->all_basin = (short)assignation_number(flog, 322, 0, keyword, num_param,
                                             num_param_components, 0., 0);

  if (par->all_basin == 1)
    {

      for (i=0; i<n; i++)
        {
          ibsn[i] = 1;
          obsn[i] = i;
        }

      nobsn = (long)ootot;

    }
  else
    {

      for (i=0; i<n; i++)
        {
          ibsn[i] = 0;
          obsn[i] = -1;
        }
      for (i=0; i<n; i++)
        {
          j = (long)assignation_number(flog, 295+i, 0, keyword, num_param,
                                       num_param_components, -1., 0);
          if (j>=1 && j<=n)
            {
              obsn[j-1] = i;
              ibsn[i] = 1;
            }
        }
      nobsn = 0;
      for (i=0; i<n; i++)
        {
          if (obsn[i] > 0) nobsn = i+1;
        }
    }

  cod = 356;
  n = 0;
  do
    {
      a = assignation_number(flog, cod, n, keyword, num_param, num_param_components,
                             (double)number_novalue, 0);
      if ((long)a != number_novalue) n++;
    }
  while ((long)a != number_novalue && n<=1000000);
  if (n==0) n=1;
  par->soil_plot_depths.reset(new Vector<double>{n});
  for (n=1; n<=par->soil_plot_depths->nh; n++)
    {
      par->soil_plot_depths->co[n] = assignation_number(flog, cod, n-1, keyword,
                                                        num_param, num_param_components, (double)number_novalue, 0);
    }

  cod = 357;
  n = 0;
  do
    {
      a = assignation_number(flog, cod, n, keyword, num_param, num_param_components,
                             (double)number_novalue, 0);
      if ((long)a != number_novalue) n++;
    }
  while ((long)a != number_novalue && n<=1000000);
  if (n==0) n=1;
  par->snow_plot_depths.reset(new Vector<double>{n});
  for (n=1; n<=par->snow_plot_depths->nh; n++)
    {
      par->snow_plot_depths->co[n] = assignation_number(flog, cod, n-1, keyword,
                                                        num_param, num_param_components, (double)number_novalue, 0);
    }

  cod = 358;
  n = 0;
  do
    {
      a = assignation_number(flog, cod, n, keyword, num_param, num_param_components,
                             (double)number_novalue, 0);
      if ((long)a != number_novalue) n++;
    }
  while ((long)a != number_novalue && n<=1000000);
  if (n==0) n=1;
  par->glac_plot_depths.reset(new Vector<double>{n});
  for (n=1; n<=par->glac_plot_depths->nh; n++)
    {
      par->glac_plot_depths->co[n] = assignation_number(flog, cod, n-1, keyword,
                                                        num_param, num_param_components, (double)number_novalue, 0);
    }

  //output snow column
  if ((long)par->snow_plot_depths->co[1] != number_novalue)
    {
      m = par->snow_plot_depths->nh;
    }
  else
    {
      m = par->max_snow_layers;
    }
  n = 6;

  osnw = (long *)malloc(n*sizeof(long));

  par->all_snow = (short)assignation_number(flog, 335, 0, keyword, num_param,
                                            num_param_components, 0., 0);

  if (par->all_snow == 1)
    {

      for (i=0; i<n; i++)
        {
          osnw[i] = i;
        }

      nosnw = n;

    }
  else
    {

      cod = 323;

      for (i=0; i<n; i++)
        {
          osnw[i] = -1;
        }
      for (i=0; i<6; i++)
        {
          j = (long)assignation_number(flog, cod+i, 0, keyword, num_param,
                                       num_param_components, -1., 0);
          if (j>=1 && j<=n) osnw[j-1] = i;
        }
      nosnw = 0;
      for (i=0; i<n; i++)
        {
          if (osnw[i] > 0) nosnw = i+1;
        }
    }


  //output glacier column
  if ((long)par->glac_plot_depths->co[1] != number_novalue)
    {
      m = par->glac_plot_depths->nh;
    }
  else
    {
      m = par->max_glac_layers;
    }
  n = 6 + 3*m + 1*par->max_glac_layers;
  oglc = (long *)malloc(n*sizeof(long));

  par->all_glac = (short)assignation_number(flog, 348, 0, keyword, num_param,
                                            num_param_components, 0., 0);

  if (par->all_glac == 1)
    {

      for (i=0; i<n; i++)
        {
          oglc[i] = i;
        }

      noglc = n;

    }
  else
    {

      cod = 335;

      for (i=0; i<n; i++)
        {
          oglc[i] = -1;
        }
      for (i=0; i<6; i++)
        {
          j = (long)assignation_number(flog, cod+i, 0, keyword, num_param,
                                       num_param_components, -1., 0);
          if (j>=1 && j<=n) oglc[j-1] = i;
        }
      for (i=6; i<9; i++)
        {
          for (k=0; k<m; k++)
            {
              j = (long)assignation_number(flog, cod+i, k, keyword, num_param,
                                           num_param_components, -1., 0);
              if (j>=1 && j<=n) oglc[j-1] = (i-6)*m + k + 6;
            }
        }
      for (i=9; i<10; i++)
        {
          for (k=0; k<par->max_glac_layers; k++)
            {
              j = (long)assignation_number(flog, cod+i, k, keyword, num_param,
                                           num_param_components, -1., 0);
              if (j>=1 && j<=n) oglc[j-1] = (i-9)*par->max_glac_layers + k + 6 + 3*m;
            }
        }
      noglc = 0;
      for (i=0; i<n; i++)
        {
          if (oglc[i] > 0) noglc = i+1;
        }
    }

  //output soil column
  n = 6;
  osl = (long *)malloc(n*sizeof(long));

  par->all_soil = (short)assignation_number(flog, 355, 0, keyword, num_param,
                                            num_param_components, 0., 0);

  if (par->all_soil == 1)
    {

      for (i=0; i<n; i++)
        {
          osl[i] = i;
        }

      nosl = n;

    }
  else
    {

      for (i=0; i<n; i++)
        {
          osl[i] = -1;
        }
      for (i=0; i<n; i++)
        {
          j = (long)assignation_number(flog, 349+i, 0, keyword, num_param,
                                       num_param_components, -1., 0);
          if (j>=1 && j<=n) osl[j-1] = i;
        }
      nosl = 0;
      for (i=0; i<n; i++)
        {
          if (osl[i] > 0) nosl = i+1;
        }
    }

  par->ric_cloud = (short)assignation_number(flog, 359, 0, keyword, num_param,
                                             num_param_components, 0., 0);
  par->vap_as_RH = (short)assignation_number(flog, 360, 0, keyword, num_param,
                                             num_param_components, 1., 0);
  par->vap_as_Td = (short)assignation_number(flog, 361, 0, keyword, num_param,
                                             num_param_components, 0., 0);
  par->ndivdaycloud = (long)assignation_number(flog, 362, 0, keyword, num_param,
                                               num_param_components, 3., 0);
  par->cast_shadow = (short)assignation_number(flog, 363, 0, keyword, num_param,
                                               num_param_components, 1., 0);
  par->wind_as_dir = (short)assignation_number(flog, 364, 0, keyword, num_param,
                                               num_param_components, 1., 0);
  par->wind_as_xy = (short)assignation_number(flog, 365, 0, keyword, num_param,
                                              num_param_components, 0., 0);

  par->snow_aging_vis = assignation_number(flog, 366, 0, keyword, num_param,
                                           num_param_components, 0.2, 0);
  par->snow_aging_nir = assignation_number(flog, 367, 0, keyword, num_param,
                                           num_param_components, 0.5, 0);

  par->DepthFreeSurface = assignation_number(flog, 368, 0, keyword, num_param,
                                             num_param_components, 0., 0);

  par->prec_as_intensity = assignation_number(flog, 369, 0, keyword, num_param,
                                              num_param_components, 0., 0);

  cod = 370;
  par->linear_interpolation_meteo = new_shortvector(nmeteo_stations);
  par->linear_interpolation_meteo->co[1] = (short)assignation_number(flog, cod,
                                           0, keyword, num_param, num_param_components, 0., 0);
  for (i=2; i<=nmeteo_stations; i++)
    {
      par->linear_interpolation_meteo->co[i] = (short)assignation_number(flog, cod,
                                               i-1, keyword, num_param, num_param_components,
                                               par->linear_interpolation_meteo->co[i-1], 0);
    }

  cod = 371;
  if (par->point_sim != 1)
    {
      par->output_vertical_distances = (short)assignation_number(flog, cod, 0,
                                                                 keyword, num_param, num_param_components, 0., 0);
      if (par->output_vertical_distances == 1)
        {
          printf("Only for point simulations the parameter %s can be assigned to 1, layers are defined vertically\n",
                 keyword[cod]);
          fprintf(flog,
                  "Only for point simulations the parameter %s can be assigned to 1, layers are defined vertically\n",
                  keyword[cod]);
        }
    }
  else
    {
      par->output_vertical_distances = (short)assignation_number(flog, cod, 0,
                                                                 keyword, num_param, num_param_components, 0., 0);
    }

  par->upwindblowingsnow = (short)assignation_number(flog, 372, 0, keyword,
                                                     num_param, num_param_components, 0., 0);

  par->UpdateK = (short)assignation_number(flog, 373, 0, keyword, num_param,
                                           num_param_components, 0., 0);

  par->ContRecovery = assignation_number(flog, 374, 0, keyword, num_param,
                                         num_param_components, 0., 0);
  par->flag1D = (short)assignation_number(flog, 375, 0, keyword, num_param,
                                          num_param_components, 0., 0);

  par->k_to_ksat = assignation_number(flog, 376, 0, keyword, num_param,
                                      num_param_components, 0., 0);
  par->RunIfAnOldRunIsPresent = (short)assignation_number(flog, 377, 0, keyword,
                                                          num_param, num_param_components, 1., 0);

  par->max_courant_land_channel = assignation_number(flog, 378, 0, keyword,
                                                     num_param, num_param_components, 0.1, 0);
  par->min_dhsup_land_channel_out = assignation_number(flog, 379, 0, keyword,
                                                       num_param, num_param_components, 1., 0);


  cod = 381;
  par->Nl_spinup = new_longvector(par->end_date->nh);
  par->Nl_spinup->co[1] = assignation_number(flog, cod, 0, keyword, num_param,
                                             num_param_components, 10000., 0);
  for (i=2; i<=par->end_date->nh; i++)
    {
      par->Nl_spinup->co[i] = assignation_number(flog, cod, i-1, keyword, num_param,
                                                 num_param_components, par->Nl_spinup->co[i-1], 0);
    }
  if (par->Nl_spinup->co[1]<10000. && par->point_sim!=1)
    {
      printf("You can use %s only if %s is set to 1\n",keyword[cod],keyword[12]);
      fprintf(flog,"You can use %s only if %s is set to 1\n",keyword[cod],
              keyword[12]);
      t_error("Not possible to continue");
    }

  par->newperiodinit = (short)assignation_number(flog, 382, 0, keyword,
                                                 num_param, num_param_components, 0., 0);
  if (par->newperiodinit != 0 && par->point_sim != 1)
    {
      printf("You can use %s only if %s is set to 1\n",keyword[382],keyword[12]);
      fprintf(flog,"You can use %s only if %s is set to 1\n",keyword[382],
              keyword[12]);
      t_error("Not possible to continue");
    }

  par->k1 = assignation_number(flog, 383, 0, keyword, num_param,
                               num_param_components, 0.484, 0);
  par->k2 = assignation_number(flog, 384, 0, keyword, num_param,
                               num_param_components, 8., 0);
  par->Lozone = assignation_number(flog, 385, 0, keyword, num_param,
                                   num_param_components, 0.3, 0);
  par->alpha_iqbal = assignation_number(flog, 386, 0, keyword, num_param,
                                        num_param_components, 1.3, 0);
  par->beta_iqbal = assignation_number(flog, 387, 0, keyword, num_param,
                                       num_param_components, 0.1, 0);

  par->albedoSWin = (short)assignation_number(flog, 388, 0, keyword, num_param,
                                              num_param_components, 0., 0);

  par->micro = (short)assignation_number(flog, 389, 0, keyword, num_param,
                                         num_param_components, 1., 0);
  par->EB = assignation_number(flog, 390, 0, keyword, num_param,
                               num_param_components, (double)number_novalue, 0);
  par->Cair = assignation_number(flog, 391, 0, keyword, num_param,
                                 num_param_components, (double)number_novalue, 0);
  par->Tsup = assignation_number(flog, 392, 0, keyword, num_param,
                                 num_param_components, (double)number_novalue, 0);

  par->Tair_default = assignation_number(flog, 393, 0, keyword, num_param,
                                         num_param_components, 5., 0);
  par->RH_default = assignation_number(flog, 394, 0, keyword, num_param,
                                       num_param_components, 70., 0)/100.;
  par->V_default = assignation_number(flog, 395, 0, keyword, num_param,
                                      num_param_components, par->Vmin, 0);
  par->Vdir_default = assignation_number(flog, 396, 0, keyword, num_param,
                                         num_param_components, 0., 0);
  par->IPrec_default = assignation_number(flog, 397, 0, keyword, num_param,
                                          num_param_components, 0., 0);

  par->soil_type_bedr_default = (long)assignation_number(flog, 399, 0, keyword,
                                                         num_param, num_param_components, 1., 0);
  if (par->soil_type_bedr_default<1
      || par->soil_type_bedr_default>par->nsoiltypes)
    {
      fprintf(flog,
              "Error:  soil_type_bedr_default lower than 0 or higher than soil types numbers");
      printf("Error:  soil_type_bedr_default lower than 0 or higher than soil types numbers");
      t_error("Fatal Error! Geotop is closed. See failing report.");
    }

  par->minP_torestore_A = assignation_number(flog, 400, 0, keyword, num_param,
                                             num_param_components, 10., 0);
  par->snow_conductivity = (short)assignation_number(flog, 401, 0, keyword,
                                                     num_param, num_param_components, 1., 0);
  par->snow_wind_compaction_1D = (short)assignation_number(flog, 402, 0,
                                                           keyword, num_param, num_param_components, 0., 0);
  if (par->snow_wind_compaction_1D == 1) par->blowing_snow = 1;

  par->DDchannel = (short)assignation_number(flog, 403, 0, keyword, num_param,
                                             num_param_components, 1., 0);
  par->DDland = (short)assignation_number(flog, 404, 0, keyword, num_param,
                                          num_param_components, 1., 0);
  par->Tbottom = assignation_number(flog, 405, 0, keyword, num_param,
                                    num_param_components, (double)number_novalue, 0);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

char **assign_string_parameter(FILE *f, long beg, long end,
                               char **string_param, char **keyword)
{

  long i;
  char **a;

  a = (char **)malloc((end-beg)*sizeof(char *));

  for (i=0; i<end-beg; i++)
    {
      a[i] = assignation_string(f, i+beg, keyword, string_param);
    }

  return (a);

}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

double assignation_number(FILE *flog, long i, long j, char **keyword,
                          double **num_param, long *num_param_components, double default_value,
                          short code_error)
{

  double a;

  if (j < num_param_components[i])
    {
      a = num_param[i][j];
    }
  else
    {
      a = (double)number_novalue;
    }

  if ((long)a == number_novalue)
    {
      if (code_error==0)
        {
          a = default_value;
          fprintf(flog,"%s[%ld] = %e (default) \n", keyword[i], j+1, a);
        }
      else
        {
          fprintf(flog, "%s[%ld] not assigned\n", keyword[i], j+1);
          fclose(flog);
          t_error("Fatal Error, See geotop.log!");
        }
    }
  else
    {
      fprintf(flog,"%s[%ld] = %e \n", keyword[i], j+1, a);
    }

  return a;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

char *assignation_string(FILE *f, long i, char **keyword, char **string_param)
{

  char *a;
  long j, dimstring = strlen(string_param[i]);

  a = (char *)malloc((dimstring+1)*sizeof(char));

  for (j=0; j<dimstring; j++)
    {
      a[j] = string_param[i][j];
    }
  a[dimstring] = 0;

  fprintf(f,"%s = %s\n", keyword[i], a);

  return (a);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

short read_soil_parameters(char *name, INIT_TOOLS *IT, SOIL *sl, long bed,
                           FILE *flog)
{

  short ok;
  long i, j, k, n, nlines, nlinesprev;
  char *temp;
  double **soildata;
  DOUBLETENSOR *old_sl_par;
  FILE *f;

  //look if there is at least 1 soil file
  i = 0;
  ok = 0;
  nlinesprev = -1;

  do
    {

      temp = namefile_i_we2(name, i+1);

      if (existing_file_wext(temp, textfile)==1)
        {
          free(temp);
          ok = 1;
          temp = namefile_i(name, i+1);
          nlines = count_lines(temp, 33, 44);
          if (nlinesprev >= 0 && nlines != nlinesprev)
            {
              f = fopen(FailedRunFile, "w");
              fprintf(f,
                      "Error:: The file %s with soil paramaters has a number of layers %ld, which different from the numbers %ld of the other soil parameter files\n",
                      temp,nlines,nlinesprev);
              fprintf(f,
                      "In Geotop it is only possible to have the same number of layers in any soil parameter files\n");
              fclose(f);
              t_error("Fatal Error! Geotop is closed. See failing report.");
            }
          nlinesprev = nlines;
        }
      else
        {
          if (i==0 && strcmp(name, string_novalue) != 0)
            {
              f = fopen(FailedRunFile, "w");
              fprintf(f,"Error:: Soil file %s not existing.\n",name);
              fclose(f);
              t_error("Fatal Error! Geotop is closed. See failing report.");
            }
        }


      free(temp);
      i++;

    }
  while (ok == 0 && i < sl->pa->ndh);

  if (ok == 1)
    {

      //save sl->pa in a new doubletensor and deallocate
      old_sl_par = new_doubletensor(sl->pa->ndh, sl->pa->nrh, sl->pa->nch);
      for (i=1; i<=sl->pa->ndh; i++)
        {
          for (n=1; n<=sl->pa->nrh; n++)
            {
              for (j=1; j<=sl->pa->nch; j++)
                {
                  old_sl_par->co[i][n][j] = sl->pa->co[i][n][j];
                }
            }
        }
      free_doubletensor(sl->pa);

      //reallocate
      sl->pa = new_doubletensor(old_sl_par->ndh, old_sl_par->nrh, nlines);

      for (i=1; i<=sl->pa->ndh; i++)
        {

          //read files
          temp = namefile_i_we2(name, i);

          if (existing_file_wext(temp, textfile)==1)
            {
              free(temp);
              temp = namefile_i(name, i);
              soildata = read_txt_matrix(temp, 33, 44, IT->soil_col_names, nsoilprop,
                                         &nlines, flog);
            }
          else
            {
              soildata = (double **)malloc(nlines*sizeof(double *));
              for (j=0; j<nlines; j++)
                {
                  k = (long)nsoilprop;
                  soildata[j] = (double *)malloc(k*sizeof(double));
                  for (n=0; n<k; n++)
                    {
                      soildata[j][n] = (double)number_absent;
                    }
                }
            }

          free(temp);

          //assign soildata to soil->pa
          for (n=1; n<=nsoilprop; n++)
            {
              for (j=1; j<=sl->pa->nch; j++)   //j is the layer index
                {
                  sl->pa->co[i][n][j] = soildata[j-1][n-1];
                }
            }

          //deallocate soildata
          for (j=0; j<nlines; j++)
            {
              free(soildata[j]);
            }
          free(soildata);

          //fix layer thickness
          n = jdz;
          for (j=1; j<=sl->pa->nch; j++)   //j is the layer index
            {
              if ((long)sl->pa->co[i][n][j] != number_novalue
                  && (long)sl->pa->co[i][n][j] != number_absent)
                {
                  if ( i > 1 && fabs( sl->pa->co[i][n][j] - sl->pa->co[i-1][n][j] ) > 1.E-5 )
                    {
                      f = fopen(FailedRunFile, "w");
                      fprintf(f,
                              "Error:: For soil type %ld it has been given a set of soil layer thicknesses different from the other ones.\n",
                              i);
                      fprintf(f,
                              "In Geotop it is only possible to have the soil layer discretization in any soil parameter files.\n");
                      fclose(f);
                      t_error("Fatal Error! Geotop is closed. See failing report.");
                    }
                }
              else if (i == 1)
                {
                  if (j <= old_sl_par->nch)
                    {
                      sl->pa->co[i][n][j] = old_sl_par->co[i][n][j];
                    }
                  else
                    {
                      sl->pa->co[i][n][j] = sl->pa->co[i][n][j-1];
                    }
                }
              else
                {
                  sl->pa->co[i][n][j] = sl->pa->co[i-1][n][j];
                }
            }

          //all other variables
          for (n=1; n<=nsoilprop; n++)
            {
              if (n != jdz)
                {
                  for (j=1; j<=sl->pa->nch; j++)   //j is the layer index
                    {
                      if ((long)sl->pa->co[i][n][j] == number_novalue
                          || (long)sl->pa->co[i][n][j] == number_absent)
                        {
                          if (j <= old_sl_par->nch)
                            {
                              sl->pa->co[i][n][j] = old_sl_par->co[i][n][j];
                            }
                          else
                            {
                              sl->pa->co[i][n][j] = sl->pa->co[i][n][j-1];
                            }
                        }
                    }
                }
            }

          //field capacity and wilting point
          for (j=1; j<=sl->pa->nch; j++)
            {
              if ( (long)sl->pa->co[i][jfc][j] == number_novalue)
                {
                  sl->pa->co[i][jfc][j] = teta_psi( (-1./3.)*1.E5/g, 0., sl->pa->co[i][jsat][j],
                                                    sl->pa->co[i][jres][j], sl->pa->co[i][ja][j],
                                                    sl->pa->co[i][jns][j], 1.-1./sl->pa->co[i][jns][j], PsiMin,
                                                    sl->pa->co[i][jss][j]);
                }

              if ( (long)sl->pa->co[i][jwp][j] == number_novalue)
                {
                  sl->pa->co[i][jwp][j] = teta_psi( -15.*1.E5/g, 0., sl->pa->co[i][jsat][j],
                                                    sl->pa->co[i][jres][j], sl->pa->co[i][ja][j],
                                                    sl->pa->co[i][jns][j], 1.-1./sl->pa->co[i][jns][j], PsiMin,
                                                    sl->pa->co[i][jss][j]);
                }
            }

          //pressure
          ok = 1;
          for (j=1; j<=sl->pa->nch; j++)
            {
              if ( (long)sl->pa->co[i][jpsi][j] == number_novalue) ok = 0;
            }
          if (ok == 1) IT->init_water_table_depth->co[i] = (double)number_novalue;
        }

      free_doubletensor(old_sl_par);

    }

  //write on the screen the soil paramater
  fprintf(flog,"\n");
  k = (long)nmet;
  fprintf(flog,"Soil Layers: %ld\n",sl->pa->nch);
  for (i=1; i<=sl->pa->ndh; i++)
    {
      fprintf(flog,"-> Soil Type: %ld\n",i);
      for (n=1; n<=nsoilprop; n++)
        {
          fprintf(flog,"%s: ",keywords_char[k+n-1]);
          for (j=1; j<=sl->pa->nch; j++)
            {
              fprintf(flog,"%f(%.2e)",sl->pa->co[i][n][j],sl->pa->co[i][n][j]);
              if (j<sl->pa->nch)fprintf(flog,", ");
            }
          fprintf(flog,"\n");
        }
    }

  //bedrock
  old_sl_par = new_doubletensor(1, IT->pa_bed->nrh, IT->pa_bed->nch);
  for (n=1; n<=IT->pa_bed->nrh; n++)
    {
      for (j=1; j<=IT->pa_bed->nch; j++)
        {
          old_sl_par->co[1][n][j] = IT->pa_bed->co[1][n][j];
        }
    }
  free_doubletensor(IT->pa_bed);

  IT->pa_bed = new_doubletensor(sl->pa->ndh, sl->pa->nrh, sl->pa->nch);
  for (i=1; i<=IT->pa_bed->ndh; i++)
    {
      for (n=1; n<=IT->pa_bed->nrh; n++)
        {
          if (n == jdz)
            {
              for (j=1; j<=IT->pa_bed->nch; j++)
                {
                  IT->pa_bed->co[i][n][j] = sl->pa->co[1][n][j];
                }
            }
          else
            {
              for (j=1; j<=IT->pa_bed->nch; j++)
                {
                  if (j <= old_sl_par->nch)
                    {
                      IT->pa_bed->co[i][n][j] = old_sl_par->co[1][n][j];
                    }
                  else
                    {
                      IT->pa_bed->co[i][n][j] = IT->pa_bed->co[i][n][j-1];
                    }
                }
              for (j=1; j<=IT->pa_bed->nch; j++)
                {
                  if ( (long)IT->pa_bed->co[i][n][j] == number_novalue ) IT->pa_bed->co[i][n][j]
                      = sl->pa->co[bed][n][j];
                }
            }
        }
    }
  free_doubletensor(old_sl_par);

  fprintf(flog,"\n");
  k = (long)nmet;
  fprintf(flog,"Soil Bedrock Layers: %ld\n",sl->pa->nch);
  for (i=1; i<=IT->pa_bed->ndh; i++)
    {
      fprintf(flog,"-> Soil Type: %ld\n",i);
      for (n=1; n<=nsoilprop; n++)
        {
          fprintf(flog,"%s: ",keywords_char[k+n-1]);
          for (j=1; j<=sl->pa->nch; j++)
            {
              fprintf(flog,"%f(%.2e)",IT->pa_bed->co[i][n][j],IT->pa_bed->co[i][n][j]);
              if (j<sl->pa->nch)fprintf(flog,", ");
            }
          fprintf(flog,"\n");
        }
    }
  fprintf(flog,"\n");

  return 1;

}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

short read_point_file(char *name, char **key_header, PAR *par, FILE *flog)
{

  DOUBLEMATRIX *chkpt2;
  double **points;
  long nlines, n, j;
  char *temp;

  if (existing_file_wext(name, textfile)==1)
    {
      temp = join_strings(name, textfile);
      printf("%s\n",temp);
      points = read_txt_matrix(temp, 33, 44, key_header, par->chkpt->nch, &nlines,
                               flog);
      free(temp);

      chkpt2 = new_doublematrix(par->chkpt->nrh, par->chkpt->nch);
      copy_doublematrix(par->chkpt, chkpt2);
      free_doublematrix(par->chkpt);

      par->chkpt = new_doublematrix(nlines, chkpt2->nch);

      for (n=1; n<=nlines; n++)
        {
          for (j=1; j<=chkpt2->nch; j++)
            {
              par->chkpt->co[n][j] = points[n-1][j-1];
              if ( (long)par->chkpt->co[n][j] == number_novalue
                   || (long)par->chkpt->co[n][j] == number_absent )
                {
                  if ( n <= chkpt2->nrh )
                    {
                      par->chkpt->co[n][j] = chkpt2->co[n][j];
                    }
                  else
                    {
                      par->chkpt->co[n][j] = chkpt2->co[chkpt2->nrh][j];
                    }
                }
            }

          if ( par->point_sim != 1 )
            {
              if ( (long)par->chkpt->co[n][ptX] == number_novalue
                   || (long)par->chkpt->co[n][ptY] == number_novalue )
                {
                  par->state_pixel = 0;
                }
            }

          free(points[n-1]);
        }

      free_doublematrix(chkpt2);
      free(points);
    }

  if (par->point_sim != 1)
    {
      for (n=1; n<=par->chkpt->nrh; n++)
        {
          if ( (long)par->chkpt->co[n][ptX] == number_novalue
               || (long)par->chkpt->co[n][ptY] == number_novalue)
            {
              fprintf(flog,
                      "Warning: The points to plot specific results are not completely specified\n");
              fprintf(flog,"Output for single point output is deactivated.\n");
              printf("Warning: The points to plot specific results are not completely specified\n");
              printf("Output for single point output is deactivated.\n");
              par->state_pixel = 0;
            }
        }
    }


  return 1;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

short read_meteostations_file(LONGVECTOR *i, METEO_STATIONS *S, char *name,
                              char **key_header, FILE *flog)
{

  double **M;
  long nlines, n, j, k;
  char *temp;

  if (existing_file_wext(name, textfile)==1)
    {
      temp = join_strings(name, textfile);
      M = read_txt_matrix(temp, 33, 44, key_header, 8, &nlines, flog);
      free(temp);

      for (j=1; j<=i->nh; j++)
        {
          for (n=1; n<=nlines; n++)
            {
              if ((long)M[n-1][0] == i->co[j])
                {
                  for (k=1; k<8; k++)
                    {
                      if ((long)M[n-1][k] != number_novalue && (long)M[n-1][k] != number_absent)
                        {
                          if (k==1)
                            {
                              S->E->co[j] = M[n-1][k];
                            }
                          else if (k==2)
                            {
                              S->N->co[j] = M[n-1][k];
                            }
                          else if (k==3)
                            {
                              S->lat->co[j] = M[n-1][k];
                            }
                          else if (k==4)
                            {
                              S->lon->co[j] = M[n-1][k];
                            }
                          else if (k==5)
                            {
                              S->Z->co[j] = M[n-1][k];
                            }
                          else if (k==6)
                            {
                              S->sky->co[j] = M[n-1][k];
                            }
                          else if (k==7)
                            {
                              S->ST->co[j] = M[n-1][k];
                            }
                        }
                    }
                }
            }
        }

      for (n=1; n<=nlines; n++)
        {
          free(M[n-1]);
        }

      free(M);

      return 1;
    }
  else
    {
      return 0;
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
