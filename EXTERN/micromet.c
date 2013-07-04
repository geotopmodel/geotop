/*

MICROMET CODE

Code written by Stefano Endrizzi by translating and adapting the idea behind the Micromet Fortran Code
by G. Liston and X. The author does not guarantee the perfect compliance of this code with the Fortran one.
However, he asks to give credit to Liston and  X, JHM YYYYY, 2006, when using it with satisfaction.

*/



#include "constant.h"
#include "struct.geotop.09375.h"
#include "micromet.h"

extern T_INIT *UV;

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void Micromet(T_INIT *UV, DOUBLEMATRIX *topo, DOUBLEMATRIX *curvature, DOUBLEMATRIX *terrain_slope, DOUBLEMATRIX *slope_az, METEO *met,
	double slopewt, double curvewt, double windspd_min, double dn, short ifill, short iobsint, long Tcode, long RHcode, long Vcode, long Vdircode,
	long Pcode, DOUBLEMATRIX *Tair_grid, DOUBLEMATRIX *RH_grid, DOUBLEMATRIX *windspd_grid, DOUBLEMATRIX *winddir_grid, DOUBLEMATRIX *sfc_pressure,
	DOUBLEMATRIX *prec_grid, double T_lapse_rate, double Td_lapse_rate, double Prec_lapse_rate){

	get_temperature(UV, met, Tcode, Tair_grid, dn, topo, ifill, iobsint, T_lapse_rate);
	get_relative_humidity(UV, met, RHcode, Tcode, RH_grid, Tair_grid, dn, topo, ifill, iobsint, Td_lapse_rate);
	get_wind(UV, met, Vcode, Vdircode, windspd_grid, winddir_grid, curvature, slope_az, terrain_slope, slopewt, curvewt, windspd_min, dn, topo, ifill, iobsint);
	get_precipitation(UV, met, Pcode, prec_grid, dn, topo, ifill, iobsint, Prec_lapse_rate);
	get_pressure(topo, sfc_pressure, UV->V->co[2]);

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************


void get_temperature(T_INIT *UV, METEO *met, long Tcode, DOUBLEMATRIX *Tair_grid, double dn, DOUBLEMATRIX *topo, short ifill, short iobsint,
	double T_lapse_rate){

	double topo_ref, delta_topo, novalue=UV->V->co[2];
	long n, col;
	long r, c;

	//Define the topographic reference surface.
	topo_ref = 0.0;

	//Convert the station data to sea level values in [K].
	for(n=1;n<=met->st->Z->nh;n++){
		delta_topo = topo_ref - met->st->Z->co[n];
		col=met->column[n-1][Tcode];
        if(col!=-1){
			if(met->var[n-1][col]!=novalue) met->var[n-1][col] += (1.E-3*T_lapse_rate * delta_topo + tk);
		}
	}

	//Use the barnes oi scheme to interpolate the station data to the grid.
	interpolate_meteo(UV, met->st, met->var, met->column, Tcode, Tair_grid, dn, ifill, iobsint);

	//Convert these grid values back to the actual gridded elevations.
    for(r=1;r<=Tair_grid->nrh;r++){
		for(c=1;c<=Tair_grid->nch;c++){
			if(topo->co[r][c]!=novalue){
				delta_topo = topo_ref - topo->co[r][c];
				Tair_grid->co[r][c] -= (1.E-3*T_lapse_rate * delta_topo + tk);
			}
		}
	}

	//Convert the station data to sea level values in [K].
	for(n=1;n<=met->st->Z->nh;n++){
		delta_topo = topo_ref - met->st->Z->co[n];
		col=met->column[n-1][Tcode];
        if(col!=-1){
			if(met->var[n-1][col]!=novalue) met->var[n-1][col] -= (1.E-3*T_lapse_rate * delta_topo + tk);
		}
	}

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void get_relative_humidity(T_INIT *UV, METEO *met, long RHcode, long Tcode, DOUBLEMATRIX *RH_grid, DOUBLEMATRIX *Tair_grid, double dn,
	DOUBLEMATRIX *topo, short ifill, short iobsint, double Td_lapse_rate){

// First convert stn relative humidity to dew-point temperature.  Use
//   the Td lapse rate to take the stn Td to sea level.  Interpolate
//   the stn Td to the grid.  Use the Td lapse rate to take the sea
//   level grid to the actual elevations.  Convert each Td to
//   relative humidity.

	double topo_ref=0; /* corrected by Emanuele Cordano on 24/9/9 */
	double delta_topo, novalue=UV->V->co[2];
	long n, colT, colRH;
	long r, c;

	// Convert the stn relative humidity to Td.
	for(n=1;n<=met->st->Z->nh;n++){
		colT=met->column[n-1][Tcode];
		colRH=met->column[n-1][RHcode];
        if(colRH!=-1){
			if(colT==-1){
				printf("RH data must be together with T data at the same station, error in station %ld\n",n);
				t_error("Correct the data and run again the code");
			}
			if(met->var[n-1][colRH]!=novalue){
				if(met->var[n-1][colT]==novalue){
					printf("RH data must be together with T data at the same station, error in station %ld\n",n);
					t_error("Correct the data and run again the code");
				}

				//convert RH in Tdew
				met->var[n-1][colRH]=Tdew(met->var[n-1][colT], Fmax(10.0,met->var[n-1][colRH])/100.0, met->st->Z->co[n]);
			}

			// Define the topographic reference surface.
			topo_ref = 0.0;

			// Convert the station data to sea level values.
			delta_topo = topo_ref - met->st->Z->co[n];
			met->var[n-1][colRH] += (1.E-3*Td_lapse_rate * delta_topo + tk);
		}
	}

	// Use the barnes oi scheme to interpolate the station data to the grid.
	interpolate_meteo(UV, met->st, met->var, met->column, RHcode, RH_grid, dn, ifill, iobsint);

	//Convert these grid values back to the actual gridded elevations.
    for(r=1;r<=RH_grid->nrh;r++){
		for(c=1;c<=RH_grid->nch;c++){
			if(topo->co[r][c]!=novalue){
				delta_topo = topo_ref - topo->co[r][c];
				RH_grid->co[r][c] -= (1.E-3*Td_lapse_rate * delta_topo + tk);
			}
		}
	}

	//Convert each Td to a gridded relative humidity.
    for(r=1;r<=RH_grid->nrh;r++){
		for(c=1;c<=RH_grid->nch;c++){
			if(topo->co[r][c]!=novalue){
				RH_grid->co[r][c] = RHfromTdew(Tair_grid->co[r][c], RH_grid->co[r][c], topo->co[r][c]);
          	}
        }
	}
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void topo_mod_winds(DOUBLEMATRIX *winddir_grid, DOUBLEMATRIX *windspd_grid, double slopewt, double curvewt, DOUBLEMATRIX *curvature,
	DOUBLEMATRIX *slope_az, DOUBLEMATRIX *terrain_slope, DOUBLEMATRIX *topo, double undef){

	long r, c, nc=topo->nch, nr=topo->nrh;
	double deg2rad=Pi/180.0,rad2deg=180.0/Pi,dirdiff,wslope_max,windwt;
	DOUBLEMATRIX *wind_slope;

	//Compute the wind modification factor which is a function of topography and wind direction following Liston and Sturm (1998).

	//Compute the slope in the direction of the wind.
	wind_slope=new_doublematrix(topo->nrh, topo->nch);
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=undef){
				wind_slope->co[r][c] = deg2rad * terrain_slope->co[r][c] * cos(deg2rad * (winddir_grid->co[r][c] - slope_az->co[r][c]));
			}
		}
	}

	//Scale the wind slope such that the max abs(wind slope) has a value
	//of abs(0.5).  Include a 1 mm slope in slope_max to prevent
	//divisions by zero in flat terrain where the slope is zero.
    wslope_max = 0.0 + 0.001;
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=undef){
				wslope_max = Fmax(wslope_max,fabs(wind_slope->co[r][c]));
			}
		}
	}
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=undef){
				wind_slope->co[r][c] /= (2.0 * wslope_max);
			}
		}
	}


	//Calculate the wind speed and direction adjustments.  The
	//curvature and wind_slope values range between -0.5 and +0.5.
	//Valid slopewt and curvewt values are between 0 and 1, with
	//values of 0.5 giving approximately equal weight to slope and
	//curvature.  I suggest that slopewt and curvewt be set such
	//that slopewt + curvewt = 1.0.  This will limit the total
	//wind weight to between 0.5 and 1.5 (but this is not required).
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=undef){

				//Compute the wind weighting factor.
				windwt = 1.0 + slopewt * wind_slope->co[r][c] + curvewt * curvature->co[r][c];

				//Generate the terrain-modified wind speed.
				windspd_grid->co[r][c] *= windwt;

				//Modify the wind direction according to Ryan (1977).  Note that it
				//is critical that "dirdiff" handles the cases where the slope
				//azimuth and the wind direction are on different sides of the
				//360-0 line.
				if (slope_az->co[r][c]>270.0 && winddir_grid->co[r][c]<90.0){
					dirdiff = slope_az->co[r][c] - winddir_grid->co[r][c] - 360.0;
				}else if(slope_az->co[r][c]<90.0 && winddir_grid->co[r][c]>270.0){
					dirdiff = slope_az->co[r][c] - winddir_grid->co[r][c] + 360.0;
				}else{
					dirdiff = slope_az->co[r][c] - winddir_grid->co[r][c];
				}
				if(fabs(dirdiff)<90.0){
					winddir_grid->co[r][c] = winddir_grid->co[r][c] - 0.5 * Fmin(wind_slope->co[r][c]*rad2deg,45.0) * sin(deg2rad * (2.0 * dirdiff));
					if(winddir_grid->co[r][c]>360.0){
						winddir_grid->co[r][c] = winddir_grid->co[r][c] - 360.0;
					}else if (winddir_grid->co[r][c]<0.0){
						winddir_grid->co[r][c] = winddir_grid->co[r][c] + 360.0;
					}
				}
			}
		}
	}

	free_doublematrix(wind_slope);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
void get_wind(T_INIT *UV, METEO *met, long Vcode, long Vdircode, DOUBLEMATRIX *windspd_grid, DOUBLEMATRIX *winddir_grid, DOUBLEMATRIX *curvature, DOUBLEMATRIX *slope_az, DOUBLEMATRIX *terrain_slope,
	double slopewt, double curvewt, double windspd_min, double dn, DOUBLEMATRIX *topo, short ifill, short iobsint){

// This program takes the station wind speed and direction, converts
//   them to u and v components, interpolates u and v to a grid,
//   converts the gridded values to speed and direction, and then
//   runs a simple wind model that adjusts those speeds and
//   directions according to topographic slope and curvature
//   relationships.  The resulting speeds and directions are
//   converted to u and v components and passed back to the main
//   program to be written to a file.  (All of the conversion
//   between u-v and speed-dir is done because of the problems
//   with interpolating over the 360/0 direction line.)

	DOUBLEMATRIX *u_grid, *v_grid;
	long ucode, vcode;
	long r, c, n, colV, colVdir, colu, colv, nc=topo->nch, nr=topo->nrh;
	double speed, u, v, novalue=UV->V->co[2];
	double deg2rad=Pi/180.0,rad2deg=180.0/Pi;

//	Filter through the original input data, and eliminate any
//  missing values (r.e., make sure each wind direction is paired
//   up with a wind speed.
	for(n=1;n<=met->st->Z->nh;n++){
		colV=met->column[n-1][Vcode];
		colVdir=met->column[n-1][Vdircode];
        if(colV!=-1){
			if(colVdir==-1){
				printf("Wind speed data must be together with Wind direction data at the same station, error in station %ld\n",n);
				t_error("Correct the data and run again the code");
			}
			if(met->var[n-1][colV]!=novalue){
				if(met->var[n-1][colVdir]==novalue){
					printf("Wind speed data must be together with Wind direction data at the same station, error in station %ld\n",n);
					t_error("Correct the data and run again the code");
				}
			}
		}
        if(colVdir!=-1){
			if(colV==-1){
				printf("Wind speed data must be together with Wind direction data at the same station, error in station %ld\n",n);
				t_error("Correct the data and run again the code");
			}
			if(met->var[n-1][colVdir]!=novalue){
				if(met->var[n-1][colV]==novalue){
					printf("Wind speed data must be together with Wind direction data at the same station, error in station %ld\n",n);
					t_error("Correct the data and run again the code");
				}
			}
		}
	}

	//Convert these station data to u and v wind components.
	ucode=Vcode;
	vcode=Vdircode;	//used to replace the meteo columns (V,Vdir) with (u,v)
	for(n=1;n<=met->st->Z->nh;n++){
		colV=met->column[n-1][Vcode];
		colVdir=met->column[n-1][Vdircode];
		if(colV!=-1){
			if(met->var[n-1][colV]!=novalue){
				speed = Fmax(windspd_min, met->var[n-1][colV]);
				u = (-speed) * sin(deg2rad * met->var[n-1][colVdir]);
				v = (-speed) * cos(deg2rad * met->var[n-1][colVdir]);
				colu=met->column[n-1][ucode];
				colv=met->column[n-1][vcode];
				met->var[n-1][colu]=u;
				met->var[n-1][colv]=v;
			}
		}
	}

	//Use the barnes oi scheme to interpolate the station data to
	//the grid.
	//U component.
	u_grid=new_doublematrix(nr,nc);
	interpolate_meteo(UV, met->st, met->var, met->column, ucode, u_grid, dn, ifill, iobsint);

	//V component.
	v_grid=new_doublematrix(nr,nc);
	interpolate_meteo(UV, met->st, met->var, met->column, vcode, v_grid, dn, ifill, iobsint);

	//Convert these u and v components to speed and directions.
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=novalue){
				winddir_grid->co[r][c] = 270.0 - rad2deg*atan2(v_grid->co[r][c],u_grid->co[r][c]);
				if (winddir_grid->co[r][c]>=360.0) winddir_grid->co[r][c] -= 360.0;
				windspd_grid->co[r][c] = pow(pow(u_grid->co[r][c], 2.0) + pow(v_grid->co[r][c], 2.0), 0.5);
			}
		}
	}

	free_doublematrix(u_grid);
	free_doublematrix(v_grid);

	//Modify the wind speed and direction according to simple wind-topography relationships.
	topo_mod_winds(winddir_grid, windspd_grid, slopewt, curvewt, curvature,  slope_az,  terrain_slope, topo, novalue);

	//Avoid problems of zero (low) winds (for example, turbulence
	//theory, log wind profile, etc., says that we must have some
	//wind.  Thus, some equations blow up when the wind speed gets
	//very small).
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=novalue){
				if (windspd_grid->co[r][c]<windspd_min) windspd_grid->co[r][c] = windspd_min;
			}
		}
	}
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void get_precipitation(T_INIT *UV, METEO *met, long Pcode, DOUBLEMATRIX *prec_grid, double dn, DOUBLEMATRIX *topo, short ifill, short iobsint,
	double Prec_lapse_rate){

// Interpolate the observed precipitation values to the grid.  Also
//   interpolate the station elevations to a reference surface.  Use
//   a precipitation "lapse rate", or adjustment factor to define
//   the precipitation on the actual elevation grid.  The reason the
//   interpolated station elevations are used as the topographic
//   reference surface (instead of something like sea level), is
//   because the precipitation adjustment factor is a non-linear
//   function of elevation difference.

// The adjustment factor that is used comes from: Thornton, P. E.,
//   S. W. Running, and M. A. White, 1997: Generating surfaces of
//   daily meteorological variables over large regions of complex
//   terrain.  J. Hydrology, 190, 214-251.

	long Zcode, n, col, r, c, nc=topo->nch, nr=topo->nrh;
	double delta_topo, alfa, novalue=UV->V->co[2];
	DOUBLEMATRIX *topo_ref_grid;

// Use the barnes oi scheme to interpolate the station data to
//   the grid.
	interpolate_meteo(UV, met->st, met->var, met->column, Pcode, prec_grid, dn, ifill, iobsint);

// Use the barnes oi scheme to interpolate the station elevation data
//   to the grid, so that it can be used as a topographic reference
//   surface.
	topo_ref_grid=new_doublematrix(nr,nc);
 	Zcode=Pcode; //used to replace the meteo column (P) with (Z)
	for(n=1;n<=met->st->Z->nh;n++){
		col=met->column[n-1][Zcode];
		if(col!=-1){
			if(met->var[n-1][col]!=novalue) met->var[n-1][col]=met->st->Z->co[n];
		}
	}
	interpolate_meteo(UV, met->st, met->var, met->column, Zcode, topo_ref_grid, dn, ifill, iobsint);

//	Convert the gridded station data to the actual gridded elevations.
	for(c=1;c<=nc;c++){
		for(r=1;r<=nr;r++){
			if(topo->co[r][c]!=novalue){
				delta_topo = topo_ref_grid->co[r][c] - topo->co[r][c];
				//Dont let the elevation difference be greater than some number
				//(like 1800 meters gives a factor of 4.4).  If it is too large
				//you get huge precipitation adjustment, a divide by zero, or
				//even negative adjustments for high elevations).
				delta_topo = Fmax(delta_topo,-1800.0);
				alfa = -1.E-3*Prec_lapse_rate * delta_topo;
				prec_grid->co[r][c] = Fmax(prec_grid->co[r][c] * (1.0 + alfa)/(1.0 - alfa), 0.0);
			}
		}
	}

	free_doublematrix(topo_ref_grid);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void get_pressure(DOUBLEMATRIX *topo, DOUBLEMATRIX *sfc_pressure, double undef){

	long r,c,nc=topo->nch,nr=topo->nrh;
	double one_atmos,scale_ht;

	one_atmos = 1013.25;
	scale_ht = 8500.0;

	//Compute the average station pressure (in bar).
	for(c=1;c<=nc;c++){
		for(r=1;r<=nr;r++){
			if(topo->co[r][c]!=undef){
				sfc_pressure->co[r][c] = one_atmos * exp((- topo->co[r][c])/scale_ht);
			}
		}
	}
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

double find_cloudfactor(double Tair, double RH, double Z, double T_lapse_rate, double Td_lapse_rate){

	double fcloud, Z_ref, press_ratio, f_max, one_minus_RHe, f_1, Td, delta_Z, Td_700, Tair_700, rh_700;

	//Assume that 700 mb is equivalent to 3000 m in a standard atmosphere.
	Z_ref = 3000.0;

	//Define the ratio of 700 mb level pressure to the surface pressure (~1000 mb).
	press_ratio = 0.7;

	//Assume dx = 80.0 km, for Walcek (1994).

	//Walcek coefficients.
	f_max = 78.0 + 80.0/15.5;
	one_minus_RHe = 0.196 + (0.76-80.0/2834.0) * (1.0 - press_ratio);
	f_1 = f_max * (press_ratio - 0.1) / 0.6 / 100.0;

	//Convert the gridded topo-surface RH to Td.
	Td = Tdew(Tair, Fmax(0.1, RH), Z);

	//Convert the topo-surface temperature values to 700 mb values.
	delta_Z = Z_ref - Z;
	Td_700 = Td - 1.E-3*Td_lapse_rate * delta_Z;
	Tair_700 = Tair - 1.E-3*T_lapse_rate * delta_Z;

	//Convert each Td to a gridded relative humidity (0-1).
    rh_700 = RHfromTdew(Tair_700, Td_700, Z);

	//Use this RH at 700 mb to define the cloud fraction (0-1).
	fcloud = f_1 * exp((rh_700 - 1.0)/one_minus_RHe);

	fcloud = Fmin(1.0, fcloud);
	fcloud = Fmax(0.0, fcloud);

	return(fcloud);

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void interpolate_meteo(T_INIT *UV, METEO_STATIONS *allmstn, double **value, long **metcol, long metcod, DOUBLEMATRIX *grid, double dn0, short ifill, short iobsint){

	long r,c,n,col;
	long nstn;// number of meteo stations
	double novalue=UV->V->co[2];
	double xmn=UV->U->co[4];
	double ymn=UV->U->co[3];
	double dx=UV->U->co[1];
	double dy=UV->U->co[2];
	double dn;
	DOUBLEVECTOR *var, *xst, *yst;

	nstn=0;
	for(n=1;n<=allmstn->Z->nh;n++){
		col=metcol[n-1][metcod];
		if(col!=-1){
			if(value[n-1][col]!=novalue) nstn++;
		}
	}
	if(nstn==0){
		printf("No data for the variable cod:%ld, num of stations=%ld\n",metcod,allmstn->Z->nh);//stop_execution();
		t_error("Recheck the meteo data and run the model again");
	}
	var=new_doublevector(nstn);
	xst=new_doublevector(nstn);
	yst=new_doublevector(nstn);
	nstn=0;
	for(n=1;n<=allmstn->Z->nh;n++){
		col=metcol[n-1][metcod];
		if(col!=-1){
			if(value[n-1][col]!=novalue){
				nstn++;
				var->co[nstn]=value[n-1][col];
				xst->co[nstn]=allmstn->E->co[n];
				yst->co[nstn]=allmstn->N->co[n];
			}
		}
	}
	if(nstn>1){
		if(iobsint==1){
			dn=dn0;
		}else{
			get_dn(grid->nch, grid->nrh, dx, dy, nstn, &dn);
		}
		barnes_oi(grid->nch, grid->nrh, dx, dy, xmn, ymn, nstn, xst, yst, var, dn, grid, novalue, ifill);
	}else{
		for(r=1;r<=grid->nrh;r++){
			for(c=1;c<=grid->nch;c++){
				grid->co[r][c]=var->co[nstn];
			}
		}
	}
	free_doublevector(var);
	free_doublevector(xst);
	free_doublevector(yst);

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void get_dn(long nc, long nr, double deltax, double deltay, long nstns, double *dn){

	//real dn_max           ! the max obs spacing, dn_r
	//real dn_min           ! dn_r, for large n

	double dn_max, dn_min;

	dn_max = pow(deltax*nc * deltay*nr, 0.5) * ((1.0 + pow((double)nstns, 0.5)) / ((double)nstns - 1.0));
	dn_min = pow((deltax*nc * deltay*nr) / (double)nstns, 0.5);

	*dn = 0.5 * (dn_min + dn_max);

}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void barnes_oi(long nc, long nr, double deltax, double deltay, double xmn, double ymn, long nstns, DOUBLEVECTOR *xstn, DOUBLEVECTOR *ystn, DOUBLEVECTOR *var,
	double dn, DOUBLEMATRIX *grid, double undef, short ifill){

	long r, c, mm, nn, nflag;
	double gamma=0.2;

	double xg,yg; //temporary x and y coords of current cell
	double w1,w2; //weights for Gauss-weighted average
	double wtot1,wtot2; //sum of weights
	double ftot1,ftot2; //accumulators for values, corrections
	double dsq;         //delx**2 + dely**2
	double xa,ya;       //x, y coords of current station
	double xb,yb;       //x, y coords of current station
	DOUBLEVECTOR *dvar;			//estimated error

	double xkappa_1;    // Gauss constant for first pass
	double xkappa_2;    // Gauss constant for second pass
	double rmax_1;      // maximum scanning radii, for first
	double rmax_2;      // and second passes
	double anum_1;      // numerator, beyond scanning radius,
	double anum_2;      // for first and second passes

	dvar=new_doublevector(nstns);

	// Compute the first and second pass values of the scaling parameter
	//   and the maximum scanning radius used by the Barnes scheme.
	//   Values above this maximum will use a 1/r**2 weighting.  Here I
	//   have assumed a gamma value of 0.2.

	// First-round values, Eqn (13).
	xkappa_1 = 5.052 * pow(2.0*dn/Pi, 2.0);

	// Define the maximum scanning radius to have weight defined by
	//   wt = 1.0 x 10**(-30) = exp(-rmax_1/xkappa_1)
	// Also scale the 1/r**2 wts so that when dsq = rmax, the wts match.
	rmax_1 = xkappa_1 * 30.0 * log(10.0);
	anum_1 = 1.E-30 * rmax_1;

	// Second-round values, Eqn (4).
	xkappa_2 = gamma * xkappa_1;
	rmax_2 = rmax_1 * gamma;
	anum_2 = 1.E-30 * rmax_2;

	// Scan each input data point and construct estimated error, dvar, at that point.
	for(nn=1;nn<=nstns;nn++){	//222

		xa = xstn->co[nn];
		ya = ystn->co[nn];
		wtot1 = 0.0;
		ftot1 = 0.0;

		for(mm=1;mm<=nstns;mm++){	//111

			xb = xstn->co[mm];
			yb = ystn->co[mm];
			dsq = pow(xb - xa, 2.0) + pow(yb - ya, 2.0);

			if(dsq<=rmax_1){
				w1 = exp((- dsq)/xkappa_1);
			}else{
				// Assume a 1/r**2 weight.
				w1 = anum_1/dsq;
			}

			wtot1 = wtot1 + w1;
			ftot1 = ftot1 + w1 * var->co[mm];
		} //end 111

		if(wtot1==0.0) printf("stn wt totals zero\n");

		dvar->co[nn]= var->co[nn] - ftot1/wtot1;

	} //end 222

	// Grid-prediction loop.  Generate the estimate using first set of
	//   weights, and correct using error estimates, dvar, and second
	//   set of weights.

	for(c=1;c<=nc;c++){ //666
		for(r=1;r<=nr;r++){	//555

			// xcoords of grid nodes at index r,c
			// ycoords of grid nodes at index r,c
			xg = xmn + deltax * (c-1);
			yg = ymn + deltay * (nr-r);

			// Scan each input data point.
			ftot1 = 0.0;
			wtot1 = 0.0;
			ftot2 = 0.0;
			wtot2 = 0.0;
			nflag = 0;

			for(nn=1;nn<=nstns;nn++){	//333

				xa = xstn->co[nn];
				ya = ystn->co[nn];
				dsq = pow(xg - xa, 2.0) + pow(yg - ya, 2.0);

				if (dsq<=rmax_2){
					w1 = exp((- dsq)/xkappa_1);
					w2 = exp((- dsq)/xkappa_2);
				}else if(dsq<=rmax_1){
					w1 = exp((- dsq)/xkappa_1);
					w2 = anum_2/dsq;
				}else{
					//Assume a 1/r**2 weight.
					w1 = anum_1/dsq;
					nflag = nflag + 1;
					//With anum_2/dsq.
					w2 = gamma * w1;
				}

				wtot1 = wtot1 + w1;
				wtot2 = wtot2 + w2;
				ftot1 = ftot1 + w1 * var->co[nn];
				ftot2 = ftot2 + w2 * dvar->co[nn];

			} //end 333

			if (wtot1==0.0 || wtot2==0.0) printf("wts total zero\n");

			if (ifill==1){
				grid->co[r][c] = ftot1/wtot1 + ftot2/wtot2;
			}else{
				if (nflag<nstns){
					grid->co[r][c] = ftot1/wtot1 + ftot2/wtot2;
				}else{
					grid->co[r][c] = undef;
				}
			}
		}//end 555
	}//end 666
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

double Tdew(double T, double RH, double Z){

	double e, Q, P, Td;
	double A, b, c;

	P=1013.25;	//pressure at sea level in [mbar]
	P*=exp(-Z*0.00013);

	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;

	e=A*exp(b*T/(c+T));
	Q=RH*0.622*e/(P-0.378*e);
	e=Q*P/(0.622+Q*0.378);
	Td=c*log(e/A)/(b-log(e/A));

	return(Td);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

double RHfromTdew(double T, double Tdew, double Z){

	double e, Q, Qs, P, RH;
	double A, b, c;

	P=1013.25;	//pressure at sea level in [mbar]
	P*=exp(-Z*0.00013);

	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;

	e=A*exp(b*Tdew/(c+Tdew));
	Q=0.622*e/(P-0.378*e);
	e=A*exp(b*T/(c+T));
	Qs=0.622*e/(P-0.378*e);
	RH=Q/Qs;

	RH=Fmin(1.0, RH);
	RH=Fmax(0.0, RH);

	return(RH);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void topo_data(double deltax, double deltay, DOUBLEMATRIX *topo, DOUBLEMATRIX *curvature, DOUBLEMATRIX *terrain_slope, DOUBLEMATRIX *slope_az,
	double curve_len_scale, double undef){

	long r,c,inc,h,k;
	long nc=topo->nch;
	long nr=topo->nrh;
	double rad2deg=180.0/Pi;
    double deltaxy,curve_max,topo1,topo2,topo3,topo4,topo5,topo6,topo7,topo8;
	DOUBLEMATRIX *dzdx, *dzdy;

	dzdx=new_doublematrix(nr,nc);
	dzdy=new_doublematrix(nr,nc);

//	Compute the topographic information required to run the wind model.
	initialize_doublematrix(curvature, 0.0);
	initialize_doublematrix(terrain_slope, 0.0);
	initialize_doublematrix(slope_az, 0.0);
	initialize_doublematrix(dzdx, 0.0);
	initialize_doublematrix(dzdy, 0.0);

// CURVATURE CALCULATIONS.

// Compute the average grid increment.
	deltaxy = 0.5 * (deltax + deltay);

// Convert the length scale to an appropriate grid increment.
	inc = Fmax(1, (long)(curve_len_scale/deltaxy));

// Compute the curvature.
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=undef){

				k=inc;
				do{
					if(r-k<1 || c-k<1) k--;
				}while(r-k<1 || c-k<1);
				do{
					if(topo->co[r-k][c-k]==undef) k--;
				}while(topo->co[r-k][c-k]==undef);
				topo1=topo->co[r-k][c-k];

				k=inc;
				do{
					if(r+k>nr || c+k>nc) k--;
				}while(r+k>nr || c+k>nc);
				do{
					if(topo->co[r+k][c+k]==undef) k--;
				}while(topo->co[r+k][c+k]==undef);
				topo2=topo->co[r+k][c+k];

				k=inc;
				do{
					if(r+k>nr || c-k<1) k--;
				}while(r+k>nr || c-k<1);
				do{
					if(topo->co[r+k][c-k]==undef) k--;
				}while(topo->co[r+k][c-k]==undef);
				topo3=topo->co[r+k][c-k];

				k=inc;
				do{
					if(r-k<1 || c+k>nc) k--;
				}while(r-k<1 || c+k>nc);
				do{
					if(topo->co[r-k][c+k]==undef) k--;
				}while(topo->co[r-k][c+k]==undef);
				topo4=topo->co[r-k][c+k];

				k=inc;
				do{
					if(r+k>nr) k--;
				}while(r+k>nr);
				do{
					if(topo->co[r+k][c]==undef) k--;
				}while(topo->co[r+k][c]==undef);
				topo5=topo->co[r+k][c];

				k=inc;
				do{
					if(r-k<1) k--;
				}while(r-k<1);
				do{
					if(topo->co[r-k][c]==undef) k--;
				}while(topo->co[r-k][c]==undef);
				topo6=topo->co[r-k][c];

				k=inc;
				do{
					if(c+k>nc) k--;
				}while(c+k>nc);
				do{
					if(topo->co[r][c+k]==undef) k--;
				}while(topo->co[r][c+k]==undef);
				topo7=topo->co[r][c+k];

				k=inc;
				do{
					if(c-k<1) k--;
				}while(c-k<1);
				do{
					if(topo->co[r][c-k]==undef) k--;
				}while(topo->co[r][c-k]==undef);
				topo8=topo->co[r][c-k];

				curvature->co[r][c] = (4.0 * topo->co[r][c] - topo1 - topo2 - topo3 - topo4) / (pow(2.0, 0.5) * 16.0 * inc * deltaxy) +
									  (4.0 * topo->co[r][c] - topo5 - topo6 - topo7 - topo8) / (16.0 * inc * deltaxy);
			}
		}
	}


// Scale the curvature such that the max abs(curvature) has a value
//   of abs(0.5).  Include a 1 mm curvature in curve_max to prevent
//   divisions by zero in flat terrain where the curvature is zero.
	curve_max = 0.0 + 0.001;
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=undef){
				curve_max = Fmax(curve_max,fabs(curvature->co[r][c]));
			}
		}
	}
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=undef){
				curvature->co[r][c] /= (2.0 * curve_max);
			}
		}
	}


// SLOPE CALCULATIONS.

// Find dzdx.
	for(r=1;r<=nr;r++){
		k=1;
		do{
			if(topo->co[r][k]==undef && k<nc) k++;
		}while(topo->co[r][k]==undef && k<nc);
		if(k<nc){
			h=nc;
            do{
				if(topo->co[r][h]==undef) h--;
			}while(topo->co[r][h]==undef);
            if(h==k){
                dzdx->co[r][h]=0.0;
            }else{
                dzdx->co[r][k] = (topo->co[r][k+1] - topo->co[r][k]) / deltax;
                if((k+1)<(h-1)){
					for(c=k+1;c<=h-1;c++){
                       dzdx->co[r][c] = (topo->co[r][c+1] - topo->co[r][c-1]) /  (2.0 * deltax);
                    }
                }
				dzdx->co[r][h] = (topo->co[r][h] - topo->co[r][h-1] ) / deltax;
			}
		}
	}

// Find dzdy.
	for(c=1;c<=nc;c++){
		k=nr;
		do{
			if(topo->co[k][c]==undef && k>1) k--;
        }while(topo->co[k][c]==undef && k>1);
		if(k>1){
            h=1;
            do{
				if(topo->co[h][c]==undef) h++;
            }while(topo->co[h][c]==undef);
            if(h==k){
                dzdy->co[h][c]=0.0;
            }else{
                dzdy->co[k][c] = (topo->co[k-1][c] - topo->co[k][c]) / deltay;
                if((k+1)>(h-1)){
					for(r=k+1;r<=h-1;r--){
                        dzdy->co[r][c] = (topo->co[r-1][c] - topo->co[r+1][c]) / (2.0 * deltay);
                    }
                }
                dzdy->co[h][c] = (topo->co[h][c] - topo->co[h+1][c]) / deltay;
			}
		}
	}

// Calculate the terrain slope and azimuth.
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(topo->co[r][c]!=undef){

				//Some compilers will not allow dzdx and dzdy to both be 0.0 in the atan2 computation.
				//if (fabs(dzdx->co[r][c])<1.E-10) dzdx->co[r][c] = 1.E-10;
				if (fabs(dzdy->co[r][c])<1.E-10) dzdy->co[r][c] = 1.E-10;

				// Compute the slope azimuth, making sure that north has zero
				//   azimuth.  Also note that for the Ryan wind rotation, the
				//   azimuth values must range from 0 to 360.
				slope_az->co[r][c] = rad2deg * (3.0 / 2.0 * Pi - atan2(dzdy->co[r][c], dzdx->co[r][c]));
				if (slope_az->co[r][c]>=360.0) slope_az->co[r][c] -= 360.0;

				// Compute the slope of the terrain.
				terrain_slope->co[r][c] = rad2deg * atan(pow(pow(dzdx->co[r][c], 2.0) + pow(dzdy->co[r][c], 2.0), 0.5));

			}
		}
	}

	free_doublematrix(dzdx);
	free_doublematrix(dzdy);
}

