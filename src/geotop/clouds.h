
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 1.225-15
 
 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 1.225-15  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */


short fill_meteo_data_with_cloudiness(double **meteo, long meteolines, double **horizon, long horizonlines, double lat, double lon, double ST, double Z, double sky,
									  double SWrefl_surr, long ndivday, double rotation, double Lozone, double alpha, double beta, double albedo);

void cloudiness(double **meteo, long meteolines, double **horizon, long horizonlines, double lat, double lon, double ST, double Z, double sky, double SWrefl_surr, 
				double *cloudtrans, long ndivday, double rotation, double Lozone, double alpha, double beta, double albedo);

double find_cloudiness(long n, double **meteo, long meteolines, double lat, double lon, double ST, double Z, double sky, double SWrefl_surr, double rotation, double Lozone, double alpha, double beta, double albedo);

double average_cloudiness(long n0, long n1, double **meteo, long meteolines, double lat, double lon, double ST, double Z, double sky, double SWrefl_surr, double rotation, double Lozone, double alpha, double beta, double albedo);

void find_sunset(long nist, long *n0, long *n1, double **meteo, long meteolines, double **horizon, long horizonlines, 
				 double lat, double lon, double ST, double rotation);


