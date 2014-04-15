void part_snow(double prec_total, double *prec_rain, double *prec_snow, double temperature, double t_rain, double t_snow)

/*	Inputs:	prec_total: 	matrix of total precipitation
			temperature:	matrix of air temperature
			t_rain:			temperature above wich all precipitation is rain
			t_snow:			temperature below wich all precipitation is snow
	Outputs:prec_snow:		matrix of solid precipitation
			prec_rain:		matrix of liquid precipitation
*/

{

	if(temperature<=t_snow){
		*prec_snow=prec_total;
		*prec_rain=0.0;
	}else if(temperature>=t_rain){
		*prec_snow=0.0;
		*prec_rain=prec_total;
	}else{
		*prec_snow=prec_total*(t_rain-temperature)/(t_rain-t_snow);
		*prec_rain=prec_total*(temperature-t_snow)/(t_rain-t_snow);
	}
	
}

