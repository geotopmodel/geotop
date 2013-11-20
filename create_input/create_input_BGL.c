#include "stdio.h"

void main() {

FILE *f;
int nc=11, nr=9, c, r, y, m, d, h, days_in_month[12];
double ds=500.0, jd=733774.0, prec=10.0, conc=0.001;
double xll=600000.0, yll=5000000.0, no_data=-9999.0;


days_in_month[0]=31; // Januar
days_in_month[1]=28; // Februar
days_in_month[2]=31; // Maerz
days_in_month[3]=30; // April
days_in_month[4]=31; // Mai
days_in_month[5]=30; // Juni
days_in_month[6]=31; // Juli
days_in_month[7]=31; // August
days_in_month[8]=30; // September
days_in_month[9]=31; // Oktober
days_in_month[10]=30; // November
days_in_month[11]=31; // Dezember


f=fopen("../input/BGL/meteo0001.txt","w");
fprintf(f,"Date,JDfrom0,Iprec,Conc\n");
for (y=2009; y<=2010; y++){
	for (m=1; m<=12; m++){
		for (d=1; d<=days_in_month[m-1]; d++){
			for (h=0; h<24; h++){
				if (m<10){
					if (d<10){
						if (h<10) {fprintf(f,"0%i/0%i/%i 0%i:00,%f,%f,%f\n",d,m,y,h,jd,prec,conc);}
						else {fprintf(f,"0%i/0%i/%i %i:00,%f,%f,%f\n",d,m,y,h,jd,prec,conc);}
					}
					else {
						if (h<10) {fprintf(f,"%i/0%i/%i 0%i:00,%f,%f,%f\n",d,m,y,h,jd,prec,conc);}
						else {fprintf(f,"%i/0%i/%i %i:00,%f,%f,%f\n",d,m,y,h,jd,prec,conc);}
					}
				}
				else {
					if (d<10){
						if (h<10) {fprintf(f,"0%i/%i/%i 0%i:00,%f,%f,%f\n",d,m,y,h,jd,prec,conc);}
						else {fprintf(f,"0%i/%i/%i %i:00,%f,%f,%f\n",d,m,y,h,jd,prec,conc);}
					}
					else {
						if (h<10) {fprintf(f,"%i/%i/%i 0%i:00,%f,%f,%f\n",d,m,y,h,jd,prec,conc);}
						else {fprintf(f,"%i/%i/%i %i:00,%f,%f,%f\n",d,m,y,h,jd,prec,conc);}
					}
				}
				if (jd>733774.000000+20.0-2.0/24.0){prec=0.0;}
				if (jd>733774.000000+40.0-2.0/24.0){prec=5.0;}
				if (jd>733774.000000+60.0-2.0/24.0){prec=0.0;}
				jd=jd+1.0/24.0;
			}
		}
	}
}
fclose(f);


}

