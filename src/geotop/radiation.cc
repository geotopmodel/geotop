
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

#include "constants.h"
#include "struct.geotop.h"
#include "radiation.h"
#include "meteo.h"
#include "tabs.h"
#include "times.h"
#include "util_math.h"
#include "geomorphology.h"
#include "math.optim.h"
#include "timer.h"

extern long number_novalue, number_absent;

extern T_INIT *UV;
extern char *FailedRunFile;
extern long Nr, Nc;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void sun(double JDfrom0, double *E0, double *Et, double *Delta)
{
  double Gamma = 2.0*GTConst::Pi*convert_JDfrom0_JD(JDfrom0)/365.25;

  /** correction sun-earth distance */
  *E0=1.00011+0.034221*cos(Gamma)+0.00128*sin(Gamma)+0.000719*cos(
          2*Gamma)+0.000077*sin(2*Gamma);

  /** correction for sideral day [rad] */
  *Et=0.000075 + 0.001868*cos(Gamma) - 0.032077*sin(Gamma) - 0.014615*cos(
          2*Gamma) - 0.04089*sin(2*Gamma);

  /** solar declination */
  *Delta=0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-0.006758*cos(
          2*Gamma)+0.000907*sin(2*Gamma)-0.002697*cos(3*Gamma)+0.00148*sin(3*Gamma);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double SolarHeight(double JD, double latitude, double Delta, double dh)
{
/**
 * angles in rad
 */
  double sine;
  /** solar hour [0.0-24.0] */
  double h = (JD-floor(JD))*24.0 + dh;

  if (h>=24)
    h-=24.0;

  if (h<0)
    h+=24.0;

  sine = sin(latitude)*sin(Delta) + cos(latitude)*cos(Delta)*cos(GTConst::omega*(12-h));

  if (sine>1)
    sine = 1.;

  if (sine<-1)
    sine = -1.;

  return std::max<double>(asin(sine), 0.0);
}

double SolarHeight_(double JD, double *others)
{
  return SolarHeight(JD, others[0], others[1], others[2]);
}

double SolarHeight__(double JD, void *others) { return SolarHeight_(JD, (double *)others); }


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double SolarAzimuth(double JD, double latitude, double Delta, double dh)
{
/**
 * angles in rad
 */
  double alpha, direction, sine, cosine;
  /** solar hour [0.0-24.0] */
  double h = (JD-floor(JD))*24.0 + dh;

  if (h>=24) h-=24.0;
  if (h<0) h+=24.0;

  /** solar height */
  sine = sin(latitude)*sin(Delta) + cos(latitude)*cos(Delta)*cos(GTConst::omega*(12-h));
  if (sine>1) sine = 1.;
  if (sine<-1) sine = -1.;
  alpha=asin(sine);

  /** solar azimuth */
  if (h<=12)
  {
    if (alpha==GTConst::Pi/2.0) /** zenith */
    {
      direction=GTConst::Pi/2.0;
    }
    else
    {
      cosine=(sin(alpha)*sin(latitude)-sin(Delta))/(cos(alpha)*cos(latitude));
      if (cosine>1) cosine = 1.;
      if (cosine<-1) cosine = -1.;
      direction=GTConst::Pi - acos(cosine);
    }
  }
  else
  {
    if (alpha==GTConst::Pi/2.0) /** zenith */
    {
      direction=3*GTConst::Pi/2.0;
    }
    else
    {
      cosine=(sin(alpha)*sin(latitude)-sin(Delta))/(cos(alpha)*cos(latitude));
      if (cosine>1) cosine = 1.;
      if (cosine<-1) cosine = -1.;
      direction=GTConst::Pi + acos(cosine);
    }
  }

  return direction;
}

double SolarAzimuth_(double JD, double *others)
{
  return SolarAzimuth(JD, others[0], others[1], others[2]);
}

double SolarAzimuth__(double JD, void *others) { return SolarAzimuth_(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double TauatmCosinc(double JD, double *others)
{
  /*double latitude = others[0];
  double Delta = others[1];
  double dh = others[2];*/
  double RH = others[3];
  double T = others[4];
  double P = others[5];
  double slope = others[6];
  double aspect = others[7];
  double Lozone = others[8];
  double alpha = others[9];
  double beta = others[10];
  double albedo = others[11];

  double height, dir;

  height = SolarHeight_(JD, others);
  if (height>0)
  {
    dir = SolarAzimuth_(JD, others);
    return atm_transmittance( std::max<double>(height,asin(0.05)), P, RH, T, Lozone, alpha, beta, albedo)
           * std::max<double>(0.0,cos(slope)*sin(height)+sin(slope)*cos(height)*cos(-aspect+dir));
  }
  else
  {
    return 0.0;
  }
}

double TauatmCosinc_(double JD, void *others) { return TauatmCosinc(JD, (double *)others); }


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double TauatmSinalpha(double JD, double *others)
{
  /*double latitude = others[0];
  double Delta = others[1];
  double dh = others[2];*/
  double RH = others[3];
  double T = others[4];
  double P = others[5];
  double Lozone = others[8];
  double alpha = others[9];
  double beta = others[10];
  double albedo = others[11];

  double height;

  height = SolarHeight_(JD, others);
  if (height>0)
  {
    return atm_transmittance(std::max<double>(height,asin(0.05)),P,RH,T,Lozone,alpha,beta,
                             albedo) * std::max<double>(sin(height),0.05);
  }
  else
  {
    return 0.0;
  }
}

double TauatmSinalpha_(double JD, void *others) { return TauatmSinalpha(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Cosinc(double JD, double *others)
{
  /*double latitude = others[0];
  double Delta = others[1];
  double dh = others[2];*/
  double slope = others[6];
  double aspect = others[7];

  double alpha, direction;

  alpha = SolarHeight_(JD, others);
  direction = SolarAzimuth_(JD, others);
  if (alpha>0)
  {
    return std::max<double>(0.0,cos(slope)*sin(alpha)+sin(slope)*cos(alpha)*cos(
            -aspect+direction));
  }
  else
  {
    return 0.0;
  }
}

double Cosinc_(double JD, void *others) { return Cosinc(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Sinalpha(double JD, double *others)
{
  /*double latitude = others[0];
  double Delta = others[1];
  double dh = others[2];*/

  double alpha;

  alpha = SolarHeight_(JD, others);
  if (alpha>0)
  {
    return std::max<double>(sin(alpha), 0.05);
  }
  else
  {
    return 0.0;
  }
}

double Sinalpha_(double JD, void *others) { return Sinalpha(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Tauatm(double JD, double *others)
{
  /*double latitude = others[0];
  double Delta = others[1];
  double dh = others[2];*/
  double RH = others[3];
  double T = others[4];
  double P = others[5];
  double Lozone = others[8];
  double alpha = others[9];
  double beta = others[10];
  double albedo = others[11];
  double height;

  height = SolarHeight_(JD, others);
  if (height < asin(0.05)) height = asin(0.05);
  return atm_transmittance(height, P, RH, T, Lozone, alpha, beta, albedo);
}

double Tauatm_(double JD, void *others) { return Tauatm(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void shortwave_radiation(double JDbeg, double JDend, double *others,
                         double sin_alpha, double E0, double sky, double SWrefl_surr,
                         double tau_cloud, short shadow, double *SWb, double *SWd, double *cos_inc_bd,
                         double *tau_atm_sin_alpha, short *SWb_yes)
{
  double kd, tau_atm, cos_inc, tau_atm_cos_inc;
  //long i;

  tau_atm = adaptiveSimpsons2(Tauatm_, others, JDbeg, JDend, 1.E-6,
                              20) / (JDend - JDbeg);
  //tau_atm = Tauatm( 0.5*(JDbeg+JDend), others);

  kd=diff2glob(tau_cloud*tau_atm);

  *tau_atm_sin_alpha = adaptiveSimpsons2(TauatmSinalpha_, others, JDbeg, JDend,
                                         1.E-6, 20) / (JDend - JDbeg);
  //*tau_atm_sin_alpha = tau_atm * sin_alpha;

  *SWd = GTConst::Isc*E0*tau_cloud*(*tau_atm_sin_alpha) * sky*kd + (1.-sky)*SWrefl_surr ;

  if (shadow == 1)
  {
    cos_inc = 0.0;
    *SWb = 0.0;
  }
  else
  {
    cos_inc = adaptiveSimpsons2(Cosinc_, others, JDbeg, JDend, 1.E-6,
                                20) / (JDend - JDbeg);
    //cos_inc = Cosinc( 0.5*(JDbeg+JDend), others);
    tau_atm_cos_inc = adaptiveSimpsons2(TauatmCosinc_, others, JDbeg, JDend,
                                        1.E-8, 20) / (JDend - JDbeg);
    //tau_atm_cos_inc = tau_atm * cos_inc;
    //cos_inc = sin_alpha;
    //tau_atm_cos_inc = *tau_atm_sin_alpha;
    *SWb = (1.-kd)*GTConst::Isc*E0*tau_cloud*tau_atm_cos_inc;
  }

  *cos_inc_bd = kd*sin_alpha + (1.-kd)*cos_inc;

  if (sin_alpha > 1.E-5)
  {
    if (shadow == 1)
    {
      *SWb_yes = 0;
    }
    else
    {
      *SWb_yes = 1;
    }
  }
  else
  {
    *SWb_yes = -1;
  }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double diff2glob(double a)
{
/**
 * Ratio between
 * - diffuse radiation
 * - global radiation
 * [Erbs et al.(1982)]
 */
  double k;
  if (a<0.22)
  {
    k = 1.0-0.09*a;
  }
  else if (a<0.80)
  {
    k = 0.9511-0.1604*a + 4.388*pow_2(a) - 16.638*pow(a,3) + 12.336*pow(a,4);
  }
  else
  {
    k = 0.165;
  }
  return (k);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double atm_transmittance(double X, double P, double RH, double T,
                         double Lozone, double a, double b, double rho_g)
{
  GEOTIMER_SECTION(__func__);
  //X = angle of the sun above the horizon [rad]
  //P = pressure [mbar]
  //RH = relative humidity [0-1]
  //T = air temperature [C]

  /*
  //from Mayers and Dale, Predicting Daily Insolation with Hourly Cloud Height and Coverage, 1983, pg 537, Journal of Climate and Applied Meteorology
  double tau_sa;//Reyleigh scattering and gas absorption transmittance
  double tau_w;//transmittance due to water vapor
  double tau_a;//transmittance due to aerosol
  double tau_atm;//global clear sky atmospheric transmittance
  double m;//optical air mass at sea level pressure
  double w;//precipitable water [cm]

  m = 35. * pow( 1224.*pow(sin(X), 2.) + 1. , -0.5 );
  tau_sa = 1.021 - 0.084 * pow(m*(0.000949*P + 0.051), 0.5);
  w = 0.493*RH*(exp(26.23-5416.0/(T+GTConst::tk)))/(T+GTConst::tk);
  tau_w = 1. - 0.077*pow(w*m, 0.3);
  tau_a = pow(0.935, m);
  tau_atm = tau_sa*tau_w*tau_a;
  return(tau_atm);*/


  /** transmissivity under cloudless sky (Iqbal par. 7.5) */
  double mr, ma, w, U1, U3, tau_r, tau_o, tau_g, tau_w, tau_a, tau_aa, rho_a;
  double tau_atm_n, tau_atm_dr, tau_atm_da, tau_atm_dm, tau_atm;
  double w0 = 0.9;
  double Fc = 0.84;

  mr = 1./(sin(X)+0.15*(power((3.885+X*180.0/GTConst::Pi),-1.253)));
  ma = mr*P/GTConst::Pa0;
  w = 0.493*RH*(exp(26.23-5416.0/(T+GTConst::tk)))/(T+GTConst::tk); /** [cm] */

  U1 = w*mr;
  U3 = Lozone*mr;

  tau_r = exp(-.0903*power(ma,.84)*(1.+ma-power(ma,1.01)));
  tau_o = 1. - (.1611*U3*power(1.+139.48*U3, -.3035) - .002715*U3/(1.+.044*U3+.0003*U3*U3));
  tau_g = exp(-.0127*power(ma,.26));
  tau_w = 1. - 2.4959*U1/(power(1.+79.034*U1,.6828)+6.385*U1);
  tau_a = .12445*a - 0.0162 + (1.003 - .125*a) * exp(-b*ma*(1.089*a + .5123)); /** from 7.4.11 Iqbal */
  tau_aa = 1. - (1. - w0)*(1. - ma + power(ma, 1.06))*(1. - tau_a);

  rho_a = .0685 + (1. - Fc)*(1. - tau_a/tau_aa);

  tau_atm_n = .9751 * tau_r * tau_o * tau_g * tau_w * tau_a;
  tau_atm_dr = .79 * tau_o * tau_g * tau_w * tau_aa * .5 * (1.-tau_r) /
               (1. - ma + power(ma, 1.02));
  tau_atm_da = .79 * tau_o * tau_g * tau_w * tau_aa * Fc * (1.-tau_a/tau_aa) /
               (1. - ma + power(ma, 1.02));
  tau_atm_dm = (tau_atm_n + tau_atm_dr + tau_atm_da) * rho_g * rho_a /
               (1. - rho_g * rho_a);
  tau_atm = tau_atm_n + tau_atm_dr + tau_atm_da + tau_atm_dm;

  return (tau_atm);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void longwave_radiation(short state, double pvap, double  /*RH*/, double T,
                        double k1, double k2, double taucloud, double *eps, double *eps_max,
                        double *eps_min)
{
  double taucloud_overcast=0.29; /** after Kimball(1928) */
  FILE *f;

  if (state==1)
  {
    *eps_min = 1.24*pow((pvap/(T+GTConst::tk)),1./7.); /** Brutsaert, 1975 */

  }
  else if (state==2)
  {
    *eps_min = 1.08*(1.0-exp(-pow(pvap,(T+GTConst::tk)/2016.0))); /** Satterlund, 1979 */

  }
  else if (state==3)
  {
    *eps_min = (0.7+5.95*0.00001*pvap*exp(1500/(T+GTConst::tk))); /** Idso(1981) */

  }
  else if (state==4)
  {
    *eps_min = (0.7+5.95*0.00001*pvap*exp(1500/(T+GTConst::tk)));
    *eps_min = -0.792 + 3.161*(*eps_min) - 1.573*(*eps_min)* (*eps_min); /** IDSO + HODGES */

  }
  else if (state==5)
  {
    *eps_min = 0.765; /** Koenig-Langlo & Augstein, 1994 */

  }
  else if (state==6)
  {
    *eps_min = (0.601+5.95*0.00001*pvap*exp(1500.0/(T+GTConst::tk))); /** Andreas and Ackley, 1982 */
  }
  else if (state==7)
  {
    *eps_min = (0.23+k1*pow((pvap*100.)/(T+GTConst::tk),1./k2)); /** Konzelmann (1994) */

  }
  else if (state==8)
  {
    *eps_min = (1.-(1.+46.5*pvap/(T+GTConst::tk))*exp(-pow(1.2+3.*46.5*pvap/(T+GTConst::tk), 0.5))); /** Prata 1996 */

  }
  else if (state==9)
  {
    *eps_min = ( 59.38 + 113.7*pow( (T+GTConst::tk)/273.16, 6.0) + 96.96*pow((465.*pvap/(T+GTConst::tk))/25., 0.5) ) /
               (5.67E-8*pow(T+GTConst::tk,4)); /** Dilley 1998 */

  }
  else
  {

    f = fopen(FailedRunFile, "w");
    fprintf(f,"Error:: Incorrect value for longwave radiation formula\n");
    fclose(f);
    t_error("Fatal Error! Geotop is closed. See failing report.");

  }

  *eps = (*eps_min) * taucloud + 1.0 * (1.-taucloud);
  *eps_max = (*eps_min) * taucloud_overcast + 1.0 * (1.-taucloud);

  /*double fcloud;
  fcloud=pow((1.0-taucloud)/0.75,1/3.4);
  *eps=(*eps_min)*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0); Pirazzini*/


}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double SB(double T)
{
  /**
   * Stefan-Boltzmann Law
   */
  double R;
  R=5.67E-8*pow(T+GTConst::tk,4);
  return (R);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double dSB_dT(double T)
{
  double dR_dT;
  dR_dT=4.0*5.67E-8*pow(T+GTConst::tk,3);
  return (dR_dT);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void rad_snow_absorption(long r, long c, Vector<double> *frac, double R,
                         STATEVAR_3D *snow)
{

  long l, m;
  double res=R, z=0.0, rho, k;

  *frac=0.;

  /** in case of snow */
  if ( (*snow->lnum)(r,c) > 1)
  {

    for (l=(*snow->lnum)(r,c); l>=1; l--)
    {
      m = (*snow->lnum)(r,c)-l+1;
      z += 0.001 * (*snow->Dzl)(l,r,c);
      rho = ((*snow->w_ice)(l,r,c) + (*snow->w_liq)(l,r,c))/ (0.001 * (*snow->Dzl)(l,r,c));
      k = rho/3.0+50.0;
      (*frac)(m) = res-R*exp(-k*z);
      res = R*exp(-k*z);
    }

    (*frac)((*snow->lnum)(r,c)+1) = res;

  }
  else
  {

    (*frac)(0) = res;

  }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double cloud_transmittance(double JDbeg, double JDend, double lat,
                           double Delta, double dh, double RH, double T,
                           double P, double SWd, double SWb, double SW, double E0, double sky,
                           double SWrefl_surr,
                           double Lozone, double alpha, double beta, double albedo)
{

  double *others;
  double tau_atm, tau_atm_sin_alpha, sin_alpha, kd, kd0;
  auto tau = (double)number_novalue;
  long j;

  others = (double *)malloc(12*sizeof(double));
  others[0] = lat;
  others[1] = Delta;
  others[2] = dh;
  others[3] = RH;
  others[4] = T;
  others[5] = P;
  others[8] = Lozone;
  others[9] = alpha;
  others[10] = beta;
  others[11] = albedo;

  tau_atm_sin_alpha = adaptiveSimpsons2(TauatmSinalpha_, others, JDbeg, JDend,
                                        1.E-6, 20) / (JDend - JDbeg);
  //tau_atm = Tauatm( 0.5*(JDbeg+JDend), others);
  //sin_alpha = Sinalpha( 0.5*(JDbeg+JDend), others);
  //tau_atm_sin_alpha = tau_atm*sin_alpha;

  if (tau_atm_sin_alpha > 0)
  {
    if ((long)SWd!=number_absent && (long)SWd!=number_novalue
        && (long)SWb!=number_absent && (long)SWb!=number_novalue)
    {
      if ( SWb+SWd > 0 && SWd > 0)
      {
        kd = SWd / (std::max<double>(0.,SWb)+SWd);
        tau = ( SWd - (1.-sky)*SWrefl_surr ) / ( GTConst::Isc*E0*tau_atm_sin_alpha*sky*kd );
      }

    }
    else if ((long)SW!=number_absent && (long)SW!=number_novalue)
    {

      kd=0.2;
      tau_atm = adaptiveSimpsons2(Tauatm_, others, JDbeg, JDend, 1.E-6,
                                  20) / (JDend - JDbeg);
      sin_alpha = adaptiveSimpsons2(Sinalpha_, others, JDbeg, JDend, 1.E-6,
                                    20) / (JDend - JDbeg);

      j=0;

      do
      {

        j++;
        kd0=kd;
        //SW = (1-kd(T))*GTConst::Isc*T*sin + sky*Kd(T)*GTConst::Isc*T*sin + (1-sky)*SWsurr
        //T=Ta*Tc
        //SW - (1-sky)*SWsurr = Tc * (GTConst::Isc*Ta*sin) * ( (1-kd) + sky*kd )
        //Tc = ( SW - (1-sky)*SWsurr ) / ( (GTConst::Isc*Ta*sin) * ( (1-kd) + sky*kd )
        tau = ( SW - (1.-sky)*SWrefl_surr ) / ( GTConst::Isc*E0*tau_atm_sin_alpha * ( (1-kd) + sky*kd ) );
        if (tau > 1) tau = 1.0;
        if (tau < 0) tau = 0.0;
        kd = diff2glob(tau * tau_atm);

      }
      while (fabs(kd0-kd)>0.005 && j<1000);

    }
  }

  if ( (long)tau != number_novalue)
  {
    if (tau<GTConst::min_tau_cloud) tau=GTConst::min_tau_cloud;
  }

  free(others);

  return tau;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double find_tau_cloud_station(double JDbeg, double JDend, long i, METEO *met,
                              double Delta, double E0, double Et,
                              double ST, double SWrefl_surr, double Lozone, double alpha, double beta,
                              double albedo)
{

  double P, RH, T, c;

  /** pressure [mbar] */
  P = pressure((*met->st->Z)(i));

  /** relative humidity */
  if ((long)met->var[i-1][iRh] != number_novalue
      && (long)met->var[i-1][iRh] != number_absent)
  {
    RH=met->var[i-1][iRh]/100.;
  }
  else
  {
    if ( (long)met->var[i-1][iT] != number_absent
         && (long)met->var[i-1][iT] != number_novalue
         && (long)met->var[i-1][iTdew] != number_absent
         && (long)met->var[i-1][iTdew] != number_novalue)
    {
      RH=RHfromTdew(met->var[i-1][iT], met->var[i-1][iTdew], (*met->st->Z)(i));
    }
    else
    {
      RH=0.4;
    }
  }
  if (RH<0.01) RH=0.01;

  T=met->var[i-1][iT];
  if ((long)T == number_novalue || (long)T == number_absent) T=0.0;

  c = cloud_transmittance(JDbeg, JDend, (*met->st->lat)(i)*GTConst::Pi/180., Delta,
                          ((*met->st->lon)(i)*GTConst::Pi/180. - ST*GTConst::Pi/12. + Et)/GTConst::omega, RH,
                          T, P, met->var[i-1][iSWd], met->var[i-1][iSWb], met->var[i-1][iSW], E0,
                          (*met->st->sky)(i), SWrefl_surr,
                          Lozone, alpha, beta, albedo);
  return c;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short shadows_point(double **hor, long n, double alpha, double azimuth,
                    double tol_mount, double tol_flat)
{
/**
 * Routine that tells you whether a point is in shadow or not, depending on the solar azimuth, elevation and horizon file at that point
 * Author: Matteo Dall'Amico, May 2011
 *
 * Inputs:
 * - Matrix<double>* hor_height: matrix of horizon_height at the point
 * - double alpha: solar altitude [degree]
 * - double azimuth: solar azimuth [degree]
 * - double tol_mount: tolerance over a mountaneaus horizon to have a reliable cloud datum [degree]
 * - double tol_flat: tolerance over a mountaneaus horizon to have a reliable cloud datum [degree]
 *
 * Output:
 * - shad: 1=the point is in shadow, 0 the point is in sun
 */
  double horiz_H; /** horizon elevation at a defined solar azimuth */
  double w; /** weight */
  long i,iend=0,ibeg=0;
  short shad;

  /** compare the current solar azimuth with the horizon matrix */
  if (azimuth>=hor[n-1][0] || azimuth<hor[0][0])
  {
    iend=0;
    ibeg=n-1;
  }
  else
  {
    for (i=1; i<=n-1; i++)
    {
      if (azimuth>=hor[i-1][0] && azimuth<hor[i][0])
      {
        iend=i;
        ibeg=i-1;
      }
    }
  }

  if (iend>ibeg)
  {
    w=(azimuth-hor[ibeg][0])/(hor[iend][0]-hor[ibeg][0]);
  }
  else if (azimuth>hor[ibeg][0])
  {
    w=(azimuth-hor[ibeg][0])/(hor[iend][0]+360.0-hor[ibeg][0]);
  }
  else
  {
    w=(azimuth-(hor[ibeg][0]-360.0))/(hor[iend][0]-(hor[ibeg][0]-360.0));
  }

  horiz_H = w * hor[iend][1] + (1.0-w) *
                               hor[ibeg][1]; /** horizon elevation at a particular time */

  if (alpha<tol_flat)
  {
    shad=1;
  }
  else if (alpha<horiz_H+tol_mount)
  {
    shad=1;
  }
  else
  {
    shad=0;
  }

  return (shad);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void shadow_haiden(Matrix<double> *Z, double alpha, double direction,
                   Matrix<short> *SH)
{
/*** Author: Thomas Haiden, Year: 16 june 2003
 * Function that calculates if each pixel (matrix of shadows) is in shadow or at sun given
 * a solar elevation angle and a solar azimuth.
 *
 * Inputs:
 * - top: elevation DTM [m]
 * - alpha: solar elevation angle [radiants]
 * - direction: solar azimuth angle (from N clockwise) [radiants]
 * - point: flag indicating whether the simulation is a point (=1) or a distributed simulation
 *
 * Outputs:
 * - shadow: shadows matrix (1 shadow, 0 sun)
 *
 * Imported by Matteo Dall'Amico on August 2009. Authorization: see email of David Whiteman on 16 June 2011
 * Revised by Stefano Endrizzi on May 2011.
 */
  long orix,oriy=0,k,l,kk,ll,r,c,rr,cc;
  double sx,sy,sz,xp,yp,x,y,zray,ztopo,q,z1,z2;

  long nk=Nr;//y
  long nl=Nc;//x
  double GDX = (*UV->U)(2);
  double GDY = (*UV->U)(1);

  sx=sin(direction)*cos(alpha);
  sy=cos(direction)*cos(alpha);
  sz=sin(alpha);

  if (fabs(sx)>fabs(sy))
  {

    if (sx>0)
    {
      orix=1;
    }
    else
    {
      orix=-1;
    }

    if (fabs(sy)>1.E-10)
    {
      if (sy>0)
      {
        oriy=1;
      }
      else
      {
        oriy=-1;
      }
    }
    else
    {
      orix=0;
    }

    for (k=0; k<nk; k++)
    {
      for (l=0; l<nl; l++)
      {
        //r=row(GDY*k+0.5*GDY, Nr, UV, number_novalue);
        //c=col(GDX*l+0.5*GDX, Nc, UV, number_novalue);
        r=Nr-k;
        c=l+1;
        kk=k;
        ll=l;
        xp=GDX*ll+0.5*GDX;
        yp=GDY*kk+0.5*GDY;
        //rr=row(GDY*kk+0.5*GDY, Nr, UV, number_novalue);
        //cc=col(GDX*ll+0.5*GDX, Nc, UV, number_novalue);
        rr=Nr-kk;
        cc=ll+1;
        zray=(*Z)(rr,cc);
        (*SH)(r,c)=0;

        while ( ((*SH)(r,c)==0) && (kk>0)&&(kk<nk-1)&&(ll>0)&&(ll<nl-1) )
        {
          q=((ll+orix)*GDX+0.5*GDX-xp)/sx;
          y=yp+q*sy;
          if (fabs(y-(GDY*kk+0.5*GDY))<GDY)
          {
            ll=ll+orix;
            xp=GDX*ll+0.5*GDX;
            yp=y;
            //rr=row(GDY*kk+0.5*GDY, Nr, UV, number_novalue);
            //cc=col(GDX*ll+0.5*GDX, Nc, UV, number_novalue);
            rr=Nr-kk;
            cc=ll+1;
            z1=(*Z)(rr,cc);
            //rr=row(GDY*kk+0.5*GDY+oriy*GDY, Nr, UV, number_novalue);
            rr=Nr-(kk+oriy);
            z2=(*Z)(rr,cc);
            ztopo=z1+(z2-z1)*(yp-(GDY*kk+0.5*GDY))/(oriy*GDY);
          }
          else
          {
            q=((kk+oriy)*GDY+0.5*GDY-yp)/sy;
            x=xp+q*sx;
            kk=kk+oriy;
            xp=x;
            yp=GDY*kk+0.5*GDY;
            //rr=row(GDY*kk+0.5*GDY, Nr, UV, number_novalue);
            //cc=col(GDX*ll+0.5*GDX, Nc, UV, number_novalue);
            rr=Nr-kk;
            cc=ll+1;
            z1=(*Z)(rr,cc);
            //cc=col(GDX*ll+0.5*GDX+orix*GDX, Nc, UV, number_novalue);
            cc=ll+orix+1;
            z2=(*Z)(rr,cc);
            ztopo=z1+(z2-z1)*(xp-(GDX*ll+0.5*GDX))/(orix*GDX);
          }
          zray=zray+q*sz;
          if (ztopo>zray) (*SH)(r,c)=1;
        }
      }
    }

  }
  else
  {

    if (sy>0)
    {
      oriy=1;
    }
    else
    {
      oriy=-1;
    }

    if (fabs(sx)>1.E-10)
    {
      if (sx>0)
      {
        orix=1;
      }
      else
      {
        orix=-1;
      }
    }
    else
    {
      orix=0;
    }

    for (k=0; k<nk; k++)
    {
      for (l=0; l<nl; l++)
      {
        //r=row(GDY*k+0.5*GDY, Nr, UV, number_novalue);
        //c=col(GDX*l+0.5*GDX, Nc, UV, number_novalue);
        r=Nr-k;
        c=l+1;
        kk=k;
        ll=l;
        xp=GDX*ll+0.5*GDX;
        yp=GDY*kk+0.5*GDY;
        //rr=row(GDY*kk+0.5*GDY, Nr, UV, number_novalue);
        //cc=col(GDX*ll+0.5*GDX, Nc, UV, number_novalue);
        rr=Nr-kk;
        cc=ll+1;
        zray=(*Z)(rr,cc);
        (*SH)(r,c)=0;

        while ( ((*SH)(r,c)==0) &&(kk>0)&&(kk<nk-1)&&(ll>0)&&(ll<nl-1))
        {
          q=((kk+oriy)*GDY+0.5*GDY-yp)/sy;
          x=xp+q*sx;
          if (fabs(x-(GDX*ll+0.5*GDX))<GDX)
          {
            kk=kk+oriy;
            yp=GDY*kk+0.5*GDY;
            xp=x;
            //rr=row(GDY*kk+0.5*GDY, Nr, UV, number_novalue);
            //cc=col(GDX*ll+0.5*GDX, Nc, UV, number_novalue);
            rr=Nr-kk;
            cc=ll+1;
            z1=(*Z)(rr,cc);
            //cc=col(GDX*ll+0.5*GDX+orix*GDX, Nc, UV, number_novalue);
            cc=ll+orix+1;
            z2=(*Z)(rr,cc);
            ztopo=z1+(z2-z1)*(xp-(GDX*ll+0.5*GDX))/(orix*GDX);
          }
          else
          {
            q=((ll+orix)*GDX+0.5*GDX-xp)/sx;
            y=yp+q*sy;
            ll=ll+orix;
            yp=y;
            xp=GDX*ll+0.5*GDX;
            //rr=row(GDY*kk+0.5*GDY, Nr, UV, number_novalue);
            //cc=col(GDX*ll+0.5*GDX, Nc, UV, number_novalue);
            rr=Nr-kk;
            cc=ll+1;
            z1=(*Z)(rr,cc);
            //rr=row(GDY*kk+0.5*GDY+oriy*GDY, Nr, UV, number_novalue);
            rr=Nr-(kk+oriy);
            z2=(*Z)(rr,cc);
            ztopo=z1+(z2-z1)*(yp-(GDY*kk+0.5*GDY))/(oriy*GDY);
          }
          zray=zray+q*sz;
          if (ztopo>zray) (*SH)(r,c)=1;
        }
      }
    }
  }
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double find_albedo(double dry_albedo, double sat_albedo, double wat_content,
                   double residual_wc, double saturated_wc)
{

  return (dry_albedo + (sat_albedo-dry_albedo) * (wat_content - residual_wc) /
                       (saturated_wc - residual_wc) );

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void find_actual_cloudiness(double *tau_cloud, double *tau_cloud_av,
                            short *tau_cloud_yes, short *tau_cloud_av_yes,
                            METEO *met, double JDb, double JDe, double Delta, double E0, double Et,
                            double ST, double SWrefl_surr,
                            double Lozone, double alpha, double beta, double albedo)
{

  /** SWdata = flag -> 0= no radiation data measured, 1=beam and diff measured, 2=global measured */
  short SWdata;
  double tc;

  if ((long)met->var[met->nstsrad-1][iSWb]!=number_absent
      && (long)met->var[met->nstsrad-1][iSWd]!=number_absent)
  {
    if ((long)met->var[met->nstsrad-1][iSWb]!=number_novalue
        && (long)met->var[met->nstsrad-1][iSWd]!=number_novalue)
    {
      SWdata=2;
    }
    else
    {
      SWdata=0;
    }
  }
  else if ((long)met->var[met->nstsrad-1][iSW]!=number_absent)
  {
    if ((long)met->var[met->nstsrad-1][iSW]!=number_novalue)
    {
      SWdata=1;
    }
    else
    {
      SWdata=0;
    }
  }
  else
  {
    SWdata=0;
  }

  if (SWdata>0)
  {
    tc = find_tau_cloud_station(JDb, JDe, met->nstsrad, met, Delta, E0, Et, ST,
                                SWrefl_surr, Lozone, alpha, beta, albedo);
    if ( (long)tc != number_novalue)
    {
      *tau_cloud_yes = 1;
      *tau_cloud = tc;
    }
    else
    {
      *tau_cloud_yes = 0;
    }
  }
  else
  {
    *tau_cloud_yes = 0;
  }

  if ( (long)met->var[met->nstcloud-1][iC]!=number_absent
       && (long)met->var[met->nstcloud-1][iC]!=number_novalue )
  {

    tc = met->var[met->nstcloud-1][iC];

    *tau_cloud_av_yes = 1;
    tc = 1. - 0.71*tc; /** from fraction of sky covered by clouds to cloud transmissivity after Kimball (1928) */
    if (tc > 1) tc = 1.;
    if (tc < 0) tc = 0.;
    *tau_cloud_av = tc;

  }
  else if ( (long)met->var[met->nstcloud-1][itauC]!=number_absent
            && (long)met->var[met->nstcloud-1][itauC]!=number_novalue )
  {

    tc = met->var[met->nstcloud-1][itauC];

    *tau_cloud_av_yes = 1;
    if (tc > 1) tc = 1.;
    if (tc < 0) tc = 0.;
    *tau_cloud_av = tc;

  }
  else
  {

    *tau_cloud_av_yes = 0;

  }
}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************



