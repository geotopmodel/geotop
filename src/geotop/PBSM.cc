/* STATEMENT:

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)

 Copyright (c), 2016 - GEOtop Foundation

 This file is part of GEOtop 2.1

 GEOtop 2.1  is a free software and is distributed under GNU General Public
 License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
 even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE

 GEOtop 2.1  is distributed as a free software in the hope to create and support
 a community of developers and users that constructively interact. If you just
 use the code, please give feedback to the authors and the community. Any way
 you use the model, may be the most trivial one, is significantly helpful for
 the future development of the GEOtop model. Any feedback will be highly
 appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.


 @brief PRAIRIE BLOWING SNOW MODEL CODE

 Code written by Stefano Endrizzi by translating and adapting the idea behind
 the PBSM Fortran Code by John Pomeroy. The author does not guarantee the
 perfect conformance of this code with the Fortran one. However, he asks to give
 credit to Pomeroy, 1993, when using it with satisfaction.

 Reference:
 John Pomeroy
 The Prairie Blowing Snow Model: characteristics, validation, operation
 Journal of Hydrology, 144 (1993) 165-192
 */

#include "PBSM.h"
#include "geotop_common.h"
#include "global_logger.h"

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

// global variables
static double T;
static double Ustar;
static double Z0;
static double Zr;

static double k_atm;
static double Diff;
static double B;
static double SvDens;
static double Sigma;

static float CC1;
static float CC2;
static float CC3;
static float CC4;
static float CC5;

/// @brief Prarie Blowing Snow Model - Pomeroy et al. (1993)

void Pbsm(long r,
          long c,
          double Fetch,
          double N,
          double dv,
          double Hv,
          double rho_sn,
          double zmeas,
          double V,
          double Ta,
          double RH,
          double *Trans,
          double *Subl,
          double *Salt,
          double Dsnow,
          double slope)
{
  //     Modified Calculations for Mean Particle Mass in this version
  //     program to calculate blowing snow horizontal flux, sublimation rate
  //     and latent heat flux due to snow sublimation for a variety of
  //     windspeeds, boundary layers and surface conditions.

  //     All variable and constants entered into the programme are in SI and
  //     use Canadian Atmospheric Environement Service Meteorological data
  //     format.  Snow transport is in kg per square meter per second
  //     from the surface to 5 metres height.  Sublimation is totaled to the top
  //     of the boundary layer for diffusion, based on the meteorological
  //     Fetch and is expressed in millimeters of blowing snow lost over
  //     a square meter of snow surface per half hour

  double A, Alpha, C, DmDt, Es, F, F1, Hsalt, Lambda, Mpm, Mpr, Nsalt, Nuss,
         TQsalt, TQsum;
  double RauTerm, Reyn, SBsalt, SBsum, SigmaZ, Ustar1, Usthr, Vs, Vsalt, Z,
         Zstb;
  long i, k;

  //  define constants for equations

  float M = 18.01;  //{molecular weight of water (kg/kmole)}
  float R = 8313;   //{universal gas constant (J/(kmole K))}
  // float ZD = 0.3;     //{height of boundary-layer at xd (m) Takeuchi (1980)}
  // float XD = 300;     //{Fetch to develop b-l to ZD (m)}
  float rho = 1.2;   //{Air density kg/m3}
  float Beta = 170;  //{Cr/Cs = 170}
  double Bound = 5.0;

  //  coefficients

  CC1 = 2.8;  //{2.3}
  CC2 = 1.6;
  CC3 = 4.2;  //{3.25} {e = 1/(CC3*Ustar)}
  CC4 = -1.55;
  CC5 = -0.544;

  // Initialization

  TQsalt = 0.0;  //{Total saltation flux}
  TQsum = 0.0;   //{Total Suspension}
  SBsalt = 0.0;
  SBsum = 0.0;
  T = Ta + GTConst::tk;  // Convert to Deg. K}

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  if (Fetch >= 300)
    {
      // Compute stubble coefficients

      Lambda = N * dv * Hv;  //{Raupach Eq. 1}, N=roughness elements per surface
      //unit, dv=diameter of roughness elements, Hv=their
      //height
      Zstb = 0.5 * Lambda;   //{Lettau, used for susp Z0''}
      RauTerm = 1. / (1. + Beta * Lambda);  //{Raupach eq. 7}

      // Check for data errors    Fluxes set to zero for the hour}

      k_atm = 0.00063 * T + 0.0673;           //{therm. cond. of atm. (J/(msK))}
      Diff = 2.06E-5 * pow(T / 273.0, 1.75);  //{diffus. of w.vap. atmos. (m2/s}
      B = GTConst::Ls * M / (R * T) - 1.0;

      // find undersaturation of w. vapour at 2 metres}

      Es = 611.15 * exp(22.452 * Ta / T);  //{sat pressure}
      SvDens = (Es * M) / (R * T);         //{sat density}
      Sigma = RH - 1.0;                    //{undersaturation at 2 m}

      // calculate Ustar

      Ustar = 0.001;  // first guess
      i = 0;
      Z0 = CC2 * 0.07519 * pow(Ustar, 2.) / (2. * GTConst::GRAVITY) +
           Zstb;  // Liston 1998
      F = (Ustar / GTConst::ka) * log(zmeas / Z0) - V;
      do
        {
          Ustar1 = Ustar;
          F1 =
            (1. / GTConst::ka) * log(zmeas / Z0) -
            (Ustar / GTConst::ka) * (CC2 * 0.07519 / GTConst::GRAVITY) * Ustar / Z0;
          k = 0;
          do
            {
              Ustar = Ustar1 - (F / F1) / pow(2., (double)k);
              Z0 = CC2 * 0.07519 * pow(Ustar, 2.) / (2. * GTConst::GRAVITY) + Zstb;
              k++;
            }
          while (fabs((Ustar / GTConst::ka) * log(zmeas / Z0) - V) > fabs(F) &&
                 k < 10);
          F = (Ustar / GTConst::ka) * log(zmeas / Z0) - V;
          i++;
        }
      while (fabs(F) > 1.E-3 && i < 500);

      if (fabs(F) > 1.E-3 || Ustar < 0 || Ustar > V)
        {
          Ustar = V / 5.;
          lg->logf(
            "Warning: Ustar does not converge r:%ld c:%ld F:%e Ustar:%f V:%f Z0:%f "
            "Zstb:%f slope:%f Dsnow:%f \n",
            r, c, F, Ustar, V, Z0, Zstb, slope, Dsnow);
        }

      // define saltation parameters and calculate saltation
      //  rate using 10/1987 MODEL OF BLOWING SNOW EQUATIONS}

      if (rho_sn <= 300)
        {
          Usthr = 0.10 * exp(0.003 * rho_sn);
        }
      else
        {
          Usthr = 0.005 * exp(0.013 * rho_sn);
        }

      // printf("%ld %ld Ustar:%f rho_sn:%f Uth:%f\n",r,c,Ustar,rho_sn,Usthr);

      if (Ustar > Usthr)
        {
          Nsalt =
            2. * rho / (CC2 * CC3 * Ustar) *
            (RauTerm - pow(Usthr, 2.0) / pow(Ustar, 2.0));  //{Eq. 4.14 updated}

          if (Nsalt > 0)
            {
              // {saltation transport}

              Hsalt = CC2 / (2. * GTConst::GRAVITY) * pow(Ustar, 2.0);  //{Eq. 4.13}
              TQsalt = CC1 * Usthr * Nsalt * Hsalt;                     //{Eq. 4.20}

              // {calculate sublimation rate in the saltation layer}

              Mpr = 100E-6;
              Alpha = 5.0;

              SigmaZ =
                Sigma *
                (1.019 + 0.027 * log(Hsalt));  //{ Eq. 6.20, Revised in May. 1997}
              if (SigmaZ > -0.01) SigmaZ = -0.01;
              Vsalt = 0.6325 * Ustar + 2.3 * Usthr;  //{Eq. 6.25}

              Reyn = (2.0 * Mpr * Vsalt) / 1.88E-5;  //{Eq. 6.22}
              Nuss = 1.79 + 0.606 * sqrt(Reyn);      //{Eq. 6.21}
              A = k_atm * T * Nuss;
              C = 1.0 / (Diff * SvDens * Nuss);
              DmDt = (2.0 * GTConst::Pi * Mpr * SigmaZ) / (GTConst::Ls * B / A + C);
              //{Eq. 6.16} {Gamma Dist. Corr.}
              Mpm = 4.0 / 3.0 * GTConst::Pi * GTConst::rho_i * Mpr * pow(Mpr, 2.0) *
                    (1.0 + 3.0 / Alpha + 2.0 / pow(Alpha, 2.0));

              Vs = DmDt / Mpm;  //{Sublimation rate coefficient Eq. 6.13}

              SBsalt = Vs * Nsalt * Hsalt;  //{Eq. 6.11}

              // calculate mass flux in the suspended layers and the sublimation
              //   rate for layers of height Inc from height r to b}

              // Loop to find the first suspended drift density level, r
              //   from the reference level Zr
              //   To preserve continuity with saltation the first suspended
              //   level drift density is less than or equal to Nsalt.}

              Zr = 0.05628 * Ustar;  //{Eq. 5.27}
              Z = Zr;

              // Nz(Z)=Nz(Zr)*exp(CC4*Zr^CC5-CC4*Z^CC5)

              i = 0;
              do
                {
                  F = 0.8 * exp(CC4 * pow(Zr, CC5) - CC4 * pow(Z, CC5)) - Nsalt;
                  F1 = 0.8 * exp(CC4 * pow(Zr, CC5) - CC4 * pow(Z, CC5)) *
                       (-CC4 * CC5) * pow(Z, CC5 - 1.);
                  Z -= F / F1;
                  i++;
                }
              while (fabs(F) > 1.E-4 && i < 10);

              // find height of fully-developed boundary layer for turbulent
              //   diffusion using a form of Pasquills plume dispersion eq.
              //   iterate towards Bound}

              /*i=0;
              Bound = 1.0;
              do{
                  F = ZD + ka*ka * (Fetch - XD) * pow(log(Bound/Z0) * log(ZD/Z0),
              -0.5) - Bound;     //{Eq. 6.9} F1 = -0.5 * ka*ka * (Fetch - XD) *
              pow(log(Bound/Z0) * log(ZD/Z0), -1.5) * (log(ZD/Z0)) / Bound - 1.0;
                  Bound -= F/F1;
                  i++;
              }while (fabs(F)>1.E-4 && i<10);
              if(Bound>5) Bound=5.;*/

              if (Bound > Z)
                {
                  TQsum = adaptiveSimpsons(suspension, Z, Bound, 1.E-8, 10);
                  SBsum = adaptiveSimpsons(sublimation, Z, Bound, 1.E-8, 10);
                }
            }
        }
    }

  *Trans = TQsum + TQsalt;          //{kg/m-width/s}
  *Subl = (SBsum + SBsalt) * (-1);  //{kg/m2/s}
  *Salt = TQsalt;                   //{kg/m-width/s}

}  //{PBSM procedure}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double suspension(double Z)
{
  double Nz, Uz, UstarZ;

  Nz = 0.8 * exp(CC4 * pow(Zr, CC5) - CC4 * pow(Z, CC5));
  UstarZ = Ustar * pow(1.2 / (1.2 + Nz), 0.5);  //{Eq. 5.17a}
  Uz = (UstarZ / GTConst::ka) * log(Z / Z0);    //{Eq. 4.17r}

  return (Nz * Uz);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double sublimation(double Z)
{
  double Nz, Uz, UstarZ, Mpr, Alpha, SigmaZ, Omega, Vsusp, Reyn, Nuss, A, C,
         DmDt, Mpm, Vs;

  Nz = 0.8 * exp(CC4 * pow(Zr, CC5) - CC4 * pow(Z, CC5));
  UstarZ = Ustar * pow(1.2 / (1.2 + Nz), 0.5);  //{Eq. 5.17a}
  Uz = (UstarZ / GTConst::ka) * log(Z / Z0);    //{Eq. 4.17r}

  Mpr = 4.6E-5 * pow(Z, -0.258);  //{Eq. 6.15}
  if (Z >= 5.0) Mpr = 30E-6;

  Alpha = 4.08 + 12.6 * Z;  //{Eq. 6.14}
  if (Z >= 1.5) Alpha = 25.0;

  SigmaZ =
    Sigma * (1.019 + 0.027 * log(Z));  //{ Eq. 6.20, Revised in May. 1997}
  if (SigmaZ > -0.01) SigmaZ = -0.01;

  Omega = 1.1E7 * pow(Mpr, 1.8);           //{Eq. 5.18}
  Vsusp = Omega + 0.0106 * pow(Uz, 1.36);  //{Eq. 6.23 - 6.24}
  Reyn = (2.0 * Mpr * Vsusp) / 1.88E-5;    //{Eq. 6.22}

  Nuss = 1.79 + 0.606 * sqrt(Reyn);  //{Eq. 6.21}
  A = k_atm * T * Nuss;
  C = 1.0 / (Diff * SvDens * Nuss);
  DmDt = (2.0 * GTConst::Pi * Mpr * SigmaZ) / (GTConst::Ls * B / A + C);
  Mpm = 1.333 * GTConst::Pi * GTConst::rho_i * pow(Mpr, 3.0) *
        (1.0 + 3.0 / Alpha +
         2.0 / pow(Alpha, 2.0));  //{Eq. 6.16} {Gamma Dist. Corr.}

  Vs = DmDt / Mpm;  //{Eq. 6.13}

  return (Vs * Nz);
}
