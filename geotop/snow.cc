
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 2.0.0
 
 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */


//Author: Stefano Endrizzi
//Date: 13 November 2005
//Contents: Snow subroutines
#include <string>
#include "snow.h"
#include "geotop_common.h"
#include "inputKeywords.h"

using namespace GTConst;

//Private Declarations
const double simpl_snow = 1E-1;
const double no_snow = 1E-6;

static void snow_compactation(double Dt, long r, long c, long l, Statevar3D *snow, double slope, Par *par);
static void snowlayer_merging(double a, long r, long c, Statevar3D *snow, long l1, long l2, long lres);
static void write_snow(long r, long c, long l, Statevar3D *snow);
static void merge_layers(double a, long r, long c, Statevar3D *snow, long l1);
static void initialize_snow(long r, long c, long l, Statevar3D *snow);
static double max_dtheta_snow(double a, double b);

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//Jordan et al., 1999
double rho_newlyfallensnow(double u, double Tatm)
{

    double rho, T;

    if(Tatm > 260.15 - tk)
    {
        T = Tatm;
        if(T >= 275.65 - tk) T = 275.65 - tk;
        rho = 500.0 * (1.0 - 0.951 * exp(-1.4 * pow(278.15 - tk - T, -1.15) - 0.008 * pow(u, 1.7)));
    }
    else
    {
        rho = 500.0 * (1.0 - 0.904 * exp(-0.008 * pow(u, 1.7)));
    }

    return(rho);

}

//Valt et al., 2014 EGU
double rho_valt(double Tatm)
{

    double rho;

    if ( Tatm >=1) rho = 188;
    else if (Tatm >= -4 && Tatm < 1) rho = 114+11*(Tatm+4);
    else if (Tatm >= -21 && Tatm < -4) rho = 73+2.41*(Tatm+21);
    else rho = 73;

    return(rho);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

static void snow_compactation(double Dt, long r, long c, long l, Statevar3D *snow, double slope, Par *par)
{

    long m;
    double theta_i, theta_w, c1, c2, c3, c4, c5, load, eta0, eta, CR1, CR2;

    theta_i = snow->w_ice[l][r][c] / (0.001 * snow->Dzl[l][r][c] * rho_i);
    theta_w = snow->w_liq[l][r][c] / (0.001 * snow->Dzl[l][r][c] * rho_w);

    if( theta_i < par->snow_maxpor)
    {

        //DESTRUCTIVE METAMORPHISM
        //gamma=theta*rho
        if(theta_i * rho_i <= par->snow_density_cutoff)
        {
            c1 = par->drysnowdef_rate;
        }
        else
        {
            c1 = par->drysnowdef_rate * exp(-0.046 * (rho_i * theta_i - par->snow_density_cutoff));
        }
        if(theta_w > 0.001) c1 *= par->wetsnowdef_rate;
        c2 = 2.777E-6;                            //[s^-1]
        c3 = 0.04;                                //[K^-1]
        CR1 = -c1 * c2 * exp(-c3 * (Tfreezing - snow->T[l][r][c]));

        //OVERBURDEN
        eta0 = par->snow_viscosity;               //[kg s m^-2]
        c4 = 0.08;                                //[K^-1]
        c5 = 0.023;                               //[m^3 kg^-1]
        load = 0.0;                               //[kg m^-2]
        for(m = l; m <= snow->lnum[r][c]; m++)
        {
            load += (snow->w_ice[m][r][c] + snow->w_liq[m][r][c]);
        }
        load *= fabs(cos(slope * Pi / 180.));
        eta = eta0 * exp(c4 * (Tfreezing - snow->T[l][r][c]) + c5 * (rho_i * theta_i));
        CR2 = -load / eta;

        // Cuell
        snow->Dzl[l][r][c] *= exp( (CR1 + CR2) * Dt );

    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void snow_layer_combination(double a, long r, long c, Statevar3D *snow, double Ta,  GeoVector<long>& inf, double SWEmax_layer, double SWEmax_tot)
{
    long l, linf, n, k, j, max = snow->Dzl.getDh() - 1;
    size_t num;
    double D, Dnew, Dmin, SWE, SWEnew, ice;
    short occurring;
    FILE *f;

    //check on SWEmax[kg/m2]
    SWE = 0.0;
    l = 1;
    while (l <= snow->lnum[r][c])
    {
        SWE += snow->w_ice[l][r][c];
        ice = snow->w_ice[l][r][c];
        snow->w_ice[l][r][c] -= Fmax(0.0, SWE - SWEmax_tot);
        if (snow->w_ice[l][r][c] < 0) snow->w_ice[l][r][c] = 0.0;
        if (ice > 0)
        {
            snow->Dzl[l][r][c] *= snow->w_ice[l][r][c] / ice;
            snow->w_liq[l][r][c] *= snow->w_ice[l][r][c] / ice;
        }
        l++;
    }

    //D=snow depth(mm)
    D = 0.0;
    SWE = 0.0;
    for(l = 1; l <= max; l++)
    {
        D += snow->Dzl[l][r][c];
        SWE += (snow->w_ice[l][r][c] + snow->w_liq[l][r][c]);
    }

    //PREPROCESSING
    //1. If the snow depth is too small, then it is reset to 0
    if(SWE < no_snow * SWEmax_layer)
    {
        snow->lnum[r][c] = 0;
        snow->type[r][c] = 0;
        for(l = 1; l <= max; l++)
        {
            snow->Dzl[l][r][c] = 0.0;
            snow->w_ice[l][r][c] = 0.0;
            snow->w_liq[l][r][c] = 0.0;
            snow->T[l][r][c] = 0.0;               //Temperatura di inizializzazione
        }

        //2. If D<Dsnow_simpl, we are in the simplified case
    }
    else if(snow->lnum[r][c] > 0 && SWE < simpl_snow * SWEmax_layer)
    {

        snow->type[r][c] = 1;
        if(snow->lnum[r][c] > 1)
        {
            for(l = snow->lnum[r][c]; l > 1; l--)
            {
                snowlayer_merging(a, r, c, snow, l, l - 1, l - 1);
            }
            for(l = 2; l <= max; l++)
            {
                snow->T[l][r][c] = 0.0;
                snow->Dzl[l][r][c] = 0.0;
                snow->w_liq[l][r][c] = 0.0;
                snow->w_ice[l][r][c] = 0.0;
            }
            snow->lnum[r][c] = 1;
        }

        //3. if z>=Dsnow_simpl, ordinary case
    }
    else if(snow->lnum[r][c] > 0 && SWE >= simpl_snow * SWEmax_layer)
    {

        snow->type[r][c] = 2;

        //4. if there is not yet a snow layer and D<Dsnow_simpl, simplified case
    }
    else if(snow->lnum[r][c] == 0 && SWE < simpl_snow * SWEmax_layer)
    {

        snow->lnum[r][c] = 1;
        snow->type[r][c] = 1;
        snow->T[1][r][c] = Fmin(Ta, -0.1);

        //5. if there is not yet a snow layer and D>=Dmin(max), simplified case
    }
    else if(snow->lnum[r][c] == 0 && SWE >= simpl_snow * SWEmax_layer)
    {

        snow->lnum[r][c] = 1;
        snow->type[r][c] = 2;
        snow->T[1][r][c] = Fmin(Ta, -0.1);

    }

    // SIMMETRICAL PARAMETERIZATION SCHEME (new)
    if(snow->type[r][c] == 2)
    {

        //remove layers < 0.01 mm
        n = 0;
        do
        {
            occurring = 0;
            if (snow->lnum[r][c] > 1)
            {
                for (l = 1; l <= snow->lnum[r][c]; l++)
                {
                    if(snow->w_ice[l][r][c] < simpl_snow * SWEmax_layer)
                    {
                        merge_layers(a, r, c, snow, l);
                        occurring = 1;
                    }
                }
            }
            n ++;
        }
        while (n <= max && occurring == 1);

        //add new layer
        if (snow->w_ice[snow->lnum[r][c]][r][c] > SWEmax_layer * (1. + simpl_snow) )
        {
            occurring = 1;
        }
        else
        {
            occurring = 0;
        }

        if (occurring == 1)
        {
            if(snow->lnum[r][c] == max)
            {
                linf = 0;
                Dmin = 9.E99;
                for (num = 1; num < inf.size(); num++)
                {
                    if (inf[num] > 0 && inf[num] < max)
                    {
                        if (Dmin > snow->w_ice[inf[num]][r][c] + snow->w_ice[inf[num] + 1][r][c])
                        {
                            Dmin = snow->w_ice[inf[num]][r][c] + snow->w_ice[inf[num] + 1][r][c];
                            linf = inf[num];
                        }
                    }
                    else if (inf[num] > 1 && inf[num] <= max)
                    {
                        if (Dmin > snow->w_ice[inf[num]][r][c] + snow->w_ice[inf[num] - 1][r][c])
                        {
                            Dmin = snow->w_ice[inf[num]][r][c] + snow->w_ice[inf[num] - 1][r][c];
                            linf = -inf[num];
                        }
                    }
                }

                if (linf > 0)
                {
                    snowlayer_merging(a, r, c, snow, linf, linf + 1, linf);
                }
                else if(linf < 0)
                {
                    snowlayer_merging(a, r, c, snow, -linf, -linf - 1, -linf - 1);
                    linf = -linf - 1;
                }
                else
                {
                    f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                    fprintf(f, "r:%ld c:%ld \n", r, c);
                    fprintf(f, "Error in snow combination - Rules to combine layers not applicable\n");
                    fclose(f);
                    t_error("Fatal Error! Geotop is closed. See failing report.");
                }

                for (l = linf + 1; l < snow->lnum[r][c]; l++)
                {
                    snow->w_ice[l][r][c] = snow->w_ice[l + 1][r][c];
                    snow->w_liq[l][r][c] = snow->w_liq[l + 1][r][c];
                    snow->T[l][r][c] = snow->T[l + 1][r][c];
                    snow->Dzl[l][r][c] = snow->Dzl[l + 1][r][c];
                }

                initialize_snow(r, c, snow->lnum[r][c], snow);

                snow->lnum[r][c] --;
            }

            snow->w_liq[snow->lnum[r][c] + 1][r][c] = snow->w_liq[snow->lnum[r][c]][r][c] * (snow->w_ice[snow->lnum[r][c]][r][c] - SWEmax_layer) / snow->w_ice[snow->lnum[r][c]][r][c];
            snow->Dzl[snow->lnum[r][c] + 1][r][c] = snow->Dzl[snow->lnum[r][c]][r][c] * (snow->w_ice[snow->lnum[r][c]][r][c] - SWEmax_layer) / snow->w_ice[snow->lnum[r][c]][r][c];
            snow->w_ice[snow->lnum[r][c] + 1][r][c] = snow->w_ice[snow->lnum[r][c]][r][c] - SWEmax_layer;
            snow->T[snow->lnum[r][c] + 1][r][c] = snow->T[snow->lnum[r][c]][r][c];

            snow->w_ice[snow->lnum[r][c]][r][c] -= snow->w_ice[snow->lnum[r][c] + 1][r][c];
            snow->w_liq[snow->lnum[r][c]][r][c] -= snow->w_liq[snow->lnum[r][c] + 1][r][c];
            snow->Dzl[snow->lnum[r][c]][r][c] -= snow->Dzl[snow->lnum[r][c] + 1][r][c];

            snow->lnum[r][c] ++;
        }

        //split layers
        if (snow->lnum[r][c] < max)
        {

            do
            {

                occurring = 0;

                //FROM UP DOWN
                k = inf.size() - 1;

                do
                {

                    l = inf[k];

                    if (snow->w_ice[l][r][c] > SWEmax_layer * 2.)
                    {

                        occurring = 1;

                        for (j = snow->lnum[r][c]; j > l; j--)
                        {
                            snow->w_ice[j + 1][r][c] = snow->w_ice[j][r][c];
                            snow->Dzl[j + 1][r][c] = snow->Dzl[j][r][c];
                            snow->w_liq[j + 1][r][c] = snow->w_liq[j][r][c];
                            snow->T[j + 1][r][c] = snow->T[j][r][c];
                        }

                        snow->Dzl[l + 1][r][c] = snow->Dzl[l][r][c] * SWEmax_layer / snow->w_ice[l][r][c];
                        snow->w_liq[l + 1][r][c] = snow->w_liq[l][r][c] * SWEmax_layer / snow->w_ice[l][r][c];
                        snow->w_ice[l + 1][r][c] = SWEmax_layer;

                        snow->T[l + 1][r][c] = snow->T[l][r][c];

                        snow->Dzl[l][r][c] -= snow->Dzl[l + 1][r][c];
                        snow->w_liq[l][r][c] -= snow->w_liq[l + 1][r][c];
                        snow->w_ice[l][r][c] -= snow->w_ice[l + 1][r][c];

                        snow->lnum[r][c] ++;
                    }

                    k--;

                }
                while (occurring == 0 && k > 0 );

            }
            while (snow->lnum[r][c] != max && occurring != 0);

        }

        //check
        Dnew = 0.0;
        SWEnew = 0.0;
        for(l = 1; l <= max; l++)
        {
            Dnew += snow->Dzl[l][r][c];
            SWEnew += (snow->w_ice[l][r][c] + snow->w_liq[l][r][c]);
        }

        if(fabs(D - Dnew) > 0.001 || fabs(SWE - SWEnew) > 0.001)
        {
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "r:%ld c:%ld Dold:%f Dnew:%f SWEold:%f SWEnew:%f\n", r, c, D, Dnew, SWE, SWEnew);
            fprintf(f, "Error in snow combination\n");
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report.");
        }

    }
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double DEPTH(long r, long c, GeoMatrix<long>& n, GeoTensor<double>& Dz)
{

    double d = 0.0;
    long l;

    for(l = 1; l <= n[r][c]; l++)
    {
        d += Dz[l][r][c];
    }

    return(d);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

static void snowlayer_merging(double a, long r, long c, Statevar3D *snow, long l1, long l2, long lres)
{

    double h;

    h = internal_energy(snow->w_ice[l1][r][c], snow->w_liq[l1][r][c], snow->T[l1][r][c]) + internal_energy(snow->w_ice[l2][r][c], snow->w_liq[l2][r][c], snow->T[l2][r][c]);
    snow->Dzl[lres][r][c] = snow->Dzl[l1][r][c] + snow->Dzl[l2][r][c];
    snow->w_ice[lres][r][c] = snow->w_ice[l1][r][c] + snow->w_ice[l2][r][c];
    snow->w_liq[lres][r][c] = snow->w_liq[l1][r][c] + snow->w_liq[l2][r][c];
    if(snow->Dzl[lres][r][c] < 0 || snow->w_ice[lres][r][c] < 0 || snow->w_liq[lres][r][c] < 0)
    {
        printf("ERROR 1 in snow layer merging r:%ld c:%ld l1:%ld l2:%ld lres:%ld\n", r, c, l1, l2, lres);
        write_snow_all(r, c, snow);
        t_error("Stop Execution");
    }
    from_internal_energy(a, h, &(snow->w_ice[lres][r][c]), &(snow->w_liq[lres][r][c]), &(snow->T[lres][r][c]));
    if(snow->T[lres][r][c] > 0)
    {
        printf("ERROR 2 in snow layer merging r:%ld c:%ld l1:%ld l2:%ld lres:%ld\n", r, c, l1, l2, lres);
        write_snow_all(r, c, snow);
        t_error("Stop Execution");
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double internal_energy(double w_ice, double w_liq, double T)
{
	
    return (c_ice * w_ice + c_liq*w_liq) * (T - Tfreezing) + Lf * w_liq;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void from_internal_energy(double a, double h, double *w_ice, double *w_liq, double *T)
{

    double SWE = (*w_ice) + (*w_liq);
    double T0;
    double A = -h / (c_ice * SWE), B = c_liq / (c_ice * a * a), C = (Lf * SWE - h) / (c_ice * SWE * a * a);
    unsigned long cont = 0;

    if(SWE > 0)
    {
        *T = 0.0;
        do
        {
            T0 = *T;
            *T = T0 - (pow(T0, 3.0) + A * pow(T0, 2.0) + B * T0 + C) / (3.0 * pow(T0, 2.0) + 2 * A * T0 + B);
            cont++;
            if(cont > 100) printf("%ld %e %e %e %e\n", cont, SWE, h, T0, *T);
        }
        while(fabs(*T - T0) > 1.E-10 && cont < 100);

        *w_liq = theta_snow(a, 1., *T) * SWE;
        *w_ice = SWE - (*w_liq);
    }
    else
    {
        *T = 0.0;
        *w_liq = 0.0;
        *w_ice = 0.0;
    }
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

static void write_snow(long r, long c, long l, Statevar3D *snow)
{

    printf("r:%ld c:%ld wice(%ld/%ld):%f wliq(%ld):%f T(%ld):%f Dz(%ld):%f\n",
           r, c, l, snow->lnum[r][c], snow->w_ice[l][r][c], l, snow->w_liq[l][r][c],
           l, snow->T[l][r][c], l, snow->Dzl[l][r][c]);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_snow_all(long r, long c, Statevar3D *snow)
{
    long l;
    for(l = 1; l <= snow->lnum[r][c]; l++)
    {
        write_snow(r, c, l, snow);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short set_snowice_min(Statevar1D *snow, long l1, long l2, double wicemin)
{

    double f, dwl, dwi, dd;

    //l1 too shallow and takes mass from l2
    if(snow->w_ice[l1] < wicemin && snow->w_ice[l2] > 1.E-6)
    {
        f = Fmin(wicemin - snow->w_ice[l1], snow->w_ice[l2]) / snow->w_ice[l2];
        dd = f * snow->Dzl[l2];
        snow->Dzl[l1] += dd;
        snow->Dzl[l2] -= dd;
        dwl = f * snow->w_liq[l2];
        snow->w_liq[l1] += dwl;
        snow->w_liq[l2] -= dwl;
        dwi = f * snow->w_ice[l2];
        snow->w_ice[l1] += dwi;
        snow->w_ice[l2] -= dwi;
        return 1;
    }
    else
    {
        return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

static void merge_layers(double a, long r, long c, Statevar3D *snow, long l1)
{

    long l;
    FILE *f;

    if(l1 > snow->lnum[r][c])
    {
        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
        fprintf(f, "Error 1 in merge_layers\n");
        fclose(f);
        t_error("Fatal Error! Geotop is closed. See failing report.");
    }

    if(l1 == snow->lnum[r][c])
    {
        snowlayer_merging(a, r, c, snow, l1, l1 - 1, l1 - 1);
        initialize_snow(r, c, l1, snow);
    }
    else
    {
        snowlayer_merging(a, r, c, snow, l1, l1 + 1, l1);
        for(l = l1 + 1; l < snow->lnum[r][c]; l++)
        {
            snow->w_ice[l][r][c] = snow->w_ice[l + 1][r][c];
            snow->w_liq[l][r][c] = snow->w_liq[l + 1][r][c];
            snow->T[l][r][c] = snow->T[l + 1][r][c];
            snow->Dzl[l][r][c] = snow->Dzl[l + 1][r][c];
        }
        initialize_snow(r, c, snow->lnum[r][c], snow);
    }
    snow->lnum[r][c] -= 1;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

static void initialize_snow(long r, long c, long l, Statevar3D *snow)
{

    snow->w_ice[l][r][c] = 0.0;
    snow->w_liq[l][r][c] = 0.0;
    snow->Dzl[l][r][c] = 0.0;
    snow->T[l][r][c] = 0.0;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_snow_age(double Psnow, double Ts, double Dt, double Prestore, double *tsnow_nondim)
{

    double r1, r2, r3;

    //effect of grain growth due to vapour diffusion
    r1 = exp(5000.0 * (1.0 / tk - 1.0 / (Ts + tk)));

    //effect melt and refreezing*/
    r2 = pow(r1, 10);
    if(r2 > 1.0) r2 = 1.0;

    //effect of dirt
    r3 = 0.3;

    //non-dimensional snow age: 10 mm of snow precipitation restore snow age Dt(s)
    *tsnow_nondim = Fmax( 0.0, (*tsnow_nondim + (r1 + r2 + r3) * Dt * 1.0E-6) * (1.0 - Psnow / Prestore) );
    if((*tsnow_nondim) != (*tsnow_nondim)) printf("tsnow no value - tausn:%f P:%f Ts:%f r1:%f r2:%f r3:%f\n", *tsnow_nondim, Psnow, Ts, r1, r2, r3);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double snow_albedo(double ground_alb, double snowD, double AEP, double freshsnow_alb, double C, double tsnow, double cosinc, double ( *F)(const double& x))
{
    double A, Fage = 1.0 - 1.0 / (1.0 + tsnow), w;
    A = freshsnow_alb * (1.0 - C * Fage);
    A += 0.4 * (1.0 - A) * (*F)(cosinc);
    if(snowD < AEP)                               //if snow is shallow (<AEP), interpolate between snow and ground albedo
    {
        w = (1.0 - snowD / AEP) * exp(-snowD * 0.5 / AEP);
        A = w * ground_alb + (1.0 - w) * A;
    }
    return(A);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Fzen(const double& cosinc)
{
    double f, b = 2.0;
    if(cosinc < 0.5)
    {
        f = (1.0 / b) * ((b + 1.0) / (1.0 + 2.0 * b * cosinc) - 1.0);
    }
    else
    {
        f = 0.0;
    }
    return(f);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void non_dimensionalize_snowage(double *snowage, double Ta)
{

    double r1, r2, r3;

    r1 = exp(5000.0 * (1.0 / 273.16 - 1.0 / (Ta + 273.16)));
    r2 = pow(r1, 10);
    if(r2 > 1.0) r2 = 1.0;
    r3 = 0.3;

    *snowage *= ((r1 + r2 + r3) * 1.0E-6);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void WBsnow(double Dt, long ns, long r, long c, Statevar3D *snow, double *Melt, double *RainOnSnow, Par *par, double slope, double Rain, Energy *E, double Evap)
{

    long l, m;
    double Se, th, thi, Wdt, Edt;

    Wdt = Rain;
    Edt = Evap;

    for (l = snow->lnum[r][c]; l >= 1; l--)
    {

        if (l > ns)
        {

            snow->w_ice[l][r][c] = 0.0;
            snow->w_liq[l][r][c] = 0.0;
            snow->Dzl[l][r][c] = 0.0;

        }
        else
        {

            m = ns - l + 1;

            //assign
            if (Edt > E->ice[m] - E->deltaw[m])
            {

                Edt -= (E->ice[m] - E->deltaw[m]);
                Wdt += (E->liq[m] + E->deltaw[m]);
                snow->w_ice[l][r][c] = 0.0;
                snow->w_liq[l][r][c] = 0.0;
                snow->Dzl[l][r][c] = 0;

            }
            else
            {

                snow->T[l][r][c] = E->Temp[m];
                snow->w_ice[l][r][c] = Fmax(0.0, E->ice[m] - E->deltaw[m] - Edt);
                snow->w_liq[l][r][c] = Fmax(0.0, E->liq[m] + E->deltaw[m] + Wdt);
                snow->Dzl[l][r][c] = 1.E3 * E->Dlayer[m];

                Edt = 0.0;
                Wdt = 0.0;

                if (snow->w_ice[l][r][c] > simpl_snow * par->max_weq_snow)
                {

                    //ACCOUNT FOR SNOW COMPACTION
                    //a)destructive metamorphism and overburden
                    snow_compactation(Dt, r, c, l, snow, slope, par);

                    //b)melting: snow depth decreases maintaining the same density
                    if(snow->w_ice[l][r][c] / E->ice[m] < 1) snow->Dzl[l][r][c] *= (snow->w_ice[l][r][c] / E->ice[m]);

                    //limit on max porosity
                    if (snow->w_ice[l][r][c] / (1.E-3 * snow->Dzl[l][r][c]*rho_w) > par->snow_maxpor)
                    {
                        snow->Dzl[l][r][c] = 1.E3 * snow->w_ice[l][r][c] / ( rho_w * par->snow_maxpor );
                    }

                    //CALCULATE LIQUID WATER GOING BELOW
                    //[-]
                    th = snow->w_liq[l][r][c] / (1.0E-3 * snow->Dzl[l][r][c] * rho_w);
                    //[-]
                    thi = snow->w_ice[l][r][c] / (1.0E-3 * snow->Dzl[l][r][c] * rho_i);
                    Se = (th - par->Sr * (1.0 - thi)) / ( (1.0 - thi) - par->Sr * (1.0 - thi));
                    if(Se < 0) Se = 0.0;
                    if(Se > 1) Se = 1.0;
                    if(th > par->Sr * (1.0 - thi)) Wdt += Fmin(5.0 * pow(Se, 3.0) * Dt , (th - par->Sr * (1.0 - thi))) * snow->Dzl[l][r][c] * 1.E-3 * rho_w;
                    snow->w_liq[l][r][c] -= Wdt;

                }
                else
                {

                    Wdt += snow->w_liq[l][r][c];
                    snow->w_liq[l][r][c] = 0.0;

                }

            }
        }
    }

    snow->lnum[r][c] = ns;

    *Melt = Wdt - Rain;

    if (Wdt < Rain)
    {
        *RainOnSnow = Rain;
    }
    else
    {
        *RainOnSnow = 0.0;
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void new_snow(double a, long r, long c, Statevar3D *snow, double P, double Dz, double T)
{

    long ns;
    double h;

    if(snow->type[r][c] == 0)
    {

        snow->Dzl[1][r][c] += Dz;
        snow->w_ice[1][r][c] += P;

    }
    else
    {

        ns = snow->lnum[r][c];
        h = internal_energy(snow->w_ice[ns][r][c], snow->w_liq[ns][r][c], snow->T[ns][r][c]);
        h += (c_ice * P) * (Fmin(T, -0.1) - Tfreezing);

        snow->Dzl[ns][r][c] += Dz;
        snow->w_ice[ns][r][c] += P;

        from_internal_energy(a, h, &(snow->w_ice[ns][r][c]), &(snow->w_liq[ns][r][c]), &(snow->T[ns][r][c]));

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void WBglacier(long ns, long ng, long r, long c, Statevar3D *glac, double *Melt, Par *par, Energy *E, double Evap)
{

    long l, m;
    double th, thi, Se, Edt;

    *Melt = 0;

    Edt = Evap;

    for (l = glac->lnum[r][c]; l >= 1; l--)
    {

        if (l > ng)
        {

            glac->w_ice[l][r][c] = 0.;
            glac->w_liq[l][r][c] = 0.;
            glac->Dzl[l][r][c] = 0.;

        }
        else
        {

            m = ns + ng - l + 1;

            //assign
            if (Edt > E->ice[m] - E->deltaw[m])
            {

                Edt -= (E->ice[m] - E->deltaw[m]);
                *Melt = *Melt +  (E->liq[m] + E->deltaw[m]);

                glac->w_ice[l][r][c] = 0.;
                glac->w_liq[l][r][c] = 0.;
                glac->Dzl[l][r][c] = 0;

            }
            else
            {

                glac->T[l][r][c] = E->Temp[m];
                glac->w_ice[l][r][c] = Fmax(0., E->ice[m] - E->deltaw[m] - Edt);
                glac->w_liq[l][r][c] = Fmax(0., E->liq[m] + E->deltaw[m]);
                glac->Dzl[l][r][c] = 1.E3 * E->Dlayer[m];

                Edt = 0.;

                if (glac->w_ice[l][r][c] > simpl_snow * par->max_weq_glac)
                {

                    //COMPACTION
                    //melting: snow depth decreases maintaining the same density
                    if(glac->w_ice[l][r][c] / E->ice[m] < 1) glac->Dzl[l][r][c] *= (glac->w_ice[l][r][c] / E->ice[m]);

                    //limit on max porosity
                    if (glac->w_ice[l][r][c] / (1.E-3 * glac->Dzl[l][r][c]*rho_w) > 0.95)
                    {
                        glac->Dzl[l][r][c] = 1.E3 * glac->w_ice[l][r][c] / ( rho_w * 0.95 );
                    }

                    //CALCULATE LIQUID WATER GOING BELOW
                    //[-]
                    th = glac->w_liq[l][r][c] / (1.0E-3 * glac->Dzl[l][r][c] * rho_w);
                    //[-]
                    thi = glac->w_ice[l][r][c] / (1.0E-3 * glac->Dzl[l][r][c] * rho_i);
                    Se = (th - par->Sr * (1.0 - thi)) / ( (1.0 - thi) - par->Sr * (1.0 - thi));
                    if(Se < 0) Se = 0.0;
                    if(Se > 1) Se = 1.0;
                    if(th > par->Sr * (1.0 - thi))
                    {
                        *Melt = *Melt + (th - par->Sr * (1.0 - thi)) * glac->Dzl[l][r][c] * 1.E-3 * rho_w;
                        glac->w_liq[l][r][c] = par->Sr * (1.0 - thi) * glac->Dzl[l][r][c] * 1.E-3 * rho_w;
                    }

                }
                else
                {

                    *Melt = *Melt + glac->w_liq[l][r][c];
                    glac->w_liq[l][r][c] = 0.;

                }
            }
        }
    }

    glac->lnum[r][c] = ns;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_SCA(Statevar3D *snow, Par *par, GeoMatrix<double>& Z, double t)
{

    long l, r, c, cont = 0, conttot = 0;
    double T, D, SWE, Tmean = 0.0, Tsmean = 0.0, Dmean = 0.0, SWEmean = 0.0, SCA;

    double JD, JDfrom0;
    long day, month, year, hour, minute;

    FILE *f;

    std::string rec = "_recNNNN" ;
    std::string crec = "_crecNNNN" ;
    std::string name, temp;

    JDfrom0 = convert_tfromstart_JDfrom0(t, par->init_date);
    convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
    convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);

    for(r = 1; r <= geotop::common::Variables::Nr; r++)
    {
        for(c = 1; c <= geotop::common::Variables::Nc; c++)
        {
            if((long)Z[r][c] != geotop::input::gDoubleNoValue)
            {
                D = 0.0;
                T = 0.0;
                SWE = 0.0;
                for(l = 1; l <= snow->lnum[r][c]; l++)
                {
                    D += snow->Dzl[l][r][c];
                    SWE += (snow->w_liq[l][r][c] + snow->w_ice[l][r][c]);
                    T += snow->T[l][r][c] * snow->Dzl[l][r][c];
                }
                if(D > 0)T /= D;

                conttot++;
                if(D > 10)
                {
                    cont++;
                    Dmean += D;
                    SWEmean += SWE;
                    Tmean += T;
                    Tsmean += snow->T[snow->lnum[r][c]][r][c];
                }

            }
        }
    }

    if(cont > 0)
    {
        Dmean /= (double)cont;
        SWEmean /= (double)cont;
        Tmean /= (double)cont;
        Tsmean /= (double)cont;
        SCA = cont / (double)conttot;
    }
    else
    {
        SCA = 0.0;
    }

    if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

    if (par->n_ContRecovery > 0)
    {
        temp = geotop::common::Variables::files[fSCA] ;
        temp += crec;
        name = temp + textfile ;
    }
    else
    {
        name = geotop::common::Variables::files[fSCA] ;
        name += textfile;
    }

    f = fopen(name.c_str() , "a");
    fprintf(f, "%ld/%ld/%ld %ld:%02.0f", day, month, year, hour, (float)minute);
    fprintf(f, ",%f,%f,%f", JDfrom0 - par->init_date, JDfrom0, JD);
    fprintf(f, ",%f,%f,%f,%f,%f,%f\n", Dmean, SWEmean, Tmean, Tsmean, (1.0 - SCA) * 100.0, SCA * 100.0);
    fclose(f);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

static double max_dtheta_snow(double a, double b)
{

    double maxT;

    maxT = -pow(b / (3 * a * a), 0.5);
    return(maxT);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double theta_snow(double a, double b, double T)
{

    double th;

    th = 1. / (b + pow(a * T, 2.));
    if(T > 0) th = 1. / (b + pow(a * 0., 2.));
    return(th);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double dtheta_snow(double a, double b, double T)
{

    double dth;

    dth = -2.*a * a * T / pow(b + pow(a * T, 2.), 2.);
    if(T > 0) dth = -2.*a * a * max_dtheta_snow(a, b) / pow(b + pow(a * max_dtheta_snow(a, b), 2.), 2.);
    return(dth);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void allocate_and_initialize_statevar_3D(Statevar3D *V, double nan, long nl, long nr, long nc)
{
    V->type.resize(nr + 1, nc + 1, 2);
    V->lnum.resize(nr + 1, nc + 1, 0);
    V->Dzl.resize(nl + 1, nr + 1, nc + 1 , 0.);
    V->T.resize(nl + 1, nr + 1, nc + 1, nan);
    V->w_ice.resize(nl + 1, nr + 1, nc + 1, 0.);
    V->w_liq.resize(nl + 1, nr + 1, nc + 1, 0.);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_statevar_3D(Statevar3D *V)
{
	free(V);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void allocate_and_initialize_statevar_1D(Statevar1D *V, double nan, long nl)
{

    V->Dzl.resize(nl + 1, 0.0);
    V->T.resize(nl + 1, nan);
    V->w_ice.resize(nl + 1, 0.0);
    V->w_liq.resize(nl + 1, 0.0);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_statevar_1D(Statevar1D *V)
{
	free(V);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short copy_statevar_from3D_to1D(long r, long c, Statevar3D *origin, Statevar1D *destination)
{

    size_t nl, l;

    nl = origin->Dzl.getDh() - 1;

    if(r < 1 || r > origin->type.getRows()) return 0;
    if(c < 1 || c > origin->type.getCols()) return 0;

    if(nl != destination->Dzl.size()) return 0;

    destination->type = origin->type[r][c];
    destination->lnum = origin->lnum[r][c];
    for (l = 1; l <= nl; l++)
    {
        destination->Dzl[l] = origin->Dzl[l][r][c];
        destination->T[l] = origin->T[l][r][c];
        destination->w_ice[l] = origin->w_ice[l][r][c];
        destination->w_liq[l] = origin->w_liq[l][r][c];
    }

    return 1;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double interpolate_snow(long r, long c, double h, long max, GeoTensor<double>& Dz, GeoTensor<double>& Q, short k)
{

    double q, z, z0 = 0.;
    long l;
    short u;

    q = geotop::input::gDoubleNoValue;

    if (h > 0)                                    //downwards
    {
        u = -1;
        l = max + 1;
    }                                             //upwards
    else
    {
        u = 1;
        l = 1;
        h = -h;
    }

    do
    {

        if (l == 1)
        {
            z = z0 + Dz[l][r][c] / 2.;
        }
        else if (l <= max)
        {
            z = z0 + Dz[l][r][c] / 2. + Dz[l - 1][r][c] / 2.;
        }
        else
        {
            z = z0 + Dz[max][r][c] / 2.;
        }

        if(fabs(h) <= z && fabs(h) > z0)
        {
            if (l == 1)
            {
                if (k == 0)
                {
                    q = Q[l][r][c];
                }
                else
                {
                    q = Q[l][r][c] / Dz[l][r][c];
                }
            }
            else if (l <= max)
            {
                if(u > 0)
                {
                    if (k == 0)
                    {
                        q = ( Q[l - 1][r][c] * (h - z0) + Q[l][r][c] * (z - h) ) / (z - z0);
                    }
                    else
                    {
                        q = ( Q[l - 1][r][c] / Dz[l - 1][r][c] * (h - z0) + Q[l][r][c] / Dz[l][r][c] * (z - h) ) / (z - z0);
                    }
                }
                else
                {
                    if (k == 0)
                    {
                        q = ( Q[l - 1][r][c] * (z - h) + Q[l][r][c] * (h - z0) ) / (z - z0);
                    }
                    else
                    {
                        q = ( Q[l - 1][r][c] / Dz[l - 1][r][c] * (z - h) + Q[l][r][c] / Dz[l][r][c] * (h - z0) ) / (z - z0);

                    }
                }
            }
            else
            {
                if (k == 0)
                {
                    q = Q[max][r][c];
                }
                else
                {
                    q = Q[max][r][c] / Dz[max][r][c];
                }

            }
        }

        z0 = z;

        l += u;

    }
    while ( (long)q == geotop::input::gDoubleNoValue && l <= max + 1 && l >= 1);

    return q;

}

//overloaded function noori
double interpolate_snow(long r, long c, double h, long max, GeoTensor<double>& Dz, GeoTensor<double>& Q)
{

    double q, z, z0 = 0.;
    long l;
    short u;

    q = geotop::input::gDoubleNoValue;

    if (h > 0)                                    //downwards
    {
        u = -1;
        l = max + 1;
    }                                             //upwards
    else
    {
        u = 1;
        l = 1;
        h = -h;
    }

    do
    {

        if (l == 1)
        {
            z = z0 + Dz[l][r][c] / 2.;
        }
        else if (l <= max)
        {
            z = z0 + Dz[l][r][c] / 2. + Dz[l - 1][r][c] / 2.;
        }
        else
        {
            z = z0 + Dz[max][r][c] / 2.;
        }

        if(fabs(h) <= z && fabs(h) > z0)
        {
            if (l == 1)
            {
                q = Q[l][r][c];
            }
            else if (l <= max)
            {
                if(u > 0)
                {
                    q = ( Q[l - 1][r][c] * (h - z0) + Q[l][r][c] * (z - h) ) / (z - z0);
                }
                else
                {
                    q = ( Q[l - 1][r][c] * (z - h) + Q[l][r][c] * (h - z0) ) / (z - z0);
                }
            }
            else
            {
                q = Q[max][r][c];
            }
        }

        z0 = z;

        l += u;

    }
    while ( (long)q == geotop::input::gDoubleNoValue && l <= max + 1 && l >= 1);

    return q;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

