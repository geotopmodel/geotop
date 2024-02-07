
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
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "util_math.h"
#include "air.energy.balance.h"
#include "water.balance.h"
#include "meteodata.h"
#include "logger.h"
#include "timer.h"
#include <ctime>
#include "math.optim.h"
#include "snow.h"

extern long number_novalue;

extern T_INIT *UV;
extern char *logfile;

extern long Nl, Nr, Nc;
extern double *odb;
extern double t_sub, t_sup;

extern char *FailedRunFile;

extern long i_sim, i_run;
extern double MM1, MM2, MMR, MMo, MS1, MS2;

//subsurface flow constants
#define tol_max_GC 1.E+5
#define tol_min_GC 1.E-13
#define max_res_adm 1.E-2
#define MM 1
#define ni 1.E-7
#define maxITER_rec_K 10
#define FixedBottomBoundary 0.0
#define FixedTopBoundary 0.0

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


short air_energy_balance(double Dt, double JD0, double JDb, double JDe,
                    SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, STATEVAR_3D *G, STATE_VEG *V,
                    Vector<double> *snowage, ALLDATA *adt, double *W)
{
   
    GEOLOG_PREFIX(__func__);
    GEOTIMER_SECTION(__func__);

    clock_t start, end;
    double Pnet, loss;
    long j;
    short a;



    start = clock();

    // GZ REVIEW-- Check if UpdateK is neccesary or not
    
    a = Energy3D(Dt, L, C, S, adt, &loss, adt->P->UpdateK);
    end=clock();
    t_sub += (end-start)/(double)CLOCKS_PER_SEC;
    if (a != 0)
    {
        return 1;
    }



    return 0;


}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/** The code solves non-linear system F(H), where H is the total head (piezometric) [mm], F is vectorial function with components equal to the number of points,
 *
 * H = P + z, where P is the pressure [mm] and z the elevation over a reference height (given by the average elevation of the domain)
 * F(H) comes out from the temporal discretization of Richards' equation PDE
 *
 * In this code: F(H) = I*f(H) + K(H)*H,
 * where I is the identity matrix, f(H) = volume storage + sink, is a vectorial function (not a matrix), and K(H) is the hydraulic conductivity matrix
 *
 * K(H) has a number of interesting properties:
 * 1) is symmetrical
 * 2) if kij is a component kii = sum(j!=i) kji
 * 3) the sum of all the components of the matrix is 0
 * 4) if T is a vector, the vector K T has the properties that the sum of its components is 0
 *
 * Here it is possible to update or not K with the new values of H.
 * In any case the derivative of K(H) is not considered in the Jacobian. So the matrix K(H) is also used in the JACOBIAN.
 *
 * K is described storing only its strict lower component (strict = without diagonal) with the 3 vectors Li, Lp, Lx (in the same way as UFMPACK) */

short Energy3D(double Dt, SOIL_STATE *L, SOIL_STATE *C,STATEVAR_3D *S, ALLDATA *adt, double *loss, short updateK)
                 
                 
{
    GEOLOG_PREFIX(__func__);
    GEOTIMER_SECTION(__func__);

    double res=0.0, res0[3], res_prev[MM], res_av, res00, lambda[3], epsilon,
            mu=0., hnew, hold=0.;
    double ds=sqrt((*UV->U)(1) * (*UV->U)(2)), area, dz, dn, dD,Tairmet,ns;
    double snowD,ThetaIce,ThetaWat,sat,ThetaAir1,ThetaAir0,totw0;
    long i, j, ch, l, r, c, m, bc, sy, cont, cont2, iter;
    long n=adt->T->lrc_cont->nrh;
    long N=adt->W->H0->nh;
    long cont_lambda_min=0;
    short out, out2;
    int sux;


    /* initialize res0 */
    res0[0]=0;
    res0[1]=0;
    res0[2]=0;
    lambda[0]=0;
    lambda[1]=0;
    lambda[2]=0;
    /* END  initialize res0 */ // ec 2012 08 24


    /*
     Layer 0 is water on the surface. This allows a more robust description of infiltration processes.
     However, surface water flow is described separately in the subroutine supflow

     Layers > 0 , represent subsurface
     */

	
     /** solution guess */
    for (i=1; i<=N; i++)
    {

        if (i<=n)
        {

            l = (*adt->T->lrc_cont)(i,1);
            r = (*adt->T->lrc_cont)(i,2);
            c = (*adt->T->lrc_cont)(i,3);
            j = adt->T->j_cont[r][c];
            Tairmet=(*adt->M->Tgrid)(r,c)+273.15;
            //ns=(*S->lnum)(r,c);
            //printf("i:%ld Layer:%ld  \n",i,l);
            
            if (l==0)
            {
                snowD=DEPTH(r, c, adt->N->S->lnum.get(),adt->N->S->Dzl.get());
				if (snowD<=adt->P->SnowDepthAirFlowLimit)
				{
                    (*adt->AE->T1)(i) = Tairmet;
                 }
                 else
                 {
                    (*adt->AE->T1)(i) =(*adt->N->S->T)(1,r,c)+273.15;

                 }
            }
            else
            {
                (*adt->AE->T1)(i) = (*adt->S->SS->Tair_internal)(l,j)+273.15;
            }
        }
        else
        {

            l = (*adt->C->lch)(i-n,1);
            ch = (*adt->C->lch)(i-n,2);
            r = (*adt->C->r)(ch);
            c = (*adt->C->c)(ch);
            Tairmet=(*adt->M->Tgrid)(r,c)+273.15;
            ns=(*S->lnum)(r,c);
            if (l==0)
            {
                snowD=DEPTH(r, c, adt->N->S->lnum.get(),adt->N->S->Dzl.get());
				if (snowD<=adt->P->SnowDepthAirFlowLimit)
                 {
                    (*adt->AE->T1)(i) = Tairmet;
                 }
                 else // IF there is snow, the boundary condition is snow temperature.
                 {
                    (*adt->AE->T1)(i) =(*adt->N->S->T)(1,r,c)+273.15;
                 }
            }
            else
            {
                (*adt->AE->T1)(i) = (*adt->C->SS->Tair_internal)(l,ch)+273.15;
            }
        }
    }

    sux = find_matrix_Kthermal_3D(Dt, L,C,S,adt->AE->Lx.get(), adt, adt->AE->T1.get());
                           
    find_fthermal_3D(Dt, adt->AE->f.get(), adt, L, C,S, adt->AE->T1.get());
              
                           
    // GZ - This function complete the term B (B= final f). B=f+LX*deltaP
    product_matrix_using_lower_part_by_vector_plus_vector2(-1., adt->AE->B.get(), adt->AE->f.get(), adt->AE->T1.get(),
                                                          adt->T->Li.get(), adt->T->Lp.get(), adt->AE->Lx.get(), adt->AF->LxJair.get(),adt->T->lrc_cont.get(),adt->C->lch.get(),n);

    res = norm_inf(adt->AE->B.get(), 1, N);

    res00 = res; /** initial norm of the residual */
    epsilon= (adt->P->AirRichardTol+ adt->P->RelTolVWb * std::min<double>( res00, sqrt((double)N) ))*100;
    //epsilon = (adt->P->tol_energy + adt->P->RelTolVWb * std::min<double>( res00, sqrt((double)N)));
    //printf("AIR ENERGY EPSILON: %e\n",epsilon);
    cont=0;
    out=0;

    //printf("res:%e\n",res);
    //The condition to go out from the Newton cycle is made on the norm of the residual < Absolute tolerance + Relative tolerance * res00
    if ( res <= std::min<double>(epsilon, max_res_adm) )
        out=1;

    //Max iteration number
    if ( cont >= adt->P->MaxiterTol )
        out=1;

    //printf("AIRENERGY out:%ld\n",out);
    
    while (out==0)
    {

        cont++;

        for (i=1; i<=N; i++)
        {
            (*adt->AE->T0)(i) = (*adt->AE->T1)(i);
            (*adt->AE->dT)(i) = 0.0;
        }

        //mu is the forcing term of Newton Raphson, defined such that ||F( x + d )|| < mu * ||F(x)|| see references
        if (cont==1)
        {
            mu = adt->P->TolCG;
        }
        else
        {
            mu *= std::min<double>( 1.0, res/res0[0] );
            if (mu < 0.5*epsilon/res) mu = 0.5*epsilon/res;
        }

        /** CALCOLATE AND STORE JACOBIAN AS SPARSE MATRIX */
        //printf("TEST AIR ENERGY BALANCE BEFORE BICGSTAB\n");
		//GZ - This term is =0 for the air. Only the terms associated to K exist.f
        sux = find_dfdHthermal_3D(Dt, adt->AE->df.get(), adt, L, C, adt->AE->T1.get(), adt->AE->Klat.get());
        
        /** it calculates only df/dH, J = I*df/dH + K, K is calculated above **/

        /** CONJUGATED GRADIENTS ALGORITHM */
        iter = BiCGSTAB_strict_lower_matrix_plus_identity_by_vector2(mu, tol_min_GC, tol_max_GC, adt->AE->dT.get(),
                                                                    adt->AE->B.get(), adt->AE->df.get(), adt->T->Li.get(),
                                                                    adt->T->Lp.get(), adt->AE->Lx.get(),adt->AF->LxJair.get(),adt->T->lrc_cont.get(),adt->C->lch.get(),n);

        //printf("AIRENERGY T1(1):%f\n",(*adt->AE->T1)(1));
        if (iter==-1){
            printf("AIR ENERGY S1\n");
            return 1; /** does not converge */}

        //non-monotonic line search (it is monotonic if M==1)
        for (m=std::min<long>(cont,MM); m>1; m--)
        {
            res_prev[m-1]=res_prev[m-2];
        }
        res_prev[0]=res;

        res_av=0.0;
        for (m=1; m<=std::min<long>(cont,MM); m++)
        {
            res_av=std::max<double>(res_prev[m-1],res_av);
        }
        cont2 = 0.0;
        res0[0] = res;

        do
        {

            cont2++;

            /** The damping factor (or path length) lambda, defined as H1 = H0 + lambda*dH is chosen by minimizing the merit function :
             * 0.5*norm(F_water(H0+lambda*dH)) interpolated with a quadratic polynome.
             * This method could not always be suited, because it has the disadvantage that it can make the solution stagnate to a point.
             * A relatively low number of MaxiterCorr (around 3-5) can prevent stagnation */

            if (cont2 == 1)
            {
                lambda[0] = 1.0;

            }
            else if (cont2 == 2)
            {
                lambda[1] = lambda[0];
                res0[1] = res;
                lambda[0] = GTConst::thmax;

            }
            else
            {
                lambda[2] = lambda[1];
                res0[2] = res0[1];
                lambda[1] = lambda[0];
                res0[1] = res;
                lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1], lambda[2], res0[2]);
            }

            for (i=1; i<=N; i++)
            {
                (*adt->AE->T1)(i) = (*adt->AE->T0)(i) + lambda[0] * (*adt->AE->dT)(i);
                if ((*adt->AE->T1)(i) != (*adt->AE->T1)(i)){
                    printf("AIR ENERGY S2\n");
                    return 1;}
            }

            if (updateK == 1 && cont <= maxITER_rec_K)
                sux = find_matrix_Kthermal_3D(Dt, L,C,S,adt->AE->Lx.get(), adt, adt->AE->T1.get());

                           
            find_fthermal_3D(Dt, adt->AE->f.get(), adt, L, C,S, adt->AE->T1.get());
                               
            // GZ - This function complete the term B (B= final f). B=f+LX*deltaP
            product_matrix_using_lower_part_by_vector_plus_vector2(-1., adt->AE->B.get(), adt->AE->f.get(), adt->AE->T1.get(),
                                                              adt->T->Li.get(), adt->T->Lp.get(), adt->AE->Lx.get(), adt->AF->LxJair.get(),adt->T->lrc_cont.get(),adt->C->lch.get(),n);
            res = norm_inf(adt->AE->B.get(), 1, N);
            //printf("..res:%e\n",res);

            out2=0;
            
            geolog << "cnt:" << cont << " res:" << res << " lambda:" << lambda[0]
                   << " Dt:" << Dt << std::endl;
            
            
            if (res <= (1.0 - ni*lambda[0]*(1.-mu))*res_av)
                out2=1;
            if (lambda[0] <= adt->P->min_lambda_wat)
                cont_lambda_min++;

            if (cont_lambda_min > adt->P->max_times_min_lambda_wat)
            {
                if (adt->P->exit_lambda_min_wat == 1)
                {
                    printf("AIR ENERGY S3, cont_lambda_min %d, %d\n",cont_lambda_min,adt->P->max_times_min_lambda_wat);
                    return 1;
                }
                else
                {
                    out2=1;
                    cont_lambda_min=0;
                    printf("AIR ENERGY S3, cont_lambda_min %d, %d\n",cont_lambda_min,adt->P->max_times_min_lambda_wat);
                }
            }
        }
        while (out2==0);

        out=0;
        /** The condition to go out from the Newton cycle is:
         * norm of the residual < Absolute tolerance + Relative tolerance * res00 **/
        if ( res <= std::min<double>( epsilon, max_res_adm ) )
            out=1;
        //Max iteration number
        if ( cont >= adt->P->MaxiterTol )
            out=1;

    }

    if ( res > epsilon )  {
    printf("AIR ENERGY S4\n"); 
    return 1;
    }

    /** it can be shown that massloss per unit pixel [mm] is the linear norm of B * Dt / total_pixel */
    *loss = norm_1(adt->AE->B.get(), 1, N)*Dt/adt->P->total_area;

    /** assign updated state variables */
    //printf("---------xxxxxxxxxxxxxxx--------AIR ENERGYY READYYYYY------------xxxxxxxxxxxxxx-----------\n");
    
    for (i=1; i<=N; i++)
    {

        if (i<=n)  //land
        {
            l = (*adt->T->lrc_cont)(i,1);
            r = (*adt->T->lrc_cont)(i,2);
            c = (*adt->T->lrc_cont)(i,3);
            j = adt->T->j_cont[r][c];
            sy = (*adt->S->type)(r,c);

            if (l>0) {
            
			ThetaIce=(*L->thi)(l,j);
			ThetaWat=(*adt->S->th)(l,j);
			sat=(*adt->S->pa)(sy,jsat,l);
			ThetaAir1=sat-ThetaIce-ThetaWat;
			
			totw0=(*L->totw0)(l,j);
			ThetaAir0=sat-totw0;
			
            //(*L->totw0)(l,j)=(*adt->S->th)(l,j)+(*L->thi)(l,j);
            (*L->Tair_internal)(l,j)=(*adt->AE->T1)(i)-273.15; //Convert from Kelvin to Celsius
            //printf("AIRENERGY i,l,ns,ThetaAir0,ThetaAir1,h,TS,TCM,Tmet,T1,T0,MaxVel,P1:, %i, %d, %d, %e, %e, %e, %f, %f, %f, %f, %f, %e, %f \n",i,l,(*S->lnum)(r,c),ThetaAir0,ThetaAir1,(*L->Hair_internal)(l,j),(*S->T)(1,r,c)+273.15,(*L->T)(l,j)+273.15,(*adt->M->Tgrid)(r,c)+273.15,(*adt->AE->T1)(i),(*L->Tair_internal)(l,j)+273.15,(*adt->AF->MaxVel)(i),(*adt->AF->P1)(i));
            if ( (r==7) and (c==7)){
            //printf("AIRENERGY i,l,ns,ThetaAir0,ThetaAir1,h,TS,TCM,Tmet,T1,T0,MaxVel,P1,T1-T0:, %i, %d, %d, %e, %e, %e, %f, %f, %f, %f, %f, %e, %f, %f \n",i,l,(*S->lnum)(r,c),ThetaAir0,ThetaAir1,(*L->Hair_internal)(l,j),(*S->T)(1,r,c)+273.15,(*L->T)(l,j)+273.15,(*adt->M->Tgrid)(r,c)+273.15,(*adt->AE->T1)(i),(*adt->S->SS->Tair_internal)(l,j)+273.15,(*adt->AF->MaxVel)(i),(*adt->AF->P1)(i),(*adt->AE->T1)(i)-((*adt->S->SS->Tair_internal)(l,j)+273.15));
            }
            }
            else{
            ThetaAir1=0;
            ThetaAir0=0;
            //printf("AIRENERGY i,l,ns,ThetaAir0,ThetaAir1,TS,Tmet,T1: %i, %d, %d, %e, %e,%f, %f, %f\n",i,l,(*S->lnum)(r,c),ThetaAir0,ThetaAir1,(*S->T)(1,r,c)+273.15,(*adt->M->Tgrid)(r,c)+273.15,(*adt->AE->T1)(i));
            
            }
            
        }
        else   /** channel */
        {
            l = (*adt->C->lch)(i-n,1);
            ch = (*adt->C->lch)(i-n,2);
            if (l>0) {
            (*C->Tair_internal)(l,ch)=(*adt->AE->T1)(i) -273.15;
            //(*C->totw0)(l,ch)=(*adt->C->th)(l,ch)+(*C->thi)(l,ch);
            } //Convert from Kelvin to Celsius}
        }
    }
    
    //PRINT SOME VARIABLES
    /*
    for (i=1; i<=Nl; i++)
    {

        if (i<=n)  //land
        {
            l = (*adt->T->lrc_cont)(i,1);
            r = (*adt->T->lrc_cont)(i,2);
            c = (*adt->T->lrc_cont)(i,3);
            j = adt->T->j_cont[r][c];
            sy = (*adt->S->type)(r,c);

            if (l>0) {

			ThetaIce=(*L->thi)(l,j);
			ThetaWat=(*adt->S->th)(l,j);
			sat=(*adt->S->pa)(sy,jsat,l);
			ThetaAir1=sat-ThetaIce-ThetaWat;
			
			totw0=(*L->totw0)(l,j);
			ThetaAir0=sat-totw0;
            //printf("AIRENERGY i,l,ns,ThetaAir0,ThetaAir1,TS,TCM,Tmet,T1,T0,Tz: %i, %d, %d, %e, %e, %f, %f, %f, %f, %f, %f\n",i,l,(*S->lnum)(r,c),ThetaAir0,ThetaAir1,(*S->T)(1,r,c)+273.15,(*L->T)(l,j)+273.15,(*adt->M->Tgrid)(r,c)+273.15,(*adt->AE->T1)(i),(*L->Tair_internal)(l,j)+273.15,(*adt->T->Z)(l,r,c)/1000);
            printf("AIRENERGY i,l,ns,ThetaAir0,ThetaAir1,h,TS,TCM,Tmet,T1,T0,MaxVel,P1:, %i, %d, %d, %e, %e, %e, %f, %f, %f, %f, %f, %e, %f \n",i,l,(*S->lnum)(r,c),ThetaAir0,ThetaAir1,(*L->Hair_internal)(l,j),(*S->T)(1,r,c)+273.15,(*L->T)(l,j)+273.15,(*adt->M->Tgrid)(r,c)+273.15,(*adt->AE->T1)(i),(*L->Tair_internal)(l,j)+273.15,(*adt->AF->MaxVel)(i),(*adt->AF->P1)(i));
            
            
            }
            else{
            ThetaAir1=0;
            ThetaAir0=0;
            printf("AIRENERGY i,l,ns,ThetaAir0,ThetaAir1,TS,Tmet,T1: %i, %d, %d, %e, %e,%f, %f, %f\n",i,l,(*S->lnum)(r,c),ThetaAir0,ThetaAir1,(*S->T)(1,r,c)+273.15,(*adt->M->Tgrid)(r,c)+273.15,(*adt->AE->T1)(i));
            
            }
            
        }
    }*/
    
    
            //printf("i:%ld Layer:%ld  \n",i,l);
            
    // GZ -- Necessary to check if we need to update variables!
     // sux = find_matrix_LxHeat_air_3D(Dt, L, C,adt->AR->LxJair.get(), adt->AF->Klat.get(), adt->AF->Kbottom.get(),
      //                       adt->C->Kbottom.get(), adt, adt->AF->P1.get());
    adt->AE->iternumber=cont;

    return 0;

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double TransferTestConstant(short HeatTransferModel, double Vel, double dp, double kthcm,
                 double Ht_a, double Ht_b, double Ht_n, double ThetaAir1,double printaux)
{

    /*
     
     Heat transfer constant (h) bewteen the air and the composite medium.
     
     if HeatTransferModel=1.0
     Equations presented by Leontini 2015 - Modelling Approaches to Natural Convection in Porouse Media (page 15)
     * 
     if HeatTransferModel=2.0 
     Equations presented by Nield 2017 - Convection in porous media - https://link.springer.com/book/10.1007/978-3-319-49562-0 (page 43)
     * 
     if HeatTransferModel=3.0 
     Equations presented by Yang, Jian 2012 - Experimental analysis of forced convective heat transfer in novel structured packed beds of particles
     * 
     * 
     * 
     if HeatTransferModel=4.0 
     Similar to Local thermal equilibirum (LTE) -> h=hmax
     * 

     */
    double h=0.;
    double Rep,Nufs,hast;
    double bet=10.0;
    double Asf,Vf,r,dh;
    double hmax=100.0;
    

    if (HeatTransferModel==1.0 )
    {
        Rep=Vel*dp/GTConst::nu_air;
		Nufs=Ht_a+Ht_b*pow(GTConst::Pr_air,0.333333)*pow(Rep,Ht_n);
		hast=Nufs*GTConst::k_air/dp;
		r=dp/2.0;
		Asf=4.0*GTConst::Pi*pow(r,2.0);
		Vf=4.0/3.0*GTConst::Pi*pow(r,3.0);
		h=Asf*(1.0-ThetaAir1)/Vf*hast;
		if (printaux==1.0)
		{
		printf("AIRENERGY, Rep,Nufs,hast,Asf,Vf,h,GTConst::Pi,r,pow(r3) =  %f, %f , %f, %f, %f , %f, %f, %f , %f  \n", Rep,Nufs,hast,Asf,Vf,h,GTConst::Pi,r,pow(r,3));
				
		}

    }
    
    if (HeatTransferModel==2.0 )
    {
		Rep=Vel*dp/GTConst::nu_air;
		Nufs=Ht_a+Ht_b*pow(GTConst::Pr_air,1.0/3.0)*pow(Rep,Ht_n);
		hast=1.0/(dp/(Nufs*GTConst::k_air)+dp/(bet*kthcm));
		h=6.0*(1.0-ThetaAir1)/dp*hast;
    }
    
    if (HeatTransferModel==3.0 )
    {
        
		r=dp/2.0;
		Asf=4.0*GTConst::Pi*pow(r,2.0);
		Vf=4.0/3.0*GTConst::Pi*pow(r,3.0);
		dh=4.0*ThetaAir1/(1.0-ThetaAir1)*Vf/Asf;
		
		Rep=Vel*dp/GTConst::nu_air;
		Nufs=Ht_a+Ht_b*pow(GTConst::Pr_air,0.333333)*pow(Rep,Ht_n)*pow(dp*ThetaAir1/dh,Ht_n);
		hast=Nufs*GTConst::k_air/dp;
		
		h=Asf*hast/Vf;
		if (printaux==1.0)
		{
		printf("HT=3 AIRENERGY, Rep,Ht_a,Ht_b,Ht_n,pow(dp*ThetaAir1/dh,Ht_n),Nufs,hast,ThetaAir1,dh,h=  %f, %f , %f, %f, %f,%f , %f, %f, %f , %f  \n", Rep,Ht_a,Ht_b,Ht_n,pow(dp*ThetaAir1/dh,Ht_n),Nufs,hast,ThetaAir1,dh,h);
				
		}

    }
    
    if (HeatTransferModel==4.0 )
    {
		h=hmax;
    }
    
    h=std::min<double>(h, hmax);


    return (h);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double k_thermalcm2b(short snow, short a, double th_liq, double th_ice,
                 double th_sat, double k_solid)
{

    /*
     K thermal for composite medium; air is not considered within the soil
     for soil:
     Quadratic parallel based on the analogy with electric lattice
     (Cosenza et al., 2003, European Journal of Soil Science, September 2003, 54, 581–587)

    if a==4 Faroku- De Vries 1981
    
    Faroku- De Vries 1981
    f_k refers to the weighting factors of de Vries (usually refered as k)
    g_k to the shape factores
     */
    double r, k=0.;
    
    double f_air=0.,f_ice=0.,f_solids=0.;
    double g_air=0.,g_ice=0.125,g_solids=0.125;
    double th_air=th_sat-th_ice-th_liq;

    if (snow==0 || (a!=2 && a!=3))
    {
        
        if (a==1){
        k = pow ( (1.-th_sat)/((1.0-th_sat)+th_liq+th_ice)*sqrt(k_solid) + 
        th_liq/((1-th_sat)+th_liq+th_ice)*sqrt(GTConst::k_liq) +
         th_ice/((1-th_sat)+th_liq+th_ice)*sqrt(
                GTConst::k_ice) , 2. ) ;
                
        }
                
       if (a==4){
        
       if (th_liq>0.09){
        g_air=0.333-(0.333-0.035)*th_air/th_sat;
        }
        else {
        g_air=0.013+0.944*th_liq;
        }       
        
        f_air=2./3.*pow(1.+(GTConst::k_air/GTConst::k_liq-1)*g_air,-1.0)+1./3.*pow(1.+(GTConst::k_air/GTConst::k_liq-1.)*(1.-2.*g_air),-1.0);
        f_ice=2./3.*pow(1.+(GTConst::k_ice/GTConst::k_liq-1)*g_ice,-1.0)+1./3.*pow(1.+(GTConst::k_ice/GTConst::k_liq-1.)*(1.-2.*g_ice),-1.0);
        f_solids=2./3.*pow(1.+(k_solid/GTConst::k_liq-1)*g_solids,-1.0)+1./3.*pow(1.+(k_solid/GTConst::k_liq-1.)*(1.-2.*g_solids),-1.0);
        
       k=(th_liq*GTConst::k_liq+f_ice*th_ice*GTConst::k_ice+f_air*th_air*GTConst::k_air +
        f_solids*(1.-th_sat)*k_solid)/(th_liq+th_ice*f_ice+th_air*f_air+(1.-th_sat)*f_solids);
        
        }

    }
    else
    {

        r = th_liq * GTConst::rho_w + th_ice * GTConst::rho_i ;

        if (a==2)
        {
            if (r < 156)
            {
                k = 0.023 + 0.234 * r*1.E-3 ;
            }
            else
            {
                k = 0.138 - 1.01 * r*1.E-3 + 3.233 * r*r*1.E-6 ;
            }
        }
        else if (a==3)
        {
            k = GTConst::k_air + (7.75E-5 * r + 1.105E-6 * r*r) * (GTConst::k_ice-GTConst::k_air) ;
        }

    }

    return (k);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double k_thermalsnow(short snow, short a, double th_liq, double th_ice,
                 double th_sat, double k_solid)
{

    /*
     for soil:
     Quadratic parallel based on the analogy with electric lattice
     (Cosenza et al., 2003, European Journal of Soil Science, September 2003, 54, 581–587)

     for snow:
     if a==1 Cosenza
     if a==2 Sturm(1997)
     if a==3 Jordan(1991)
     */

    double r, k=0.;

    if (snow==0 || (a!=2 && a!=3))
    {
        k = pow ( (1.-th_sat)*sqrt(k_solid) + th_liq*sqrt(GTConst::k_liq) + th_ice*sqrt(
                GTConst::k_ice) + (th_sat-th_liq-th_ice)*sqrt(GTConst::k_air), 2. ) ;

    }
    else
    {

        r = th_liq * GTConst::rho_w + th_ice * GTConst::rho_i ;

        if (a==2)
        {
            if (r < 156)
            {
                k = 0.023 + 0.234 * r*1.E-3 ;
            }
            else
            {
                k = 0.138 - 1.01 * r*1.E-3 + 3.233 * r*r*1.E-6 ;
            }
        }
        else if (a==3)
        {
            k = GTConst::k_air + (7.75E-5 * r + 1.105E-6 * r*r) * (GTConst::k_ice-GTConst::k_air) ;
        }

    }

    return (k);

}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_matrix_Kthermal_3D(double  /*Dt*/, SOIL_STATE *SL, SOIL_STATE *SC,STATEVAR_3D *S,
                     Vector<double> *Lx, ALLDATA *adt, Vector<double> *T)
{
    GEOTIMER_SECTION(__func__);

    long i, l, r, c, j, I, R, C, J, sy, syn, ch, cnt=0;
    long n=(Nl+1)*adt->P->total_pixel;
    double dz=0.0, dzn=0.0, dD=0.0, k=0.0, kn=0.0, kmax=0.0, kmaxn=0.0;
    double area, ds=sqrt((*UV->U)(1)*(*UV->U)(2)), dn,ns;
    double ThetaIce,ThetaWat,sat,ThetaAir;
    double ThetaIce2,ThetaWat2,sat2,ThetaAir2;// GZ - This variables are calculated in the cell where where there is flux exchange

    //double psi, ice, a, ns, res, sat, ss, Temp;


    /**GZ REVIEW - General comments. Necessary to check :
    -Check Units
    -Check border conditions
    **/
    for (i=1; i<=T->nh; i++)
    {
        /** VERTICAL FLUXES */
        if ( i<=n) //land
        {

            l=(*adt->T->lrc_cont)(i,1);
            r=(*adt->T->lrc_cont)(i,2);
            c=(*adt->T->lrc_cont)(i,3);
            j=adt->T->j_cont[r][c];
            sy=(*adt->S->type)(r,c);
			k = GTConst::k_air;
            
            ch=(*adt->C->ch)(r,c);
            area=ds*ds/cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0);
            if (ch>0) area-=(*adt->C->length)(ch) * adt->P->w_dx *
                            ds; //area of the pixel[m2]

            /** vertical flux */
            if (l>0)
            {
                            
                // GZ - Evaluate if create variable ThetaAir inside AirEnergy
                ThetaIce=(*adt->S->SS->thi)(l,j);
                ThetaWat=(*adt->S->th)(l,j); // This variable is updated at the end of the water Balance line 500
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;
                dz = (*adt->S->pa)(sy,jdz,l)/1000.0;
            }

            /** flux from cell below */
            if (l<Nl)
            {

                I = i+1;

                if (l==0) /** overland flow */
                {
                    /*
                    ThetaIce2=(*adt->S->SS->thi)(l+1,j);
                    ThetaWat2=(*adt->S->th)(l+1,j);
                    sat2=(*adt->S->pa)(sy,jsat,l+1);
                    ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                    
                    dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
                    dD = 0.5*dzn;
                    kn = GTConst::k_air;
                    cnt++;
                    if ((ThetaAir2>GTConst::ThetaAirMin)){
                    (*Lx)(cnt) = (-kn*(ThetaAir2)/1.0)*area/dD*1000.0;
                     }
					else{
                    (*Lx)(cnt) =0.0;
					}*/
					
					cnt++;
					(*Lx)(cnt) =0.0;
					
                }
                else    /** subsurface flow */
                {

                    ThetaIce2=(*adt->S->SS->thi)(l+1,j);
                    ThetaWat2=(*adt->S->th)(l+1,j);
                    sat2=(*adt->S->pa)(sy,jsat,l+1);
                    ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                    
                    dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
                    dD = 0.5*dz + 0.5*dzn;
                    kn = GTConst::k_air;
                    cnt++;
                    if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                    (*Lx)(cnt) = (-kn*(ThetaAir+ThetaAir2)/2.0)*area/dD*1000.0; 
                    //printf("AIR flux\n");
                    }
                    else{
                        (*Lx)(cnt) =0.0;
                        //printf("AIR NOT flux\n");
                        }
                    //printf("AIR ENERGY LXKAIR %e \n",(*adt->AF->LxJair)(cnt));
                    //printf("AIR ENERGY VERTICAL i,I,l,Kn,ThetaAir,ThetaAir2,LX: %i,%i, %d ,%f %e %e %e\n",i,I,l,kn,ThetaAir,ThetaAir2,(*Lx)(cnt));
                    //printf("AIR ENERGY VERT i,l,area,dz,ThetaAir,ThetaAir2,LX: %i, %d ,%f ,%f , %f %f %e\n",i,l,area,dz,ThetaAir,ThetaAir2,(*Lx)(cnt));
                }
                
            }  

        }
        else  /** channel */
        {

            l=(*adt->C->lch)(i-n,1);
            ch=(*adt->C->lch)(i-n,2);
            r=(*adt->C->r)(ch);
            c=(*adt->C->c)(ch);
            sy=(*adt->C->soil_type)(ch);
            

            area=(*adt->C->length)(ch) * adt->P->w_dx * ds;

            /** vertical hydraulic conductivity */
            if (l>0)
            {
               ThetaIce=(*adt->C->SS->thi)(l,ch);
               ThetaWat=(*adt->C->th)(l,ch);
               sat=(*adt->S->pa)(sy,jsat,l);
               ThetaAir=sat-ThetaIce-ThetaWat;
            
                dz = (*adt->S->pa)(sy,jdz,l)/1000.0;
            }

            /** flux from cell below */
            if (l<Nl)
            {

                I = i+1;

                if (l==0) /** OVERLAND flow */
                {

                    cnt++;
                    (*Lx)(cnt) =0.0;


                }
                else    /** SUBSURFACE flow */
                {

                    ThetaIce2=(*adt->C->SS->thi)(l+1,ch);
                    ThetaWat2=(*adt->C->th)(l+1,ch);
                    sat2=(*adt->S->pa)(sy,jsat,l+1);
                    ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                    
                    dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
                    dD = 0.5*dz + 0.5*dzn;
                    kn = GTConst::k_air;
                    cnt++;
                    if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                    (*Lx)(cnt) = (-kn*(ThetaAir+ThetaAir2)/2.0)*area/dD*1000.0; 
                    }
                    else{
                        (*Lx)(cnt) =0.0;
                        }
                }

            }  
        }

        /** LATERAL FLUXES */
        if (i<=n)
        {


            
            
            /** lateral hydraulic conductivity */
            if (l>0)
            {
                
                ThetaIce=(*adt->S->SS->thi)(l,j);
                ThetaWat=(*adt->S->th)(l,j);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;
                
                k = GTConst::k_air;
                kmax = GTConst::k_air;

            }

            /** 4 neighbouring cells */
            R = r-1;
            C = c;
            /** -------------------- (1) -------------------- */
            if (R>=1 && R<=Nr && C>=1 && C<=Nc)
            {
                if ((long)(*adt->L->LC)(R,C)!=number_novalue && adt->T->i_cont[l][R][C]>i)
                {

                    I = adt->T->i_cont[l][R][C];
                    syn = (*adt->S->type)(R,C);
                    J = adt->T->j_cont[R][C];
                    
                    dD = find_3Ddistance(ds, (*adt->T->Z0)(r,c)-(*adt->T->Z0)(R,C)) ;//[m]
                    dn = ds/cos(0.5*atan((*adt->T->dzdE)(r,c))+0.5*atan((*adt->T->dzdE)(R,C)));//[m]

                    if (l>0)
                    {

                        
                        ThetaIce2=(*adt->S->SS->thi)(l,J);
                        ThetaWat2=(*adt->S->th)(l,J);
                        sat2=(*adt->S->pa)(syn,jsat,l);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                        /** SUBSURFACE Flow */
                        kn = GTConst::k_air;
                        cnt++;
                        if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                        (*Lx)(cnt) = (-kn*(ThetaAir+ThetaAir2)/2.0)*(dn*dz/dD)*1000.0; 
                        }
                        else{
                        (*Lx)(cnt) =0.0;
                        }
                        //printf("AIR ENERGY HORIZ1 i,I,l,dn,dz,dD,ThetaAir,ThetaAir2,LX: %i,%i, %d ,%f,%f ,%f , %f %f %e\n",i,I,l,dn,dz,dD,ThetaAir,ThetaAir2,(*Lx)(cnt));

                    }
                    else
                    {

                        /** SURFACE Flow */
                        kn = 0.;
                        cnt++;
                        (*Lx)(cnt) =0.0;
                        //(*Lx)(cnt) = (-kn*ThetaAir+(*adt->AF->Jair)(i)*GTConst::rho_air*GTConst::c_air)*(dn*1.E-3*dz/dD); //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
                    }
                }
            }

            R = r+1;
            C = c;

            /** -------------------- (2) -------------------- */
            if (R>=1 && R<=Nr && C>=1 && C<=Nc)
            {
                if ((long)(*adt->L->LC)(R,C)!=number_novalue && adt->T->i_cont[l][R][C]>i)
                {

                    I = adt->T->i_cont[l][R][C];
                    syn = (*adt->S->type)(R,C);
                    J = adt->T->j_cont[R][C];
                    
                    dD = find_3Ddistance(ds,(*adt->T->Z0)(r,c)-(*adt->T->Z0)(R,C)) ;//[m]
                    dn = ds/cos(0.5*atan((*adt->T->dzdE)(r,c))+0.5*atan((*adt->T->dzdE)(R,C)));//[m]

                    if (l>0)
                    {

                        ThetaIce2=(*adt->S->SS->thi)(l,J);
                        ThetaWat2=(*adt->S->th)(l,J);
                        sat2=(*adt->S->pa)(syn,jsat,l);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                        kn = GTConst::k_air;
                        cnt++;
                        if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                        (*Lx)(cnt) = (-kn*(ThetaAir+ThetaAir2)/2.0)*(dn*dz/dD)*1000.0; 
                        }
                        else{
                        (*Lx)(cnt) =0.0;
                        }
                        //printf("AIR ENERGY HORIZ2 i,I,l,dn,dz,dD,ThetaAir,ThetaAir2,LX: %i, %i,%d ,%f,%f ,%f , %f %f %e\n",i,I,l,dn,dz,dD,ThetaAir,ThetaAir2,(*Lx)(cnt));

                    }
                    else
                    {

                        /** SURFACE Flow */
                        kn = 0.;
                        cnt++;
                        (*Lx)(cnt) =0.;
                        //(*Lx)(cnt) = (-kn*ThetaAir+(*adt->AF->Jair)(i)*GTConst::rho_air*GTConst::c_air)*(dn*1.E-3*dz/dD); //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
                   }
                }
            }

            R = r;
            C = c-1;

            /** -------------------- (3) -------------------- */
            if (R>=1 && R<=Nr && C>=1 && C<=Nc)
            {
                if ((long)(*adt->L->LC)(R,C)!=number_novalue && adt->T->i_cont[l][R][C]>i)
                {

                    I = adt->T->i_cont[l][R][C];
                    syn = (*adt->S->type)(R,C);
                    J = adt->T->j_cont[R][C];
                    
                    dD = find_3Ddistance(ds,(*adt->T->Z0)(r,c)-(*adt->T->Z0)(R,C));  /// [mm]
                    dn = ds/cos(0.5*atan((*adt->T->dzdE)(r,c))+0.5*atan((*adt->T->dzdE)(R,C))); /// [m]

                    if (l>0)
                    {
                        ThetaIce2=(*adt->S->SS->thi)(l,J);
                        ThetaWat2=(*adt->S->th)(l,J);
                        sat2=(*adt->S->pa)(syn,jsat,l);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                        kn = GTConst::k_air;

                        cnt++;
                        if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                        (*Lx)(cnt) = (-kn*(ThetaAir+ThetaAir2)/2.0)*(dn*dz/dD)*1000.0; 
                        }
                        else{
                        (*Lx)(cnt) =0.0;
                        }
                        //printf("AIR ENERGY HORIZ3 i,I,l,dn,dz,dD,ThetaAir,ThetaAir2,LX: %i, %i,%d ,%f,%f ,%f , %f %f %e\n",i,I,l,dn,dz,dD,ThetaAir,ThetaAir2,(*Lx)(cnt));

                    }
                    else
                    {

                        /** SURFACE Flow */
                        kn = 0.;

                        cnt++;
                        (*Lx)(cnt) =0.;
                        //(*Lx)(cnt) = (-kn*ThetaAir+(*adt->AF->Jair)(i)*GTConst::rho_air*GTConst::c_air)*(dn*1.E-3*dz/dD); //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]

                    }
                }
            }

            R = r;
            C = c+1;

            /** -------------------- (4) -------------------- */
            if (R>=1 && R<=Nr && C>=1 && C<=Nc)
            {
                if ((long)(*adt->L->LC)(R,C)!=number_novalue && adt->T->i_cont[l][R][C]>i)
                {

                    I = adt->T->i_cont[l][R][C];
                    syn = (*adt->S->type)(R,C);
                    J = adt->T->j_cont[R][C];
                    
                    dD = find_3Ddistance(ds,(*adt->T->Z0)(r,c)-(*adt->T->Z0)(R,C)) ;//[m]
                    dn = ds/cos(0.5*atan((*adt->T->dzdE)(r,c))+0.5*atan((*adt->T->dzdE)(R,C)));//[m]

                    if (l>0)
                    {

                        ThetaIce2=(*adt->S->SS->thi)(l,J);
                        ThetaWat2=(*adt->S->th)(l,J);
                        sat2=(*adt->S->pa)(syn,jsat,l);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                        
                        kn = GTConst::k_air;

                        cnt++;
                        if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                        (*Lx)(cnt) = (-kn*(ThetaAir+ThetaAir2)/2.0)*(dn*dz/dD)*1000.0; 
                        }
                        else{
                        (*Lx)(cnt) =0.0;
                        }
                        //printf("AIR ENERGY HORIZ4 i,I,l,dn,dz,dD,ThetaAir,ThetaAir2,LX: %i, %i,%d ,%f,%f ,%f , %f %f %e\n",i,I,l,dn,dz,dD,ThetaAir,ThetaAir2,(*Lx)(cnt));

                    }
                    else
                    {

                        /** SURFACE Flow */
                        kn = 0.;

                        cnt++;
                        (*Lx)(cnt) =0.;
                        //(*Lx)(cnt) = (-kn*ThetaAir+(*adt->AF->Jair)(i)*GTConst::rho_air*GTConst::c_air)*(dn*1.E-3*dz/dD); //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]

                    }
                }
            }

            /** exchange with channels */
            if (l>0 && (*adt->C->ch)(r,c) > 0)
            {

                cnt++;
                (*Lx)(cnt) =0.;
                //(*Lx)(cnt) = (2. * (*adt->C->length)((*adt->C->ch)(r,c)) * 1.E-3*dz)*(-kn+(*adt->AF->Jair)(i)*GTConst::rho_air*GTConst::c_air)/dD;
            }
        }
    }
    return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_dfdHthermal_3D(double Dt, Vector<double> *df, ALLDATA *adt, SOIL_STATE *L,
                 SOIL_STATE *C, Vector<double> *T, Matrix<double> *Klat)
{
    long i, l, r, c, j, sy, ch, bc;
    long n=(Nl+1)*adt->P->total_pixel;
    double dz, dn, dD, T1,T0,V1,TempCM,Tairmet;
    double area, ds=sqrt((*UV->U)(1)*(*UV->U)(2));
    double ThetaIce,ThetaWat,sat,ThetaAir1,ThetaAir0,totw0;
    double h;

    for (i=1; i<=T->nh; i++)
    {

        (*df)(i) = 0.;

        if (i<=n)
        {
            l = (*adt->T->lrc_cont)(i,1);
            r = (*adt->T->lrc_cont)(i,2);
            c = (*adt->T->lrc_cont)(i,3);
            j = adt->T->j_cont[r][c];
            sy = (*adt->S->type)(r,c);
            bc = (*adt->T->BC_counter)(r,c);
            ch = (*adt->C->ch)(r,c);
            area = ds*ds/cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0);
            if (ch>0)
                area-= (*adt->C->length)(ch) * adt->P->w_dx * ds; //area of the pixel[m2]
            
            T1 = (*T)(i);

           
            Tairmet=(*adt->M->Tgrid)(r,c)+273.15;
            

            if (l>0)
            {
                h=(*L->Hair_internal)(l,j);
                T0 = (*adt->S->SS->Tair_internal)(l,j)+273.15;

                TempCM=(*adt->S->SS->T)(l,j)+273.15;
                ThetaIce=(*adt->S->SS->thi)(l,j);
                ThetaWat=(*adt->S->th)(l,j);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir1=sat-ThetaIce-ThetaWat;
                
                totw0=(*adt->S->SS->totw0)(l,j);
                ThetaAir0=sat-totw0;
            }

        }
        else
        {
            l = (*adt->C->lch)(i-n,1);
            ch = (*adt->C->lch)(i-n,2);
            r = (*adt->C->r)(ch);
            c = (*adt->C->c)(ch);
            sy = (*adt->C->soil_type)(ch);
            bc=0;
            area=(*adt->C->length)(ch) * adt->P->w_dx * ds;
            
            T1 = (*T)(i);
            
            
            Tairmet=(*adt->M->Tgrid)(r,c)+273.15;

            if (l>0)
            {

                h=(*C->Hair_internal)(l,ch);
                T0 = (*adt->C->SS->Tair_internal)(l,ch)+273.15;
                TempCM=(*adt->C->SS->T)(l,ch)+273.15;
                ThetaIce=(*adt->C->SS->thi)(l,ch);
                ThetaWat=(*adt->C->th)(l,ch);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir1=sat-ThetaIce-ThetaWat;
                
                totw0=(*adt->C->SS->totw0)(l,ch);
                ThetaAir0=sat-totw0;
            }

        }

        /** hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt) */
        if (l==0) 
        {
            (*df)(i) += 1.0;

        }
        else
        {
            dz = (*adt->S->pa)(sy,jdz,l)/1000.0;

            //V1 = area*dz *GTConst::rho_ref*GTConst::c_air*(1-GTConst::beta*(T0-GTConst::Tref))*(ThetaAir1);
            
            //dTheta/dt=!0
            
           
           if (ThetaAir1>GTConst::ThetaAirMin){
           
				
				V1 = area*dz *GTConst::rho_ref*GTConst::c_air*(ThetaAir1);
				
				//GZ For now dtheta/dt neglected
                //V1 = area*dz *GTConst::rho_ref*GTConst::c_air*( 2*ThetaAir1-ThetaAir0);
                (*df)(i) += (V1)/Dt*1000.0;
                
                // GZ - Include heat transference with the composite medium.
                // GZ -- Decide what to do at l=0

                (*df)(i) += area*dz*h*1000.0;
            }
            else{
                (*df)(i) += 1.0;
            }
            /*
            if ((FixedTopBoundary==1) and (l==1))
			{
			  if (ThetaAir1>GTConst::ThetaAirMin){
			  (*df)(i)+=1.0;
			  }
			}
					   
			// Heat Flux at the bottom
			// GZ - Fixed value for the heat flux from below - CHECK the sign
			
			
			if ((FixedBottomBoundary==1) and (l==Nl))
			{
				if (ThetaAir1>GTConst::ThetaAirMin){
				(*df)(i)+=1.0;

				}
				
			}*/
           
            
            
        }  

        
    }
    return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_fthermal_3D(double Dt, Vector<double> *f, ALLDATA *adt, SOIL_STATE *SL,
              SOIL_STATE *SC,STATEVAR_3D *S, Vector<double> *T)
{
  GEOTIMER_SECTION(__func__);
  
    long i;
    long n=(Nl+1)*adt->P->total_pixel;
    //double k=0.0, kn=0.0, kmax=0.0, kmaxn=0.0,krel=0.0,kp=0.0,kheat2=0.0;
    

#pragma omp parallel for
    for (i=1; i<=T->nh; i++)
    {
        long  l, r, c, j, sy, ch, bc,ns;
        double dz, dn, dD, V0, V1, psi1, psi0, ice=0.0;
        double area,Ttop,Tairmet,TempCM,TempCM0,T1,T0, ds=sqrt((*UV->U)(1)*(*UV->U)(2));
        double ThetaIce,ThetaWat,sat,ThetaAir1,ThetaAir0,totw0,cons,Res,kmax,krel,kp;
		double kt=0.0,kthcm=0.0;
        double dp=0.05, bet=10.0,h=0.0,h0=0.0;
        double snowD=0.0;
        double H_Explicit=0.0;
        
        if (i<=n)
        {
            l = (*adt->T->lrc_cont)(i,1);
            r = (*adt->T->lrc_cont)(i,2);
            c = (*adt->T->lrc_cont)(i,3);
            j = adt->T->j_cont[r][c];
            sy = (*adt->S->type)(r,c);
            
            
            bc = (*adt->T->BC_counter)(r,c);
            ch = (*adt->C->ch)(r,c);
            area = ds*ds/cos((*adt->T->slope)(r,c)*GTConst::Pi/180.);
           
            if (ch>0)
                area -= (*adt->C->length)(ch) * adt->P->w_dx * ds; //area of the pixel[m2]
                
            //std::cout << "HERE: i is equal to " << i; 
            T1 = (*T)(i);

            Tairmet=(*adt->M->Tgrid)(r,c)+273.15;
            ns=(*adt->N->S->lnum)(r,c);
            
            snowD=DEPTH(r, c, S->lnum.get(), S->Dzl.get());
			if (snowD<=adt->P->SnowDepthAirFlowLimit)
            {
                Ttop= Tairmet;
            }
            else
            {
                Ttop=(*adt->N->S->T)(1,r,c)+273.15;
            }
            
            if (l>0)
            {
                TempCM=(*adt->S->SS->T)(l,j)+273.15;
                TempCM0=(*adt->S->SS->T0)(l,j)+273.15;
                T0 = (*adt->S->SS->Tair_internal)(l,j)+273.15;
                
                if ((FixedTopBoundary==1) and (l==1))
				{
				 T0=Ttop;
				}
						   
				// Heat Flux at the bottom
				// GZ - Fixed value for the heat flux from below - CHECK the sign
				
				
				if ((FixedBottomBoundary==1) and (l==Nl))
				{
					T0=adt->P->Tboundary+273.15;
				}
                
                
                ice = (*adt->S->SS->thi)(l,j);
                ThetaIce=ice;
                ThetaWat=(*adt->S->th)(l,j);
                sat=(*adt->S->pa)(sy,jsat,l);
                dp=(*adt->S->pa)(sy,jdp,l);
                //printf("AIRENERGY l,dpaux: %i, %f,\n",l,dpaux);
				
                ThetaAir1=sat-ThetaIce-ThetaWat;
                
                totw0=(*adt->S->SS->totw0)(l,j);
                ThetaAir0=sat-totw0;
                
                kmax=(*adt->S->pa)(sy,jKn,l)/1000.0;
				kp=kmax*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
				krel=ThetaAir1/sat; // krel=relative permeability (0-1) We need to improve this value
            }

        }
        else
        {
            l = (*adt->C->lch)(i-n,1);
            ch = (*adt->C->lch)(i-n,2);
            r = (*adt->C->r)(ch);
            c = (*adt->C->c)(ch);
            sy = (*adt->C->soil_type)(ch);
            
            bc = 0;
            area = (*adt->C->length)(ch) * adt->P->w_dx * ds;
            
            T1 = (*T)(i);
            
            Tairmet=(*adt->M->Tgrid)(r,c)+273.15;
            ns=(*S->lnum)(r,c);

            snowD=DEPTH(r, c, S->lnum.get(), S->Dzl.get());
			if (snowD<=adt->P->SnowDepthAirFlowLimit)
            {
                Ttop= Tairmet;
            }
            else
            {
                Ttop=(*S->T)(1,r,c)+273.15;
            }
            if (l>0)
            {
                TempCM=(*adt->C->SS->T)(l,ch)+273.15;
                TempCM0=(*adt->C->SS->T0)(l,ch)+273.15;
                T0 = (*adt->C->SS->Tair_internal)(l,ch)+273.15;
 
                ice = (*adt->C->SS->thi)(l,ch);
                ThetaIce=ice;
                ThetaWat=(*adt->C->th)(l,j);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir1=sat-ThetaIce-ThetaWat;
                dp=(*adt->S->pa)(sy,jdp,l);
                totw0=(*adt->C->SS->totw0)(l,ch);
                ThetaAir0=sat-totw0;
            }

        }

        /** hydraulic capacity (diagonal term) */
        // GZ -- Decide what to do at l=0
        if (l==0) 
        {
            V1 = 0.;
            V0 = 0.;
            //(*f)(i) = (V1-V0)/Dt;
            //T1(l==0)=Tairmet
            //printf("AIRENERGY i,l,ns,TS,Tmet,Ttop: %i, %d, %d, %f, %f , %f\n",i,l,(*S->lnum)(r,c),(*S->T)(1,r,c),(*adt->M->Tgrid)(r,c),Ttop);
            
            (*f)(i)=T1;
            (*f)(i)-=Ttop;
        }
        else
        {
            dz = (*adt->S->pa)(sy,jdz,l)/1000.0;
            cons=  area*dz*GTConst::rho_ref*GTConst::c_air/Dt;
			//printf("AIRENERGY FFF i,l,ns,TS,Tmet,Ttop: %i, %d, %f, %f , %f\n",i,l,area,dz,GTConst::rho_ref*(1-GTConst::beta*(T0-GTConst::Tref)));
            //printf("AIR ENERGY FF i,l,totw0,totw1: %i, %d, %f, %f \n",i,l,totw0,ThetaWat);
            //(*f)(i) = cons*(ThetaAir1*(T1-T0));
            
            //dTheta/dt=!0
            /*
             if ((r==10) and (c==10)){
                //printf("AIRENERGY l,dpaux: %i, %f,\n",l,dp);
				//printf("AIRENERGY l,T0,TempCM,TempCM0,dp,h,Vel,Res,Fi = %i,%i, %f, %f , %f, %f, %f, %f , %f, %f, %f  \n",l,j,T0,TempCM,TempCM0,dp,h,h0,(*adt->AF->MaxVel)(i),Res*area*dz,cons*(ThetaAir1*(T1-T0)));
				//printf("%i, %f, %f , %f, %f, %f , %f, %f ,%i, %f , %f, %f  \n",l,T0,TempCM,TempCM0,h,fabs((*adt->AF->MaxVel)(i)),kthcm,GTConst::k_air,adt->P->HeatTransferModel,adt->P->Ht_a,adt->P->Ht_b,adt->P->Ht_n);
				//printf("%i, %f, %f, %f \n",l,dz,area,GTConst::rho_ref*GTConst::c_air);
				
				printf("AIRENERGY l,T0,TempCM,TempCM0,h,Vel,kthcm,kari,transfermodel,a,b,n,dp,sat,ThetaAir1 = %i, %f, %f , %f, %f, %f , %f, %f ,%i, %f , %f, %f, %f, %f, %f  \n",l,T0,TempCM,TempCM0,h,fabs((*adt->AF->MaxVel)(i)),kthcm,GTConst::k_air,adt->P->HeatTransferModel,adt->P->Ht_a,adt->P->Ht_b,adt->P->Ht_n,dp,sat,ThetaAir1);
				}
			*/
            if (ThetaAir1>GTConst::ThetaAirMin){
                
                (*f)(i) = cons*(ThetaAir1*(T1-T0))*1000.0;
                
                //GZ For now dtheta/dt neglected
                //(*f)(i) = cons*(ThetaAir1*(T1-T0)+(T1-GTConst::Tref)*(ThetaAir1-ThetaAir0));
                
                
                // GZ - Include heat transference with the composite medium.
				kthcm =  k_thermalcm2b(0, 4, ThetaWat, ice, sat,(*adt->S->pa)(sy,jkt,l));
				h=TransferTestConstant(adt->P->HeatTransferModel, fabs((*adt->AF->MaxVel)(i)), dp, kthcm,
                adt->P->Ht_a, adt->P->Ht_b, adt->P->Ht_n, ThetaAir1,0);
                
				 //Calculate residue of source term
                if (i<=n)
                {
                    h0=(*adt->S->SS->Hair_internal)(l,j);
                    
                }
                else{
                    h0=(*adt->C->SS->Hair_internal)(l,ch);
                }
                
                
                //if ((l==1) and (r==7) and (c==7)){
                /*
                if ((r==10) and (c==10)){
                //printf("AIRENERGY l,dpaux: %i, %f,\n",l,dp);
				//printf("AIRENERGY l,T0,TempCM,TempCM0,dp,h,Vel,Res,Fi = %i,%i, %f, %f , %f, %f, %f, %f , %f, %f, %f  \n",l,j,T0,TempCM,TempCM0,dp,h,h0,(*adt->AF->MaxVel)(i),Res*area*dz,cons*(ThetaAir1*(T1-T0)));
				//printf("%i, %f, %f , %f, %f, %f , %f, %f ,%i, %f , %f, %f  \n",l,T0,TempCM,TempCM0,h,fabs((*adt->AF->MaxVel)(i)),kthcm,GTConst::k_air,adt->P->HeatTransferModel,adt->P->Ht_a,adt->P->Ht_b,adt->P->Ht_n);
				//printf("%i, %f, %f, %f \n",l,dz,area,GTConst::rho_ref*GTConst::c_air);
				
				printf("AIRENERGY l,T0,TempCM,TempCM0,h,Vel,kthcm,kari,transfermodel,a,b,n,dp = %i, %f, %f , %f, %f, %f , %f, %f ,%i, %f , %f, %f, %f  \n",l,T0,TempCM,TempCM0,h,fabs((*adt->AF->MaxVel)(i)),kthcm,GTConst::k_air,adt->P->HeatTransferModel,adt->P->Ht_a,adt->P->Ht_b,adt->P->Ht_n,dp);
				}*/
				

               
                if (H_Explicit==1.0){
                Res=0.0;}
                else{Res=h0*(TempCM-TempCM0);}
                
                (*f)(i) += area*dz*(T1*h-TempCM*h-Res)*1000.0;
                
                //Viscous dissipation
                /*
                if ((l==1) and (r==7) and (c==7)){
                printf("AIR BALANCE TRANSFERENCE %f, %f, %f, %f,  %f, %f , %f  \n", area*dz*(T1*h-TempCM*h-Res)*1000.0, area*dz*Res*1000.0,T1-TempCM,TempCM-TempCM0,h,h0,Dt);
                }*/
               
				
                //(*f)(i) -= ThetaAir1*area*dz*pow(fabs((*adt->AF->MaxVel)(i)),2.0)*GTConst::mu_air/(kp*krel);
                if (i<=n)
                {
                    (*SL->Hair_internal)(l,j)=h;
                }
                else{
                    (*SC->Hair_internal)(l,ch)=h;
                }

                //printf("AIRENERGY i,l,ns,TS,Tmet,Ttop,TCM,T1,T0,h: %i, %d, %d, %f, %f , %f,%f, %f , %f, %f\n",i,l,(*S->lnum)(r,c),(*S->T)(1,r,c)+273.15,Tairmet,Ttop,TempCM,T1,T0,h);
            }
            else{
                (*f)(i)=T1;
                (*f)(i)-=TempCM;
                h=0;
                
                if (i<=n)
                {
                     (*SL->Hair_internal)(l,j)=h;
                }
                else{
                    (*SC->Hair_internal)(l,ch)=h;
                }
                
            }

            
        }

		/*
		if ((FixedTopBoundary==1) and (l==1))
		{
		  if (ThetaAir1>GTConst::ThetaAirMin){
		  (*f)(i)+=T1;
          (*f)(i)-=Ttop;
          }
		}
                   
        // Heat Flux at the bottom
        // GZ - Fixed value for the heat flux from below - CHECK the sign
        
        
        if ((FixedBottomBoundary==1) and (l==Nl))
        {
            if (ThetaAir1>GTConst::ThetaAirMin){
            (*f)(i)+=T1;
            (*f)(i)-=adt->P->Tboundary+273.15;
            //(*f)(i) -= ThetaAir1*area*adt->P->Fboundary;
            //printf("AIRENERGY i,l.Fboundary,%i, %d, %e\n",i,l,adt->P->Tboundary);
            }
            
        }*/
		
        /** lateral drainage at the border */
        // GZ -- No Heat flux in the lateral boundaries

    }
    


    return 0;
}
