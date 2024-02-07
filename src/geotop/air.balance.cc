
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
#include "water.balance.h"
#include "air.balance.h"
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
#define closedtop 0
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
short air_balance(double Dt, double JD0, double JD1, double JD2,
                    SOIL_STATE *L, SOIL_STATE *C,STATEVAR_3D *S, ALLDATA *adt, Vector<double> *Vsub,
                    Vector<double> *Vsup,
                    double *Voutnet, double *Voutlandsub, double *Voutlandsup,
                    double *Voutlandbottom)
{
    GEOLOG_PREFIX(__func__);
    GEOTIMER_SECTION(__func__);

    clock_t start, end;
    double Pnet, loss;
    long j;
    short a;
    

	start = clock();

	a = AirRichards3D(Dt, L, C,S, adt, &loss, Vsub, Voutlandbottom, Voutlandsub, &Pnet, adt->P->UpdateK);
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

short AirRichards3D(double Dt, SOIL_STATE *L, SOIL_STATE *C,STATEVAR_3D *S, ALLDATA *adt, double *loss, Vector<double> *Vsub,
                 double *Vbottom, double *Vlatsub, double *Total_Pnet, short updateK)
{
    GEOLOG_PREFIX(__func__);
    GEOTIMER_SECTION(__func__);

    double res=0.0, res0[3], res_prev[MM], res_av, res00, lambda[3], epsilon,
            mu=0., hnew, hold=0.;
    double ds=sqrt((*UV->U)(1) * (*UV->U)(2)), area, dz, dn, dD;
    long i, j, ch, l, r, c, m, bc, sy, cont, cont2, iter;
    long n=adt->T->lrc_cont->nrh;
    long N=adt->W->H0->nh;
    long cont_lambda_min=0;
    short out, out2;
    int sux;
    double ThetaIce,ThetaWat,sat,ThetaAir1,ThetaAir0,totw0,ns,snowD,nsaux,Tairmet;

    *Total_Pnet = 0.;

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
     
     /** INITIALIZE TEMP and water content  */
      
	for (i=1; i<=N; i++)
	{
		if (i<=n)
		{

			l = (*adt->T->lrc_cont)(i,1);
			r = (*adt->T->lrc_cont)(i,2);
			c = (*adt->T->lrc_cont)(i,3);
			j = adt->T->j_cont[r][c];
			Tairmet=(*adt->M->Tgrid)(r,c)+273.15;
			sy = (*adt->S->type)(r,c);
			ns=(*adt->N->S->lnum)(r,c);
			
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
			sy = (*adt->C->soil_type)(ch);
			ns=(*S->lnum)(r,c);
			if (l==0)
			{
				snowD=DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get());
				
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
	
    /** solution guess */

    for (i=1; i<=N; i++)
    {
		(*adt->AF->MaxVel)(i)=0.0;
		(*adt->AF->VxW)(i)=0.0;
		(*adt->AF->VxE)(i)=0.0;
		
		(*adt->AF->VyN)(i)=0.0;
		(*adt->AF->VyS)(i)=0.0;
		
		(*adt->AF->VzT)(i)=0.0;
		(*adt->AF->VzB)(i)=0.0;

        if (i<=n)
        {

            l = (*adt->T->lrc_cont)(i,1);
            r = (*adt->T->lrc_cont)(i,2);
            c = (*adt->T->lrc_cont)(i,3);
            j = adt->T->j_cont[r][c];
			(*adt->AF->P1)(i) = (*adt->AF->P0)(i);
			
        }
        else
        {

            l = (*adt->C->lch)(i-n,1);
            ch = (*adt->C->lch)(i-n,2);
            r = (*adt->C->r)(ch);
            c = (*adt->C->c)(ch);
			(*adt->AF->P1)(i) = (*adt->AF->P0)(i);
        }
    }
	
    //GZ - Build Lx (K*area/d)
    sux = find_matrix_Kair_3D(Dt, L, C,S,adt->AF->Lx.get(),adt->AF->LxB.get(),adt->AF->FluxDir.get(), adt, adt->AF->P1.get());

    //GZ - Fill f with the terms asociated to the water and ice content, the flux term associated to the buoyancy and the border conditions.
    find_fair_3D(Dt, adt->AF->f.get(), adt, L, C,S, adt->AF->P1.get());

    // GZ - This function complete the term B (B= final f). B=f+LX*deltaP
    product_matrix_using_lower_part_by_vector_plus_vector3(-1., adt->AF->B.get(), adt->AF->f.get(), adt->AF->P1.get(),
                                                          adt->T->Li.get(), adt->T->Lp.get(), adt->AF->Lx.get(),adt->AF->LxB.get(),adt->T->lrc_cont.get(),adt->C->lch.get(),n);

    res = norm_inf(adt->AF->B.get(), 1, N);

    res00 = res; /** initial norm of the residual */
    epsilon= (adt->P->AirRichardTol+ adt->P->RelTolVWb * std::min<double>( res00, sqrt((double)N) ));
    
    //printf("AIR BALANCE EPSILON,AirRichardTol : %e ,%e\n",epsilon,adt->P->AirRichardTol);
    cont=0;
    out=0;

    //printf("res:%e\n",res);
    //The condition to go out from the Newton cycle is made on the norm of the residual < Absolute tolerance + Relative tolerance * res00
    if ( res <= std::min<double>(epsilon, max_res_adm) )
        out=1;

    //Max iteration number
    if ( cont >= adt->P->MaxiterTol )
        out=1;
    
    //printf("AIR BALANCE INI out:%ld\n",out);

    while (out==0)
    {

        cont++;

        for (i=1; i<=N; i++)
        {
            (*adt->AF->P0)(i) = (*adt->AF->P1)(i);
            (*adt->AF->dP)(i) = 0.0;
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

		//GZ - This term is = 0 for the air. Only the terms associated to K exist. 
        sux = find_dfdPair_3D(Dt, adt->AF->df.get(), adt, L, C,S, adt->AF->P1.get());
        
        /** it calculates only df/dH, J = I*df/dH + K, K is calculated above **/

        /** CONJUGATED GRADIENTS ALGORITHM */
        //printf("AIR BALANCE BEFORE BICGSTAB ITE:%ld \n",cont);
        iter = BiCGSTAB_strict_lower_matrix_plus_identity_by_vector3(mu, tol_min_GC, tol_max_GC, adt->AF->dP.get(),
                                                                    adt->AF->B.get(), adt->AF->df.get(), adt->T->Li.get(),
                                                                    adt->T->Lp.get(), adt->AF->Lx.get(),adt->T->lrc_cont.get(),adt->C->lch.get(),n);
        if (iter==-1){
            printf("AIR BALANCE S1\n");
            return 1; /** does not converge */
            }

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
                
   
                (*adt->AF->P1)(i) = (*adt->AF->P0)(i) + lambda[0] * (*adt->AF->dP)(i);
                
                if ((*adt->AF->P1)(i) != (*adt->AF->P1)(i)){
                    
                    //printf("AIR BALANCE S2 i:%i - l:%d - r:%d c:%d P:%e F:%e DF:%e P0:%e DP:%e \n : ",i,(*adt->T->lrc_cont)(i,1),(*adt->T->lrc_cont)(i,2),(*adt->T->lrc_cont)(i,3),
                    //(*adt->AF->P1)(i),(*adt->AF->B)(i),(*adt->AF->df)(i),(*adt->AF->P0)(i),(*adt->AF->dP)(i));
                    return 1;
                    }
            }
            
            if ( cont <= maxITER_rec_K && (adt->AF->Courant)>0){
                     sux = find_matrix_Kair_3D(Dt, L, C,S,adt->AF->Lx.get(),adt->AF->LxB.get(),adt->AF->FluxDir.get(), adt, adt->AF->P1.get());}
			/*
            if (updateK == 1 && cont <= maxITER_rec_K){
                printf("UPDATEKKK-----------------------------------------------------------------------------------------*************************----------------");
                sux = find_matrix_Kair_3D(Dt, L, C,S,adt->AF->Lx.get(),adt->AF->LxB.get(), adt, adt->AF->P1.get());}
            */
            /*
            if ( cont <= maxITER_rec_K){
                sux = find_matrix_Kair_3D(Dt, L, C,S,adt->AF->Lx.get(),adt->AF->LxB.get(), adt, adt->AF->P1.get());}
                */
            

            //GZ - Fill f with the terms asociated to the water and ice content, the flux term associated to the buoyancy and the border conditions.
            find_fair_3D(Dt, adt->AF->f.get(), adt, L, C,S, adt->AF->P1.get());
            sux = find_dfdPair_3D(Dt, adt->AF->df.get(), adt, L, C,S, adt->AF->P1.get());

            product_matrix_using_lower_part_by_vector_plus_vector3(-1., adt->AF->B.get(), adt->AF->f.get(), adt->AF->P1.get(),
                                                          adt->T->Li.get(), adt->T->Lp.get(), adt->AF->Lx.get(),adt->AF->LxB.get(),adt->T->lrc_cont.get(),adt->C->lch.get(),n);
            res = norm_inf(adt->AF->B.get(), 1, N);
            //printf("..res:%e\n",res);
            
            /*
            for (i=1; i<=N; i++)
            {
                if (fabs((*adt->AF->B)(i))>epsilon){
                
                printf("AIR BALANCE RESS i,l,r,c,res: %i, %d, %d, %d, %f \n",i,(*adt->T->lrc_cont)(i,1),(*adt->T->lrc_cont)(i,2),(*adt->T->lrc_cont)(i,3),(*adt->AF->B)(i));
                
                }
            }*/

            out2=0;

            geolog << "cnt:" << cont << " res:" << res << " lambda:" << lambda[0]
                   << " Dt:" << Dt << " P:" << *Total_Pnet << std::endl;

            if (res <= (1.0 - ni*lambda[0]*(1.-mu))*res_av)
                out2=1;
            if (lambda[0] <= adt->P->min_lambda_wat)
                cont_lambda_min++;

            if (cont_lambda_min > adt->P->max_times_min_lambda_wat)
            {
                if (adt->P->exit_lambda_min_wat == 1)
                {
                    printf("AIR BALANCE S3\n");
                    /*
                    for (i=1; i<=N; i++)
                    {
                        
                        if (i<=n)
                        {
                            l = (*adt->T->lrc_cont)(i,1);
                            r = (*adt->T->lrc_cont)(i,2);
                            c = (*adt->T->lrc_cont)(i,3);
                            j = adt->T->j_cont[r][c];
                            
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
                            }
                            else{
                            ThetaAir1=0.0;
                            ThetaAir0=0.0;
                            }
                            
                            if (fabs((*adt->AF->B)(i))>epsilon){
                            //printf("AIR BALANCE S3 i,l,r,c,ThetaAir0,ThetaAir1,P,T,RES,ns, %i, %d, %d, %d, %e, %e, %e, %f,%e,%f\n",i,l,r,c,ThetaAir0,ThetaAir1,(*adt->AF->P1)(i),(*adt->AE->T1)(i),(*adt->AF->B)(i),nsaux);
                            //printf("AIRENERGY i,l,ns,ThetaAir0,ThetaAir1,TS,TCM,Tmet,T1,T0: %i, %d, %d, %e, %e, %f,%f, %f , %f,%f\n",i,l,(*S->lnum)(r,c),ThetaAir0,ThetaAir1,(*S->T)(1,r,c)+273.15,(*L->T)(l,c)+273.15,(*adt->M->Tgrid)(r,c)+273.15,(*adt->AE->T1)(i),(*L->Tair_internal)(l,j)+273.15);
							}
                        }
 
                    }*/
                    return 1;
                }
                else
                {
                    out2=1;
                    cont_lambda_min=0;
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

    if ( res > epsilon ) return 1;

    /** it can be shown that massloss per unit pixel [mm] is the linear norm of B * Dt / total_pixel */
    *loss = norm_1(adt->AF->B.get(), 1, N)*Dt/adt->P->total_area;

    /** assign updated state variables */
    // GZ -- Necessary to check if we need to update variables!
    sux = find_matrix_LxHeat_air_3D(Dt, L, C,adt->AF->LxJair.get(),adt->AF->FluxDir.get(), adt->AF->Klat.get(), adt->AF->Kbottom.get(),
                           adt->C->Kbottom.get(), adt, adt->AF->P1.get());
    //printf("----------------------------------AIR BALANCE READYYYYY-------------------------------------------\n");
    /*
    for (i=1; i<=N; i++)
    {
        
        if (i<=n)
        {
            l = (*adt->T->lrc_cont)(i,1);
            r = (*adt->T->lrc_cont)(i,2);
            c = (*adt->T->lrc_cont)(i,3);
            j = adt->T->j_cont[r][c];
            
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
            }
            else{
            ThetaAir1=0.0;
            ThetaAir0=0.0;
            printf("AIR BALANCE FINAL i,l,r,c,P,T %i, %d, %d, %d, %f, %f, %f, \n",i,l,r,c,(*adt->T->Z)(l,r,c)/1000,(*adt->AF->P1)(i),(*adt->AE->T1)(i));
            }
            //printf("AIR BALANCE FINAL i,l,r,c,ThetaAir0,ThetaAir1,P,T,B-RES %i, %d, %d, %d, %e, %e, %e, %f,%e\n",i,l,r,c,ThetaAir0,ThetaAir1,(*adt->AF->P1)(i),(*adt->AE->T1)(i),(*adt->AF->B)(i));
            //printf("AIRENERGY i,l,ns,ThetaAir0,ThetaAir1,TS,TCM,Tmet,T1,T0: %i, %d, %d, %e, %e, %f,%f, %f , %f,%f\n",i,l,(*S->lnum)(r,c),ThetaAir0,ThetaAir1,(*S->T)(1,r,c)+273.15,(*L->T)(l,c)+273.15,(*adt->M->Tgrid)(r,c)+273.15,(*adt->AE->T1)(i),(*L->Tair_internal)(l,j)+273.15);

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

        }
    }*/

    return 0;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_matrix_Kair_3D(double  /*Dt*/, SOIL_STATE *SL, SOIL_STATE *SC,STATEVAR_3D *S,
                     Vector<double> *Lx,Vector<double> *LxB,Vector<short>*FluxDir, ALLDATA *adt, Vector<double> *P)
{
    GEOTIMER_SECTION(__func__);

    long i, l, r, c, j, I, R, C, J, sy, syn, ch, cnt=0;
    long n=(Nl+1)*adt->P->total_pixel;
    double dz=0.0, dzn=0.0, dD=0.0, k=0.0, kn=0.0,kheat=0.0,kheatn=0.0, kmax=0.0, kmaxn=0.0,krel=0.0,kp=0.0;
    double kn2=0.0,krel2=0.0,kp2=0.0;
    double area, ds=sqrt((*UV->U)(1)*(*UV->U)(2)), dn,ns;
    double ThetaIce,ThetaWat,sat,ThetaAir,Ha;
    double ThetaIce2,ThetaWat2,sat2,ThetaAir2;// GZ - This variables are calculated in the cell where where there is flux exchange
	double Pmi,PmI; //Presion motriz i=current cel, I= cell with flux exchange
	double fe,knHar,kheatHar; //Harmonic mean. 
	short FluxDirAux=0.0;
	double snowD=0.0;
	short Upwind=0.0;
	double SurfaceKincrease=2.0;
	double Klateral=100.0;
	double dzT=0.5;
	double Klatmax=200;

    //double psi, ice, a, ns, res, sat, ss, Temp;
    //equation written in [m2*mm/s]
	// In this function Lx is filled with flow term associated to the presuure only. The term associated to the buoyancy are included in find_f_air_3D. 
	//In fin_f_air_3D is also included the border conditions
    for (i=1; i<=P->nh; i++)
    {
        /** VERTICAL FLUXES */
        if ( i<=n) //land
        {

            l=(*adt->T->lrc_cont)(i,1);
            r=(*adt->T->lrc_cont)(i,2);
            c=(*adt->T->lrc_cont)(i,3);
            j=adt->T->j_cont[r][c];
            sy=(*adt->S->type)(r,c);

            ch=(*adt->C->ch)(r,c);
            area=ds*ds/cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0);
            if (ch>0) area-=(*adt->C->length)(ch) * adt->P->w_dx *
                            ds; //area of the pixel[m2]

            /** vertical hydraulic conductivity */
            if (l>0)
            {
                ThetaIce=(*adt->S->SS->thi)(l,j);
                ThetaWat=(*adt->S->th)(l,j);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;
                //ThetaAir=std::max<double>(ThetaAir, GTConst::ThetaAirMin);
                
                kmax=(*adt->S->pa)(sy,jKn,l)/1000.0; //k[mm/s]/1000=k[m/s] 
                dz = (*adt->S->pa)(sy,jdz,l)/1000.0;
            }

            /** flux from cell below */
            if (l<Nl)
            {

                I = i+1;

                if (l==0) /** overland flow */
                {
                   /* /// (1) ----- Top cell flux -Air flux in the top -- Boundary condition P(l==0)=Patm (if Ns==0). If Ns>0 Air flux top=0*/
					snowD=DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get());
					//printf(" SNOOOOW DEPTHH, %f, \n",snowD);adt->P->HeatTransferModel
					if ((snowD<=adt->P->SnowDepthAirFlowLimit) and (closedtop==0))
					{
						dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
						dD = 0.5*dzn;
						
						ThetaIce2=(*adt->S->SS->thi)(l+1,j);
						ThetaWat2=(*adt->S->th)(l+1,j);
						sat2=(*adt->S->pa)(sy,jsat,l+1);
						ThetaAir2=sat2-ThetaIce2-ThetaWat2;

						krel2=kair_from_psi(jKn, (*SL->P)(l+1,j),(*SL->thi)(l+1,j),l+1,adt->S->pa->matrix(sy), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air);
							
						kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
					
						krel=1.0;// krel=relative permeability (0-1) We need to improve this value
						kmaxn=(*adt->S->pa)(sy,jKn,l+1)/1000.0;
						kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn=SurfaceKincrease*(kp2*krel/GTConst::mu_air);
						kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
						
						knHar=(kn+kn2)/2.0;
						kheatHar=(kheat+kheatn)/2.0;
						
						if (Upwind==1){
						if ((((*P)(I)-(*P)(i))-kheatHar*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))>0)
						{
						/** upward flux */
							knHar=kn2;
							kheatHar=kheatn;
							FluxDirAux=1.0;
						}
						else
						{
						/** downward flow */
							knHar=kn; 
							kheatHar=kheat;
							FluxDirAux=-1.0;
						}
						}
                    
						cnt++;

						(*Lx)(cnt) = -(area*knHar/dD)*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
						(*LxB)(cnt) = -(area*knHar/dD*(kheatHar*(-dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))))*1000.0  ;
						(*FluxDir)(cnt)=FluxDirAux;
                            
                            
					}
					else
					{
 
                        cnt++;

                        (*Lx)(cnt) = 0.0;
                        (*LxB)(cnt) = 0.0;
                        (*FluxDir)(cnt)=0.0;
                        //printf("AIR BALANCE i,l, (*Lx)(cnt),(*Lx)(cnt)/area,%i,%i,%e,%e \n",i,l, (*Lx)(cnt),(*Lx)(cnt)/area/1000);

					}
                    
                }
                else    /** subsurface flow */
                {
                    dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
                    dD = 0.5*dz + 0.5*dzn;
                    fe=0.5*dzn/dD;

                    ThetaIce2=(*adt->S->SS->thi)(l+1,j);
                    ThetaWat2=(*adt->S->th)(l+1,j);
                    sat2=(*adt->S->pa)(sy,jsat,l+1);
                    ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                    //ThetaAir2=std::max<double>(ThetaAir2, GTConst::ThetaAirMin);
                    
                    krel2=kair_from_psi(jKn, (*SL->P)(l+1,j),(*SL->thi)(l+1,j),l+1,adt->S->pa->matrix(sy), adt->P->k_to_ksat,ThetaAir2)/1000.0;
                    kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
                    kn2=(kp2/GTConst::mu_air);
                    
                    kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                    
                    krel=kair_from_psi(jKn, (*SL->P)(l,j),(*SL->thi)(l,j),l,adt->S->pa->matrix(sy), adt->P->k_to_ksat,ThetaAir)/1000.0;
                    kp=krel*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
                    kn=(kp/GTConst::mu_air);
                    
                    kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                    
                    //printf("AIR BALANCE FF i,l,krel,kn: %i, %d, %f, %e \n",i,l,krel,kp);
                    knHar=pow(((1.0-fe)/kn+fe/kn2),-1.0);
                    //kheatHar=pow(((1.0-fe)/kheat+fe/kheatn),-1.0);
                    
                    //knHar=(kn+kn2)/2.0;
                    kheatHar=(kheat+kheatn)/2.0;
                    
                    if (Upwind==1){
                    if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)>0)
					{
					/** upward flux */
						knHar=kn2;
						kheatHar=kheatn;
						FluxDirAux=1.0;
					}
					else
					{
					/** downward flow */
						knHar=kn;
						kheatHar=kheat;
						FluxDirAux=-1.0;
					}}

                    cnt++;

                    (*Lx)(cnt) = -area/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
                    (*LxB)(cnt) = -area/dD*knHar*(kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)*1000.0  ;
                    (*FluxDir)(cnt)=FluxDirAux;
                    
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
                //ThetaAir=sat-ThetaIce-ThetaWat;
                //ThetaAir=std::max<double>(ThetaAir, GTConst::ThetaAirMin);
                
                kmax=(*adt->S->pa)(sy,jKn,l)/1000.0;
                dz = (*adt->S->pa)(sy,jdz,l)/1000.0;
                
            }

            /** flux from cell below */
            if (l<Nl)
            {

                I = i+1;

                if (l==0) /** OVERLAND flow */
                {

					/* /// (1) ----- Top cell flux -Air flux in the top -- Boundary condition P(l==0)=Patm (if Ns==0). If Ns>0 Air flux top=0*/
				
					ns=(*adt->N->S->lnum)(r,c);
					if ((ns==0) and (closedtop==0))
					{
						/// (1) ----- Top cell flux
						dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
						dD = 0.5*dzn;

						ThetaIce2=(*adt->C->SS->thi)(l+1,ch);
						ThetaWat2=(*adt->C->th)(l+1,ch);
						sat2=(*adt->S->pa)(sy,jsat,l+1);
						//ThetaAir2=sat2-ThetaIce2-ThetaWat2;
						//ThetaAir2=std::max<double>(ThetaAir2, GTConst::ThetaAirMin);

						kmaxn=(*adt->S->pa)(sy,jKn,l+1)/1000.0;
                        
                        kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
                        krel2=ThetaAir2/sat2; // krel=relative permeability (0-1) We need to improve this value
                        kn2=(kp2*krel2/GTConst::mu_air);
                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                    
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));

						cnt++;
                        
                        (*Lx)(cnt) = -area*kn/dD*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m3/s]
                        (*LxB)(cnt) = -area*kn/dD*((kheat+kheatn)/2.0*(-dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0)))*1000.0  ;
          
				
					}
					else
					{
						
                        /// (1) ----- Top cell flux
						dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
						dD = 0.5*dzn;

						ThetaIce2=(*adt->C->SS->thi)(l+1,ch);
						ThetaWat2=(*adt->C->th)(l+1,ch);
						sat2=(*adt->S->pa)(sy,jsat,l+1);
						ThetaAir2=sat2-ThetaIce2-ThetaWat2;
						ThetaAir2=std::max<double>(ThetaAir2, GTConst::ThetaAirMin);

						kmaxn=(*adt->S->pa)(sy,jKn,l+1)/1000.0;
                        
                        kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
                        krel2=ThetaAir2/sat2; // krel=relative permeability (0-1) We need to improve this value
                        kn2=(kp2*krel2/GTConst::mu_air);
                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                    
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));

						cnt++;
                        
						(*Lx)(cnt) =0.0;
						(*LxB)(cnt) =0.0;
					}

                }
                else    /** SUBSURFACE flow */
                {

                    dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
                    dD = 0.5*dz + 0.5*dzn;
					fe=0.5*dzn/dD;
					
                    ThetaIce2=(*adt->C->SS->thi)(l+1,ch);
                    ThetaWat2=(*adt->C->th)(l+1,ch);
                    sat2=(*adt->S->pa)(sy,jsat,l+1);
                    ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                    //ThetaAir2=std::max<double>(ThetaAir2, GTConst::ThetaAirMin);
                    
                    kmaxn=(*adt->S->pa)(sy,jKn,l+1)/1000.0;
                    kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
                    krel2=ThetaAir2/sat2; // krel=relative permeability (0-1) We need to improve this value
                    kn2=(kp2*krel2/GTConst::mu_air);
                    kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                
                    kp=kmax*GTConst::mu_l/(GTConst::rho_w*GTConst::g);//  k=pkn*mu/(rho*g), kp=permeability
                    krel=ThetaAir/sat;// krel=relative permeability (0-1) We need to improve this value
                    kn=(kp*krel/GTConst::mu_air);
                    kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                    
                    knHar=pow(((1.0-fe)/kn+fe/kn2),-1.0);
                    kheatHar=pow(((1.0-fe)/kheat+fe/kheatn),-1.0);
                    
                    cnt++;

                    (*Lx)(cnt) = -(area/dD*knHar)*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
                    (*LxB)(cnt) = -(area/dD*knHar*(kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0))*1000.0  ;

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
                //ThetaAir=std::max<double>(ThetaAir, GTConst::ThetaAirMin);
                               
                kmax=(*adt->S->pa)(sy,jKl,l)/1000.0;
				
				krel=kair_from_psi(jKl, (*SL->P)(l,j),(*SL->thi)(l,j),l,adt->S->pa->matrix(sy), adt->P->k_to_ksat,ThetaAir)/1000.0;
				kp=krel*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
				k=(kp/GTConst::mu_air);
                
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
                    fe=0.5*dzn/dD;

                    if (l>0)
                    {

                        ThetaIce2=(*adt->S->SS->thi)(l,J);
                        ThetaWat2=(*adt->S->th)(l,J);
                        sat2=(*adt->S->pa)(syn,jsat,l);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                        //ThetaAir2=std::max<double>(ThetaAir2, GTConst::ThetaAirMin);
                        
                        krel2=kair_from_psi(jKl, (*SL->P)(l,J),(*SL->thi)(l,J),l,adt->S->pa->matrix(syn), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air);

                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                    
                        kn=k;
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                        
                        knHar=2.0*kn*kn2/(kn+kn2);
						//kheatHar=2.0*kheat*kheatn/(kheat+kheatn);
						
						//knHar=(kn+kn2)/2.0;
						kheatHar=(kheat+kheatn)/2.0;
						
						if (Upwind==1){
						if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0)
						{
						/** Inbound flux */
							knHar=kn2;
							kheatHar=kheatn;
							FluxDirAux=1.0;
						}
						else
						{
						/** Outbound flux */
							knHar=kn;
							kheatHar=kheat;
							FluxDirAux=-1.0;
						}}
                    
                        cnt++;

                        (*Lx)(cnt) = -(dn*dz)/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
                        (*LxB)(cnt) = -(dn*dz)/dD*knHar*kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))*1000.0;
                        (*FluxDir)(cnt)=FluxDirAux;

                    }
                    else
                    {

                        snowD=std::max<double>(DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get()),DEPTH(R,C, adt->N->S->lnum.get(), adt->N->S->Dzl.get()));
						//printf(" SNOOOOW DEPTHH, %f, \n",snowD);adt->P->HeatTransferModel
						if ((snowD<=adt->P->SnowDepthAirFlowLimit) and (closedtop==0))
						{
							kmaxn=Klateral*(*adt->S->pa)(syn,jKn,l+1)/1000.0;
							//printf("AIR BALANCE 1 - i,l,kmaxn,%i,%i,%f \n",i,l, kmaxn);
							kmaxn=std::min<double>(kmaxn,Klatmax);
							//printf("AIR BALANCE 1 - i,l,kmaxn,%i,%i,%f \n",i,l, kmaxn);
							kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
							krel2=1.0; // krel=relative permeability (0-1) We need to improve this value
							kn2=(kp2*krel2/GTConst::mu_air);
							kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
						
							krel=1.0;// krel=relative permeability (0-1) We need to improve this value
							kn=(kp2*krel/GTConst::mu_air);
							kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
							
							knHar=(kn+kn2)/2.0;
							kheatHar=(kheat+kheatn)/2.0;
							
							if (Upwind==1){
							if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0)
							{
							/** Inbound flux */
								knHar=kn2;
								kheatHar=kheatn;
								FluxDirAux=1.0;
							}
							else
							{
							/** Outbound flux */
								knHar=kn;
								kheatHar=kheat;
								FluxDirAux=-1.0;
							}}
						
						
							cnt++;

							(*Lx)(cnt) = -(dn*dzT)/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
							(*LxB)(cnt) = -(dn*dzT)/dD*knHar*kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))*1000.0;
							(*FluxDir)(cnt)=FluxDirAux;
								
								
						}
						else
						{
	 
							cnt++;

							(*Lx)(cnt) = 0.0;
							(*LxB)(cnt) = 0.0;
							(*FluxDir)(cnt)=FluxDirAux;
							//printf("AIR BALANCE i,l, (*Lx)(cnt),(*Lx)(cnt)/area,%i,%i,%e,%e \n",i,l, (*Lx)(cnt),(*Lx)(cnt)/area/1000);

						}

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
                        //ThetaAir2=std::max<double>(ThetaAir2, GTConst::ThetaAirMin);
                        
                        krel2=kair_from_psi(jKl, (*SL->P)(l,J),(*SL->thi)(l,J),l,adt->S->pa->matrix(syn), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air);


                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                    
                        kn=k;
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                        
                        knHar=2.0*kn*kn2/(kn+kn2);
						//kheatHar=2.0*kheat*kheatn/(kheat+kheatn);
						
						//knHar=(kn+kn2)/2.0;
						kheatHar=(kheat+kheatn)/2.0;
						
						if (Upwind==1){
							if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0)
							{
							/** Inbound flux */
								knHar=kn2;
								kheatHar=kheatn;
								FluxDirAux=1.0;
							}
							else
							{
							/** Outbound flux */
								knHar=kn;
								kheatHar=kheat;
								FluxDirAux=-1.0;
							}}
                    
                        cnt++;

                        (*Lx)(cnt) = -(dn*dz)/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
                        (*LxB)(cnt) = -(dn*dz)/dD*knHar*kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))*1000.0;
                        (*FluxDir)(cnt)=FluxDirAux;

                    }
                    else
                    {

                        snowD=std::max<double>(DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get()),DEPTH(R,C, adt->N->S->lnum.get(), adt->N->S->Dzl.get()));
						//printf(" SNOOOOW DEPTHH, %f, \n",snowD);adt->P->HeatTransferModel
						if ((snowD<=adt->P->SnowDepthAirFlowLimit) and (closedtop==0))
						{
							kmaxn=Klateral*(*adt->S->pa)(syn,jKn,l+1)/1000.0;
							//printf("AIR BALANCE 2 - i,l,kmaxn,%i,%i,%f \n",i,l, kmaxn);
							kmaxn=std::min<double>(kmaxn,Klatmax);
							//printf("AIR BALANCE 2 - i,l,kmaxn,%i,%i,%f \n",i,l, kmaxn);
							kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
							krel2=1.0; // krel=relative permeability (0-1) We need to improve this value
							kn2=(kp2*krel2/GTConst::mu_air);
							kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
						
							krel=1.0;// krel=relative permeability (0-1) We need to improve this value
							kn=(kp2*krel/GTConst::mu_air);
							kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
							
							knHar=(kn+kn2)/2.0;
							kheatHar=(kheat+kheatn)/2.0;
							
							if (Upwind==1){
							if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0)
							{
							/** Inbound flux */
								knHar=kn2;
								kheatHar=kheatn;
								FluxDirAux=1.0;
							}
							else
							{
							/** Outbound flux */
								knHar=kn;
								kheatHar=kheat;
								FluxDirAux=-1.0;
							}}
						
							cnt++;

							(*Lx)(cnt) = -(dn*dzT)/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
							(*LxB)(cnt) = -(dn*dzT)/dD*knHar*kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))*1000.0;
							(*FluxDir)(cnt)=FluxDirAux;
								
								
						}
						else
						{
	 
							cnt++;

							(*Lx)(cnt) = 0.0;
							(*LxB)(cnt) = 0.0;
							(*FluxDir)(cnt)=00;
							//printf("AIR BALANCE i,l, (*Lx)(cnt),(*Lx)(cnt)/area,%i,%i,%e,%e \n",i,l, (*Lx)(cnt),(*Lx)(cnt)/area/1000);

						}

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

                    dD = find_3Ddistance(ds,(*adt->T->Z0)(r,c)-(*adt->T->Z0)(R,C)) ; /// [m]
                    dn = ds/cos(0.5*atan((*adt->T->dzdE)(r,c))+0.5*atan((*adt->T->dzdE)(R,C))); /// [m]

                    if (l>0)
                    {
                        
                        ThetaIce2=(*adt->S->SS->thi)(l,J);
                        ThetaWat2=(*adt->S->th)(l,J);
                        sat2=(*adt->S->pa)(syn,jsat,l);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                        //ThetaAir2=std::max<double>(ThetaAir2, GTConst::ThetaAirMin);
                        
                        krel2=kair_from_psi(jKl, (*SL->P)(l,J),(*SL->thi)(l,J),l,adt->S->pa->matrix(syn), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air);

                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                    
                        kn=k;
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                        
                        knHar=2.0*kn*kn2/(kn+kn2);
						//kheatHar=2.0*kheat*kheatn/(kheat+kheatn);
						
						//knHar=(kn+kn2)/2.0;
						kheatHar=(kheat+kheatn)/2.0;
						
						if (Upwind==1){
							if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0)
							{
							/** Inbound flux */
								knHar=kn2;
								kheatHar=kheatn;
								FluxDirAux=1.0;
							}
							else
							{
							/** Outbound flux */
								knHar=kn;
								kheatHar=kheat;
								FluxDirAux=-1.0;
							}}
							
                    
                        cnt++;

                        (*Lx)(cnt) = -(dn*dz)/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
                        (*LxB)(cnt) = -(dn*dz)/dD*knHar*kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))*1000.0;
                        (*FluxDir)(cnt)=FluxDirAux;


                    }
                    else
                    {

                        snowD=std::max<double>(DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get()),DEPTH(R,C, adt->N->S->lnum.get(), adt->N->S->Dzl.get()));
						//printf(" SNOOOOW DEPTHH, %f, \n",snowD);adt->P->HeatTransferModel
						if ((snowD<=adt->P->SnowDepthAirFlowLimit) and (closedtop==0))
						{
							kmaxn=Klateral*(*adt->S->pa)(syn,jKn,l+1)/1000.0;
							//printf("AIR BALANCE 3 - i,l,kmaxn,%i,%i,%f \n",i,l, kmaxn);
							kmaxn=std::min<double>(kmaxn,Klatmax);
							//printf("AIR BALANCE 3 - i,l,kmaxn,%i,%i,%f \n",i,l, kmaxn);
							kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
							krel2=1.0; // krel=relative permeability (0-1) We need to improve this value
							kn2=(kp2*krel2/GTConst::mu_air);
							kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
						
							krel=1.0;// krel=relative permeability (0-1) We need to improve this value
							kn=(kp2*krel/GTConst::mu_air);
							kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
							
							knHar=(kn+kn2)/2.0;
							kheatHar=(kheat+kheatn)/2.0;
							
							if (Upwind==1){
							if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0)
							{
							/** Inbound flux */
								knHar=kn2;
								kheatHar=kheatn;
								FluxDirAux=1.0;
							}
							else
							{
							/** Outbound flux */
								knHar=kn;
								kheatHar=kheat;
								FluxDirAux=-1.0;
							}}
							
							cnt++;

							(*Lx)(cnt) = -(dn*dzT)/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
							(*LxB)(cnt) = -(dn*dzT)/dD*knHar*kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))*1000.0;
							(*FluxDir)(cnt)=FluxDirAux;
								
						}
						else
						{
							cnt++;

							(*Lx)(cnt) = 0.0;
							(*LxB)(cnt) = 0.0;
							(*FluxDir)(cnt)=0.0;
							//printf("AIR BALANCE i,l, (*Lx)(cnt),(*Lx)(cnt)/area,%i,%i,%e,%e \n",i,l, (*Lx)(cnt),(*Lx)(cnt)/area/1000);

						}

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
                        //ThetaAir2=std::max<double>(ThetaAir2, GTConst::ThetaAirMin);

                        krel2=kair_from_psi(jKl, (*SL->P)(l,J),(*SL->thi)(l,J),l,adt->S->pa->matrix(syn), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air);

                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                    
                        kn=k;
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                        
                        knHar=2.0*kn*kn2/(kn+kn2);
						//kheatHar=2.0*kheat*kheatn/(kheat+kheatn);
						
						//knHar=(kn+kn2)/2.0;
						kheatHar=(kheat+kheatn)/2.0;
                        //knHar=2.0*kn*kn2/(kn+kn2);
						//kheatHar=2.0*kheat*kheatn/(kheat+kheatn);
						if (Upwind==1){
							if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0)
							{
							/** Inbound flux */
								knHar=kn2;
								kheatHar=kheatn;
								FluxDirAux=1.0;
							}
							else
							{
							/** Outbound flux */
								knHar=kn;
								kheatHar=kheat;
								FluxDirAux=-1.0;
							}}
                    
                        cnt++;

                        (*Lx)(cnt) = -(dn*dz)/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
                        (*LxB)(cnt) = -(dn*dz)/dD*knHar*kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))*1000.0;
                        (*FluxDir)(cnt)=FluxDirAux;
 
                    }
                    else
                    {

                        snowD=std::max<double>(DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get()),DEPTH(R,C, adt->N->S->lnum.get(), adt->N->S->Dzl.get()));
						//printf(" SNOOOOW DEPTHH, %f, \n",snowD);adt->P->HeatTransferModel
						if ((snowD<=adt->P->SnowDepthAirFlowLimit) and (closedtop==0))
						{
							kmaxn=Klateral*(*adt->S->pa)(syn,jKn,l+1)/1000.0;
							//printf("AIR BALANCE 4 - i,l,kmaxn,%i,%i,%f \n",i,l, kmaxn);
							kmaxn=std::min<double>(kmaxn,Klatmax);
							//printf("AIR BALANCE 4 - i,l,kmaxn,%i,%i,%f \n",i,l, kmaxn);
							kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
							krel2=1.0; // krel=relative permeability (0-1) We need to improve this value
							kn2=(kp2*krel2/GTConst::mu_air);
							kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));

							krel=1.0;// krel=relative permeability (0-1) We need to improve this value
							kn=(kp2*krel/GTConst::mu_air);
							kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
							
							knHar=(kn+kn2)/2.0;
							kheatHar=(kheat+kheatn)/2.0;
							
							if (Upwind==1){
							if ((((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0)
							{
							/** Inbound flux */
								knHar=kn2;
								kheatHar=kheatn;
								FluxDirAux=1.0;
							}
							else
							{
							/** Outbound flux */
								knHar=kn;
								kheatHar=kheat;
								FluxDirAux=-1.0;
							}}
						
							cnt++;
							(*Lx)(cnt) = -(dn*dzT)/dD*knHar*1000.0 ;  //Area[m2] * kn[m2/Pa] * dP[Pa]/dD[m], equation written in [m2*mm/s]
							(*LxB)(cnt) = -(dn*dzT)/dD*knHar*kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))*1000.0;
							(*FluxDir)(cnt)=FluxDirAux;
								
								
						}
						else
						{
							cnt++;

							(*Lx)(cnt) = 0.0;
							(*LxB)(cnt) = 0.0;
							(*FluxDir)(cnt)=0.0;
							//printf("AIR BALANCE i,l, (*Lx)(cnt),(*Lx)(cnt)/area,%i,%i,%e,%e \n",i,l, (*Lx)(cnt),(*Lx)(cnt)/area/1000);

						}
                    }
                }
            }

            /** exchange with channels */
            if (l>0 && (*adt->C->ch)(r,c) > 0)
            {

                ch = (*adt->C->ch)(r,c);
                syn = (*adt->C->soil_type)(ch);
                I = n + adt->C->ch3[l][ch];

                kn = 0.0;

                /** Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
                 * Area[m2] = 2 * channel length * layer thickness = 2 * length[m] * dz[mm] * 1.E-3[m/mm]
                 * dD[mm] = 0.25 * ds * (1+w_dx) */
                dD = find_3Ddistance(ds * (1.0 + adt->P->w_dx) / 4.0,
                                     1.E-3*adt->P->depr_channel) * 1.E3;//[mm]


                cnt++;
                (*Lx)(cnt) =0.0;
                (*LxB)(cnt) =0.0;
                (*FluxDir)(cnt)=0.0;
            }
        }
    }
    return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_matrix_LxHeat_air_3D(double  Dt, SOIL_STATE *SL, SOIL_STATE *SC,
                     Vector<double> *LxJair,Vector<short> *FluxDir, Matrix<double> *Klat, Matrix<double> *Kbottom_l,
                     Vector<double> *Kbottom_ch, ALLDATA *adt, Vector<double> *P)
{
    GEOTIMER_SECTION(__func__);

    long i, l, r, c, j, I, R, C, J, sy, syn, ch, cnt=0;
    long n=(Nl+1)*adt->P->total_pixel;
    double dz=0.0, dzn=0.0, dD=0.0, k=0.0, kn=0.0,knaux=0.0,kheat=0.0,kheatn=0.0, kmax=0.0, kmaxn=0.0,krel=0.0,kp=0.0;
    double kn2=0.0,krel2=0.0,kp2=0.0;
    double area, ds=sqrt((*UV->U)(1)*(*UV->U)(2)), dn;
    double ThetaIce,ThetaWat,sat,ThetaAir,Ha;
    double ThetaIce2,ThetaWat2,sat2,ThetaAir2,ns;// GZ - This variables are calculated in the cell where where there is flux exchange
	double Pmi,PmI; //Presion motriz i=current cel, I= cell with flux exchange
    int show=0;
    double fe,knHar,knHarvel,kheatHar,knvel,kn2vel; //Harmonic mean. 
    double FluxTest=0.0;
    double FluxTest2=0.0;
    double KairMult=1000;
    double Cu=0.0;
    short Upwind=0.0;
    double snowD=0.0;
    double SurfaceKincrease=2.0;
/* 

This function fills  LxJair. LxJair is used in energy balance function to calculate the heat fluxes across each boundary. LxJiar*Temp is the total convective heat flux
*/ 
    for (i=1; i<=P->nh; i++)
    {
        /** VERTICAL FLUXES */
        if ( i<=n) //land
        {
            l=(*adt->T->lrc_cont)(i,1);
            r=(*adt->T->lrc_cont)(i,2);
            c=(*adt->T->lrc_cont)(i,3);
            j=adt->T->j_cont[r][c];
            sy=(*adt->S->type)(r,c);

            ch=(*adt->C->ch)(r,c);
            area=ds*ds/cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0);
            if (ch>0) area-=(*adt->C->length)(ch) * adt->P->w_dx *
                            ds; //area of the pixel[m2]

            /** vertical hydraulic conductivity */
            if (l>0)
            {
                ThetaIce=(*adt->S->SS->thi)(l,j);
                ThetaWat=(*adt->S->th)(l,j);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;
                
                kmax=(*adt->S->pa)(sy,jKn,l)/1000.0; //k[mm/s]/1000=k[m/s]

                dz = (*adt->S->pa)(sy,jdz,l)/1000.0;
                
            }

            /** flux from cell below */
            if (l<Nl)
            {
                I = i+1;

                if (l==0) /** overland flow */
                {
                    ns=(*adt->N->S->lnum)(r,c);
                    snowD=DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get());
					//printf(" SNOOOOW DEPTHH, %f, \n",snowD);
					if ((snowD<=adt->P->SnowDepthAirFlowLimit) and (closedtop==0))
                    {
						dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
						dD =  0.5*dzn;
					
                        PmI=((*P)(I)-(*P)(i))+GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref))*(-dD);
                        ThetaIce2=(*adt->S->SS->thi)(l+1,j);
                        ThetaWat2=(*adt->S->th)(l+1,j);
                        sat2=(*adt->S->pa)(sy,jsat,l+1);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                        
                        krel2=kair_from_psi(jKn, (*SL->P)(l+1,j),(*SL->thi)(l+1,j),l+1,adt->S->pa->matrix(sy), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;// units J m/(K N s)
						
						kmaxn=(*adt->S->pa)(sy,jKn,l+1)/1000.0;
						krel=1.0;// krel=relative permeability (0-1) We need to improve this value
                        kp=krel*SurfaceKincrease*kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g);//  k=pkn*mu/(rho*g), kp=permeability
						kn=(kp/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;
						
						kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
						kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
							
						knvel=kp/GTConst::mu_air;
						kn2vel=kp2/GTConst::mu_air;

						knHar=(kn+kn2)/2.0;
						knHarvel=(knvel+kn2vel)/2.0;
						kheatHar=(kheat+kheatn)/2.0;

						cnt++;
						if (Upwind==1){
						if ((*FluxDir)(cnt)>0)
						{
						// upward flux 
							knHar=kn2;
							knHarvel=kn2vel;
							kheatHar=kheatn;
							/*
							if ((((*P)(I)-(*P)(i))-kheatn*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))<0 )
							{
							//printf("AIR DIF SIGN V1 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						}
						else
						{
						// downward flow 
							knHar=kn;
							knHarvel=knvel;
							kheatHar=kheat;
							/*
							if ((((*P)(I)-(*P)(i))-kheatn*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))>0 )
							{
							//printf("AIR DIF SIGN V2 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						}}
						
	
                        if (ThetaAir2>GTConst::ThetaAirMin){
                        (*LxJair)(cnt) = -area*knHar/dD*(((*P)(I)-(*P)(i))-kheatHar*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))*1000.0 ; 
						
						if (fabs(-(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-kheatHar*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0)))>fabs((*adt->AF->MaxVel)(i))){
							(*adt->AF->MaxVel)(i)=-knHarvel/dD*(((*P)(I)-(*P)(i))-kheatHar*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0));
						}
						
						if (fabs((kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-kheatHar*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0)))>fabs((*adt->AF->MaxVel)(I))){
							(*adt->AF->MaxVel)(I)=-knHarvel/dD*(((*P)(I)-(*P)(i))-kheatHar*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0));
						}
						
						
						(*adt->AF->VzB)(i)=knHarvel/dD*(((*P)(I)-(*P)(i))-kheatHar*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0));
						(*adt->AF->VzT)(I)=knHarvel/dD*(((*P)(I)-(*P)(i))-kheatHar*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0));
						Cu=std::max<double>(Cu, (*adt->AF->VzB)(i)*Dt/dD);
						Cu=std::max<double>(Cu, (*adt->AF->VzT)(I)*Dt/dD);
						
						/*
						if ((i==3693)and (show==2)){
					
						FluxTest=FluxTest-area*(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))*1000.0;
						printf("AIR BALANCE VERT cnt,i,I,l,Flux,FluxTOT,B: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,-area*(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0)),FluxTest,(*adt->AF->B)(i));
						
						}
						
						if ((I==3693)and (show==2)){
					
						FluxTest=FluxTest+area*(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))*1000.0;
						printf("AIR BALANCE VERT cnt,i,I,l,Flux,FluxTOT,B: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,+area*(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0)),FluxTest,(*adt->AF->B)(i));
						
						}
						
						if ((i==2)and (show==2)){
					
						FluxTest2=FluxTest2-area*(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))*1000.0;
						printf("AIR BALANCE VERT2 cnt,i,I,l,Flux,FluxTOT,B: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,-area*(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0)),FluxTest2,(*adt->AF->B)(i));
						
						}
						
						if ((I==2)and (show==2)){
						printf( "AIR BALLL VERT2 AUX %f \n",FluxTest2);
						FluxTest2=FluxTest2+area*(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))*1000.0;
						printf("AIR BALANCE VERT2b cnt,i,I,l,Flux,FluxTOT,B: %i, %i, %i, %i,%e,%e,%e  \n",cnt,i,I,l,area*(kp2*krel2/GTConst::mu_air)/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0)),FluxTest2,(*adt->AF->B)(i));
						}*/
						
						}
                        else{
                        (*LxJair)(cnt) =0.0;
                        }
                        /*
                        if (show==1){
                        //printf("AIR BALANCE VERT cnt,i,I,l,LxJair,pmI,pmi,Z: %i, %i, %i,%i ,%e,%f,%f,%f \n",cnt,i,I,l,(*LxJair)(cnt),PmI,((*P)(I)-(*P)(i))+GTConst::rho_ref*GTConst::g*(1-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref))*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000,(*adt->T->Z)(l,r,c)/1000);
						}*/
                        //printf("AIR VERT 111  PI,pi,TI,ti,dD,ThetaAir2,area,kn  %f,%f,%f,%f,%f,%f,%f,%f\n",(*P)(I),(*P)(i),(*adt->AE->T1)(I),(*adt->AE->T1)(i),dD,ThetaAir2,area,kn);
                        //printf("AIR BALANCE  cnt,l,area,kn,dD,LxJair,Vel: %i, %i ,%e,%e,%e,%e,%e \n",cnt,l,area,kn,dD,(*LxJair)(cnt),knaux/dD*(((*P)(I)-Patm)+kheat*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000));

                        //This includes the real direction of the flux.   >0 flow out from the current cell. 
                    }
                    else
                    {

						cnt++;

						(*LxJair)(cnt) = 0.;  //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						/*
                        if (show==1){
                        printf("AIR BALANCE VERT cnt,i,I,l,LxJair,pmI,pmi,Z: %i, %i, %i,%i ,%e,%f,%f,%f \n",cnt,i,I,l,(*LxJair)(cnt),PmI,((*P)(I)-(*P)(i))+GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref))*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0,(*adt->T->Z)(l,r,c)/1000.0);
						}*/
                    }
                    

                }
                else    /** subsurface flow */
                {

                    dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
                    dD = 0.5*dz + 0.5*dzn;
                    fe=0.5*dzn/dD;
                    
                    ThetaIce2=(*adt->S->SS->thi)(l+1,j);
                    ThetaWat2=(*adt->S->th)(l+1,j);
                    sat2=(*adt->S->pa)(sy,jsat,l+1);
                    ThetaAir2=sat2-ThetaIce2-ThetaWat2;
 
                    krel2=kair_from_psi(jKn, (*SL->P)(l+1,j),(*SL->thi)(l+1,j),l,adt->S->pa->matrix(sy), adt->P->k_to_ksat,ThetaAir2)/1000.0;
					kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
					kn2=(kp2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;// units J m/(K N s)

                    kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));

                    krel=kair_from_psi(jKn, (*SL->P)(l,j),(*SL->thi)(l,j),l,adt->S->pa->matrix(sy), adt->P->k_to_ksat,ThetaAir)/1000.0;
                    kp=krel*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
                    kn=(kp/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;// units J m/(K N s)
                    
                    kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                        
                    knvel=kp/GTConst::mu_air;
                    kn2vel=kp2/GTConst::mu_air;
 
					//knHar=(kn+kn2)/2.0;
					knHar=pow(((1.0-fe)/kn+fe/kn2),-1.0);
					//knHarvel=(knvel+kn2vel)/2.0;
					knHarvel=pow(((1.0-fe)/knvel+fe/kn2vel),-1.0);
					kheatHar=(kheat+kheatn)/2.0;
                    
                    cnt++;
                    if (Upwind==1){
						if ((*FluxDir)(cnt)>0)
						{
						// upward flux 
							knHar=kn2;
							knHarvel=kn2vel;
							kheatHar=kheatn;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)<0 )
							{
							//printf("AIR DIF SIGN 1 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						
						}
						else
						{
						// downward flow 
							knHar=kn;
							knHarvel=knvel;
							kheatHar=kheat;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)>0 )
							{
							//printf("AIR DIF SIGN 2 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
					}}
					
                    if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                    (*LxJair)(cnt) = -area/dD*knHar*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)*1000.0 ;  //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
                    
                    if (fabs(-1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0))>fabs((*adt->AF->MaxVel)(i))){
                    (*adt->AF->MaxVel)(i)=-1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0);
					}
					if (fabs(1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0))>fabs((*adt->AF->MaxVel)(I))){
                    (*adt->AF->MaxVel)(I)=1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0);
					}
					
					(*adt->AF->VzB)(i)=1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0);
					(*adt->AF->VzT)(I)=1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0);
					Cu=std::max<double>(Cu, (*adt->AF->VzB)(i)*Dt/dD);
					Cu=std::max<double>(Cu, (*adt->AF->VzT)(I)*Dt/dD);
					
					/*
					if ((i==3693)and (show==2)){
					
						FluxTest=FluxTest-area/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)*1000.0;
						printf("AIR BALANCE VERT cnt,i,I,l,knHarvel,Flux,FluxTOT,B: %i, %i, %i, %i,%e,%e,%e ,%e \n",cnt,i,I,l,knHarvel,-area/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0),FluxTest,(*adt->AF->B)(i));
						
					}
						
					if ((I==3693)and (show==2)){
					
						FluxTest=FluxTest+area/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)*1000.0;
						printf("AIR BALANCE VERT cnt,i,I,l,knHarvel,Flux,FluxTOT,B: %i, %i, %i, %i,%e,%e,%e,%e  \n",cnt,i,I,l,knHarvel,area/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0),FluxTest,(*adt->AF->B)(I));
						
					}
					
					if ((i==2)and (show==2)){
					
						FluxTest2=FluxTest2-area/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)*1000.0;
						printf("AIR BALANCE VERT2 cnt,i,I,l,knHarvel,Flux,FluxTOT,B: %i, %i, %i, %i,%e,%e,%e,%e  \n",cnt,i,I,l,knHarvel,-area/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0),FluxTest2,(*adt->AF->B)(i));
						
					}
						
					if ((I==21)and (show==2)){
					
						FluxTest2=FluxTest2+area/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0)*1000.0;
						printf("AIR BALANCE VERT2 cnt,i,I,l,knHarvel,Flux,FluxTOT,B: %i, %i, %i, %i,%e,%e,%e,%e  \n",cnt,i,I,l,knHarvel,area/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0),FluxTest2,(*adt->AF->B)(I));
						
					}*/
					
					}
                    else{
                        (*LxJair)(cnt) =0.0;
                    }
                    
                    //printf("AIR BALANCE  cnt,l,area,kn,dD,LxJair,Vel: %i, %i ,%e,%e,%e,%e,%e \n",cnt,l,area,kn,dD,(*LxJair)(cnt),knaux/dD*(((*P)(I)-(*P)(i))+kheat*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000));
					/*
					if (show==1){
                    printf("AIR BALANCE VERT cnt,i,I,l,LxJair,pmI,pmi,Z: %i, %i, %i,%i ,%e,%f,%f,%f \n",cnt,i,I,l,(*LxJair)(cnt),PmI,((*P)(I)-(*P)(i))+GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref))*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0,(*adt->T->Z)(l,r,c)/1000.0);
					}*/
                    //printf("AIR BALANCE VERT cnt,i,I,l,PI,Pi,RGTI,TI,RGTi,Ti: %i, %i, %i,%i ,%f, %f, %f, %f , %f, %f \n",cnt,i,I,l,(*P)(I),(*P)(i),GTConst::rho_ref*GTConst::g*(1-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref)),(*adt->AE->T1)(I),GTConst::rho_ref*GTConst::g*(1-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref)),(*adt->AE->T1)(i));
					//printf("AIR BALANCE VERT  DT:%f  dD :%f\n", ((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000, dD);
                    //This includes the real direction of the flux.   >0 flow out from the current cell. 
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
            printf("AIR CHANNEL i,l %i, %i\n",i,l);
            area=(*adt->C->length)(ch) * adt->P->w_dx * ds;

            /** vertical hydraulic conductivity */
            if (l>0)
            {
                
                ThetaIce=(*adt->C->SS->thi)(l,ch);
                ThetaWat=(*adt->C->th)(l,ch);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;
                kmax=(*adt->S->pa)(sy,jKn,l)/1000.0;
               
                dz = (*adt->S->pa)(sy,jdz,l)/1000.0;
            }

            /** flux from cell below */
            if (l<Nl)
            {

                I = i+1;

                if (l==0) /** OVERLAND flow */
                {
                    ns=(*adt->N->S->lnum)(r,c);
					if (ns==0)
					{
						dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
						dD = 0.5*dzn;

                        ThetaIce2=(*adt->C->SS->thi)(l+1,ch);
                        ThetaWat2=(*adt->C->th)(l+1,ch);
                        sat2=(*adt->S->pa)(sy,jsat,l+1);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                        
                        
                        kmaxn=(*adt->S->pa)(sy,jKn,l+1)/1000.0;
                        kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
                        krel2=ThetaAir2/sat2; // krel=relative permeability (0-1) We need to improve this value
                        kn2=(kp2*krel2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;// units J m/(K N s)
                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                    
                        kn=(kp2*krel2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;// units J m/(K N s)
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));

						cnt++;

                        if (ThetaAir2>GTConst::ThetaAirMin){
                        (*LxJair)(cnt) = -area*(kn2+kn)/2.0/dD*(((*P)(I)-(*P)(i))-(kheat+kheatn)/2.0*dD*cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0))*1000.0 ;  //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
                        }
                        else{
                        (*LxJair)(cnt) =0.0;
                        }
					}
					else
					{
						kn=0.;
						cnt++;
						(*LxJair)(cnt) = 0.0;  //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
					}

                }
                else    /** SUBSURFACE flow */
                {

                    dzn = (*adt->S->pa)(sy,jdz,l+1)/1000.0;
                    dD = 0.5*dz + 0.5*dzn;
                    
                    ThetaIce2=(*adt->C->SS->thi)(l+1,ch);
                    ThetaWat2=(*adt->C->th)(l+1,ch);
                    sat2=(*adt->S->pa)(sy,jsat,l+1);
                    ThetaAir2=sat2-ThetaIce2-ThetaWat2;
                    fe=0.5*dzn/dD;
                    
                    kmaxn=(*adt->S->pa)(sy,jKn,l+1)/1000.0;
                    kp2=kmaxn*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
                    krel2=ThetaAir2/sat2; // krel=relative permeability (0-1) We need to improve this value
                    //kn2=(kp2*krel2/GTConst::mu_air)*GTConst::rho_ref*(1-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref))*GTConst::c_air;
                    kn2=(kp2*krel2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;
                    kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));
                
                    kp=kmax*GTConst::mu_l/(GTConst::rho_w*GTConst::g);//  k=pkn*mu/(rho*g), kp=permeability
                    krel=ThetaAir/sat;// krel=relative permeability (0-1) We need to improve this value
                    //kn=(kp*krel/GTConst::mu_air)*GTConst::rho_ref*(1-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref))*GTConst::c_air;
                    kn=(kp*krel/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;
                    kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                        
                    knvel=kp*krel/GTConst::mu_air;
					kn2vel=kp2*krel2/GTConst::mu_air;
										
					knHar=pow(((1.0-fe)/kn+fe/kn2),-1.0);
					knHarvel=pow(((1.0-fe)/knvel+fe/kn2vel),-1.0);
					kheatHar=pow(((1.0-fe)/kheat+fe/kheatn),-1.0);

					cnt++;
					if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
					(*LxJair)(cnt) = -area/dD*knHar*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0) ;  //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]


					if (fabs(-1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0))>fabs((*adt->AF->MaxVel)(i))){
					(*adt->AF->MaxVel)(i)=-1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0);
					}
					if (fabs(1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0))>fabs((*adt->AF->MaxVel)(I))){
					(*adt->AF->MaxVel)(I)=1.0/dD*knHarvel*(((*P)(I)-(*P)(i))+kheatHar*((*adt->T->Z)(l+1,r,c)-(*adt->T->Z)(l,r,c))/1000.0);
					}
					}
					else{
						(*LxJair)(cnt) =0.0;
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
 
				krel=kair_from_psi(jKl, (*SL->P)(l,j),(*SL->thi)(l,j),l,adt->S->pa->matrix(sy), adt->P->k_to_ksat,ThetaAir)/1000.0;
				kp=krel*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
				k=(kp/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;
				
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
    
                        krel2=kair_from_psi(jKl, (*SL->P)(l,J),(*SL->thi)(l,J),l,adt->S->pa->matrix(syn), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;
						
                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));

                        kn = k;
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                            
                        knvel=kp/GTConst::mu_air;
						kn2vel=kp2/GTConst::mu_air;
						
						//knHar=(kn+kn2)/2.0;
						knHar=2.0*kn*kn2/(kn+kn2);
						//knHarvel=(knvel+kn2vel)/2.0;
						knHarvel=2.0*knvel*kn2vel/(knvel+kn2vel);
						kheatHar=(kheat+kheatn)/2.0;
						
                        cnt++;
                        if (Upwind==1){
						if ((*FluxDir)(cnt)>0)
						{
						// upward flux 
							knHar=kn2;
							knHarvel=kn2vel;
							kheatHar=kheatn;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))<0 )
							{
							//printf("AIR DIF SIGN 1 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						
						}
						else
						{
						// downward flow 
							knHar=kn;
							knHarvel=knvel;
							kheatHar=kheat;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0 )
							{
							//printf("AIR DIF SIGN 2 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						}}
                        
                        if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                        
                           (*LxJair)(cnt) = -(dn*dz)/dD *knHar*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0 ; //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
							
							if (fabs(-1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))))>fabs((*adt->AF->MaxVel)(i))){
							(*adt->AF->MaxVel)(i)= -1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							}
							if (fabs(1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))))>fabs((*adt->AF->MaxVel)(I))){
							(*adt->AF->MaxVel)(I)= 1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							}

							(*adt->AF->VyS)(i)=1./dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							(*adt->AF->VyN)(I)=1./dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							Cu=std::max<double>(Cu, (*adt->AF->VyS)(i)*Dt/ds);
							Cu=std::max<double>(Cu, (*adt->AF->VyN)(I)*Dt/ds);

							/*
							if ((i==3693)and (show==2)){
						
								FluxTest=FluxTest-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat1 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest);
							
							}
								
							if ((I==3693)and (show==2)){
							
								FluxTest=FluxTest+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat1 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest);
							
							}
							
							if ((i==2) and (show==2)){
						
								FluxTest2=FluxTest2-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat1 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest2);
							
							}
								
							if ((I==2) and (show==2)){
							
								FluxTest2=FluxTest2+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat1 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest2);
							
							}*/
                        }
                        else{
                                (*LxJair)(cnt) =0.0;
                        }

                    }
                    else
                    {

                        kn = 0.0;
                        cnt++;
                        (*LxJair)(cnt) =0.0;

                    }
                    /*
                    if (show==1){
                    printf("AIR BALANCE LAT1 cnt,i,I,l,FLUX,pmI,pmi,Z: %i,%i, %i,%i ,%e,%e,%e %f \n",cnt,i,I,l,-(dn*dz)/dD *(kp*krel/GTConst::mu_air)*((*P)(I)-(*P)(i)),(*P)(I),(*P)(i),(*adt->T->Z)(l,r,c)/1000.0);
					}*/
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

                        krel2=kair_from_psi(jKl, (*SL->P)(l,J),(*SL->thi)(l,J),l,adt->S->pa->matrix(syn), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;
						
                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));

                        kn = k;
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                            
                        knvel=kp/GTConst::mu_air;
						kn2vel=kp2/GTConst::mu_air;
						
						//knHar=(kn+kn2)/2.0;
						knHar=2.0*kn*kn2/(kn+kn2);
						//knHarvel=(knvel+kn2vel)/2.0;
						knHarvel=2.0*knvel*kn2vel/(knvel+kn2vel);
						kheatHar=(kheat+kheatn)/2.0;
						
						
                        cnt++;
                        if (Upwind==1){
						if ((*FluxDir)(cnt)>0)
						{
						// upward flux 
							knHar=kn2;
							knHarvel=kn2vel;
							kheatHar=kheatn;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))<0 )
							{
							//printf("AIR DIF SIGN 1 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						
						}
						else
						{
						// downward flow 
							knHar=kn;
							knHarvel=knvel;
							kheatHar=kheat;
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0 )
							{
							//printf("AIR DIF SIGN 2 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}
						}}

                        if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                        
                           (*LxJair)(cnt) =- (dn*dz)/dD *knHar*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0 ; //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
							
							if (fabs(-1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))))>fabs((*adt->AF->MaxVel)(i))){
							(*adt->AF->MaxVel)(i)= -1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							}
							if (fabs(1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))))>fabs((*adt->AF->MaxVel)(I))){
							(*adt->AF->MaxVel)(I)= 1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							}
							

							(*adt->AF->VyN)(i)=-1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							(*adt->AF->VyS)(I)=-1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
						
							Cu=std::max<double>(Cu, (*adt->AF->VyN)(i)*Dt/ds);
							Cu=std::max<double>(Cu, (*adt->AF->VyS)(I)*Dt/ds);
							/*
							if ((i==3693)and (show==2)){
						
								FluxTest=FluxTest-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat2 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest);
							
							}
								
							if ((I==3693)and (show==2)){
							
								FluxTest=FluxTest+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat2 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest);
							
							}
							
							if ((i==2) and (show==2)){
						
								FluxTest2=FluxTest2-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat2 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest2);
							
							}
								
							if ((I==2) and (show==2)){
							
								FluxTest2=FluxTest2+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat2 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest2);
							
							}*/

                        }
                        else{
                                (*LxJair)(cnt) =0.0;
                        }

                    }
                    else
                    {

                        kn = 0.0;
                        cnt++;
                        (*LxJair)(cnt) =0.0;

                    }
                    /*
                    if (show==1){
                    printf("AIR BALANCE LAT2 cnt,i,I,l,LxJair,pi,pI,Z: %i,%i, %i,%i ,%e,%e,%e,%f \n",cnt,i,I,l,(*LxJair)(cnt),PmI,((*P)(I)-(*P)(i))+GTConst::rho_ref*GTConst::g*(1-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref))*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)),(*adt->T->Z)(l,r,c)/1000);
					}*/
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

                    dD = find_3Ddistance(ds,(*adt->T->Z0)(r,c)-(*adt->T->Z0)(R,C)) ; /// [m]
                    dn = ds/cos(0.5*atan((*adt->T->dzdE)(r,c))+0.5*atan((*adt->T->dzdE)(R,C))); /// [m]

                    if (l>0)
                    {
                        ThetaIce2=(*adt->S->SS->thi)(l,J);
                        ThetaWat2=(*adt->S->th)(l,J);
                        sat2=(*adt->S->pa)(syn,jsat,l);
                        ThetaAir2=sat2-ThetaIce2-ThetaWat2;

                        krel2=kair_from_psi(jKl, (*SL->P)(l,J),(*SL->thi)(l,J),l,adt->S->pa->matrix(syn), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;
						
                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));

                        kn = k;
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                            
                        knvel=kp/GTConst::mu_air;
						kn2vel=kp2/GTConst::mu_air;
						
						//knHar=(kn+kn2)/2.0;
						knHar=2.0*kn*kn2/(kn+kn2);
						//knHarvel=(knvel+kn2vel)/2.0;
						knHarvel=2.0*knvel*kn2vel/(knvel+kn2vel);
						kheatHar=(kheat+kheatn)/2.0;
						
						
                        cnt++;
                        if (Upwind==1){
						if ((*FluxDir)(cnt)>0)
						{
						// upward flux 
							knHar=kn2;
							knHarvel=kn2vel;
							kheatHar=kheatn;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))<0 )
							{
							//printf("AIR DIF SIGN 1 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						
						}
						else
						{
						// downward flow 
							knHar=kn;
							knHarvel=knvel;
							kheatHar=kheat;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0 )
							{
							//printf("AIR DIF SIGN 2 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						}}

                        if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                        
                           (*LxJair)(cnt) = -(dn*dz)/dD *knHar*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0 ; //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
							
							if (fabs(-1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))))>fabs((*adt->AF->MaxVel)(i))){
							(*adt->AF->MaxVel)(i)= -1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							}
							if (fabs(1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))))>fabs((*adt->AF->MaxVel)(I))){
							(*adt->AF->MaxVel)(I)= 1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							}
							
							(*adt->AF->VxW)(i)=1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							(*adt->AF->VxE)(I)=1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							
							Cu=std::max<double>(Cu, (*adt->AF->VxW)(i)*Dt/ds);
							Cu=std::max<double>(Cu, (*adt->AF->VxE)(I)*Dt/ds);
			
							/*
							if ((i==3693)and (show==2)){
						
								FluxTest=FluxTest-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat3cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest);
							
							}
								
							if ((I==3693)and (show==2)){
							
								FluxTest=FluxTest+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat3 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest);
							
							}
							
							if ((i==2) and (show==2)){
						
								FluxTest2=FluxTest2-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat3 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest2);
							
							}
								
							if ((I==2) and (show==2)){
							
								FluxTest2=FluxTest2+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat3 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest2);
							
							}*/
                        }
                        else{
                                (*LxJair)(cnt) =0.0;
                        }

                    }
                    else
                    {

                        kn = 0.;
                        cnt++;
                        (*LxJair)(cnt) =0.;

                    }
                    /*
                    if (show==1){
                    printf("AIR BALANCE LAT3 cnt,i,I,l,LxJair,pi,pI,Z: %i,%i, %i,%i ,%e,%e,%e,%f \n",cnt,i,I,l,(*LxJair)(cnt),PmI,((*P)(I)-(*P)(i))+GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref))*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)),(*adt->T->Z)(l,r,c)/1000.0);
					}*/
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

                        krel2=kair_from_psi(jKl, (*SL->P)(l,J),(*SL->thi)(l,J),l,adt->S->pa->matrix(syn), adt->P->k_to_ksat,ThetaAir2)/1000.0;
						kp2=krel2*GTConst::mu_l/(GTConst::rho_w*GTConst::g); //  k=kn*mu/(rho*g), k=permeability
						kn2=(kp2/GTConst::mu_air)*GTConst::rho_ref*GTConst::c_air;
						
                        kheatn=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(I)-GTConst::Tref));

                        kn = k;
                        kheat=GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref));
                            
                        knvel=kp/GTConst::mu_air;
						kn2vel=kp2/GTConst::mu_air;
						
						//knHar=(kn+kn2)/2.0;
						knHar=2.0*kn*kn2/(kn+kn2);
						//knHarvel=(knvel+kn2vel)/2.0;
						knHarvel=2.0*knvel*kn2vel/(knvel+kn2vel);
						kheatHar=(kheat+kheatn)/2.0;
						
						
                        cnt++;
                        if (Upwind==1){
						if ((*FluxDir)(cnt)>0)
						{
						// upward flux 
							knHar=kn2;
							knHarvel=kn2vel;
							kheatHar=kheatn;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))<0 )
							{
							//printf("AIR DIF SIGN 1 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						
						}
						else
						{
						// downward flow 
							knHar=kn;
							knHarvel=knvel;
							kheatHar=kheat;
							/*
							if ((((*P)(I)-(*P)(i))+kheatn*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))>0 )
							{
							//printf("AIR DIF SIGN 2 cnt,i,I,l,: %i, %i, %i, %i, \n",cnt,i,I,l);
							}*/
						}}
                        
                        if ((ThetaAir2>GTConst::ThetaAirMin) and (ThetaAir>GTConst::ThetaAirMin)){
                        
                           (*LxJair)(cnt) = -(dn*dz)/dD *knHar*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0 ; //Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
							
							if (fabs(-1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))))>fabs((*adt->AF->MaxVel)(i))){
							(*adt->AF->MaxVel)(i)= -1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							}
							if (fabs(1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))))>fabs((*adt->AF->MaxVel)(I))){
							(*adt->AF->MaxVel)(I)= 1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							}
							
							
							(*adt->AF->VxE)(i)=-1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							(*adt->AF->VxW)(I)=-1.0/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)));
							
							Cu=std::max<double>(Cu, (*adt->AF->VxE)(i)*Dt/ds);
							Cu=std::max<double>(Cu, (*adt->AF->VxW)(I)*Dt/ds);
							/*				
							if ((i==3693)and (show==2)){
						
								FluxTest=FluxTest-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat4 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest);
							
							}
								
							if ((I==3693)and (show==2)){
							
								FluxTest=FluxTest+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat4 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest);
							
							}
							
							if ((i==2) and (show==2)){
						
								FluxTest2=FluxTest2-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat4 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,-(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest2);
							
							}
								
							if ((I==2) and (show==2)){
							
								FluxTest2=FluxTest2+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)))*1000.0;
								printf("AIR BALANCE lat4 cnt,i,I,l,knHarvel,Flux,FluxTOT: %i, %i, %i, %i,%e,%e,%e \n",cnt,i,I,l,knHarvel,+(dn*dz)/dD *knHarvel*((*P)(I)-(*P)(i)+kheatHar*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c))),FluxTest2);
							
							}*/
                        }
                        else{
                                (*LxJair)(cnt) =0.0;
                        }

                    }
                    else
                    {

                        kn = 0.0;
                        cnt++;
                        (*LxJair)(cnt) =0.0;
                    }
                    /*
                    if (show==1){
                    printf("AIR BALANCE LAT4 cnt,i,I,l,LxJair,pi,pI,Z: %i,%i, %i,%i ,%e,%e,%e%f \n",cnt,i,I,l,(*LxJair)(cnt),PmI,((*P)(I)-(*P)(i))+GTConst::rho_ref*GTConst::g*(1.0-GTConst::beta*( (*adt->AE->T1)(i)-GTConst::Tref))*((*adt->T->Z0)(R,C)-(*adt->T->Z0)(r,c)),(*adt->T->Z)(l,r,c)/1000.0);
					}*/
                }
            }

            /** exchange with channels */
            if (l>0 && (*adt->C->ch)(r,c) > 0)
            {

                ch = (*adt->C->ch)(r,c);
                syn = (*adt->C->soil_type)(ch);
                I = n + adt->C->ch3[l][ch];

                kn = 0.;

                /** Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
                 * Area[m2] = 2 * channel length * layer thickness = 2 * length[m] * dz[mm] * 1.E-3[m/mm]
                 * dD[mm] = 0.25 * ds * (1+w_dx) */
                dD = find_3Ddistance(ds * (1.0 + adt->P->w_dx) / 4.0,
                                     1.E-3*adt->P->depr_channel) * 1.E3;//[mm]


                cnt++;
                (*LxJair)(cnt) =0.;
            }
        }
    }
    (adt->AF->Courant)=Cu;
    return 0;
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_fair_3D(double Dt, Vector<double> *f, ALLDATA *adt, SOIL_STATE *L,
              SOIL_STATE *C,STATEVAR_3D *S, Vector<double> *P)
{
  GEOTIMER_SECTION(__func__);
  
  
	// In this funcion is included ther term asociated to the water and ice content, the flux term associated to the buoyancy and the border conditions.
    long i,I;
    long n=(Nl+1)*adt->P->total_pixel;

    
    //double ThetaAirMin=0.01;

#pragma omp parallel for
    for (i=1; i<=P->nh; i++)
    {
        long  l, r, c, j, sy, ch, bc;
        double dz, dn,dzn, dD, V0, V1, totw1=0.0, totw0=0.0, ice=0.0;
        double area, ds=sqrt((*UV->U)(1)*(*UV->U)(2));
        double ThetaIce,ThetaWat,sat,ThetaAir,ns,kmax,krel,kp,kn,snowD=0.0;
        const double PI = 3.141592653589793238463;

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
            {
                area -= (*adt->C->length)(ch) * adt->P->w_dx * ds; //area of the pixel[m2]
            }

            if (l>0)
            {
                
                ThetaIce=(*adt->S->SS->thi)(l,j);
                ThetaWat=(*adt->S->th)(l,j);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;
                //ThetaAir=std::max<double>(ThetaAir, ThetaAirMin);
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
            ns=(*S->lnum)(r,c);
            if (l>0)
            {
                //totw1=(*adt->C->th)(l,ch)+(*C->thi)(l,ch);
                //totw0=(*C->totw0)(l,ch);
                ThetaIce=(*adt->C->SS->thi)(l,ch);
                ThetaWat=(*adt->C->th)(l,ch);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;

			}
        }

        /** hydraulic capacity (diagonal term) */
        if (l==0)
        {
 
            //(*f)(i)-=GTConst::Pa0*100.0;
            snowD=DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get());
				
			if ((snowD>adt->P->SnowDepthAirFlowLimit) and (closedtop==0))
			{
				(*f)(i)=(*P)(i);
				//(*f)(i)-=(GTConst::Pa0*100.0-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000.0);
				(*f)(i)-=(*adt->AF->P0)(i);
			}
			
			if (closedtop==1){
				(*f)(i)=(*P)(i);
				//(*f)(i)-=(GTConst::Pa0*100.0-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000.0);
				(*f)(i)-=(*adt->AF->P0)(i);
			}
			
			
		}
        /*
        else
        {

            
            dz = (*adt->S->pa)(sy,jdz,l);
            V1 = area*dz * totw1;
            V0 = area*dz * totw0;
            //printf("AIR BALANCE  i,l,totw0,totw1: %i, %d ,%f,%f\n",i,l,totw0,totw1);
           
            if (ThetaAir>GTConst::ThetaAirMin){
            //(*f)(i) = 0;
           
            //(*f)(i) = -(V1-V0)/Dt;
            
            if ((ns>0) and (i==Nl+1)){

            //(*f)(i) -= (GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000)*area*kn/dz*10000;
            
            //(*f)(i)=(*P)(i)*(area*kn/dz)/1e3;
            //(*f)(i)-=(*adt->AF->P0)(i)*(area*kn/dz)/1e3;
            

            (*f)(i)=(*P)(i)/1e7;
            (*f)(i)-=(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000)/1e7;
            
            //(*f)(i)= erf(((*P)(i)-(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000))/1e7)*(adt->P->TolVWb*500);
            //(*f)(i)=pow((*P)(i)-(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000),2)*pow(adt->P->TolVWb,2)/100;
            //(*f)(i)=pow(((*P)(i)-(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000))/100000,4)*adt->P->TolVWb*100;
            //printf("AIR BALANCE  i,l,ns,(*adt->AF->P0)(i): %i, %d ,%f,%f\n",i,l,ns,(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000));
            }
            }
            else
            {
            (*f)(i)=(*P)(i);
            (*f)(i)-=(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000);
            //(*f)(i)-=(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*((*adt->T->Z)(l,r,c)));
            }
            
            
 
        }*/


    } 
	
    return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_dfdPair_3D(double Dt, Vector<double> *df, ALLDATA *adt, SOIL_STATE *L,
                 SOIL_STATE *C,STATEVAR_3D *S, Vector<double> *P)
{
    long i,I, l, r, c, j, sy, ch, bc;
    long n=(Nl+1)*adt->P->total_pixel;
    double dz, dn, dD, psi1, ice=0.0;
    double area, ds=sqrt((*UV->U)(1)*(*UV->U)(2));
    double ThetaIce,ThetaWat,sat,ThetaAir,ns,snowD=0.0;
    const double PI = 3.141592653589793238463;
   

    for (i=1; i<=P->nh; i++)
    {

        (*df)(i) = 0.;
        
        if (i<=n)
        {
            l = (*adt->T->lrc_cont)(i,1);
            r = (*adt->T->lrc_cont)(i,2);
            c = (*adt->T->lrc_cont)(i,3);
            j = adt->T->j_cont[r][c];
            sy = (*adt->S->type)(r,c);
            ch = (*adt->C->ch)(r,c);
            area = ds*ds/cos((*adt->T->slope)(r,c)*GTConst::Pi/180.0);
            ns=(*S->lnum)(r,c);
            if (ch>0)
                area-= (*adt->C->length)(ch) * adt->P->w_dx * ds; //area of the pixel[m2]
                
            if (l>0)
            {
                ThetaIce=(*adt->S->SS->thi)(l,j);
                ThetaWat=(*adt->S->th)(l,j);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;
             }
        }
        else
        {
            l = (*adt->C->lch)(i-n,1);
            ch = (*adt->C->lch)(i-n,2);
            r = (*adt->C->r)(ch);
            c = (*adt->C->c)(ch);
            sy = (*adt->C->soil_type)(ch);
			ns=(*S->lnum)(r,c);
            area=(*adt->C->length)(ch) * adt->P->w_dx * ds;
            if (l>0)
            {
                ThetaIce=(*adt->C->SS->thi)(l,ch);
                ThetaWat=(*adt->C->th)(l,ch);
                sat=(*adt->S->pa)(sy,jsat,l);
                ThetaAir=sat-ThetaIce-ThetaWat;
                //ThetaAir=std::max<double>(ThetaAir, ThetaAirMin);
			}
 
        }
        
        if (l==0)
        {
            snowD=DEPTH(r, c, adt->N->S->lnum.get(), adt->N->S->Dzl.get());
				
			if ((snowD>adt->P->SnowDepthAirFlowLimit) and (closedtop==0))
			{
				(*df)(i) += 1.0;
			}
			if (closedtop==1){
				(*df)(i) += 1.0;
			}

        }
        /*
        else
        {
            
			dz = (*adt->S->pa)(sy,jdz,l);
            if (ThetaAir<=GTConst::ThetaAirMin){
            (*df)(i) += 1.0;
            //printf("AIR BALANCE  i,l,ThetaAir1: %i, %d ,%f\n",i,l,ThetaAir);
            }
            else{
            
            if ((ns>0) and (i==Nl+1)){
                            
            //(*f)(i) -= (GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000)*area*kn/dz*10000;
            
			
 
            (*df)(i)=1.0/1e7;
            //(*df)(i)+=2* adt->P->TolVWb*500 *exp(-pow((*P)(i)-(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000),2)/pow(1e7,2))/(1e7* sqrt(PI));
            
            //(*df)(i)+=2*pow(adt->P->TolVWb,2)/100*((*P)(i)-(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000));
            //(*df)(i)+=4*adt->P->TolVWb*100*pow((*P)(i)-(GTConst::Pa0*100-GTConst::rho_ref*GTConst::g*(*adt->T->Z)(l,r,c)/1000),3)/pow(100000,4);
            
            //printf("AIR BALANCE222  i,l,ns: %i, %d ,%f\n",i,l,ns);
            //printf("HEREEE2\n");
            }
            }
            
            
            
            
        }*/
        
			

    }
    return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
