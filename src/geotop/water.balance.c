
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


#include "constants.h"
#include "struct.geotop.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "util_math.h"
#include "water.balance.h"
#include "meteodata.h"

#include <time.h>

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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short water_balance(double Dt, double JD0, double JD1, double JD2, SOIL_STATE *L, SOIL_STATE *C, ALLDATA *adt, DOUBLEVECTOR *Vsub, DOUBLEVECTOR *Vsup, 
					double *Voutnet, double *Voutlandsub, double *Voutlandsup, double *Voutlandbottom){
	
	clock_t start, end;
	FILE *flog;
	double Pnet, loss;
	long j;
	short a;
	
	double mm1, mm2, mmo;

	//double m1=0., m2=0., mo=0.;
	//double ds, area, dz;
	//long r, c, l, sy;
	
	flog = fopen(logfile, "a");
	
	if(adt->P->qin==1){
		time_interp_linear(JD0, JD1, JD2, adt->M->qinv, adt->M->qins, adt->M->qinsnr, 2, 0, 0, &(adt->M->qinline));		
	}
			
	if (adt->P->point_sim != 1) {//distributed simulations
		
		//surface flow: 1st half of time step
		start=clock();
		supflow(adt->P->DDland, adt->P->DDchannel, Dt/2., adt->I->time, L->P->co[0], adt->W->h_sup->co, C->P->co[0], adt->C->h_sup->co, adt->T, adt->L, adt->W, adt->C, adt->P, adt->M, Vsup, Voutnet, Voutlandsup, flog, &mm1, &mm2, &mmo);
		end=clock();
		
		/*MMo += mmo;
		MS1 = mm1;
		MS2 = mm2;	
		MM1 = 0.;
		MM2 = 0.;
		printf("%e %e %e %e %e %e\n",MM1,MM2,MMR,MS1,MS2,MMo);*/

		t_sup += (end-start)/(double)CLOCKS_PER_SEC;	
		
		//subsurface flow with time step Dt0 (decreasing if not converging)
		start = clock();
		
		/*ds=sqrt(UV->U->co[1]*UV->U->co[2]);
		for (j=1; j<=adt->W->H1->nh; j++) {
			l=adt->T->lrc_cont->co[j][1];
			r=adt->T->lrc_cont->co[j][2];
			c=adt->T->lrc_cont->co[j][3];
			sy=adt->S->type->co[r][c];
			area=ds*ds/cos(adt->T->slope->co[r][c]*Pi/180.);
			if(l==0){
				m1 += area * 1.E-3*Fmax(0.0, L->P->co[l][adt->T->j_cont[r][c]]) / cos(adt->T->slope->co[r][c]*Pi/180.);
			}else {
				dz = adt->S->pa->co[sy][jdz][l];		
				m1 += area*1.E-3*dz * theta_from_psi(L->P->co[l][adt->T->j_cont[r][c]], 0, l, adt->S->pa->co[sy], PsiMin);
				
			}
			if(l==0) mo += area * 1.E-3 * adt->W->Pnet->co[r][c];
		}*/
			
		
		a = Richards3D(Dt, L, C, adt, flog, &loss, Vsub, Voutlandbottom, Voutlandsub, &Pnet, adt->P->UpdateK);
		end=clock();
		t_sub += (end-start)/(double)CLOCKS_PER_SEC;
		if (a != 0){
			fclose(flog);
			return 1;
		}
				
		/*ds=sqrt(UV->U->co[1]*UV->U->co[2]);
		for (j=1; j<=adt->W->H1->nh; j++) {
			l=adt->T->lrc_cont->co[j][1];
			r=adt->T->lrc_cont->co[j][2];
			c=adt->T->lrc_cont->co[j][3];
			sy=adt->S->type->co[r][c];
			area=ds*ds/cos(adt->T->slope->co[r][c]*Pi/180.);
			if(l==0){
				m2 += area * 1.E-3*Fmax(0.0, L->P->co[l][adt->T->j_cont[r][c]]) / cos(adt->T->slope->co[r][c]*Pi/180.);
			}else {
				dz = adt->S->pa->co[sy][jdz][l];		
				m2 += area*1.E-3*dz * theta_from_psi(L->P->co[l][adt->T->j_cont[r][c]], 0, l, adt->S->pa->co[sy], PsiMin);
			}
		}
		
		printf("SUB: m1:%e m2:%e mo:%e dm:%e\n",m1,m2,mo,fabs(m1+mo-m2));
		
		MM1 = m1;
		MM2 = m2;
		MMR += mo;*/
		
		//surface flow: 2nd half of time step
		start=clock();
		supflow(adt->P->DDland, adt->P->DDchannel, Dt/2., adt->I->time, L->P->co[0], adt->W->h_sup->co, C->P->co[0], adt->C->h_sup->co, adt->T, adt->L, adt->W, adt->C, adt->P, adt->M, Vsup, Voutnet, Voutlandsup, flog, &mm1, &mm2, &mmo);
		end=clock();
		
		/*MMo += mmo;
		MS1 = mm1;
		MS2 = mm2;
		printf("%e %e %e %e %e %e\n",MM1,MM2,MMR,MS1,MS2,MMo);*/

		t_sup += (end-start)/(double)CLOCKS_PER_SEC;	

	}else {//point simulations
		
		start = clock();
		for (j=1; j<=adt->P->total_pixel; j++) {	
			a = Richards1D(j, Dt, L, adt, flog, &loss, Voutlandbottom, Voutlandsub, &Pnet, adt->P->UpdateK);
			if (a != 0){
				fclose(flog);
				return 1;
			}			
			if( L->P->co[0][j] > 0 ) L->P->co[0][j] = Fmin( L->P->co[0][j], Fmax(0.,-adt->T->BC_DepthFreeSurface->co[j])*cos(adt->T->slope->co[1][j]*Pi/180.) );
		}
		end=clock();
		t_sub += (end-start)/(double)CLOCKS_PER_SEC;
		
	}
	
	
	odb[oopnet] = Pnet;
	odb[oomasserror] = loss;	
	
	fclose(flog);
	
	return 0;
	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/* The code solves non-linear system F(H), where H is the total head (piezometric) [mm], F is vectorial function with components equal to the number of points,
 
 H = P + z, where P is the pressure [mm] and z the elevation over a reference height (given by the average elevation of the domain)
 
 F(H) comes out from the temporal discretization of Richards' equation PDE
 
 In this code: F(H) = I*f(H) + K(H)*H,
 
 where I is the identity matrix, f(H) = volume storage + sink, is a vectorial function (not a matrix), and K(H) is the hydraulic conductivity matrix
 
 K(H) has a number of interesting properties:
 
 1) is symmetrical 
 
 2) if kij is a component kii = sum(j!=i) kji
 
 3) the sum of all the components of the matrix is 0
 
 4) if T is a vector, the vector K T has the properties that the sum of its components is 0
 
 Here it is possible to update or not K with the new values of H, In any case the derivative of K(H) is not considered in the Jacobian. So the matrix K(H) is also
 
 used in the JACOBIAN.
 
 K is described storing only its strict lower component (strict = without diagonal) with the 3 vectors Li, Lp, Lx (in the same way as UFMPACK) */


short Richards3D(double Dt, SOIL_STATE *L, SOIL_STATE *C, ALLDATA *adt, FILE *flog, double *loss, DOUBLEVECTOR *Vsub, double *Vbottom, double *Vlatsub, double *Total_Pnet, short updateK){
	
	double res=0.0, res0[3], res_prev[MM], res_av, res00, lambda[3], epsilon, mu=0., hnew, hold=0.;
	double ds=sqrt(UV->U->co[1]*UV->U->co[2]), area, dz, dn, dD;	
	long i, j, ch, l, r, c, m, bc, sy, cont, cont2, iter;
	long n=adt->T->lrc_cont->nrh;	
	long N=adt->W->H0->nh;
	long cont_lambda_min=0;
	short out, out2;	
	int sux;
	
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
	
	for(i=1; i<=N; i++){
		
		if (i<=n) {
			
			l = adt->T->lrc_cont->co[i][1];
			r = adt->T->lrc_cont->co[i][2];
			c = adt->T->lrc_cont->co[i][3];
			j = adt->T->j_cont[r][c];
			
			//precipitation
			if (l == 0 && adt->W->Pnet->co[r][c] > 0) {
				*Total_Pnet = *Total_Pnet + (adt->W->Pnet->co[r][c]/cos(adt->T->slope->co[r][c]*Pi/180.)) / (double)adt->P->total_pixel;
			}
			
			//solution guess
			if (adt->W->Pnet->co[r][c] > 0 && l == 0) {
				adt->W->H1->co[i] = Fmax(0., L->P->co[l][j]) + (adt->W->Pnet->co[r][c]/cos(adt->T->slope->co[r][c]*Pi/180.)) + adt->T->Z->co[l][r][c];
			}else {
				adt->W->H1->co[i] = L->P->co[l][j] + adt->T->Z->co[l][r][c];
			}

		}else {
			
			l = adt->C->lch->co[i-n][1];
			ch = adt->C->lch->co[i-n][2];
			r = adt->C->r->co[ch];
			c = adt->C->c->co[ch];
			
			//solution guess
			if (adt->W->Pnet->co[r][c] > 0 && l == 0) {
				adt->W->H1->co[i] = Fmax(0., C->P->co[l][ch]) + (adt->W->Pnet->co[r][c]/cos(adt->T->slope->co[r][c]*Pi/180.)) + ( adt->T->Z->co[l][r][c] - adt->P->depr_channel );
			}else {
				adt->W->H1->co[i] = C->P->co[l][ch] + ( adt->T->Z->co[l][r][c] - adt->P->depr_channel );	
			}
			
		}
	}
	
	sux = find_matrix_K_3D(Dt, L, C, adt->W->Lx, adt->W->Klat, adt->W->Kbottom, adt->C->Kbottom, adt, adt->W->H1);
	
	find_f_3D(Dt, adt->W->f, adt, L, C, adt->W->H1, adt->W->Klat, adt->W->Kbottom, adt->C->Kbottom);
	
	product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
	
	res = norm_inf(adt->W->B, 1, N);
	
	res00 = res; //initial norm of the residual
	epsilon = adt->P->TolVWb + adt->P->RelTolVWb * Fmin( res00 , sqrt((double)N) );
	
	cont=0;
	out=0;
	
	//The condition to go out from the Newton cycle is made on the norm of the residual < Absolute tolerance + Relative tolerance * res00
	if( res <= Fmin( epsilon , max_res_adm ) ) out=1;
	
	//Max iteration number
	if( cont >= adt->P->MaxiterTol ) out=1;	
	
	while (out==0) {
		
		cont++;
		
		for(i=1; i<=N; i++){
			adt->W->H0->co[i] = adt->W->H1->co[i];
			adt->W->dH->co[i] = 0.;
		}
				
		//mu is the forcing term of Newton Raphson, defined such that ||F( x + d )|| < mu * ||F(x)|| see references
		if (cont==1) {
			mu = adt->P->TolCG;
		}else{
			mu *= Fmin( 1.0 , res/res0[0] );
			if(mu < 0.5*epsilon/res) mu = 0.5*epsilon/res;
		}
		
		//CALCOLATE AND STORE JACOBIAN AS SPARSE MATRIX
		
		sux = find_dfdH_3D(Dt, adt->W->df, adt, L, C, adt->W->H1, adt->W->Klat);	//it calcolates only df/dH, J = I*df/dH + K, K is calculated above
		
		//CONJUGATED GRADIENTS ALGORITHM
		iter = BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->W->df, adt->T->Li, adt->T->Lp, adt->W->Lx);
		if(iter==-1) return 1;	//does not converge 

		//non-monotonic line search (it is monotonic if M==1)
		for(m=Fminlong(cont,MM); m>1; m--){
			res_prev[m-1]=res_prev[m-2];
		}
		res_prev[0]=res;
		
		res_av=0.0;
		for(m=1; m<=Fminlong(cont,MM); m++){
			res_av=Fmax(res_prev[m-1],res_av);
		}				
		cont2 = 0.0;
		res0[0] = res;
		
		do{
			
			cont2++;
			
			/*The damping factor (or path length) lambda, defined as H1 = H0 + lambda*dH is chosen by minimizing the merit function :
			 0.5*norm(F_water(H0+lambda*dH)) interpolated with a quadratic polynome.
			 This method could not always be suited, because it has the disadvantage that it can make the solution stagnate to a point.
			 A relatively low number of MaxiterCorr (around 3-5) can prevent stagnation*/
			
			if(cont2 == 1){
				lambda[0] = 1.0;
				
			}else if(cont2 == 2){
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = thmax;
				
			}else{
				lambda[2] = lambda[1];
				res0[2] = res0[1];
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1], lambda[2], res0[2]);
				
			}
			
			for(i=1; i<=N; i++){
				adt->W->H1->co[i] = adt->W->H0->co[i] + lambda[0] * adt->W->dH->co[i];
				if(adt->W->H1->co[i] != adt->W->H1->co[i]) return 1;
			}
			
			if(updateK == 1 && cont <= maxITER_rec_K) sux = find_matrix_K_3D(Dt, L, C, adt->W->Lx, adt->W->Klat, adt->W->Kbottom, adt->C->Kbottom, adt, adt->W->H1);		
			
			find_f_3D(Dt, adt->W->f, adt, L, C, adt->W->H1, adt->W->Klat, adt->W->Kbottom, adt->C->Kbottom);
			
			product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
			
			res = norm_inf(adt->W->B, 1, N);		
			
			out2=0;
			
			//printf("cnt:%ld res:%e lambda:%e Dt:%f P:%f\n",cont,res,lambda[0],Dt,*Total_Pnet);

			if(res <= (1.0 - ni*lambda[0]*(1.-mu))*res_av) out2=1;
			if(lambda[0] <= adt->P->min_lambda_wat) cont_lambda_min++;
			
			if(cont_lambda_min > adt->P->max_times_min_lambda_wat){
				if(adt->P->exit_lambda_min_wat == 1){
					return 1;
				}else {
					out2=1;
					cont_lambda_min=0;
				}
			}
						
		}while(out2==0);	
		
		out=0;
		//The condition to go out from the Newton cycle is made on the norm of the residual < Absolute tolerance + Relative tolerance * res00
		if( res <= Fmin( epsilon , max_res_adm ) ) out=1;
		//Max iteration number
		if( cont >= adt->P->MaxiterTol ) out=1;	
				
	}
	
	if( res > epsilon ) return 1;
	
	//it can be shown that massloss per unit pixel [mm] is the linear norm of B * Dt / total_pixel
	*loss = norm_1(adt->W->B, 1, N)*Dt/adt->P->total_area;

	//assign updated state variables	
	for(i=1; i<=N; i++){
		
		if (i<=n) {//land
			
			l = adt->T->lrc_cont->co[i][1];
			r = adt->T->lrc_cont->co[i][2];
			c = adt->T->lrc_cont->co[i][3];
			j = adt->T->j_cont[r][c];
			sy = adt->S->type->co[r][c];			
			ch = adt->C->ch->co[r][c];
			bc = adt->T->BC_counter->co[r][c];
			
			L->P->co[l][j] = adt->W->H1->co[i] - adt->T->Z->co[l][r][c];
						
			//update variables
			if(l>0){
				adt->S->th->co[l][j] = theta_from_psi(L->P->co[l][j], L->thi->co[l][j], l, adt->S->pa->co[sy], PsiMin);
				adt->S->Ptot->co[l][j] = psi_from_theta(adt->S->th->co[l][j]+L->thi->co[l][j], 0., l, adt->S->pa->co[sy], PsiMin);									 
				adt->S->th->co[l][j] = Fmin( adt->S->th->co[l][j] , adt->S->pa->co[sy][jsat][l]-L->thi->co[l][j] );	
			}
			
			//volume lost at the bottom
			if(l==Nl){
				area = ds*ds/cos(adt->T->slope->co[r][c]*Pi/180.);
				if (ch>0) area -= adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
				*Vbottom = *Vbottom + area * adt->W->Kbottom->co[r][c] * 1.E-3 * Dt;
			}
			
			//lateral drainage at the border
			if (bc>0) {
				
				if (l>0) {
					
					if (adt->T->pixel_type->co[r][c] == 1 || adt->T->pixel_type->co[r][c] == 11) {
						//The depth of the free surface is multiplied by cosine since Z's are the layer depths in vertical direction
						if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.) && adt->W->H1->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
							
							dz = adt->S->pa->co[sy][jdz][l];//[mm]						
							
							if ((long)adt->L->LC->co[r+1][c]==number_novalue || (long)adt->L->LC->co[r-1][c]==number_novalue) {
								dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdN->co[r][c]));//mm
								dn = ds / cos(atan(adt->T->dzdE->co[r][c]));//m
							}else {
								dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdE->co[r][c]));
								dn = ds / cos(atan(adt->T->dzdN->co[r][c]));//m
							}
							
							*Vlatsub = *Vlatsub + Dt * (dn*dz*1.E-3) * 1.E-3*adt->W->Klat->co[bc][l]*adt->P->free_drainage_lateral*(adt->W->H1->co[i] - adt->T->Z->co[l][r][c]) / dD;
							
						}
						
					}else if (adt->T->pixel_type->co[r][c] == 2 || adt->T->pixel_type->co[r][c] == 12) {
						
						if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.) ) {
							
							dz = adt->S->pa->co[sy][jdz][l];//[mm]						
							
							if ((long)adt->L->LC->co[r+1][c]==number_novalue || (long)adt->L->LC->co[r-1][c]==number_novalue) {
								dn = ds / cos(atan(adt->T->dzdE->co[r][c]));//m
							}else {
								dn = ds / cos(atan(adt->T->dzdN->co[r][c]));//m
							}
							
							*Vlatsub = *Vlatsub + Dt * (dn*dz*1.E-3) * 1.E-3*adt->W->Klat->co[bc][l]*adt->P->free_drainage_lateral;
							
						}
					}
				}
			}
			
			
		}else {//channel
			
			l = adt->C->lch->co[i-n][1];
			ch = adt->C->lch->co[i-n][2];
			r = adt->C->r->co[ch];
			c = adt->C->c->co[ch];
			sy = adt->C->soil_type->co[ch];
			
			if (l==0){
				//hold and hnew are normal
				hold = Fmax(0., C->P->co[l][ch]) / cos(adt->T->slope->co[r][c] * Pi/180.);
			}
			
			//depr channel is defined vertical
			C->P->co[l][ch] = adt->W->H1->co[i] - ( adt->T->Z->co[l][r][c] - adt->P->depr_channel );
			
			if(l>0){
				adt->C->th->co[l][ch] = theta_from_psi(C->P->co[l][ch], C->thi->co[l][ch], l, adt->S->pa->co[sy], PsiMin);
				adt->C->th->co[l][ch] = Fmin( adt->C->th->co[l][ch] , adt->S->pa->co[sy][jsat][l]-C->thi->co[l][ch] );														 
			}
			
			if (l==0){
				//hold and hnew are normal
				hnew = Fmax(0., C->P->co[l][ch]) / cos(adt->T->slope->co[r][c] * Pi/180.);	
				Vsub->co[ch] += 1.E-3 * ( hnew - hold ) * adt->C->length->co[ch] * adt->P->w_dx * ds;
			}
			
			if(l==Nl){
				area = adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
				*Vbottom = *Vbottom + area * adt->C->Kbottom->co[ch] * 1.E-3 * Dt;			
			}
			
		}
	}
	
	return 0;			
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short Richards1D(long c, double Dt, SOIL_STATE *L, ALLDATA *adt, FILE *flog, double *loss, double *Vbottom, double *Vlat, double *Total_Pnet, short updateK){
	
	double res=0.0, res0[3], res_prev[MM], res_av, res00, lambda[3], epsilon, mu=0.;
	double ds=sqrt(UV->U->co[1]*UV->U->co[2]), area, dz, dn, dD;
	
	long i, l, r=1, m, bc, sy, cont, cont2, iter;
	long N=adt->W->H0->nh;
	long cont_lambda_min=0;
	short out, out2;	
	int sux;
	FILE *f;
	
	*Total_Pnet = 0.;
	
	res0[0]=0.;
	res0[1]=0.;
	res0[2]=0.;
	lambda[0]=0.;
	lambda[1]=0.;
	lambda[2]=0.;	
	
	for(i=1; i<=N; i++){//layers
		
		l = i-1;
		
		if (l == 0 && adt->W->Pnet->co[r][c] > 0) {
			*Total_Pnet = *Total_Pnet + (adt->W->Pnet->co[r][c]/cos(Fmin(max_slope,adt->T->slope->co[r][c])*Pi/180.)) / (double)adt->P->total_pixel;
		}
		
		//solution guess		
		if (adt->W->Pnet->co[r][c] > 0 && l == 0) {
			adt->W->H1->co[i] = Fmax(0., L->P->co[l][c]) + (adt->W->Pnet->co[r][c]/cos(Fmin(max_slope,adt->T->slope->co[r][c])*Pi/180.)) + adt->T->Z->co[l][r][c];
		}else {
			adt->W->H1->co[i] = L->P->co[l][c] + adt->T->Z->co[l][r][c];
		}
		
	}
	
	sux = find_matrix_K_1D(c, Dt, L, adt->W->Lx, adt->W->Klat, adt->W->Kbottom, adt, adt->W->H1);
		
	find_f_1D(c, Dt, L, adt->W->f, adt, adt->W->H1, adt->W->Klat, adt->W->Kbottom);
	
	product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
	
	res = norm_inf(adt->W->B, 1, N);
	
	res00 = res; //initial norm of the residual
	epsilon = adt->P->TolVWb + adt->P->RelTolVWb * Fmin( res00 , sqrt((double)N) );
	
	cont=0;
	
	out=0;
	//The condition to go out from the Newton cycle is made on the norm of the residual < Absolute tolerance + Relative tolerance * res00
	if( res <= Fmin( epsilon , max_res_adm ) ) out=1;
	//Max iteration number
	if( cont >= adt->P->MaxiterTol ) out=1;	
	
	while (out==0) {
		
		cont++;
		
		for(i=1; i<=N; i++){
			adt->W->H0->co[i] = adt->W->H1->co[i];
			adt->W->dH->co[i] = 0.;	
		}
				
		//mu is the forcing term of Newton Raphson, defined such that ||F( x + d )|| < mu * ||F(x)|| see references
		if (cont==1) {
			mu = adt->P->TolCG;
		}else{
			mu *= Fmin( 1.0 , res/res0[0] );
			if(mu < 0.5*epsilon/res) mu = 0.5*epsilon/res;
		}
		
		//CALCOLATE AND STORE JACOBIAN AS SPARSE MATRIX
		sux = find_dfdH_1D(c, Dt, L, adt->W->df, adt, adt->W->H1, adt->W->Klat);	//it calcolates only df/dH, J = I*df/dH + K, K is calculated above
		
		//CONJUGATED GRADIENTS ALGORITHM
		iter = BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->W->df, adt->T->Li, adt->T->Lp, adt->W->Lx);
		if(iter==-1){
			return 1;	//does not converge 
		}
		
		//non-monotonic line search (it is monotonic if M==1)
		for(m=Fminlong(cont,MM); m>1; m--){
			res_prev[m-1]=res_prev[m-2];
		}
		res_prev[0]=res;
		
		res_av=0.0;
		for(m=1; m<=Fminlong(cont,MM); m++){
			res_av=Fmax(res_prev[m-1],res_av);
		}				
		cont2 = 0.0;
		res0[0] = res;
		
		do{
			
			cont2++;
			
			/*The damping factor (or path length) lambda, defined as H1 = H0 + lambda*dH is chosen by minimizing the merit function :
			 0.5*norm(F_water(H0+lambda*dH)) interpolated with a quadratic polynome.
			 This method could not always be suited, because it has the disadvantage that it can make the solution stagnate to a point.
			 A relatively low number of MaxiterCorr (around 3-5) can prevent stagnation*/
			
			if(cont2 == 1){
				lambda[0] = 1.0;
				
			}else if(cont2 == 2){
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = thmax;
				
			}else{
				lambda[2] = lambda[1];
				res0[2] = res0[1];
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1], lambda[2], res0[2]);
				
			}
			
			for(i=1; i<=N; i++){
				
				adt->W->H1->co[i] = adt->W->H0->co[i] + lambda[0] * adt->W->dH->co[i];
				
				if(adt->W->H1->co[i] != adt->W->H1->co[i]) {
					f = fopen(FailedRunFile, "w");
					fprintf(f, "Simulation Period:%ld\n",i_sim);
					fprintf(f, "Run Time:%ld\n",i_run);
					fprintf(f, "Number of days after start:%f\n",adt->I->time/86400.);					
					fprintf(f, "Error: no value psi Richards3D l:%ld point:%ld\n",i-1,c);
					fclose(f);
					t_error("Fatal Error! Geotop is closed. See failing report.");	
				}
				
			}
			
			if(updateK == 1 && cont <= maxITER_rec_K) sux = find_matrix_K_1D(c, Dt, L, adt->W->Lx, adt->W->Klat, adt->W->Kbottom, adt, adt->W->H1);		
			
			find_f_1D(c, Dt, L, adt->W->f, adt, adt->W->H1, adt->W->Klat, adt->W->Kbottom);
			
			product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
			
			res = norm_inf(adt->W->B, 1, N);	
			
			out2=0;
			
			if(res <= (1.0 - ni*lambda[0]*(1.-mu))*res_av) out2=1;
			if(lambda[0] <= adt->P->min_lambda_wat) cont_lambda_min++;
			
			if(cont_lambda_min > adt->P->max_times_min_lambda_wat){
				if(adt->P->exit_lambda_min_wat == 1){
					return 1;
				}else {
					out2=1;
					cont_lambda_min=0;
				}
			}
			
			
		}while(out2==0);	
		
		out=0;
		//The condition to go out from the Newton cycle is made on the norm of the residual < Absolute tolerance + Relative tolerance * res00
		if( res <= Fmin( epsilon , max_res_adm ) ) out=1;
		//Max iteration number
		if( cont >= adt->P->MaxiterTol ) out=1;	
		
		//printf("Dt:%f %ld %ld res:%e(%e) lambda:%e cont_lambda_min:%ld i:%ld\n",Dt,cont,cont2,res,epsilon,lambda[0],cont_lambda_min,c);
		
	}
		
	if( res > epsilon ){
		return 1;
	}
	
	//it can be shown that massloss per unit pixel [mm] is the linear norm of B * Dt / total_pixel
	*loss = norm_1(adt->W->B, 1, N)*Dt/adt->P->total_area;
	
	//assign updated state variables	
	for(i=1; i<=N; i++){
		
		l = i-1;
		sy = adt->S->type->co[r][c];
		bc = adt->T->BC_counter->co[r][c];
		
		L->P->co[l][c] = adt->W->H1->co[i] - adt->T->Z->co[l][r][c];
		
		//update variables
		if(l>0){
			adt->S->th->co[l][c] = theta_from_psi(L->P->co[l][c], L->thi->co[l][c], l, adt->S->pa->co[sy], PsiMin);
			adt->S->Ptot->co[l][c] = psi_from_theta(adt->S->th->co[l][c]+L->thi->co[l][c], 0., l, adt->S->pa->co[sy], PsiMin);									 
			adt->S->th->co[l][c] = Fmin( adt->S->th->co[l][c] , adt->S->pa->co[sy][jsat][l]-L->thi->co[l][c] );														 
		}
		
		//volume lost at the bottom
		if(l==Fminlong(adt->P->Nl_spinup->co[i_sim],Nl)){
			area = ds*ds;
			*Vbottom = *Vbottom + area * adt->W->Kbottom->co[r][c] * 1.E-3 * Dt;
		}
		
		//lateral drainage at the border
		if (bc>0 && l>0) {
			if (adt->T->pixel_type->co[r][c] == 1) {
				if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.) && adt->W->H1->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
					dz = adt->S->pa->co[sy][jdz][l];//[mm]
					dn = ds;
					dD = 0.5 * 1.E3*ds;
					*Vlat = *Vlat + Dt * (dn*dz*1.E-3) * 1.E-3*adt->W->Klat->co[bc][l]*adt->P->free_drainage_lateral*(adt->W->H1->co[i] - adt->T->Z->co[l][r][c]) / dD;
				}
			}
		}
	}
	
	return 0;			
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double cm_h(double cm0, double h, double h_thres1, double h_thres2){
	
	double cm;
	
	if(h > h_thres2){
		cm = cm0;
	}else if(h > h_thres1){
		cm = ( 0.0 * (h_thres2 - h) + cm0*(h) ) / h_thres2;
	}else{
		cm = 0.0;
	}
	
	return(cm);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//cnt is the counter of Li Lp Lx (lower diagonal without diagonal)

int find_matrix_K_3D(double Dt, SOIL_STATE *SL, SOIL_STATE *SC, DOUBLEVECTOR *Lx, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom_l, DOUBLEVECTOR *Kbottom_ch, ALLDATA *adt, DOUBLEVECTOR *H){
	
	long i, l, r, c, j, I, R, C, J, sy, syn, ch, cnt=0;
	long n=(Nl+1)*adt->P->total_pixel;
	double dz=0.0, dzn=0.0, dD=0.0, k=0.0, kn=0.0, kmax=0.0, kmaxn=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]), dn;
	//double psi, ice, a, ns, res, sat, ss, Temp;
	
	for(i=1;i<=H->nh;i++){
		
		//VERTICAL FLUXES
		if( i<=n){//land
			
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			j=adt->T->j_cont[r][c];
			sy=adt->S->type->co[r][c];
			
			ch=adt->C->ch->co[r][c];
			area=ds*ds/cos(adt->T->slope->co[r][c]*Pi/180.);
			if (ch>0) area-=adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			
			//vertical hydraulic conductivity
			if(l>0){
				dz = adt->S->pa->co[sy][jdz][l];
				if (l==Nl && adt->P->free_drainage_bottom>0) Kbottom_l->co[r][c] = k_from_psi(jKn, H->co[i] - adt->T->Z->co[l][r][c], SL->thi->co[l][c], SL->T->co[l][c], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat	);		
			}
			
			//flux from cell below
			if (l<Nl) {
				
				I = i+1;
				
				if(l==0){	//overland flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dzn;
					
					if( H->co[i] < H->co[I] ){	
						//upward flux
						kn = k_from_psi(jKn, H->co[I] - adt->T->Z->co[l+1][r][c], SL->thi->co[l+1][j], SL->T->co[l+1][j], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}else{	
						//downward flow
						kn = k_from_psi(jKn, psisat_from(SL->thi->co[l+1][j], l+1, adt->S->pa->co[sy]), SL->thi->co[l+1][j], SL->T->co[l+1][j], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}
					
				}else{	//subsurface flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dz + 0.5*dzn;
					
					if( H->co[i] < H->co[I] ){	
						//upward flux
						kn = k_from_psi(jKn, H->co[I] - adt->T->Z->co[l+1][r][c], SL->thi->co[l+1][j], SL->T->co[l+1][j], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}else{	
						//downward flow
						kn = k_from_psi(jKn, H->co[i] - adt->T->Z->co[l][r][c], SL->thi->co[l][j], SL->T->co[l][j], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}
					
					/*psi = Arithmetic_Mean(dz, dzn, H->co[i] - adt->T->Z->co[l][r][c], H->co[I] - adt->T->Z->co[l+1][r][c]);
					ice = Arithmetic_Mean(dz, dzn, SL->thi->co[l][j], SL->thi->co[l+1][j]);
					Temp = Arithmetic_Mean(dz, dzn, SL->T->co[l][j], SL->T->co[l+1][j]);
					a = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][ja][l+1]);
					ns = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jns][l+1]);
					res = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][jres][l+1]);
					sat = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jsat][l+1]);
					ss = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jss][l], adt->S->pa->co[sy][jss][l+1]);
					k = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jKn][l], adt->S->pa->co[sy][jKn][l+1]);
					kn = k_hydr_soil(psi, k, adt->P->imp, ice, sat, res, a, ns, 1.-1./ns, 0.5, Temp, adt->P->k_to_ksat);*/
					
					kmax = k_from_psi( jKn,  psisat_from( SL->thi->co[l][j], l, adt->S->pa->co[sy]), SL->thi->co[l][j], SL->T->co[l][j], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					kmaxn = k_from_psi( jKn,  psisat_from( SL->thi->co[l+1][j], l+1, adt->S->pa->co[sy]), SL->thi->co[l+1][j], SL->T->co[l+1][j], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					kmaxn = Fmin(kmax, kmaxn);
					
					kn = Fmin(kn, kmaxn);
										
				}
				
				cnt++;
				Lx->co[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
				
			}
			
		}else{//channel
			
			l=adt->C->lch->co[i-n][1];
			ch=adt->C->lch->co[i-n][2];
			r=adt->C->r->co[ch];
			c=adt->C->c->co[ch];
			sy=adt->C->soil_type->co[ch];
			
			area=adt->C->length->co[ch] * adt->P->w_dx * ds;
			
			//vertical hydraulic conductivity
			if(l>0){
				dz = adt->S->pa->co[sy][jdz][l];
				if (l==Nl && adt->P->free_drainage_bottom>0) Kbottom_ch->co[ch] = k_from_psi(jKn, H->co[i] - (adt->T->Z->co[l][r][c]-adt->P->depr_channel), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
			}
			
			//flux from cell below
			if (l<Nl) {
				
				I = i+1;
				
				if(l==0){	//overland flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dzn;
					
					if( H->co[i] < H->co[I] ){	
						//upward flux
						kn = k_from_psi(jKn, H->co[I] - (adt->T->Z->co[l+1][r][c]-adt->P->depr_channel), SC->thi->co[l+1][ch], SC->T->co[l+1][ch], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}else{	
						//downward flow
						kn = k_from_psi(jKn, psisat_from(SC->thi->co[l+1][ch], l+1, adt->S->pa->co[sy]), SC->thi->co[l+1][ch], SC->T->co[l+1][ch], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}
					
				}else{	//subsurface flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dz + 0.5*dzn;
					
					if( H->co[i] < H->co[I] ){	
						//upward flux
						kn = k_from_psi(jKn, H->co[I] - (adt->T->Z->co[l+1][r][c]-adt->P->depr_channel), SC->thi->co[l+1][ch], SC->T->co[l+1][ch], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}else{	
						//downward flow
						kn = k_from_psi(jKn, H->co[i] - (adt->T->Z->co[l][r][c]-adt->P->depr_channel), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					}			
					
					/*psi = Arithmetic_Mean(dz, dzn, H->co[i] - (adt->T->Z->co[l][r][c]-adt->P->depr_channel), H->co[I] - (adt->T->Z->co[l+1][r][c]-adt->P->depr_channel));
					ice = Arithmetic_Mean(dz, dzn, SC->thi->co[l][ch], SC->thi->co[l+1][ch]);
					Temp = Arithmetic_Mean(dz, dzn, SC->T->co[l][ch], SC->T->co[l+1][ch]);
					a = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][ja][l+1]);
					ns = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jns][l+1]);
					res = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][jres][l+1]);
					sat = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jsat][l+1]);
					ss = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jss][l], adt->S->pa->co[sy][jss][l+1]);
					k = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jKn][l], adt->S->pa->co[sy][jKn][l+1]);
					kn = k_hydr_soil(psi, k, adt->P->imp, ice, sat, res, a, ns, 1.-1./ns, 0.5, Temp, adt->P->k_to_ksat);*/
					
					kmax = k_from_psi(jKn, psisat_from(SC->thi->co[l][ch], l, adt->S->pa->co[sy]), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					kmaxn = k_from_psi(jKn, psisat_from(SC->thi->co[l+1][ch], l+1, adt->S->pa->co[sy]), SC->thi->co[l+1][ch], SC->T->co[l+1][ch], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					kmaxn = Fmin(kmax, kmaxn);
					
					kn = Fmin(kn, kmaxn);
					
				}
				
				cnt++;
				Lx->co[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
			}
		}
		
		//LATERAL FLUXES
		if (i<=n){
			
			//lateral hydraulic conductivity
			if(l>0){
				k = k_from_psi(jKl, H->co[i] - adt->T->Z->co[l][r][c], SL->thi->co[l][j], SL->T->co[l][j], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
				kmax = k_from_psi(jKl, psisat_from(SL->thi->co[l][j], l, adt->S->pa->co[sy]), SL->thi->co[l][l], SL->T->co[l][j], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
				if(adt->T->BC_counter->co[r][c]>0){
					if(adt->T->pixel_type->co[r][c] == 1 || adt->T->pixel_type->co[r][c] == 2 || adt->T->pixel_type->co[r][c] == 11 || adt->T->pixel_type->co[r][c] == 12) Klat->co[adt->T->BC_counter->co[r][c]][l] = k;
				}
			}
			
			//4 neighbouring cells
			R = r-1;
			C = c;
			//1.
			if(R>=1 && R<=Nr && C>=1 && C<=Nc){
				if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type->co[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));//[m]
					
					if(l>0){
						
						//Subsurface Flow
						if (H->co[I] > H->co[i]) {
							kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi->co[l][J], l, adt->S->pa->co[syn]), SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
												
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						//Surface Flow
						if (H->co[I] > H->co[i]) {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[I] - adt->T->Z->co[l][R][C]) / cos(adt->T->slope->co[R][C]*Pi/180.);
						}else {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[i] - adt->T->Z->co[l][r][c]) / cos(adt->T->slope->co[r][c]*Pi/180.);
						}
						
						if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			R = r+1;
			C = c;
			dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));
			//2.
			if(R>=1 && R<=Nr && C>=1 && C<=Nc){
				if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type->co[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));//[m]
					
					if(l>0){
						if (H->co[I] > H->co[i]) {
							kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi->co[l][J], l, adt->S->pa->co[syn]), SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						if (H->co[I] > H->co[i]) {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[I] - adt->T->Z->co[l][R][C]) / cos(adt->T->slope->co[R][C]*Pi/180.);
						}else {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[i] - adt->T->Z->co[l][r][c]) / cos(adt->T->slope->co[r][c]*Pi/180.);
						}
						
						if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
					
				}
			}
			
			R = r;
			C = c-1;
			dn = ds/cos(0.5*atan(adt->T->dzdN->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));
			//3.
			if(R>=1 && R<=Nr && C>=1 && C<=Nc){
				if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type->co[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));//[m]
					
					if(l>0){
						if (H->co[I] > H->co[i]) {
							kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi->co[l][J], l, adt->S->pa->co[syn]), SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						if (H->co[I] > H->co[i]) {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[I] - adt->T->Z->co[l][R][C]) / cos(adt->T->slope->co[R][C]*Pi/180.);
						}else {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[i] - adt->T->Z->co[l][r][c]) / cos(adt->T->slope->co[r][c]*Pi/180.);
						}
						
						if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			R = r;
			C = c+1;
			dn = ds/cos(0.5*atan(adt->T->dzdN->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));
			//4.
			if(R>=1 && R<=Nr && C>=1 && C<=Nc){
				if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type->co[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));//[m]
					
					if(l>0){
						if (H->co[I] > H->co[i]) {
							kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi->co[l][J], l, adt->S->pa->co[syn]), SL->thi->co[l][J], SL->T->co[l][J], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						if (H->co[I] > H->co[i]) {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[I] - adt->T->Z->co[l][R][C]) / cos(adt->T->slope->co[R][C]*Pi/180.);
						}else {
							kn = adt->L->ty->co[(long)adt->L->LC->co[R][C]][jcm];
							dz = Fmax(0., H->co[i] - adt->T->Z->co[l][r][c]) / cos(adt->T->slope->co[r][c]*Pi/180.);
						}
						
						if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			//exchange with channels
			if (l>0 && adt->C->ch->co[r][c] > 0) {
				
				ch = adt->C->ch->co[r][c];
				syn = adt->C->soil_type->co[ch];
				I = n + adt->C->ch3[l][ch];
				
				if (H->co[I] > H->co[i]) {
					kn = k_from_psi(jKl, H->co[I] - (adt->T->Z->co[l][r][c]-adt->P->depr_channel), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
				}else {
					kn = k;
				}
				
				kmaxn = k_from_psi(jKl, psisat_from(SC->thi->co[l][ch], l, adt->S->pa->co[syn]), SC->thi->co[l][ch], SC->T->co[l][ch], l, adt->S->pa->co[syn], adt->P->imp, adt->P->k_to_ksat);	
				kmaxn = Fmin(kmax, kmaxn);
				kn = Fmin(kn, kmaxn);
				
				//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]	
				//Area[m2] = 2 * channel length * layer thickness = 2 * length[m] * dz[mm] * 1.E-3[m/mm]
				//dD[mm] = 0.25 * ds * (1+w_dx)
				dD = find_3Ddistance(ds * (1.0 + adt->P->w_dx) / 4.0, 1.E-3*adt->P->depr_channel) * 1.E3;//[mm]
				
				cnt++;
				Lx->co[cnt] = -(2.*adt->C->length->co[adt->C->ch->co[r][c]]*1.E-3*dz)*kn/dD;						
				
			}
		}	
	}										
	
	
	return 0;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_matrix_K_1D(long c, double Dt, SOIL_STATE *L, DOUBLEVECTOR *Lx, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom, ALLDATA *adt, DOUBLEVECTOR *H){
	
	long i, l, r=1, I, sy, cnt=0;
	double dz=0.0, dzn=0.0, dD=0.0, kn=0.0, kmax=0.0, kmaxn=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	//double psi, ice, a, ns, res, sat, ss, Temp, k;
	
	for(i=1;i<=H->nh;i++){
		
		//VERTICAL FLUXES
		l = i-1;	
		sy=adt->S->type->co[r][c];
		area=ds*ds;
		
		//vertical hydraulic conductivity
		if(l>0){
			dz = adt->S->pa->co[sy][jdz][l];
			if (l==Fminlong(adt->P->Nl_spinup->co[i_sim],Nl) && adt->P->free_drainage_bottom>0) Kbottom->co[r][c] = k_from_psi(jKn, H->co[i] - adt->T->Z->co[l][r][c], L->thi->co[l][c], L->T->co[l][c], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);		
		}
		
		//flux from cell below
		if (l<Fminlong(adt->P->Nl_spinup->co[i_sim],Nl)) {
			
			I = i+1;
			
			if(l==0){	//overland flow
				
				dzn = adt->S->pa->co[sy][jdz][l+1];
				dD = 0.5*dzn;
				
				if( H->co[i] < H->co[I] ){	
					//upward flux
					kn = k_from_psi(jKn, H->co[I] - adt->T->Z->co[l+1][r][c], L->thi->co[l+1][c], L->T->co[l+1][c], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					
				}else{	
					//downward flow
					kn = k_from_psi(jKn, psisat_from(L->thi->co[l+1][c], l+1, adt->S->pa->co[sy]), L->thi->co[l+1][c], L->T->co[l+1][c], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
					
				}				
				
			}else{	//subsurface flow
				
				dzn = adt->S->pa->co[sy][jdz][l+1];
				dD = 0.5*dz + 0.5*dzn;
				
				if( H->co[i] < H->co[I] ){	
					//upward flux
					kn = k_from_psi(jKn, H->co[I] - adt->T->Z->co[l+1][r][c], L->thi->co[l+1][c], L->T->co[l+1][c], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
				}else{	
					//downward flow
					kn = k_from_psi(jKn, H->co[i] - adt->T->Z->co[l][r][c], L->thi->co[l][c], L->T->co[l][c], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
				}
				
				/*psi = Arithmetic_Mean(dz, dzn, H->co[i] - adt->T->Z->co[l][r][c], H->co[I] - adt->T->Z->co[l+1][r][c]);
				ice = Arithmetic_Mean(dz, dzn, L->thi->co[l][c], L->thi->co[l+1][c]);
				Temp = Arithmetic_Mean(dz, dzn, L->T->co[l][c], L->T->co[l+1][c]);
				a = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][ja][l+1]);
				ns = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jns][l+1]);
				res = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][jres][l+1]);
				sat = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jsat][l+1]);
				ss = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jss][l], adt->S->pa->co[sy][jss][l+1]);
				k = Arithmetic_Mean(dz, dzn, adt->S->pa->co[sy][jKn][l], adt->S->pa->co[sy][jKn][l+1]);
				kn = k_hydr_soil(psi, k, adt->P->imp, ice, sat, res, a, ns, 1.-1./ns, 0.5, Temp, adt->P->k_to_ksat);*/
				
				kmax = k_from_psi( jKn,  psisat_from( L->thi->co[l][c], l, adt->S->pa->co[sy]), L->thi->co[l][c], L->T->co[l][c], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
				kmaxn = k_from_psi( jKn,  psisat_from( L->thi->co[l+1][c], l+1, adt->S->pa->co[sy]), L->thi->co[l+1][c], L->T->co[l+1][c], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
				kmaxn = Fmin(kmax, kmaxn);
				
				kn = Fmin(kn, kmaxn);
				
			}
			
			cnt++;
			Lx->co[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
		}
		
		//LATERAL FLUXES
		if (adt->T->pixel_type->co[r][c] == 1){
			
			//lateral hydraulic conductivity
			if (l>0) Klat->co[adt->T->BC_counter->co[r][c]][l] = k_from_psi(jKl, H->co[i] - adt->T->Z->co[l][r][c], L->thi->co[l][c], L->T->co[l][c], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);	
			
		}
		
	}
	
	return 0;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_dfdH_3D(double Dt, DOUBLEVECTOR *df, ALLDATA *adt, SOIL_STATE *L, SOIL_STATE *C, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat){
	
	long i, l, r, c, j, sy, ch, bc;
	long n=(Nl+1)*adt->P->total_pixel;
	double dz, dn, dD, psi1, ice=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	
	for(i=1;i<=H->nh;i++){
		
		df->co[i] = 0.;
		
		if (i<=n) {
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			j=adt->T->j_cont[r][c];
			sy=adt->S->type->co[r][c];
			bc=adt->T->BC_counter->co[r][c];
			ch=adt->C->ch->co[r][c];
			area=ds*ds/cos(adt->T->slope->co[r][c]*Pi/180.);
			if (ch>0) area-=adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			psi1 = H->co[i] - adt->T->Z->co[l][r][c];
			if(l>0) ice = L->thi->co[l][j];
			
		}else {
			l=adt->C->lch->co[i-n][1];
			ch=adt->C->lch->co[i-n][2];
			r=adt->C->r->co[ch];
			c=adt->C->c->co[ch];
			sy=adt->C->soil_type->co[ch];
			bc=0;
			area=adt->C->length->co[ch] * adt->P->w_dx * ds;
			psi1 = H->co[i] - (adt->T->Z->co[l][r][c] - adt->P->depr_channel);
			if(l>0) ice = C->thi->co[l][ch];
		}
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			if(psi1>0) df->co[i] += ( area / cos(adt->T->slope->co[r][c]*Pi/180.) ) / Dt;
			
		}else{
			dz = adt->S->pa->co[sy][jdz][l];
			df->co[i] += dtheta_dpsi_from_psi(psi1, ice, l, adt->S->pa->co[sy], PsiMin) * area * dz / Dt;
			
		}
		
		//lateral drainage at the border
		if (bc>0){
			
			if (l>0) {
				if (adt->T->pixel_type->co[r][c] == 1 || adt->T->pixel_type->co[r][c] == 11) {
					if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.) && H->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
						if ((long)adt->L->LC->co[r+1][c]==number_novalue || (long)adt->L->LC->co[r-1][c]==number_novalue) {
							dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdN->co[r][c]));//mm
							dn = ds / cos(atan(adt->T->dzdE->co[r][c]));//m
						}else {
							dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdE->co[r][c]));
							dn = ds / cos(atan(adt->T->dzdN->co[r][c]));//m
						}
						
						dz = adt->S->pa->co[sy][jdz][l];//[mm]
						df->co[i] += (dn*dz*1.E-3) * Klat->co[bc][l]*adt->P->free_drainage_lateral / dD;
						
					}
				}
				
			}
		}
		
	}
	
	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_dfdH_1D(long c, double Dt, SOIL_STATE *L, DOUBLEVECTOR *df, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat){
	
	long i, r=1, l, sy, bc;
	double dz, dn, dD, psi1, ice=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	
	for(i=1;i<=H->nh;i++){
		
		df->co[i] = 0.;
		
		l = i-1;	
		sy=adt->S->type->co[r][c];
		bc=adt->T->BC_counter->co[r][c];
		area=ds*ds;
		psi1 = H->co[i] - adt->T->Z->co[l][r][c];
		if(l>0) ice = L->thi->co[l][c];
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			if(psi1>0) df->co[i] += ( area / cos( Fmin(max_slope,adt->T->slope->co[r][c])*Pi/180.) ) / Dt;
		}else{
			dz = adt->S->pa->co[sy][jdz][l];
			df->co[i] += dteta_dpsi(psi1, ice, adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], 
									adt->S->pa->co[sy][jns][l], 1.-1./adt->S->pa->co[sy][jns][l], PsiMin, adt->S->pa->co[sy][jss][l])*
			area * dz / Dt;
		}
		
		//lateral drainage at the border
		if (bc>0 && l>0) {
			if (adt->T->pixel_type->co[r][c] == 1) {
				if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.) && H->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
					dz = adt->S->pa->co[sy][jdz][l];//[mm]
					dn = ds;
					dD = 0.5 * 1.E3*ds;
					df->co[i] += (dn*dz*1.E-3) * Klat->co[bc][l]*adt->P->free_drainage_lateral / dD;
				}
			}
		}
	}
	
	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_f_3D(double Dt, DOUBLEVECTOR *f, ALLDATA *adt, SOIL_STATE *L, SOIL_STATE *C, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom_l, DOUBLEVECTOR *Kbottom_ch){
	
	long i, l, r, c, j, sy, ch, bc;
	long n=(Nl+1)*adt->P->total_pixel;
	double dz, dn, dD, V0, V1, psi1, psi0, ice=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	
	for(i=1;i<=H->nh;i++){
		
		if (i<=n) {
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			j=adt->T->j_cont[r][c];
			sy=adt->S->type->co[r][c];
			bc=adt->T->BC_counter->co[r][c];
			ch=adt->C->ch->co[r][c];
			area=ds*ds/cos(adt->T->slope->co[r][c]*Pi/180.);
			if (ch>0) area-=adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			psi0 = L->P->co[l][j];
			psi1 = H->co[i] - adt->T->Z->co[l][r][c];
			if(l>0) ice = L->thi->co[l][j];
			
		}else {
			l=adt->C->lch->co[i-n][1];
			ch=adt->C->lch->co[i-n][2];
			r=adt->C->r->co[ch];
			c=adt->C->c->co[ch];
			sy=adt->C->soil_type->co[ch];
			bc=0;
			area=adt->C->length->co[ch] * adt->P->w_dx * ds;
			psi0 = C->P->co[l][ch];
			psi1 = H->co[i] - (adt->T->Z->co[l][r][c] - adt->P->depr_channel);
			if(l>0) ice = C->thi->co[l][ch];
			
		}
		
		//hydraulic capacity (diagonal term)
		if(l==0){
			V1 = area * Fmax(0.0, psi1) / cos(adt->T->slope->co[r][c]*Pi/180.);
			V0 = area * Fmax(0.0, psi0) / cos(adt->T->slope->co[r][c]*Pi/180.);
		}else{
			dz = adt->S->pa->co[sy][jdz][l];		
			V1 = area*dz * theta_from_psi(psi1, ice, l, adt->S->pa->co[sy], PsiMin);
			V0 = area*dz * theta_from_psi(psi0, ice, l, adt->S->pa->co[sy], PsiMin);
		}
		
		f->co[i] = (V1-V0)/Dt;
		
		//drainage at the bottom
		if (l==Nl){
			if (i<=n) {
				f->co[i] += area*Kbottom_l->co[r][c];
			}else {
				f->co[i] += area*Kbottom_ch->co[ch];
			}
		}
		
		//lateral drainage at the border
		if (bc>0) {
			
			if (l>0) {
				if (adt->T->pixel_type->co[r][c] == 1 || adt->T->pixel_type->co[r][c] == 11) {
					if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.) && H->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
						if ((long)adt->L->LC->co[r+1][c]==number_novalue || (long)adt->L->LC->co[r-1][c]==number_novalue) {
							dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdN->co[r][c]));//mm
							dn = ds / cos(atan(adt->T->dzdE->co[r][c]));//m
						}else {
							dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdE->co[r][c]));
							dn = ds / cos(atan(adt->T->dzdN->co[r][c]));//m
						}
						
						dz = adt->S->pa->co[sy][jdz][l];//[mm]
						f->co[i] += (dn*dz*1.E-3) * Klat->co[bc][l]*adt->P->free_drainage_lateral*(H->co[i] - adt->T->Z->co[l][r][c]) / dD;
					}
					
				}else if (adt->T->pixel_type->co[r][c] == 2 || adt->T->pixel_type->co[r][c] == 12) {
					if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.) ) {
						if ((long)adt->L->LC->co[r+1][c]==number_novalue || (long)adt->L->LC->co[r-1][c]==number_novalue) {
							dn = ds / cos(atan(adt->T->dzdE->co[r][c]));//m
						}else {
							dn = ds / cos(atan(adt->T->dzdN->co[r][c]));//m
						}
						
						dz = adt->S->pa->co[sy][jdz][l];//[mm]
						f->co[i] += (dn*dz*1.E-3) * Klat->co[bc][l]*adt->P->free_drainage_lateral;
					}
				}
				
			}
		}
		
		//evaporation and precipitation
		if(l>0){
			if(i<=n){
				f->co[i] += area*adt->S->ET->co[l][r][c]/Dt;
			}else {
				ch=adt->C->ch->co[r][c];
				f->co[i] += area*adt->C->ET->co[l][ch]/Dt;
			}
		}else {
			f->co[i] -= area*adt->W->Pnet->co[r][c]/Dt;
		}
		
		
	}
	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_f_1D(long c, double Dt, SOIL_STATE *L, DOUBLEVECTOR *f, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom){
	
	long i, l, r=1, sy, bc;
	double dz, dn, dD, V0, V1, psi1, psi0, ice=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	
	for(i=1;i<=H->nh;i++){
		
		l = i-1;	
		sy=adt->S->type->co[r][c];
		bc=adt->T->BC_counter->co[r][c];
		area=ds*ds;
		psi0 = L->P->co[l][c];
		psi1 = H->co[i] - adt->T->Z->co[l][r][c];
		if(l>0) ice = L->thi->co[l][c];
		
		//hydraulic capacity (diagonal term)
		if(l==0){
			V1 = area * Fmax(0.0, psi1) / cos( Fmin(max_slope,adt->T->slope->co[r][c])*Pi/180.);
			V0 = area * Fmax(0.0, psi0) / cos( Fmin(max_slope,adt->T->slope->co[r][c])*Pi/180.);
			
		}else{
			
			dz = adt->S->pa->co[sy][jdz][l];					
			V1 = area*dz * theta_from_psi(psi1, ice, l, adt->S->pa->co[sy], PsiMin);
			V0 = area*dz * theta_from_psi(psi0, ice, l, adt->S->pa->co[sy], PsiMin);
			
		}
		
		f->co[i] = (V1-V0)/Dt;
				
		//drainage at the bottom
		if (l==Fminlong(adt->P->Nl_spinup->co[i_sim],Nl)){
			f->co[i] += area*Kbottom->co[r][c];
		}		
		
		//lateral drainage at the border
		if (bc>0 && l>0) {
			if (adt->T->pixel_type->co[r][c] == 1) {
				if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.) && H->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
					dz = adt->S->pa->co[sy][jdz][l];//[mm]
					dn = ds;
					dD = 0.5 * 1.E3*ds;
					f->co[i] += (dn*dz*1.E-3) * Klat->co[bc][l]*adt->P->free_drainage_lateral*(H->co[i] - adt->T->Z->co[l][r][c]) / dD;
				}
			}
		}
		
		
		//evaporation 
		if(l>0){
			f->co[i] += area*adt->S->ET->co[l][r][c]/Dt;
		}else {
			f->co[i] -= area*adt->W->Pnet->co[r][c]/Dt;
		}
		
		
	}
	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double find_3Ddistance(double horizontal_distance, double vertical_distance){
	
	return sqrt(pow(horizontal_distance,2.)+pow(vertical_distance,2.));
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_dt_max(short DD, double Courant, double *h, LAND *land, TOPO *top, CHANNEL *cnet, PAR *par, METEO *met, double t, double *dt){
	
	double q, ds=sqrt(UV->U->co[1]*UV->U->co[2]), area, Vmax, H;
	short d;
	long r, c, j, ch;
	
	for (j=1; j<=par->total_pixel; j++) {
		
		r = top->rc_cont->co[j][1];
		c = top->rc_cont->co[j][2];
		
		H = Fmax(0.0 , h[j]) / cos(top->slope->co[r][c]*Pi/180.);//h[i] is the pressure at the surface, H is the depth of water normal to the surface
		
		if(H > par->min_hsup_land){
			
			if(DD==1){
				draining_land(1., j, top, land, par, cnet, h, top->Jdown->co[j], top->Qdown->co[j]);
			}else {
				draining_land(0., j, top, land, par, cnet, h, top->Jdown->co[j], top->Qdown->co[j]);
			}
			
			area = ds*ds;
			area /= cos(top->slope->co[r][c]*Pi/180.);
			ch = cnet->ch->co[r][c];
			if (ch>0) area -= cnet->length->co[ch] * par->w_dx * ds; //area of the pixel[m2]
			
			Vmax = area * 1.E-3 * H;
			
			q = 0.;
			for (d=1; d<=4; d++) {
				q += top->Qdown->co[j][d];	//outgoing discharge
			}
			
			Vmax = Fmax(Vmax, 1.E-10);
			q = Fmax(q, 1.E-30);
			
			if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
			
		}
	}
	
	if(*dt < par->dtmin_sup) *dt = par->dtmin_sup;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void supflow(short DDland, short DDch, double Dt, double t, double *h, double *dV, double *hch, double *dhch, TOPO *top, LAND *land, 
			 WATER *wat, CHANNEL *cnet, PAR *par, METEO *met, DOUBLEVECTOR *Vsup, double *Voutnet, double *Voutland, FILE *flog,
			 double *mm1, double *mm2, double *mmo) {
	
	
	long d, r, c, j, R, C, ch;                                    
	double ds=sqrt(UV->U->co[1]*UV->U->co[2]), area, Vmax, H;							     
	double q, q0, tb, te=0.0, dt;
	long cnt=0,cnt2=0,cnt3=0;
	
	double m1=0.,m2=0.,mo=0.;
	
	for (j=1; j<=par->total_pixel; j++) {
		r = top->rc_cont->co[j][1];
		c = top->rc_cont->co[j][2];
		H = Fmax(0., h[j]) / cos(top->slope->co[r][c]*Pi/180.);
		area = ds*ds;
		area /= cos(top->slope->co[r][c]*Pi/180.);
		m1 += H*1.E-3*area;
	}
	
	do{
		
		tb=te;
		dt=Dt;
		
		find_dt_max(DDland, par->max_courant_land, h, land, top, cnet, par, met, t, &dt);		
		cnt++;
		
		te=tb+dt;
		if(te>Dt){
			te=Dt;
			dt=te-tb;
		}		
		
		for (j=1; j<=par->total_pixel; j++) {
			r = top->rc_cont->co[j][1];
			c = top->rc_cont->co[j][2];
			
			H = Fmax(0., h[j]) / cos(top->slope->co[r][c]*Pi/180.);
			
			dV[j] = 0.0;
			
			if(H > par->min_hsup_land){
				
				area = ds*ds;
				area /= cos(top->slope->co[r][c]*Pi/180.);
				ch = cnet->ch->co[r][c];
				if (ch>0) area -= cnet->length->co[ch] * par->w_dx * ds; //area of the pixel[m2]
				
				Vmax = area * 1.E-3 * H;
				
				q = 0.;
				for (d=1; d<=4; d++) {
					q += top->Qdown->co[j][d];
				}
				
				if (q*dt > Vmax){
					q0 = q;
					q = Vmax/dt;
					for (d=1; d<=4; d++) {
						top->Qdown->co[j][d] *= (q/q0);
					}
				}
				
				for (d=1; d<=4; d++) {
					if(top->Jdown->co[j][d]==0 ){
						if (top->BC_counter->co[r][c] > 0 && (top->pixel_type->co[r][c] == 1 || top->pixel_type->co[r][c] == 2 || top->pixel_type->co[r][c] == 11 || top->pixel_type->co[r][c] == 12)){
							*Voutland = *Voutland + top->Qdown->co[j][d]*dt;
							mo += top->Qdown->co[j][d]*dt;
						}
					}
					if (q > 0) top->Qdown->co[j][d] /= q;
				}
				
				dV[j] = q*dt;
				
			}
		}
		
		for (j=1; j<=par->total_pixel; j++) {
			r = top->rc_cont->co[j][1];
			c = top->rc_cont->co[j][2];
			
			area = ds*ds;
			area /= cos(top->slope->co[r][c]*Pi/180.);
			ch = cnet->ch->co[r][c];
			if (ch>0) area -= cnet->length->co[ch] * par->w_dx * ds; //area of the pixel[m2]
			
			h[j] -= (1.E3*dV[j]/area) * cos(top->slope->co[r][c]*Pi/180.);
			
			if (top->BC_counter->co[r][c] > 0 && top->pixel_type->co[r][c] == -1) {
				if(h[j] > 0){
					h[j] += (1.E3*met->qinv[1]*ds*Dt/area) * cos(top->slope->co[r][c]*Pi/180.);
				}else{
					if( met->qinv[1]*ds > 0) h[j] = (1.E3*met->qinv[1]*ds*Dt/area) * cos(top->slope->co[r][c]*Pi/180.);
				}
			}
			
			for (d=1; d<=4; d++) {
				if(top->Jdown->co[j][d]>0){
					
					R = top->rc_cont->co[top->Jdown->co[j][d]][1];
					C = top->rc_cont->co[top->Jdown->co[j][d]][2];
					
					area = ds*ds;
					area /= cos(top->slope->co[R][C]*Pi/180.);
					ch = cnet->ch->co[R][C];
					if (ch>0) area -= cnet->length->co[ch] * par->w_dx * ds; //area of the pixel[m2]

					if(h[top->Jdown->co[j][d]]>0){
						h[top->Jdown->co[j][d]] += (1.E3*dV[j]*top->Qdown->co[j][d]/area) * cos(top->slope->co[R][C]*Pi/180.);
					}else{
						if( dV[j]*top->Qdown->co[j][d] > 0) h[top->Jdown->co[j][d]] = (1.E3*dV[j]*top->Qdown->co[j][d]/area) * cos(top->slope->co[R][C]*Pi/180.);
					}
					
				}
			}
			
		}			
		
		
		supflow_chla(dt, t, h, hch, top, wat, cnet, par, Vsup, flog, &cnt2);
		channel_flow(dt, t, DDch, hch, dhch, top, cnet, par, land, Voutnet, flog, &cnt3);
		
	}while(te<Dt);
	
	//printf("%f %f %f\n",Dt/cnt,Dt/cnt2,Dt/cnt3);
	for (j=1; j<=par->total_pixel; j++) {
		r = top->rc_cont->co[j][1];
		c = top->rc_cont->co[j][2];
		H = Fmax(0., h[j]) / cos(top->slope->co[r][c]*Pi/180.);
		area = ds*ds;
		area /= cos(top->slope->co[r][c]*Pi/180.);
		m2 += H*1.E-3*area;
	}	
	
	*mm1 = m1;
	*mm2 = m2;
	*mmo = mo;
	
	//printf("SUP: m1:%e m2:%e mo:%e dm:%e\n",m1,m2,mo,fabs(m1-mo-m2));

}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_dt_max_chla(double Courant, double *h, double *hch, TOPO *top, CHANNEL *cnet, PAR *par, double t, double *dt){
	
	double q, ds=sqrt(UV->U->co[1]*UV->U->co[2]), area, areach, Vmax, H, Hch, DH;
	long r, c, ch;
	
	for (ch=1; ch<=par->total_channel; ch++) {
		
		r = cnet->r->co[ch];
		c = cnet->c->co[ch];
		
		H = Fmax(0.0 , h[ch]) / cos(top->slope->co[r][c]*Pi/180.);//h[i] is the pressure at the surface, H is the depth of water normal to the surface
		area = ds*ds/cos(top->slope->co[r][c]*Pi/180.) - cnet->length->co[ch] * par->w_dx * ds;				
		
		Hch = Fmax(0., hch[ch] ) / cos(top->slope->co[r][c]*Pi/180.) - par->depr_channel * cos(top->slope->co[r][c]*Pi/180.);
		areach = cnet->length->co[ch] * par->w_dx * ds;
		
		if(H > par->min_hsup_land){
			
			if( Hch < -par->min_dhsup_land_channel_in ){	//free flow
				
				DH = H;
				q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
				Vmax = Fmin( areach*1.E-3*(-Hch) , area*1.E-3*H );
				
				Vmax = Fmax(Vmax, 1.E-10);
				q = Fmax(q, 1.E-30);
				
				if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
				
			}else if( H - Hch > par->min_dhsup_land_channel_in ){//submerged flow towards channel
				
				DH = H - Hch;
				q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
				Vmax = 1.E-3*DH / (1./area + 1./areach);
				
				Vmax = Fmax(Vmax, 1.E-10);
				q = Fmax(q, 1.E-30);
				
				if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
				
			}
		}
		
		if ( Hch - H > par->min_dhsup_land_channel_out ) {
			
			DH = Hch - H;
			q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*Hch;//m3/s
			
			Vmax = 1.E-3*DH / (1./area + 1./areach);
			
			Vmax = Fmax(Vmax, 1.E-10);
			q = Fmax(q, 1.E-30);
			
			if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
			
		}
	}
	
	if(*dt < par->dtmin_sup) *dt = par->dtmin_sup;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void supflow_chla(double Dt, double t, double *h, double *hch, TOPO *top, WATER *wat, CHANNEL *cnet, PAR *par, DOUBLEVECTOR *Vsup, FILE *flog, long *cnt){
	
	long ch, r, c;
	double ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	double H, Hch, DH, area, areach, q, tb, te=0., dt, Vmax;
	
	do{
		
		tb=te;
		dt=Dt;
		
		find_dt_max_chla(par->max_courant_land_channel, h, hch, top, cnet, par, t, &dt);		
		*cnt = *cnt + 1;
		
		te=tb+dt;
		if(te>Dt){
			te=Dt;
			dt=te-tb;
		}		
		
		//Superficial flow to the channels
		for(ch=1;ch<=par->total_channel;ch++){
			
			r = cnet->r->co[ch];
			c = cnet->c->co[ch];
			
			H = Fmax(0., h[top->j_cont[r][c]]) / cos(top->slope->co[r][c]*Pi/180.);
			Hch = Fmax(0., hch[ch] ) / cos(top->slope->co[r][c]*Pi/180.) - par->depr_channel * cos(top->slope->co[r][c]*Pi/180.);
			
			area = ds*ds/cos(top->slope->co[r][c]*Pi/180.) - cnet->length->co[ch] * par->w_dx * ds;
			areach = cnet->length->co[ch] * par->w_dx * ds;
			
			if(H > par->min_hsup_land){
				
				if( Hch < -par->min_dhsup_land_channel_in ){	//free flow
					
					DH = H;
					q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
					
					Vmax = Fmin( areach*1.E-3*(-Hch) , area*1.E-3*H );
					if (q*dt > Vmax) q = Vmax/dt;
					
					Vsup->co[ch] += q*dt;
					
					h[top->j_cont[r][c]] -= (1.E3 * q*dt/area) * cos(top->slope->co[r][c]*Pi/180.);
					
					if(hch[ch]>0){
						hch[ch] += (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*Pi/180.);	//mm;
					}else{
						if( q > 0 ) hch[ch] = (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*Pi/180.);	//mm;
					}				
					
				}else if( H - Hch > par->min_dhsup_land_channel_in ){//submerged flow towards channel
					
					DH = H - Hch;
					q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
					
					Vmax = 1.E-3*DH / (1./area + 1./areach);
					if (q*dt > Vmax) q = Vmax/dt;
					
					Vsup->co[ch] += q*dt;
					
					h[top->j_cont[r][c]] -= (1.E3 * q*dt/area) * cos(top->slope->co[r][c]*Pi/180.);
					
					if(hch[ch]>0){
						hch[ch] += (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*Pi/180.);	//mm;
					}else{
						if( q > 0 ) hch[ch] = (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*Pi/180.);	//mm;
					}		
				}
			}
			
			if ( Hch - H > par->min_dhsup_land_channel_out ) {
				
				DH = Hch - H;
				q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*Hch;//m3/s
				
				Vmax = 1.E-3*DH / (1./area + 1./areach);
				if (q*dt > Vmax) q = Vmax/dt;
				
				Vsup->co[ch] -= q*dt;
				
				hch[ch] -= (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*Pi/180.);
				
				if(h[top->j_cont[r][c]]>0){
					h[top->j_cont[r][c]] += (1.E3 * q*dt/area) * cos(top->slope->co[r][c]*Pi/180.);	//mm;
				}else{
					if( q > 0 ) h[top->j_cont[r][c]] = (1.E3 * q*dt/area) * cos(top->slope->co[r][c]*Pi/180.);	//mm;
				}
			}
		}
		
	}while(te<Dt);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_dt_max_channel(short DDcomplex, double Courant, double *h, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double t, double *dt){
	
	long r, c, ch, R, C;		
	double Ks, q, Vmax, i, H, dn, dD, ds;
	
	ds = sqrt(UV->U->co[1]*UV->U->co[2]);
	dn = par->w_dx*ds;
	
	for(ch=1;ch<=par->total_channel;ch++){
		
		//if(DDcomplex!=1 && t==0) draining_channel(0., ch, top->Z0, h, cnet, &(cnet->ch_down->co[ch]));
		
		r = cnet->r->co[ch];
		c = cnet->c->co[ch];
		H = Fmax(0., h[ch]) / cos(top->slope->co[r][c]*Pi/180.);
		
		if(H > par->min_hsup_channel){
			
			if(DDcomplex==1){
				draining_channel(1., ch, top->Z0, h, cnet, &(cnet->ch_down->co[ch]));
			}else {
				draining_channel(0., ch, top->Z0, h, cnet, &(cnet->ch_down->co[ch]));
			}
			
			Vmax = 1.E-3*H*dn*cnet->length->co[ch];	//m3
			
			if(top->is_on_border->co[r][c] == 1 && cnet->ch_down->co[ch]==ch){//outlet section
				
				q = Cd*(2./3.)*sqrt(2.*g*1.E-3*H)*(1.E-3*H)*dn;	//[m3/s]
				
			}else{
				
				R = cnet->r->co[cnet->ch_down->co[ch]];
				C = cnet->c->co[cnet->ch_down->co[ch]];
				
				if( (R-r==1 || R-r==-1) && (C-c==1 || C-c==-1) ){
					dD = find_3Ddistance(ds*sqrt(2.), top->Z0->co[r][c] - top->Z0->co[R][C]);
				}else {
					dD = find_3Ddistance(ds, top->Z0->co[r][c] - top->Z0->co[R][C]);
				}
				
				Ks = cm_h(par->Ks_channel, H, 1., par->thres_hchannel);
				
				i = ( (top->Z0->co[r][c] - top->Z0->co[R][C] ) + 1.E-3*(Fmax(0.0, h[ch]) - Fmax(0.0, h[cnet->ch_down->co[ch]])) ) / dD;	
				
				if(i<0) i=0.;
				
				q = dn * Ks * pow( 1.E-3 * H , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
				
			}
			
			Vmax = Fmax(Vmax, 1.E-10);
			q = Fmax(q, 1.E-10);
			if(Courant*Vmax/q<(*dt)) *dt=Courant*Vmax/q; 
			if(*dt<par->dtmin_sup) *dt=par->dtmin_sup;
			
		}else{
			
			cnet->ch_down->co[ch] = ch;
			
		}	
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void channel_flow(double Dt, double t, short DDcomplex, double *h, double *dV, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double *Vout, FILE *f, long *cnt)

{
	long r,c,ch,R,C;                                    
	double ds, dn, dD;
	double Ks;											  // the Strickler's coefficent
	double i;											  // hydraulic gradient
	double q,tb,te,dt,H;
	
	ds = sqrt(UV->U->co[1]*UV->U->co[2]);
	dn = par->w_dx*ds;
	
	if( par->point_sim==0 && cnet->r->co[1]!=0 ){	//if it is not point simulation and there are channels
		
		dn = par->w_dx * UV->U->co[1];		//transversal length [m]
		
		te=0.0;
		
		do{
			
			tb=te;
			dt=Dt;
			
			find_dt_max_channel(DDcomplex, par->max_courant_channel, h, top, cnet, par, land, t, &dt);
			*cnt = *cnt + 1;
			
			te=tb+dt;
			if(te>Dt){
				te=Dt;
				dt=te-tb;
			}		
			
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				
				dV[ch] = 0.0;
				
				H = Fmax(0., h[ch]) / cos(top->slope->co[r][c]*Pi/180.);
				
				if(H > par->min_hsup_channel){
					
					R = cnet->r->co[cnet->ch_down->co[ch]];
					C = cnet->c->co[cnet->ch_down->co[ch]];
					
					if(top->is_on_border->co[r][c] == 1 && cnet->ch_down->co[ch]==ch){//outlet section
						
						q = Cd*(2./3.)*sqrt(2.*g*1.E-3*H)*(1.E-3*H)*dn;	//[m3/s]
						
					}else{
						
						if( (R-r==1 || R-r==-1) && (C-c==1 || C-c==-1) ){
							dD = find_3Ddistance(ds*sqrt(2.), top->Z0->co[r][c] - top->Z0->co[R][C]);
						}else {
							dD = find_3Ddistance(ds, top->Z0->co[r][c] - top->Z0->co[R][C]);
						}
						
						Ks = cm_h(par->Ks_channel, H, 1., par->thres_hchannel);
						
						i= ( (top->Z0->co[r][c] - top->Z0->co[R][C] ) + 1.E-3*(Fmax(0.0, h[ch]) - Fmax(0.0, h[cnet->ch_down->co[ch]])) ) / dD;		
						
						if(i<0) i=0.;
						
						q = dn * Ks * pow( 1.E-3 * H , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
						
					}
					
					dV[ch] = Fmin( q*dt , 1.E-3*H*dn*cnet->length->co[ch] );	//m3
					
				}				
			}
			
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				
				h[ch] -= (1.E3*dV[ch]/(dn*cnet->length->co[ch])) * cos(top->slope->co[r][c]*Pi/180.);
				
				if(top->is_on_border->co[r][c] == 1 && cnet->ch_down->co[ch]==ch){//outlet section
					*Vout = *Vout + dV[ch];	//m3
				}else {
					R = cnet->r->co[cnet->ch_down->co[ch]];
					C = cnet->c->co[cnet->ch_down->co[ch]];					
					if(h[cnet->ch_down->co[ch]]>0){
						h[cnet->ch_down->co[ch]] += (1.E3*dV[ch]/(dn*cnet->length->co[cnet->ch_down->co[ch]])) * cos(top->slope->co[R][C]*Pi/180.);	//mm;
					}else{
						if( dV[ch] > 0) h[cnet->ch_down->co[ch]] = (1.E3*dV[ch]/(dn*cnet->length->co[cnet->ch_down->co[ch]])) * cos(top->slope->co[R][C]*Pi/180.);	//mm;
					}
				}							
			}
			
		}while(te<Dt);
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void draining_land(double alpha, long i, TOPO *T, LAND *L, PAR *P, CHANNEL *cnet, double *h, long *I, double *Q){
	
	double H, p, pn, dD, dn, Ks;
	double ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	
	long d, r, c;
	long ir[5] = {0, -1, 1, 0,  0};
	long ic[5] = {0,  0, 0, -1, 1};
	
	if (h[i] > 0) {
		
		r = T->rc_cont->co[i][1];
		c = T->rc_cont->co[i][2];
		H = Fmax(h[i], 0.)/cos(T->slope->co[r][c]*Pi/180.);
		p = T->Z0->co[r][c] + alpha*1.E-3*Fmax(h[i], 0.);
		Ks = cm_h(L->ty->co[(short)L->LC->co[r][c]][jcm], H, P->thres_hsup_1,  P->thres_hsup_2);
		
		for (d=1; d<=4; d++) {
			
			if (r+ir[d]>=1 && r+ir[d]<=Nr && c+ic[d]>=1 && c+ic[d]<=Nc) {
				I[d] = T->j_cont[r+ir[d]][c+ic[d]];
			}else {
				I[d] = 0;
			}
						
			if (I[d]>0) {		
				
				dD = find_3Ddistance(ds, T->Z0->co[r][c] - T->Z0->co[r+ir[d]][c+ic[d]]);
				dn = ds;
				
				if(ir[d]==1 || ir[d]==-1){
					dn /= cos(0.5*atan(T->dzdE->co[r][c])+0.5*atan(T->dzdE->co[r+ir[d]][c+ic[d]]));
				}else {
					dn /= cos(0.5*atan(T->dzdN->co[r][c])+0.5*atan(T->dzdN->co[r+ir[d]][c+ic[d]]));
				}		
				
				pn = T->Z0->co[r+ir[d]][c+ic[d]] + alpha*1.E-3*Fmax(h[I[d]], 0.);
				
				if (pn < p) {
					Q[d] = Ks*dn*pow(1.E-3*H, 1.0+P->gamma_m)*sqrt((p-pn)/dD);					
				}else {
					Q[d] = 0.;
				}
				
			}else {
				
				if(T->BC_counter->co[r][c] > 0){
					if (H >= -T->BC_DepthFreeSurface->co[T->BC_counter->co[r][c]] ){
						Q[d] = Cd*(2./3.)*sqrt(2.*g*1.E-3*H)*(1.E-3*H)*ds;
						if(cnet->ch->co[r][c]>0) Q[d] = Q[d] * (1.-P->w_dx);
					}else {
						I[d] = i;
						Q[d] = 0.;
					}
				}else {
					I[d] = i;
					Q[d] = 0.;
				}
				
			}
			
		}
		
	}else {
		
		for (d=1; d<=4; d++) {
			I[d] = i;
			Q[d] = 0.;
		}
		
	}
}



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void draining_channel(double alpha, long ch, DOUBLEMATRIX *Z, double *h, CHANNEL *cnet, long *CH){
	
	long d, r, c;
	double elev, elev1;
	long ir[9] = {0, -1, -1, 0, 1, 1, 1, 0, -1};
	long ic[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
	
	*CH = ch;
	
	r = cnet->r->co[ch];
	c = cnet->c->co[ch];
	elev = Z->co[r][c] + alpha*1.E-3*Fmax(h[ch], 0.);	
	
	for (d=1; d<=8; d++) {
		if (r+ir[d]>=1 && r+ir[d]<=Nr && c+ic[d]>=1 && c+ic[d]<=Nc) {
			if(cnet->ch->co[r+ir[d]][c+ic[d]] > 0){
				elev1 = Z->co[r+ir[d]][c+ic[d]] + alpha*1.E-3*Fmax(h[cnet->ch->co[r+ir[d]][c+ic[d]]], 0.);
				if( elev1 < elev){
					elev = elev1;
					*CH = cnet->ch->co[r+ir[d]][c+ic[d]];
				}
			}
		}
		
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


