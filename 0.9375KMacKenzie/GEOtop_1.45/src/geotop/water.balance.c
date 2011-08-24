
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "constants.h"
#include "struct.geotop.h"
#include "pedo.funct.h"
#include "../libraries/math/sparse_matrix.h"
#include "../libraries/math/util_math.h"
#include "water.balance.h"

#include <time.h>

extern long number_novalue, number_absent;

extern T_INIT *UV;
extern char **files;
extern char *logfile;

extern long Nl, Nr, Nc;
extern double *outdata_basin;
extern double t_sub, t_sup;

//subsurface flow constants
#define tol_max_GC 1.E+5
#define tol_min_GC 1.E-13
#define max_res_adm 1.E-2
#define M 1
#define ni 1.E-7
#define max_slope 89.999

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void water_balance(ALLDATA *adt){

	short out, sux;
	long r, c, l, k=0, ch, iter, sy;
	double Dt0, Dt, tb, te, te0, Loss;
	clock_t start, end;
	FILE *flog;		
	
	flog = fopen(logfile, "a");
		
	initialize_doublevector(adt->C->Vsub, 0.0);
	initialize_doublevector(adt->C->Vsup, 0.0);
	
	adt->C->Vout = 0.;
	adt->W->Voutland = 0.;

	te0=0.0;
		
	do{

		//find Dt0 to solve water balance
		Dt0=adt->P->Dt;
		if(te0+Dt0 > adt->P->Dt) Dt0=adt->P->Dt-te0;
		
		//surface flow: 1st half of time step
		start=clock();
		supflow(Dt0/2., adt->I->time, 1, 1, adt->S->P->co[0], adt->W->h_sup->co, adt->C->P->co[0], adt->C->h_sup->co, 
				adt->T, adt->L, adt->W, adt->C, adt->P, &(adt->C->Vout), &(adt->W->Voutland), flog);

		end=clock();
		t_sup += (end-start)/(double)CLOCKS_PER_SEC;	
		
		//subsurface flow with time step Dt0 (decreasing if not converging)
		start = clock();
		te=te0;

		do{
			tb=te;	
			Dt=Dt0; 
			
			do{
				if (te+Dt>Dt0) Dt=Dt0-te;
				if (Dt<adt->P->Dt){
					printf("Water Balance Time step:%f\n",Dt);	
					fprintf(flog,"Water Balance Time step:%f\n",Dt);	
				}
				sux = Richards(Dt, &Loss, &iter, adt, flog);
				out=1;
				//if not converging, it reduces the time step				
				if(sux==0){
					Dt/=2.;
					k++;
					out=0;
				}
			}while( out==0 && k<=adt->P->max_times_halving_time_step_wat); 
									
			if(sux==0){
				printf("ERROR:Water balance does not converge\n");
				printf("It is not possible to continue, Please check the parameters in the block 2 of the parameter file\n");
				printf("or reduce the time step or review the initial conditions\n\n");
				printf("If you think that everything is right, please report this error to the GEOtop community\n\n");
				
				fprintf(flog,"ERROR:Water balance does not converge\n");
				fprintf(flog,"It is not possible to continue, Please check the parameters in the block 2 of the parameter file\n");
				fprintf(flog,"or reduce the time step or review the initial conditions\n\n");
				fprintf(flog,"If you think that everything is right, please report this error to the GEOtop community\n\n");
				fclose(flog);
				
				t_error("Code execution interrupted");
			}
			
			if (Dt<Dt0){
				Dt*=2.;
				k--;
			}
			
			te=tb+Dt;   
						
			outdata_basin[oomasserror] += fabs(Loss);
		
		}while(te<Dt0);
		end=clock();
		t_sub += (end-start)/(double)CLOCKS_PER_SEC;

		//surface flow: 2nd half of time step
		start=clock();
		supflow(Dt0/2., adt->I->time, 1, 1, adt->S->P->co[0], adt->W->h_sup->co, adt->C->P->co[0], adt->C->h_sup->co, 
				adt->T, adt->L, adt->W, adt->C, adt->P, &(adt->C->Vout), &(adt->W->Voutland), flog);
		end=clock();
		t_sup += (end-start)/(double)CLOCKS_PER_SEC;
		
		te0+=Dt0;	
		
	}while(te0<adt->P->Dt);
	
	//Final stuff
		
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if((long)adt->L->LC->co[r][c]!=number_novalue){
				for(l=1;l<=Nl;l++){
					sy = adt->S->type->co[r][c];
					
					adt->S->th->co[l][r][c] = theta_from_psi( adt->S->P->co[l][r][c], l, r, c, adt->S, PsiMin );

					//total water pressure (liq + ice)
					adt->S->Ptot->co[l][r][c] = psi_teta(adt->S->th->co[l][r][c]+adt->S->thice->co[l][r][c], 0.0, adt->S->pa->co[sy][jsat][l],
														 adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 
														 1.-1./adt->S->pa->co[sy][jns][l], PsiMin, adt->S->pa->co[sy][jss][l]);
					
					//th is the water content (always <= saturation)
					adt->S->th->co[l][r][c] = Fmin( adt->S->th->co[l][r][c] , adt->S->pa->co[sy][jsat][l]-adt->S->thice->co[l][r][c] );
					
				}
			}
		}
	}
	
	
	for (ch=1; ch<=adt->P->total_channel; ch++) {
		sy = adt->C->soil_type->co[ch];
		for(l=1;l<=Nl;l++){
			adt->C->th->co[l][ch] = teta_psi(adt->C->P->co[l][ch], adt->C->thice->co[l][ch], adt->S->pa->co[sy][jsat][l],
										 adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l],
										 1.-1./adt->S->pa->co[sy][jns][l], PsiMin, adt->S->pa->co[sy][jss][l]);
			adt->C->th->co[l][ch] = Fmin(adt->C->th->co[l][ch], adt->S->pa->co[sy][jsat][l]-adt->C->thice->co[l][ch]);
		}
	}
	
	fclose(flog);
				
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
 
 
short Richards(double Dt, double *loss, long *iter, ALLDATA *adt, FILE *f){
	
	double res=0.0, res0[3], res_prev[M], res_av, res00, lambda[3], epsilon, mu, hold, hnew;
	double ds=sqrt(UV->U->co[1]*UV->U->co[2]);

	long i, ch, l, r, c, m, cont, cont2, iter_tot=0;
	long n=(Nl+1)*adt->P->total_pixel;	
	long N=adt->W->H0->nh;
	long cont_lambda_min=0;
	static long cum_iter;	
	short out, out2;	
	int sux;
					
	if(adt->I->time == 0.0) cum_iter = 0;
		
	*iter = 0;	//conjugated gradient iteration number	
		
	for(i=1; i<=N; i++){
		
		if (i<=n) {
			l = adt->T->lrc_cont->co[i][1];
			r = adt->T->lrc_cont->co[i][2];
			c = adt->T->lrc_cont->co[i][3];
			adt->W->P0->co[i] = adt->S->P->co[l][r][c];
			
			//add precipitation first
			if (adt->W->Pnet->co[r][c] > 0 && l==0){
				if(adt->W->P0->co[i]<0) {
					adt->W->P0->co[i] = find_sup_pressure(adt->W->Pnet->co[r][c]*(Dt/adt->P->Dt), adt->T->slope->co[r][c]);
				}else {
					adt->W->P0->co[i] += find_sup_pressure(adt->W->Pnet->co[r][c]*(Dt/adt->P->Dt), adt->T->slope->co[r][c]);
				}
			}

			//calculates piezometric head (=pressure head + elevation)
			adt->W->H0->co[i] = adt->W->P0->co[i] + adt->T->Z->co[l][r][c];
			
		}else {
			l = adt->C->lch->co[i-n][1];
			ch = adt->C->lch->co[i-n][2];
			r = adt->C->r->co[ch];
			c = adt->C->c->co[ch];
			adt->W->P0->co[i] = adt->C->P->co[l][ch];
			
			//add precipitation first
			if (adt->W->Pnet->co[r][c] > 0 && l==0){
				if(adt->W->P0->co[i]<0) {
					adt->W->P0->co[i] = find_sup_pressure(adt->W->Pnet->co[r][c]*(Dt/adt->P->Dt), adt->T->slope->co[r][c]);
				}else {
					adt->W->P0->co[i] += find_sup_pressure(adt->W->Pnet->co[r][c]*(Dt/adt->P->Dt), adt->T->slope->co[r][c]);
				}
			}
			
			//calculates piezometric head (=pressure head + elevation)
			adt->W->H0->co[i] = adt->W->P0->co[i] + ( adt->T->Z->co[l][r][c] - adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.) );	
		}
		
		
		/*
		Layer 0 is water on the surface. This allows a more robust description of infiltration processes.
		However, surface water flow is described separately in the subroutine supflow
		
		Layers > 0 , represent subsurface
		*/
				
		adt->W->H1->co[i] = adt->W->H0->co[i];
				
		if(adt->W->H1->co[i] != adt->W->H1->co[i]) printf("no value in r:%ld c:%ld l:%ld\n",r,c,l);
		
	}
		
	sux = find_matrix_K(adt->W->Lx, adt->W->Klat, adt, adt->W->H1, Dt);
			
	find_f(adt->W->f, adt, adt->W->H1, adt->W->Klat, Dt);
	
	product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
		
	res = norm_inf(adt->W->B, 1, N);
	
	res00 = res; //initial norm of the residual
	epsilon = adt->P->TolVWb + adt->P->RelTolVWb * Fmin( res00 , sqrt((double)N) );

	if(res!=res) printf("res in no value\n");

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
		}

		initialize_doublevector(adt->W->dH, 0.);
		
		//mu is the forcing term of Newton Raphson, defined such that ||F( x + d )|| < mu * ||F(x)|| see references
		if (cont==1) {
			mu = adt->P->TolCG;
		}else{
			mu *= Fmin( 1.0 , res/res0[0] );
			if(mu < 0.5*epsilon/res) mu = 0.5*epsilon/res;
		}
		
		//CALCOLATE AND STORE JACOBIAN AS SPARSE MATRIX
		sux = find_dfdH(adt->W->df, adt, adt->W->H1, adt->W->Klat, Dt);	//it calcolates only df/dH, J = I*df/dH + K, K is calculated above
		
		//CONJUGATED GRADIENTS ALGORITHM
		*iter = BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->W->df, adt->T->Li, adt->T->Lp, adt->W->Lx);
		if(*iter==-1) return 0;	//does not converge 
		iter_tot += (*iter);//The number of the cumulated GC iterations is stored
		
		//non-monotonic line search (it is monotonic if M==1)
		for(m=Fminlong(cont,M); m>1; m--){
			res_prev[m-1]=res_prev[m-2];
		}
		res_prev[0]=res;
		
		res_av=0.0;
		for(m=1; m<=Fminlong(cont,M); m++){
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
					printf("no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);
					stop_execution();
				}
				
			}
			
			if(adt->P->UpdateK == 1) sux = find_matrix_K(adt->W->Lx, adt->W->Klat, adt, adt->W->H1, Dt);		
									
			find_f(adt->W->f, adt, adt->W->H1, adt->W->Klat, Dt);

			product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
			
			res = norm_inf(adt->W->B, 1, N);		
									
			out2=0;
			
			if(res <= (1.0 - ni*lambda[0]*(1.-mu))*res_av) out2=1;
			if(lambda[0] <= adt->P->min_lambda_wat) cont_lambda_min++;
			//printf("res:%e lambda:%e cont:%ld Dt:%f %ld\n",res,lambda[0],cont,Dt,cont_lambda_min);
			if(cont_lambda_min > adt->P->max_times_min_lambda_wat){
				if(adt->P->exit_lambda_min_wat == 1){
					return 0;
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
		
	cum_iter += iter_tot;
	
	if( res > epsilon ) return 0;

	//it can be shown that massloss per unit pixel [mm] is the linear norm of B * Dt / total_pixel
	*loss = norm_1(adt->W->B, 1, N)*Dt/adt->P->total_area;
		
	//assign updated state variables	
	for(i=1; i<=N; i++){
		
		if (i<=n) {
			l = adt->T->lrc_cont->co[i][1];
			r = adt->T->lrc_cont->co[i][2];
			c = adt->T->lrc_cont->co[i][3];
			adt->S->P->co[l][r][c] = adt->W->H1->co[i] - adt->T->Z->co[l][r][c];
			
		}else {
			l = adt->C->lch->co[i-n][1];
			ch = adt->C->lch->co[i-n][2];
			r = adt->C->r->co[ch];
			c = adt->C->c->co[ch];
			
			adt->W->P0->co[i] = adt->W->H1->co[i] - ( adt->T->Z->co[l][r][c] - adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.) );

			if (l==0){
				hold = find_hsup(Fmax(0., adt->C->P->co[l][ch]), adt->T->slope->co[r][c]);
				hnew = find_hsup(Fmax(0., adt->W->P0->co[i]), adt->T->slope->co[r][c]);			
				adt->C->Vsub->co[ch] += 1.E-3 * ( hnew - hold ) * adt->C->length->co[ch] * adt->P->w_dx * ds;
			}
			
			adt->C->P->co[l][ch] = adt->W->P0->co[i];

		}
	}
	
	return 1;
			
	
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
int find_matrix_K(DOUBLEVECTOR *Lx, DOUBLEMATRIX *Klat, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i, l, r, c, I, R, C, sy, syn, ch, cnt=0;
	long n=(Nl+1)*adt->P->total_pixel;
	double dz=0.0, dzn=0.0, dD=0.0, k=0.0, kn=0.0, kmax=0.0, kmaxn=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]), dn;
	
	for(i=1;i<=H->nh;i++){
		
		//VERTICAL FLUXES
		if( i<=n){//land

			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			sy=adt->S->type->co[r][c];

			ch=adt->C->ch->co[r][c];
			
			area=ds*ds;
			if (adt->P->point_sim!=1) area/=cos(adt->T->slope->co[r][c]*Pi/180.);
			if (ch>0) area-=adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
				
			//vertical hydraulic conductivity
			if(l>0){
				dz = adt->S->pa->co[sy][jdz][l];
				k = k_from_psi(jKn, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
						
			//flux from cell below
			if (l<Nl) {
				
				I = i+1;
				
				if(l==0){	//overland flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dzn;

					if( H->co[i] < H->co[I] ){	
						//upward flux
						kn = k_from_psi(jKn,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
					}else{	
						//downward flow
						kn = k_from_psi(jKn,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
					}
					
					
				}else{	//subsurface flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dz + 0.5*dzn;
					
					kn = k_from_psi(jKn,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
					kn = Mean(adt->P->harm_or_arit_mean_normal, dz, dzn, k, kn);
					kmax = k_from_psi(jKn,  psisat_from(l,   r, c, adt->S), l,   r, c, adt->S, adt->P->imp);
					kmaxn = k_from_psi(jKn,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp);
					kmaxn = Harmonic_Mean(dz, dzn, kmax, kmaxn);
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
				k = K( H->co[i] - (adt->T->Z->co[l][r][c]-adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.)), 
					  adt->S->pa->co[sy][jKn][l], adt->P->imp, adt->C->thice->co[l][ch], adt->S->pa->co[sy][jsat][l], 
					  adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 
					  1.-1./adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], PsiMin, adt->C->T->co[l][ch] );
			}
			
			//flux from cell below
			if (l<Nl) {
				
				I = i+1;
				
				if(l==0){	//overland flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dzn;
					
					if( H->co[i] < H->co[I] ){	
						//upward flux
						kn = K( H->co[I] - (adt->T->Z->co[l+1][r][c]-adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.)), 
							   adt->S->pa->co[sy][jKn][l+1], adt->P->imp, adt->C->thice->co[l+1][ch], adt->S->pa->co[sy][jsat][l+1], 
							   adt->S->pa->co[sy][jres][l+1], adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 
							   1.-1./adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], PsiMin, adt->C->T->co[l+1][ch] );
					}else{	
						//downward flow
						kn = K( psi_saturation(adt->C->thice->co[l+1][ch], adt->S->pa->co[sy][jsat][l+1], adt->S->pa->co[sy][jres][l+1],
											   adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 1.-1./adt->S->pa->co[sy][jns][l+1]),
							   adt->S->pa->co[sy][jKn][l+1], adt->P->imp, adt->C->thice->co[l+1][ch], adt->S->pa->co[sy][jsat][l+1], 
							   adt->S->pa->co[sy][jres][l+1], adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 
							   1.-1./adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], PsiMin, adt->C->T->co[l+1][ch] );
					}
										
				}else{	//subsurface flow
					
					dzn = adt->S->pa->co[sy][jdz][l+1];
					dD = 0.5*dz + 0.5*dzn;
					
					kn = K( H->co[I] - (adt->T->Z->co[l+1][r][c]-adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.)), 
						   adt->S->pa->co[sy][jKn][l+1], adt->P->imp, adt->C->thice->co[l+1][ch], adt->S->pa->co[sy][jsat][l+1], 
						   adt->S->pa->co[sy][jres][l+1], adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 
						   1.-1./adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], PsiMin, adt->C->T->co[l+1][ch] );
					kn = Mean(adt->P->harm_or_arit_mean_normal, dz, dzn, k, kn);
					
					kmax = K( psi_saturation(adt->C->thice->co[l][ch], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
											 adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1.-1./adt->S->pa->co[sy][jns][l]),
						   adt->S->pa->co[sy][jKn][l], adt->P->imp, adt->C->thice->co[l][ch], adt->S->pa->co[sy][jsat][l], 
						   adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 
						   1.-1./adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], PsiMin, adt->C->T->co[l][ch] );
					
					kmaxn = K( psi_saturation(adt->C->thice->co[l+1][ch], adt->S->pa->co[sy][jsat][l+1], adt->S->pa->co[sy][jres][l+1],
											  adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 1.-1./adt->S->pa->co[sy][jns][l+1]),
							  adt->S->pa->co[sy][jKn][l+1], adt->P->imp, adt->C->thice->co[l+1][ch], adt->S->pa->co[sy][jsat][l+1], 
							  adt->S->pa->co[sy][jres][l+1], adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 
							  1.-1./adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], PsiMin, adt->C->T->co[l+1][ch] );
					
					kmaxn = Harmonic_Mean(dz, dzn, kmax, kmaxn);
					kn = Fmin(kn, kmaxn);
					
				}
				
				cnt++;
				Lx->co[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
				
			}
		}
			
		//LATERAL FLUXES
		if (i<=n && l>0){
				
			//lateral hydraulic conductivity
			k = k_from_psi(jKl, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );			
			if (adt->T->pixel_type->co[r][c] == 1) Klat->co[adt->T->BC_counter->co[r][c]][l] = k;
			
			if (adt->P->point_sim!=1){
				
				//4 neighbouring cells
				R = r-1;
				C = c;
				dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));
				//1.
				if(R>=1 && R<=Nr && C>=1 && C<=Nc){
					if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
						
						I=adt->T->i_cont[l][R][C];	
						
						kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
						kn = Mean(adt->P->harm_or_arit_mean_parallel, 1., 1., k, kn);
						kmax = k_from_psi(jKl, psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
						kmaxn = k_from_psi(jKl,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp );
						kmaxn = Harmonic_Mean(1., 1., kmax, kmaxn);
						kn = Fmin(kn, kmaxn);

						dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
				
				R = r+1;
				C = c;
				dn = ds/cos(0.5*atan(adt->T->dzdE->co[r][c])+0.5*atan(adt->T->dzdE->co[R][C]));
				//2.
				if(R>=1 && R<=Nr && C>=1 && C<=Nc){
					if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
						
						I=adt->T->i_cont[l][R][C];	
						
						kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
						kn = Mean(adt->P->harm_or_arit_mean_parallel, 1., 1., k, kn);
						kmax = k_from_psi(jKl, psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
						kmaxn = k_from_psi(jKl,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp );
						kmaxn = Harmonic_Mean(1., 1., kmax, kmaxn);
						kn = Fmin(kn, kmaxn);

						dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]

					}
				}
				
				R = r;
				C = c-1;
				dn = ds/cos(0.5*atan(adt->T->dzdN->co[r][c])+0.5*atan(adt->T->dzdN->co[R][C]));
				//3.
				if(R>=1 && R<=Nr && C>=1 && C<=Nc){
					if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
						
						I=adt->T->i_cont[l][R][C];	
						
						kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
						kn = Mean(adt->P->harm_or_arit_mean_parallel, 1., 1., k, kn);
						kmax = k_from_psi(jKl, psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
						kmaxn = k_from_psi(jKl,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp );
						kmaxn = Harmonic_Mean(1., 1., kmax, kmaxn);
						kn = Fmin(kn, kmaxn);

						dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]

					}
				}
				
				
				R = r;
				C = c+1;
				dn = ds/cos(0.5*atan(adt->T->dzdN->co[r][c])+0.5*atan(adt->T->dzdN->co[R][C]));
				//4.
				if(R>=1 && R<=Nr && C>=1 && C<=Nc){
					if((long)adt->L->LC->co[R][C]!=number_novalue && adt->T->i_cont[l][R][C]>i){ 
						
						I=adt->T->i_cont[l][R][C];	
						
						kn = k_from_psi(jKl, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
						kn = Mean(adt->P->harm_or_arit_mean_parallel, 1., 1., k, kn);
						kmax = k_from_psi(jKl, psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
						kmaxn = k_from_psi(jKl,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp );
						kmaxn = Harmonic_Mean(1., 1., kmax, kmaxn);
						kn = Fmin(kn, kmaxn);

						dD = find_3Ddistance(ds, adt->T->Z0->co[r][c]-adt->T->Z0->co[R][C]) * 1.E3;//[mm]
						
						cnt++;
						Lx->co[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]

					}
				}
				
				//exchange with channels
				if (adt->C->ch->co[r][c] > 0) {
					
					ch = adt->C->ch->co[r][c];
					syn = adt->C->soil_type->co[ch];
					I = n + adt->C->ch3[l][ch];
					
					kn = K( H->co[I] - (adt->T->Z->co[l][r][c]-adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.)),
						   adt->S->pa->co[syn][jKl][l], adt->P->imp, adt->C->thice->co[l][ch], adt->S->pa->co[syn][jsat][l], 
						   adt->S->pa->co[syn][jres][l], adt->S->pa->co[syn][ja][l], adt->S->pa->co[syn][jns][l], 
						   1.-1./adt->S->pa->co[syn][jns][l], adt->S->pa->co[syn][jv][l], PsiMin, adt->C->T->co[l][ch] );
										
					kn = Mean(adt->P->harm_or_arit_mean_parallel, sqrt(area), sqrt(adt->C->length->co[ch] * adt->P->w_dx * ds), k, kn);
					
					kmax = k_from_psi(jKl, psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
					
					kmaxn = K( psi_saturation(adt->C->thice->co[l][ch], adt->S->pa->co[syn][jsat][l], adt->S->pa->co[syn][jres][l],
											  adt->S->pa->co[syn][ja][l], adt->S->pa->co[syn][jns][l], 1.-1./adt->S->pa->co[syn][jns][l]),
							  adt->S->pa->co[syn][jKn][l], adt->P->imp, adt->C->thice->co[l][ch], adt->S->pa->co[syn][jsat][l], 
							  adt->S->pa->co[syn][jres][l], adt->S->pa->co[syn][ja][l], adt->S->pa->co[syn][jns][l], 
							  1.-1./adt->S->pa->co[syn][jns][l], adt->S->pa->co[syn][jv][l], PsiMin, adt->C->T->co[l][ch] );
					
					kmaxn = Harmonic_Mean(sqrt(area), sqrt(adt->C->length->co[ch] * adt->P->w_dx * ds), kmax, kmaxn);
					kn = Fmin(kn, kmaxn);

					//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]	
					//Area[m2] = 2 * channel length * layer thickness = 2 * length[m] * dz[mm] * 1.E-3[m/mm]
					//dD[mm] = 0.25 * ds * (1+w_dx)
					dD = find_3Ddistance(ds * (1.0 + adt->P->w_dx) / 4.0, 1.E-3*adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.)) * 1.E3;//[mm]
					
					cnt++;
					Lx->co[cnt] = -(2.*adt->C->length->co[adt->C->ch->co[r][c]]*1.E-3*dz)*kn/dD;						
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

int find_dfdH(DOUBLEVECTOR *df, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *K, double Dt){
	
	long i, r, l, c, sy, ch, bc;
	long n=(Nl+1)*adt->P->total_pixel;
	double dz, Hboundary, psi1, ice=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]);
		
	for(i=1;i<=H->nh;i++){
		
		if (i<=n) {
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			sy=adt->S->type->co[r][c];
			bc=adt->T->BC_counter->co[r][c];
			ch=adt->C->ch->co[r][c];
			area=ds*ds;
			if (adt->P->point_sim!=1) area/=cos(adt->T->slope->co[r][c]*Pi/180.);
			if (ch>0) area-=adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			psi1 = H->co[i] - adt->T->Z->co[l][r][c];
			if(l>0) ice = adt->S->thice->co[l][r][c];

		}else {
			l=adt->C->lch->co[i-n][1];
			ch=adt->C->lch->co[i-n][2];
			r=adt->C->r->co[ch];
			c=adt->C->c->co[ch];
			sy=adt->C->soil_type->co[ch];
			bc=0;
			area=adt->C->length->co[ch] * adt->P->w_dx * ds;
			psi1 = H->co[i] - (adt->T->Z->co[l][r][c] - adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.));
			if(l>0) ice = adt->C->thice->co[l][ch];
		}
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			if(psi1>0){
				df->co[i] = area*find_dhsup(adt->T->slope->co[r][c]) /Dt;
			}else{
				df->co[i] = 0.;
			}
		}else{
			dz = adt->S->pa->co[sy][jdz][l];
			
			df->co[i] = dteta_dpsi(psi1, ice, adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], 
								   adt->S->pa->co[sy][jns][l], 1.-1./adt->S->pa->co[sy][jns][l], PsiMin, adt->S->pa->co[sy][jss][l])*
						area * dz / Dt;
		
		}
		
		//lateral drainage at the border
		if (bc>0 && l>0) {
			//Head at the boundary: constant in the vertical (hydrostatic)
			Hboundary = adt->T->Z->co[0][r][c] - 1.E3*adt->T->BC_LatDistance->co[bc]*sin(adt->T->slope->co[r][c]*Pi/180.)
						- adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.);
			//Not admitted flow towards the boundary pixel from the border (the border is always a sink)
			if (H->co[i] - Hboundary > 0) {
				dz = adt->S->pa->co[sy][jdz][l];//[mm]
				df->co[i] += (ds*dz*1.E-3) * K->co[bc][l]/(1.E3*adt->T->BC_LatDistance->co[bc]);
			}
		}
	}

	return 0;
}
		
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_f(DOUBLEVECTOR *f, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *K, double Dt){

	long i, l, r, c, sy, ch, bc;
	long n=(Nl+1)*adt->P->total_pixel;
	double dz, V0, V1, Hboundary, psi1, psi0, ice=0.0;
	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	
	for(i=1;i<=H->nh;i++){
		
		psi0 = adt->W->P0->co[i];

		if (i<=n) {
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			sy=adt->S->type->co[r][c];
			bc=adt->T->BC_counter->co[r][c];
			ch=adt->C->ch->co[r][c];
			area=ds*ds;
			if (adt->P->point_sim!=1) area/=cos(adt->T->slope->co[r][c]*Pi/180.);
			if (ch>0) area-=adt->C->length->co[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			psi1 = H->co[i] - adt->T->Z->co[l][r][c];
			if(l>0) ice = adt->S->thice->co[l][r][c];
			
		}else {
			l=adt->C->lch->co[i-n][1];
			ch=adt->C->lch->co[i-n][2];
			r=adt->C->r->co[ch];
			c=adt->C->c->co[ch];
			sy=adt->C->soil_type->co[ch];
			bc=0;
			area=adt->C->length->co[ch] * adt->P->w_dx * ds;
			psi1 = H->co[i] - (adt->T->Z->co[l][r][c] - adt->P->depr_channel*cos(adt->T->slope->co[r][c]*Pi/180.));
			if(l>0) ice = adt->C->thice->co[l][ch];

		}

		
		//hydraulic capacity (diagonal term)
		if(l==0){
			V1 = area * find_hsup(Fmax(0.0, psi1), adt->T->slope->co[r][c]);
			V0 = area * find_hsup(Fmax(0.0, psi0), adt->T->slope->co[r][c]);
		}else{
			dz = adt->S->pa->co[sy][jdz][l];		
			
			V1 = area*dz * teta_psi(psi1, ice, adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
									adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 
									1.-1./adt->S->pa->co[sy][jns][l], PsiMin, adt->S->pa->co[sy][jss][l]);
			V0 = area*dz * teta_psi(psi0, ice, adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
									adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 
									1.-1./adt->S->pa->co[sy][jns][l], PsiMin, adt->S->pa->co[sy][jss][l]);
			
		}
		
		f->co[i] = (V1-V0)/Dt;
		
		//lateral drainage at the border
		if (bc>0 && l>0) {
			//Head at the boundary: constant in the vertical (hydrostatic)
			Hboundary = adt->T->Z->co[0][r][c] - 1.E3*adt->T->BC_LatDistance->co[bc]*sin(adt->T->slope->co[r][c]*Pi/180.)
						- adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*Pi/180.);
			//Not admitted flow towards the boundary pixel from the border (the border is always a sink)
			if (H->co[i] - Hboundary > 0) {
				dz = adt->S->pa->co[sy][jdz][l];//[mm]
				f->co[i] += (ds*dz*1.E-3) * K->co[bc][l]*(H->co[i] - Hboundary)/(1.E3*adt->T->BC_LatDistance->co[bc]);
			}
		}
					
		//evaporation 
		if(l>0){
			if(i<=n){
				f->co[i] += area*adt->S->ET->co[l][r][c]/adt->P->Dt;
			}else {
				ch=adt->C->ch->co[r][c];
				f->co[i] += area*adt->C->ET->co[l][ch]/adt->P->Dt;
			}
		}
						
	}
	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double find_hsup(double Psurface, double slope){

	return Fmax(0.,Psurface)/cos(Fmin(max_slope,slope)*Pi/180.);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double find_dhsup(double slope){
	
	return 1./cos(Fmin(max_slope,slope)*Pi/180.);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double find_sup_pressure(double hsup, double slope){
	
	return Fmax(0.,hsup)*cos(Fmin(max_slope,slope)*Pi/180.);
	
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

void find_dt_max(short DDcomplex, double Courant, double **h, double *hch, LAND *land, TOPO *top, CHANNEL *cnet, PAR *par, double t, double *dt){
	
	double dn, dD;
	double i, q, ds=sqrt(UV->U->co[1]*UV->U->co[2]), Ks, area, AREA, areach, Vmax, H, Hch, DH;
	double alpha = 0.0;
	short lu;
	long r, c, R, C, ch, CH;
	double min_dhsup_land_channel=20.;
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			
			if( (long)land->LC->co[r][c]!=number_novalue){
				
				if(DDcomplex!=1 && t==0) draining_land(0., r, c, top->Z0, h, land->LC, &(top->Rdown->co[r][c]), &(top->Cdown->co[r][c]));
				
				H = find_hsup(h[r][c], top->slope->co[r][c]);
				
				if(H > par->min_hsup_land){
					
					if(DDcomplex==1) draining_land(1., r, c, top->Z0, h, land->LC, &(top->Rdown->co[r][c]), &(top->Cdown->co[r][c]));	
					R=top->Rdown->co[r][c];
					C=top->Cdown->co[r][c];
					
					area = ds*ds/cos(top->slope->co[r][c]*Pi/180.);
					ch = cnet->ch->co[r][c];
					if (ch>0) area -= cnet->length->co[ch] * par->w_dx * ds; //area of the pixel[m2]
					
					Vmax = area * 1.E-3 * H;
					
					if(top->is_on_border->co[r][c]==1 && R==r && C==c && cnet->ch->co[r][c]==0){
						
						q = Cd*(2./3.)*sqrt(2.*g*1.E-3*H)*1.E-3*H*ds;//m3/s
						
						Vmax = Fmax(Vmax, 1.E-10);
						q = Fmax(q, 1.E-10);
						
						if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
						
					}else{
						
						AREA = ds*ds/cos(top->slope->co[R][C]*Pi/180.);
						CH = cnet->ch->co[R][C];
						if (CH>0) area -= cnet->length->co[CH] * par->w_dx * ds; //area of the pixel[m2]
						
						if( (R-r==1 || R-r==-1) && (C-c==1 || C-c==-1) ){
							dn = ds*( 0.5*sqrt(pow(1./cos(atan(top->dzdE->co[r][c])),2.)+pow(1./cos(atan(top->dzdN->co[R][C])),2.) ) +
									 0.5*sqrt(pow(1./cos(atan(top->dzdN->co[r][c])),2.)+pow(1./cos(atan(top->dzdE->co[R][C])),2.) ) );
							dD = find_3Ddistance(ds*sqrt(2.), top->Z0->co[r][c] - top->Z0->co[R][C]);
						}else if(R-r==1 || R-r==-1){
							dn = ds/cos(atan(top->dzdE->co[r][c]));
							dD = find_3Ddistance(ds, top->Z0->co[r][c] - top->Z0->co[R][C]);
						}else {
							dn = ds/cos(atan(top->dzdN->co[r][c]));
							dD = find_3Ddistance(ds, top->Z0->co[r][c] - top->Z0->co[R][C]);
						}
						
						lu = (short)land->LC->co[r][c];
						Ks = cm_h(land->ty->co[lu][jcm], H, par->thres_hsup_1,  par->thres_hsup_2);
						
						i = ( (top->Z0->co[r][c] - top->Z0->co[R][C] ) + 1.E-3*(Fmax(0.0, h[r][c]) - Fmax(0.0, h[R][C])) ) / dD;		
						if(i<0) i=0.;
						
						q = dn*Ks*pow(1.E-3*H,1.0+par->gamma_m)*sqrt(i);//m3/s
						
						if(i>1.E-3) Vmax = Fmin(Vmax, (1.-alpha)*i*dD / ( find_sup_pressure(1./area, top->slope->co[r][c]) + find_sup_pressure(1./AREA, top->slope->co[R][C])));
						
						Vmax = Fmax(Vmax, 1.E-10);
						q = Fmax(q, 1.E-10);
						
						if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
						
					}																
					
					if(*dt < par->dtmin_sup) *dt = par->dtmin_sup;
					
					
					if(ch > 0){	//channel pixel
						
						Hch = find_hsup(hch[ch], top->slope->co[r][c]) - par->depr_channel;
						areach = cnet->length->co[ch] * par->w_dx * ds;
						
						if( Hch < -min_dhsup_land_channel ){	//free flow
							
							DH = H;
							q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
							Vmax = Fmin( areach*1.E-3*(-Hch) , area*1.E-3*H );
							
							Vmax = Fmax(Vmax, 1.E-10);
							q = Fmax(q, 1.E-10);
							
							if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
							
						}else if( H - Hch > min_dhsup_land_channel ){//submerged flow towards channel
							
							DH = H - Hch;
							q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
							Vmax = 1.E-3*DH / (1./area + 1./areach);
							
							Vmax = Fmax(Vmax, 1.E-10);
							q = Fmax(q, 1.E-10);
							
							if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
							
						}
						
						if(*dt < par->dtmin_sup) *dt = par->dtmin_sup;
						
					}			
					
				}else{
					
					if(DDcomplex==1){
						top->Rdown->co[r][c] = r;
						top->Cdown->co[r][c] = c;
					}
					
				}
			}
			
		}
	}
	
	//submerged flow from the channel
	for(ch=1;ch<=par->total_channel;ch++){
		
		r = cnet->r->co[ch];
		c = cnet->c->co[ch];
		
		H = find_hsup(h[r][c], top->slope->co[r][c]);
		Hch = find_hsup(hch[ch], top->slope->co[r][c]) - par->depr_channel;
		
		area = ds*ds/cos(top->slope->co[r][c]*Pi/180.) - cnet->length->co[ch] * par->w_dx * ds;				
		areach = cnet->length->co[ch] * par->w_dx * ds;
		
		if ( Hch - H > min_dhsup_land_channel ) {
			
			DH = Hch - H;
			q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*Hch;//m3/s
			
			Vmax = 1.E-3*DH / (1./area + 1./areach);
			
			Vmax = Fmax(Vmax, 1.E-10);
			q = Fmax(q, 1.E-10);
			
			if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
			
			if(*dt < par->dtmin_sup) *dt = par->dtmin_sup;
			
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void supflow(double Dt, double t, short DDland, short DDchannel, double **h, double **dV, double *hch, double *dhch, 
			 TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, double *Voutnet, double *Voutland, FILE *flog)

{
	long r, c, R, C, ch;                                    
	double ds=sqrt(UV->U->co[1]*UV->U->co[2]), area, areach, Vmax, H, Hch, DH;							     
	double Ks;											  // the Strickler's coefficent calculated with a standard deviation
	double dn, dD;											
	double i;											  // hydraulic gradient
	double q, tb, te, dt;
	short lu;
	
	if(par->point_sim!=1){	//distributed simulations
		
		te=0.0;
				
		//fprintf(flog, "Surface flow: ");
		
		do{
			
			tb=te;
			dt=Dt;
			
			find_dt_max(DDland, par->max_courant_land, h, hch, land, top, cnet, par, t, &dt);
						
			te=tb+dt;
			if(te>Dt){
				te=Dt;
				dt=te-tb;
			}		
						
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						
						H = find_hsup(h[r][c], top->slope->co[r][c]);
						
						if(H > par->min_hsup_land){
							
							area = ds*ds/cos(top->slope->co[r][c]*Pi/180.);
							ch = cnet->ch->co[r][c];
							if (ch>0) area -= cnet->length->co[ch] * par->w_dx * ds; //area of the pixel[m2]
							
							lu=(short)land->LC->co[r][c];
							Ks=cm_h(land->ty->co[lu][jcm], H, par->thres_hsup_1,  par->thres_hsup_2);
							
							R = top->Rdown->co[r][c];
							C = top->Cdown->co[r][c];
							
							Vmax = area * 1.E-3 * H;
							
							if(top->is_on_border->co[r][c]==1 && R==r && C==c && cnet->ch->co[r][c]==0){
								
								q = Cd*(2./3.)*sqrt(2.*g*1.E-3*H)*(1.E-3*H)*ds;//m3/s
								if (q*dt > Vmax) q = Vmax/dt;
								*Voutland = *Voutland + q*dt;	//m3
								
							}else{
								
								if( (R-r==1 || R-r==-1) && (C-c==1 || C-c==-1) ){
									dn = ds*( 0.5*sqrt(pow(1./cos(atan(top->dzdE->co[r][c])),2.)+pow(1./cos(atan(top->dzdN->co[R][C])),2.) ) +
											  0.5*sqrt(pow(1./cos(atan(top->dzdN->co[r][c])),2.)+pow(1./cos(atan(top->dzdE->co[R][C])),2.) ) );
									dD = find_3Ddistance(ds*sqrt(2.), top->Z0->co[r][c] - top->Z0->co[R][C]);
								}else if(R-r==1 || R-r==-1){
									dn = ds/cos(atan(top->dzdE->co[r][c]));
									dD = find_3Ddistance(ds, top->Z0->co[r][c] - top->Z0->co[R][C]);
								}else {
									dn = ds/cos(atan(top->dzdN->co[r][c]));
									dD = find_3Ddistance(ds, top->Z0->co[r][c] - top->Z0->co[R][C]);
								}
								
								i = ( (top->Z0->co[r][c] - top->Z0->co[R][C] ) + 1.E-3*(Fmax(0.0, h[r][c]) - Fmax(0.0, h[R][C])) ) / dD;		
								if(i<0) i=0.;

								q = dn*Ks*pow(1.E-3*H,1.0+par->gamma_m)*sqrt(i);//m3/s
								if (q*dt > Vmax) q = Vmax/dt;

							}
							
						}else{
							
							q = 0.0;
							
						}
						
						dV[r][c] = q*dt;
						
					}
				}
			}
			
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						
						area = ds*ds/cos(top->slope->co[r][c]*Pi/180.);
						ch = cnet->ch->co[r][c];
						if (ch>0) area -= cnet->length->co[ch] * par->w_dx * ds; //area of the pixel[m2]
						
						h[r][c] -= find_sup_pressure(1.E3*dV[r][c]/area, top->slope->co[r][c]);
						
						R = top->Rdown->co[r][c];
						C = top->Cdown->co[r][c];
						
						area = ds*ds/cos(top->slope->co[R][C]*Pi/180.);
						ch = cnet->ch->co[R][C];
						if (ch>0) area -= cnet->length->co[ch] * par->w_dx * ds; //area of the pixel[m2]						
						
						if(r!=R || c!=C){
							if(h[R][C]>0){
								h[R][C] += find_sup_pressure(1.E3*dV[r][c]/area, top->slope->co[R][C]);
							}else{
								if( dV[r][c] > 0) h[R][C] = find_sup_pressure(1.E3*dV[r][c]/area, top->slope->co[R][C]);
							}
						}
					}
				}
			}
			
			
			//Superficial flow to the channels
			for(ch=1;ch<=par->total_channel;ch++){
				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				
				H = find_hsup(h[r][c], top->slope->co[r][c]);
				Hch = find_hsup(hch[ch], top->slope->co[r][c]);

				area = ds*ds/cos(top->slope->co[r][c]*Pi/180.) - cnet->length->co[ch] * par->w_dx * ds;
				areach = cnet->length->co[ch] * par->w_dx * ds;
				
				if( Hch <= par->depr_channel ){	//free flow
					
					DH = H;
					q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
					
					Vmax = Fmin( areach*1.E-3*(par->depr_channel - Hch) , area*1.E-3*H );
					if (q*dt > Vmax) q = Vmax/dt;

				}else if( H > Hch - par->depr_channel ){//submerged flow towards channel

					DH = H - ( Hch-par->depr_channel );
					q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
					
					Vmax = 1.E-3*DH / (1./area + 1./areach);
					if (q*dt > Vmax) q = Vmax/dt;

				}else {	
					q = 0.;
				}
								
				cnet->Vsup->co[ch] += q*dt;
				
				h[r][c] -= find_sup_pressure(1.E3 * q*dt/area, top->slope->co[r][c]);
								
				if(hch[ch]>0){
					hch[ch] += find_sup_pressure(1.E3 * q*dt/areach, top->slope->co[r][c]);	//mm;
				}else{
					if( q > 0 ) hch[ch] = find_sup_pressure(1.E3 * q*dt/areach, top->slope->co[r][c]);	//mm;
				}
								
			}
			
			//Superficial flow from the channels
			for(ch=1;ch<=par->total_channel;ch++){
				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];

				H = find_hsup(h[r][c], top->slope->co[r][c]);
				Hch = find_hsup(hch[ch], top->slope->co[r][c]);

				area = ds*ds/cos(top->slope->co[r][c]*Pi/180.) - cnet->length->co[ch] * par->w_dx * ds;				
				areach = cnet->length->co[ch] * par->w_dx * ds;
				
				if ( Hch - par->depr_channel > H ) {
					DH = ( Hch-par->depr_channel ) - H;
					q = Cd*(2./3.)*sqrt(2.*g*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*(Hch - par->depr_channel);//m3/s
					
					cnet->Vsup->co[ch] -= q*dt;

					Vmax = 1.E-3*DH / (1./area + 1./areach);
					if (q*dt > Vmax) q = Vmax/dt;
					
					hch[ch] -= find_sup_pressure(1.E3 * q*dt/areach, top->slope->co[r][c]);

					if(h[r][c]>0){
						h[r][c] += find_sup_pressure(1.E3 * q*dt/area, top->slope->co[r][c]);	//mm;
					}else{
						if( q > 0 ) h[r][c] = find_sup_pressure(1.E3 * q*dt/area, top->slope->co[r][c]);	//mm;
					}
				}
			}
							
			channel_flow(dt, t, DDchannel, hch, dhch, top, cnet, par, land, Voutnet, flog);
						
		}while(te<Dt);
		
		//fprintf(flog, "\n\n");
		
	}else{	//point simulation  
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if( (long)land->LC->co[r][c]!=number_novalue){
					h[r][c]=Fmin(h[r][c], 0.0);
				}
			}
		}
	}
		
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
				
		if(DDcomplex!=1 && t==0) draining_channel(0., ch, top->Z0, h, cnet, &(cnet->ch_down->co[ch]));
		
		r = cnet->r->co[ch];
		c = cnet->c->co[ch];
		H = find_hsup(h[ch], top->slope->co[r][c]);
		
		if(H > par->min_hsup_channel){
			
			if(DDcomplex==1) draining_channel(1., ch, top->Z0, h, cnet, &(cnet->ch_down->co[ch]));
			
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
				
				Ks=cm_h(par->Ks_channel, H, 1., par->thres_hchannel);
				
				i = ( (top->Z0->co[r][c] - top->Z0->co[R][C] ) + 1.E-3*(Fmax(0.0, h[ch]) - Fmax(0.0, h[cnet->ch_down->co[ch]])) ) / dD;	
				
				if(i<0) i=0.;
				
				q = dn * Ks * pow( 1.E-3 * H , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
				
			}
						
			Vmax = Fmax(Vmax, 1.E-10);
			q = Fmax(q, 1.E-10);
			if(Courant*Vmax/q<(*dt)) *dt=Courant*Vmax/q; 
			if(*dt<par->dtmin_sup) *dt=par->dtmin_sup;
			
		}else{
			
			if(DDcomplex==1) cnet->ch_down->co[ch] = ch;
			
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void channel_flow(double Dt, double t, short DDcomplex, double *h, double *dV, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double *Vout, FILE *f)

{
	long r,c,ch,R,C;                                    
	double ds, dn, dD;
	double Ks;											  // the Strickler's coefficent
	double i;											  // hydraulic gradient
	double q,tb,te,dt,H;
	
	ds = sqrt(UV->U->co[1]*UV->U->co[2]);
	dn = par->w_dx*ds;
	
	if( par->point_sim!=1 && cnet->r->co[1]!=0 ){	//if it is not point simulation and there are channels
		
		dn = par->w_dx * UV->U->co[1];		//transversal length [m]
		
		te=0.0;
		
		do{
			
			tb=te;
			dt=Dt;
			
			find_dt_max_channel(DDcomplex, par->max_courant_channel, h, top, cnet, par, land, t, &dt);
						
			te=tb+dt;
			if(te>Dt){
				te=Dt;
				dt=te-tb;
			}		
			
			//fprintf(f,"%f/%f ",dt,Dt);						
			
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];

				H = find_hsup(h[ch], top->slope->co[r][c]);

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
						
						Ks=cm_h(par->Ks_channel, H, 1., par->thres_hchannel);
						
						i= ( (top->Z0->co[r][c] - top->Z0->co[R][C] ) + 1.E-3*(Fmax(0.0, h[ch]) - Fmax(0.0, h[cnet->ch_down->co[ch]])) ) / dD;		
						
						if(i<0) i=0.;
						
						q = dn * Ks * pow( 1.E-3 * H , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
						
					}
					
					dV[ch] = Fmin( q*dt , 1.E-3*H*dn*cnet->length->co[ch] );	//m3
					
				}else{
					
					dV[ch] = 0.0;
					
				}
				
			}
			
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
								
				h[ch] -= find_sup_pressure(1.E3*dV[ch]/(dn*cnet->length->co[ch]), top->slope->co[r][c]);
				
				if(top->is_on_border->co[r][c] == 1 && cnet->ch_down->co[ch]==ch){//outlet section
					*Vout = *Vout + dV[ch];	//m3
				}else {
					R = cnet->r->co[cnet->ch_down->co[ch]];
					C = cnet->c->co[cnet->ch_down->co[ch]];					
					if(h[cnet->ch_down->co[ch]]>0){
						h[cnet->ch_down->co[ch]] += find_sup_pressure(1.E3*dV[ch]/(dn*cnet->length->co[cnet->ch_down->co[ch]]), top->slope->co[R][C]);	//mm;
					}else{
						if( dV[ch] > 0) h[cnet->ch_down->co[ch]] = find_sup_pressure(1.E3*dV[ch]/(dn*cnet->length->co[cnet->ch_down->co[ch]]), top->slope->co[R][C]);	//mm;
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

void draining_land(double alpha, long r, long c, DOUBLEMATRIX *Z, double **h, DOUBLEMATRIX *LC, long *R, long *C){
	
	long d;
	double elev, elev1;
	long ir[9] = {0, -1, -1, 0, 1, 1, 1, 0, -1};
	long ic[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
	
	*R = r;
	*C = c;
	
	elev = Z->co[r][c] + alpha*1.E-3*Fmax(h[r][c], 0.);
	
	for (d=1; d<=8; d++) {
		if (r+ir[d]>=1 && r+ir[d]<=Nr && c+ic[d]>=1 && c+ic[d]<=Nc) {
			if ((long)LC->co[r+ir[d]][c+ic[d]] != number_novalue) {
				elev1 = Z->co[r+ir[d]][c+ic[d]] + alpha*1.E-3*Fmax(h[r+ir[d]][c+ic[d]], 0.);
				if( elev1 < elev){
					elev = elev1;
					*R = r+ir[d];
					*C = c+ic[d];
				}
			}
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

