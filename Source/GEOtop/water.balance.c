
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */
  
    
#include "constant.h"
#include "keywords_file.h"
#include "struct.geotop.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "util_math.h"
#include "water.balance.h"

extern T_INIT *UV;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;
extern DOUBLEVECTOR *outdata_basin;

#define tol_max_GC 1.E+5
#define tol_min_GC 1.E-13
#define max_res_adm 1.E-2
#define min_lambda 1.E-6
#define max_Courant_sup 0.25
#define max_Courant_channel 0.25
#define M 1
#define ni 1.E-4

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void water_balance(ALLDATA *adt){

	DOUBLETENSOR *P;
	
	short out, sux;
	long r, c, l, iter, sy;
	double n, Dt0, Dt, tb, te, te0, Loss, Vout=0.;
	
	FILE *f;		

	P=new_doubletensor0(Nl, Nr, Nc);

	te0=0.0;
	
	do{

		Dt0=adt->P->Dt;
		if(te0+Dt0 > adt->P->Dt) Dt0=adt->P->Dt-te0;
		//find_dt_max(adt->P->max_Courant_sup_sub, adt->S->P->co[0], adt->L, adt->T, adt->C, adt->P, &Dt0);								
		te0 += Dt0;	
			
		n=1.;
		te=0;
	
		do{
	
			tb=te;				
		
			do{
		
				Dt=Dt0/n;

				f=fopen(files->co[ferr]+1,"a");
				fprintf(f,"\nn:%f Dt:%f Dt0:%f te0:%f tb:%f \n",n,Dt,Dt0,te0,tb);	
				fclose(f);
				
				if(Dt<adt->P->Dt) printf("Time step:%f/%f\n",Dt,te0);
				
				sux = Richards(Dt, P, &Loss, &iter, adt);
								
				out=1;
								
				if(sux==0 && Dt>adt->P->DtminWb){
				
					n*=adt->P->nredDtWb;
					out=0;
					
					f=fopen(files->co[ferr]+1,"a");
					fprintf(f,"Reduced time step: dt:%f Loss:%f \n\n",Dt0/n,Loss);
					fclose(f);					
				
				}else{
					
					f=fopen(files->co[ferr]+1,"a");
					fprintf(f,"Time step ok: dt:%f Loss:%f \n\n",Dt0/n,Loss);
					fclose(f);	
				
				}
									
			}while( out==0 ); 
			
			if(sux==0){
				printf("ERROR:Water balance does not converge\n");
				printf("It is not possible to continue, Please check the parameters in the block 2 of the parameter file\n");
				printf("or reduce the time step or review the initial conditions\n\n");
				printf("If you think that everything is right, please report this error to the GEOtop community\n\n");
				t_error("Code execution interrupted");
			}
			
			te=tb+Dt;                                 

			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(adt->L->LC->co[r][c]!=NoV){
						for(l=0;l<=Nl;l++){
							adt->S->P->co[l][r][c] = P->co[l][r][c];
						}
					}
				}
			}
			
			supflow(Dt, adt->S->P->co[0], adt->W->dh_sup->co, adt->T, adt->L, adt->W, adt->C, adt->P, &Vout);
						
			outdata_basin->co[oomasserror] += fabs(Loss);
		
		}while(te<Dt0);
		
	}while(te0<adt->P->Dt);
	
	adt->C->Q_out = Vout/adt->P->Dt;
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(adt->L->LC->co[r][c]!=NoV){
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
	
	free_doubletensor(P);
				
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
 
 1) is SPD (symmetrical positive defined)
 
 2) if kij is a component kii = sum(j!=i) kji
 
 3) the sum of all the components of the matrix is 0
 
 4) if T is a vector, the vector K T has the properties that the sum of its components is 0
 
 Here it is possible to update or not K with the new values of H, In any case the derivative of K(H) is not considered in the Jacobian. So the matrix K(H) is also
 
 used in the JACOBIAN.
 
 K is described storing only its strict lower component (strict = without diagonal) with the 3 vectors Li, Lp, Lx (in the same way as UFMPACK) */
 
 
short Richards(double Dt, DOUBLETENSOR *P, double *loss, long *iter, ALLDATA *adt){
	
	//DOUBLEVECTOR *diag, *udiag;

	double res=0.0, res0[3], res_prev[M], res_av, res00, lambda[3], epsilon, mu;

	long i, l, r, c, m, cont, cont2, iter_tot=0, n=(Nl+1)*adt->P->total_pixel;	
	static long cum_iter;	
	short out, out2;	
	int sux;
		
	FILE *f;
			
	if(adt->I->time == 0.0) cum_iter = 0;
		
	//diag=new_doublevector(n);
	//udiag=new_doublevector(n-1);
	
	*iter = 0;	//conjugated gradient iteration number
	
		
	for(i=1;i<=n;i++){
		
		l = adt->T->lrc_cont->co[i][1];
		r = adt->T->lrc_cont->co[i][2];
		c = adt->T->lrc_cont->co[i][3];
		
		/*
		Layer 0 is water on the surface. This allows a more robust description of infiltration processes.
		However, surface water flow is described separately in the subroutine supflow
		
		Layers > 0 , represent subsurface
		*/
		
		adt->W->H0->co[i] = adt->S->P->co[l][r][c] + adt->T->Z->co[l][r][c];		
		
		adt->W->H1->co[i] = adt->W->H0->co[i];
				
		if(adt->W->H1->co[i] != adt->W->H1->co[i]) printf("no value in r:%ld c:%ld l:%ld\n",r,c,l);
		
	}
	
	sux = find_matrix_K(adt->W->Lx, adt, adt->W->H1, Dt);
	
	
	find_f(adt->W->f, adt, adt->W->H1, Dt);
	product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
	res = norm_2(adt->W->B, n);
	
	res00 = res; //initial norm of the residual
	epsilon = adt->P->TolVWb + adt->P->RelTolVWb * Fmin( res00 , sqrt((double)n) );

	if(res!=res) printf("res in no value\n");

	cont=0;

	do{
		
		cont++;
		
		for(i=1;i<=n;i++){	
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
		//sux = J_water(adt->W->Lx, adt, adt->W->H1, Dt);
		//sux = J_water_with_transposed(adt->W->Lx, adt->W->Ux, adt, adt->W->H1, Dt);
		sux = find_dfdH(adt->W->df, adt, adt->W->H1, Dt);	//it calcolates only df/dH, J = I*df/dH + K, K is calculated above
		
		//NEGLECTING THE LATERAL FLUXES (neglecting the extra-tridiagonal part)
		/*get_diag_lower_matrix(diag, udiag, adt->T->Ai, adt->T->Ap, adt->W->Ax);
		 tridiag(0, 0, 0, n, udiag, diag, udiag, adt->W->B, adt->W->dH);*/
		
		//VARIOUS CONJUGATED GRADIENTS ALGORITHM
		//*iter = BiCGSTAB( mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->H0, Dt, adt->W->B, J_d_water_, adt);	
		//*iter = BiCGSTAB_lower(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->T->Li, adt->T->Lp, adt->W->Lx);
		//*iter = BiCGSTAB_upper(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->T->Ui, adt->T->Up, adt->W->Ux);
		//*iter = BiCGSTAB_LU_SSOR(1.0, mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->T->Li, adt->T->Lp, adt->W->Lx, adt->T->Ui, adt->T->Up, adt->W->Ux);
		//*iter = BiCGSTAB_diag(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->T->Li, adt->T->Lp, adt->W->Lx);
		//*iter = BiCGSTAB_unpreconditioned(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->T->Li, adt->T->Lp, adt->W->Lx);
		*iter = BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(mu, tol_min_GC, tol_max_GC, adt->W->dH, adt->W->B, adt->W->df, adt->T->Li, adt->T->Lp, adt->W->Lx);
		iter_tot += (*iter);//The number of the cumulated GC iterations is stored
		
		//using the UFMPACK
		/*sux = J_water_triplet(adt->W->Jtriplet, adt, adt->W->H1, Dt);
		sux = update_umfpack_matrix_from_UMFPACK_REAL_TRIPLET(adt->W->Jtriplet, adt->W->Jmatrix);
		sux = solve_linear_system_with_UMFPACK_REAL_MATRIX(adt->W->Jmatrix, adt->W->dH, adt->W->B);*/
				
		//non-monotonic line search (it is monotonic if M==1)
		for(m=Fminlong(cont,M);m>1;m--){
			res_prev[m-1]=res_prev[m-2];
		}
		res_prev[0]=res;
		
		res_av=0.0;
		for(m=1;m<=Fminlong(cont,M);m++){
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
			
			for(i=1;i<=n;i++){
				
				adt->W->H1->co[i] = adt->W->H0->co[i] + lambda[0] * adt->W->dH->co[i];
				
				if(adt->W->H1->co[i] != adt->W->H1->co[i]) {
					printf("no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);
					stop_execution();
				}
				
			}
			
			if(adt->P->UpdateK == 1) sux = find_matrix_K(adt->W->Lx, adt, adt->W->H1, Dt);		
			//sux = find_matrix_K(adt->W->Lx, adt, adt->W->H1, Dt);
									
			find_f(adt->W->f, adt, adt->W->H1, Dt);
			product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
			res = norm_2(adt->W->B, n);		
			
			out2=0;
			
			if(res <= (1.0 - ni*lambda[0]*(1.-mu))*res_av) out2=1;
			if(lambda[0] <= min_lambda) out2=1;
			if(cont==1) out2=1;
					
		}while(out2==0);	
		
		out=0;
		//The condition to go out from the Newton cycle is made on the norm of the residual < Absolute tolerance + Relative tolerance * res00
		if( res <= Fmin( epsilon , max_res_adm ) ) out=1;
		//Max iteration number
		if( cont >= adt->P->MaxiterTol ) out=1;	
		
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"iter:%ld/%ld cont:%ld cont2:%ld res:%e res_av:%e epsilon:%e mu:%e lambda:%e\n",*iter,iter_tot,cont,cont2,res,res_av,epsilon,mu,lambda[0]);
		fclose(f);

	}while(out==0);
					
	cum_iter += iter_tot;
	
	//it can be shown that massloss per unit pixel [mm] is the linear norm of B * Dt / total_pixel
	*loss = sum(adt->W->B, n)*Dt/adt->P->total_pixel;
	f=fopen(files->co[ferr]+1,"a");
	fprintf(f,"cumiter:%ld massloss:%e\n",cum_iter,*loss);
	fclose(f);
	
	//assign updated state variables
	for(i=1;i<=n;i++){
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		P->co[l][r][c] = adt->W->H1->co[i] - adt->T->Z->co[l][r][c];
	}	
	
	//calculate Qsub
	for(i=1;i<=adt->C->r->nh;i++){
		r = adt->C->r->co[i];
		c = adt->C->c->co[i];
		//there are channel pixels, if there are not channels r vector has 1 component and is equal to 0
		if(r>0) adt->C->Qsub->co[i] = adt->P->Kch_b * (adt->P->w_dx*UV->U->co[1]*adt->C->length->co[i]) * 1.E-3*Fmax( P->co[0][r][c] - adt->C->h_sup->co[i], 0.0 );//m3/s
	}
			
	//free_doublevector(diag);
	//free_doublevector(udiag);
	
	//end subroutine
	if( res <= epsilon ){
		return 1;//converges
	}else{
		return 0;//does not converge
	}
		
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double J_d_water(long i, DOUBLEVECTOR *d, DOUBLEVECTOR *H, double Dt, ALLDATA *adt){

	/*Calculates the i-th component of the Jacobian(H) multiplied by d (=Newton direction vector)
	
	F(H) = (V(H) - V0)/Dt - Ah*(d/dz * ( K(H) * dH/dz ) ) - Av*(d/dx * ( K(H) * dH/dx ) ) - Av*(d/dy * ( K(H) * dH/dy ) ) + Ah*dz*sink(H)
	
	   Ah = horizontal area (dx*dy)
	   Av = vertical area (dz*dx = dz*dy)
	 
	   J(H) = dF/dH
	 
	   K(H) is not considered in the J 
	 
	   F and J are here considered divided by Ah */
	

	long l, r, c, I, R, C;
	long sy;
	double a = 0.0;
	double dz=0.0, k=0.0, kn, dzn, ds, h, klim;
	
	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];
	
	//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
	if(l==0){
		h = H->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
		if(h>0) a += d->co[i]*(1./Dt);
	}else{
		dz = adt->S->pa->co[sy][jdz][l];
		a += d->co[i]*(dtheta_dpsi_from_psi(H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin)*dz/Dt);
	}
	
	//vertical hydraulic conductivity
	if(l>0){
		if(adt->P->UpdateK==1){
			k = k_from_psi(jKv, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
		}else{
			k = k_from_psi(jKv, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
		}
	}		
	
	//Vertical fluxes (diagonal and tridiagonal terms)
	if(l<Nl){
		I = adt->T->i_cont[l+1][r][c];	
		
		if(l==0){	//overland flow
			
			if( H->co[i]< H->co[I] ){	//upward flux
				kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
			}else{	//downward flow
				kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
			}
			
			dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];
			
		}else{	//subsurface flow
			
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKv,  adt->S->P->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
			}
			
			dzn = adt->S->pa->co[sy][jdz][l+1];
			kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			dzn = 0.5*dz + 0.5*dzn;
			
		}
		
		a += (kn/dzn)*(d->co[i]-d->co[I]);
		
	}
	
	if(l>0){
		I=adt->T->i_cont[l-1][r][c];
		if(l==1){	
			
			if( H->co[I] < H->co[i] ){	//upward flux
				kn = k_from_psi(jKv,  H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{	//downward flow
				kn = k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
			}
			
			dzn = 0.5*adt->S->pa->co[sy][jdz][l];
			
		}else{
			
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKv,  adt->S->P->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
			}

			dzn = adt->S->pa->co[sy][jdz][l-1];
			kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l-1, r, c, adt->S), l-1, r, c, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			dzn = 0.5*dz + 0.5*dzn;
			
		}
		
		a+=(kn/dzn)*(d->co[i]-d->co[I]);
	}
	
	//lateral hydraulic conductivity
	if(l>0){
		if(adt->P->UpdateK==1){
			k = k_from_psi(jKh, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
		}else{
			k = k_from_psi(jKh, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
		}
	}		
	
	//lateral fluxes
	//4 neighbouring cells
	R = r+1;
	C = c;
	ds = 1.E3*UV->U->co[2];	//[mm]
	//1.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];				
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}
			kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			a += kn*((d->co[i]-d->co[I])/ds)*(dz/ds);
		}
	}
	
	R = r-1;
	C = c;
	ds = 1.E3*UV->U->co[2];	//[mm]
	//2.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];				
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}
			kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			a += kn*((d->co[i]-d->co[I])/ds)*(dz/ds);
		}
	}
	
	
	R = r;
	C = c+1;
	ds = 1.E3*UV->U->co[1];	//[mm]
	//3.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];				
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}
			kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			a += kn*((d->co[i]-d->co[I])/ds)*(dz/ds);
		}
	}	
	
	
	R = r;
	C = c-1;
	ds = 1.E3*UV->U->co[1];	//[mm]
	//4.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];				
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}
			kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			a += kn*((d->co[i]-d->co[I])/ds)*(dz/ds);
		}
	}
	
	
	//subsuperficial flow to the channels (sink depending on H)
	if(adt->T->pixel_type->co[r][c]>=10 && l==0){ 
		h = H->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
		if( h > 0 ) a += adt->P->Kch_b * adt->P->w_dx * d->co[i];
	}

	return(a);
	
}

double J_d_water_(long i, DOUBLEVECTOR *d, DOUBLEVECTOR *H, double Dt, void *adt) { 
	return J_d_water(i, d, H, Dt, (ALLDATA *)adt); 
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double F_water(long i, DOUBLEVECTOR *H, double Dt, ALLDATA *adt){
	
	long l, r, c, I, R, C;
	short sy;
	double a = 0.0;
	double dz=0.0, k=0.0, kn, dzn, ds, h, klim, V0, V1;
	
	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];
		
	//hydraulic capacity (diagonal term)
	if(l==0){
		V1 = Fmax(0.0, H->co[i] - adt->T->Z->co[l][r][c]);
		V0 = Fmax(0.0, adt->S->P->co[l][r][c]);
	}else{
		dz = adt->S->pa->co[sy][jdz][l];		
		V1 = dz * theta_from_psi( H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin );
		V0 = dz * theta_from_psi( adt->S->P->co[l][r][c], l, r, c, adt->S, PsiMin );
	}
	a += (V1-V0)/Dt;
			
	//vertical hydraulic conductivity
	if(l>0){
		if(adt->P->UpdateK==1){
			k = k_from_psi(jKv, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
		}else{
			k = k_from_psi(jKv, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
		}
	}		
	
	//Vertical fluxes (diagonal and tridiagonal terms)
	if(l<Nl){
		I = adt->T->i_cont[l+1][r][c];	
		
		if(l==0){	//overland flow
			
			if( H->co[i] < H->co[I] ){	//upward flux
				kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
			}else{	//downward flow
				kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
			}
			dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];
			
		}else{	//subsurface flow
			
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKv,  adt->S->P->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
			}
			
			dzn = adt->S->pa->co[sy][jdz][l+1];
			kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			dzn = 0.5*dz + 0.5*dzn;
			
		}
		
		a += (kn/dzn)*(H->co[i]-H->co[I]);

	}
	
	
	
	if(l>0){
		I=adt->T->i_cont[l-1][r][c];
		if(l==1){	
			
			if( H->co[I] < H->co[i] ){	//upward flux
				kn = k_from_psi(jKv,  H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{	//downward flow
				kn = k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
			}
			dzn = 0.5*adt->S->pa->co[sy][jdz][l];
			
		}else{
			
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKv,  adt->S->P->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
			}
			
			dzn = adt->S->pa->co[sy][jdz][l-1];
			kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l-1, r, c, adt->S), l-1, r, c, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			dzn = 0.5*dz + 0.5*dzn;
			
		}
		
		a += (kn/dzn)*(H->co[i]-H->co[I]);
	}

	

	//lateral hydraulic conductivity
	if(l>0){
		if(adt->P->UpdateK==1){
			k = k_from_psi(jKh, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
		}else{
			k = k_from_psi(jKh, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
		}
	}		
	
	//lateral fluxes
	//4 neighbouring cells
	R = r+1;
	C = c;
	ds = 1.E3*UV->U->co[2];	//[mm]
	//1.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];				
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}
			kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			a += kn*((H->co[i]-H->co[I])/ds)*(dz/ds);
		}
	}
	
	R = r-1;
	C = c;
	ds = 1.E3*UV->U->co[2];	//[mm]
	//2.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];				
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}
			kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			a += kn*((H->co[i]-H->co[I])/ds)*(dz/ds);
			if( kn*((H->co[i]-H->co[I])/ds)*(dz/ds) != kn*((H->co[i]-H->co[I])/ds)*(dz/ds) ) printf("5\n");
		}
	}
	
	
	R = r;
	C = c+1;
	ds = 1.E3*UV->U->co[1];	//[mm]
	//3.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];				
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}
			kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			a += kn*((H->co[i]-H->co[I])/ds)*(dz/ds);
			if( kn*((H->co[i]-H->co[I])/ds)*(dz/ds) != kn*((H->co[i]-H->co[I])/ds)*(dz/ds) ) printf("6\n");
		}
	}	
	
	
	R = r;
	C = c-1;
	ds = 1.E3*UV->U->co[1];	//[mm]
	//4.
	if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
		if(adt->L->LC->co[R][C]!=NoV){
			I=adt->T->i_cont[l][R][C];				
			if(adt->P->UpdateK==1){
				kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}else{
				kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
			}
			kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
			klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
								 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
			if(kn > klim) kn=klim;
			a += kn*((H->co[i]-H->co[I])/ds)*(dz/ds);
			if( kn*((H->co[i]-H->co[I])/ds)*(dz/ds) != kn*((H->co[i]-H->co[I])/ds)*(dz/ds) ) printf("7\n");
		}
	}
	
	
	//subsuperficial flow to the channels
	if(adt->T->pixel_type->co[r][c]>=10 && l==0){ 
		h = H->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
		if( h > 0 ) a += adt->P->Kch_b * adt->P->w_dx * h;
	}
	
	
	//evaporation and precipitation
	if(l>0){
		a += adt->S->ET->co[l][r][c];
	}else{
		a -= adt->W->Pnet->co[r][c];		
	}
		
	return(a);
	
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double cm_h(double cm0, double h, double h_thres){
	
	double cm;
	
	if(h > h_thres){
		cm = cm0;
	}else if(h > 0){
		cm = ( 0.0 * (h_thres - h) + cm0*(h) ) / h_thres;
	}else{
		cm = 0.0;
	}
	
	return(cm);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void supflow(double Dt, double **h, double **dh, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, double *Vout)

{
	long r,c,R,C,ch;                                    
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	double dx,dy;                                         // the two dimensions of a pixel
	double Ks;											  // the Strickler's coefficent calculated with a standard deviation
	double b[10];                                         // area perpendicular to the superficial flow divided by hsup
	double i;											  // hydraulic gradient
	double q,tb,te,dt;
	short lu;
	
	if(par->point_sim!=1){	//distributed simulations
		
		dx=UV->U->co[1];                                    
		dy=UV->U->co[2];                                    
		
		b[1]=0.0;  
		b[2]=dy;             
		b[3]=sqrt(dx*dx+dy*dy); 
		b[4]=dx;            
		b[5]=sqrt(dx*dx+dy*dy); 
		b[6]=dy;   
		b[7]=sqrt(dx*dx+dy*dy); 
		b[8]=dx;             
		b[9]=sqrt(dx*dx+dy*dy); 
				
		te=0.0;
		
		do{
			
			tb=te;
			dt=Dt;
			
			find_dt_max(max_Courant_sup, h, land, top, cnet, par, &dt);
			
			te=tb+dt;
			if(te>Dt){
				te=Dt;
				dt=te-tb;
			}			
						
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(land->LC->co[r][c]!=NoV){
						
						if(h[r][c] > 0){
							lu=(short)land->LC->co[r][c];
							Ks=cm_h(land->ty->co[lu][jcm], h[r][c], par->thres_hsup);
							
							R=r+r_DD[top->DD->co[r][c]];
							C=c+c_DD[top->DD->co[r][c]];
							if(land->LC->co[R][C]!=NoV){
							
								i=0.001*( h[r][c] - h[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];				
								if(i<0) i=0.0;							
								q = b[top->DD->co[r][c]+1]*Ks*pow(h[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;
								
							}else{//pixel draining out of the basin
								
								i=top->i_DD->co[r][c];				
								if(i<0) i=0.0;							
								q = b[top->DD->co[r][c]+1]*Ks*pow(h[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;
								
								*Vout = *Vout + q*dx*dy*dt/1000.;	//m3
								
							}
								
						}else{
							q = 0.0;
						}
						
						dh[r][c] = Fmin( q*dt , Fmax( h[r][c] , 0.0 ) );
												
					}
				}
			}
			
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(land->LC->co[r][c]!=NoV){
						
						h[r][c] -= dh[r][c];
						
						R=r+r_DD[top->DD->co[r][c]];
						C=c+c_DD[top->DD->co[r][c]];
						
						if(land->LC->co[R][C]!=NoV){
			
							if(h[R][C]>0){
								h[R][C] += dh[r][c];
							}else{
								if( dh[r][c] > 0) h[R][C] = dh[r][c];
							}
						
						}
					}
				}
			}
			
			
			//Superficial flow to the channels
			for(ch=1;ch<=cnet->r->nh;ch++){
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				
				if(r>0){
				
					if(cnet->h_sup->co[ch] <= par->depr_channel){	//free flow
									
						if(h[r][c] > 0){
							q = Cd*(2./3.)*sqrt(2.*g*1.E-3*h[r][c])*((2.*cnet->length->co[ch])/(dx*dy))*h[r][c];//mm/s
							if(q*dt > h[r][c]) q = h[r][c]/dt;
						}else{
							q = 0.0;
						}
						
					}else if(h[r][c] > cnet->h_sup->co[ch] + par->depr_channel){//submerged flow towards channel
					
						q = Cd*(2./3.)*sqrt(2.*g*1.E-3*(h[r][c]-(cnet->h_sup->co[ch] + par->depr_channel)))*((2.*cnet->length->co[ch])/(dx*dy))*h[r][c];//mm/s
						if(q*dt > h[r][c]) q = h[r][c]/dt;					
					
					}else{	//submerged flow towards land
						
						q = 0.0;
						
					}
				
					h[r][c] -= q*dt;
					
					cnet->Qsup->co[ch] = 1.E-3*q*dx*dy;	//m3/s
					
				}
				
			}
						
			channel_flow(dt, cnet->h_sup, cnet->dh_sup, top, cnet, par, Vout);
			
			
		}while(te<Dt);
		
		
	}else{	//point simulation  
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(land->LC->co[r][c]!=NoV){
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

void find_dt_max(double Courant, double **h, LAND *land, TOPO *top, CHANNEL *cnet, PAR *par, double *dt){
	
	double b[10];    //area perpendicular to the superficial flow divided by hsup
	double i, q, dx, dy, Ks;											  
	short lu;
	long r, c, R, C, ch;
	
	short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	
	dx=UV->U->co[1];                                    
	dy=UV->U->co[2];                                    
	
	b[1]=0.0;  
	b[2]=dy;             
	b[3]=sqrt(dx*dx+dy*dy); 
	b[4]=dx;            
	b[5]=sqrt(dx*dx+dy*dy); 
	b[6]=dy;   
	b[7]=sqrt(dx*dx+dy*dy); 
	b[8]=dx;             
	b[9]=sqrt(dx*dx+dy*dy); 
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(h[r][c]>0 && land->LC->co[r][c]!=NoV){
				
				lu=(short)land->LC->co[r][c];
				Ks=cm_h(land->ty->co[lu][jcm], h[r][c], par->thres_hsup);
				
				R=r+r_DD[top->DD->co[r][c]];
				C=c+c_DD[top->DD->co[r][c]];
				
				if(land->LC->co[R][C]!=NoV){
					i=0.001*( h[r][c] - h[R][C] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];				
				}else{
					i=top->i_DD->co[r][c];
				}
				if(i<0) i=0;  //correct false slopes
				
				q=b[top->DD->co[r][c]+1]*Ks*pow(h[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;	//mm/s
				if(q>0){
					if(Courant*h[r][c]/q<(*dt)) *dt=Courant*h[r][c]/q; 
				}
				
				if(top->pixel_type->co[r][c]>=10){	//channel pixel
					
					ch = cnet->ch->co[r][c];
					q = 0.;
					
					if(cnet->h_sup->co[ch] <= par->depr_channel){	//free flow
						
						q = Cd*(2./3.)*sqrt(2.*g*1.E-3*h[r][c])*((2.*cnet->length->co[ch])/(dx*dy))*h[r][c];//mm/s
						
					}else if(h[r][c] > cnet->h_sup->co[ch] + par->depr_channel){//submerged flow towards channel
						
						q = Cd*(2./3.)*sqrt(2.*g*1.E-3*(h[r][c]-(cnet->h_sup->co[ch] + par->depr_channel)))*((2.*cnet->length->co[ch])/(dx*dy))*h[r][c];//mm/s
					
					}
						
					if(q>0){
						if(Courant*h[r][c]/q<(*dt)) *dt=Courant*h[r][c]/q; 
					}
				}
				
					
			}
		}
	}
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/*void routing(CHANNEL *cnet){
	
	long s, ch;
	
	initialize_doublevector(cnet->Qsup_spread, 0.0);
	initialize_doublevector(cnet->Qsub_spread, 0.0);
	
	for(ch=1;ch<=cnet->r->nh;ch++){
		for(s=1;s<=cnet->Qsup_spread->nh;s++){
			cnet->Qsup_spread->co[s] += (cnet->Qsup->co[ch])*(cnet->fraction_spread->co[ch][s]);//m3/s
			cnet->Qsub_spread->co[s] += (cnet->Qsub->co[ch])*(cnet->fraction_spread->co[ch][s]);//m3/s
		}
	} 
	
	for(s=1;s<=cnet->Qsup_spread->nh-1;s++){
		cnet->Q_sub_s->co[s]=cnet->Q_sub_s->co[s+1]+cnet->Qsub_spread->co[s];
		cnet->Q_sup_s->co[s]=cnet->Q_sup_s->co[s+1]+cnet->Qsup_spread->co[s];
	}
	
	cnet->Q_sub_s->co[cnet->Qsub_spread->nh]=cnet->Qsub_spread->co[cnet->Qsub_spread->nh];
	cnet->Q_sup_s->co[cnet->Qsup_spread->nh]=cnet->Qsup_spread->co[cnet->Qsup_spread->nh];
}*/


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//for UMFPACK triplet
/*int J_water_triplet(UMFPACK_REAL_TRIPLET *triplet, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i,l,r,c,I,R,C,sy,cnt=0,cntdiag;
	double h=0.0, dz=0.0, dzn=0.0, k=0.0, kn=0.0, klim=0.0, ds;
			
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
				
		cnt++;//the cell itself
		cntdiag=cnt;
		//triplet->Ti[cnt]=(UF_long)(i-1);
		//triplet->Tj[cnt]=(UF_long)(i-1);
		triplet->Tx[cnt]=0.0;
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			h = H->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
			if(h>0){
				triplet->Tx[cnt] += 1./Dt;
				//add subsuperficial flow to the channels (sink depending on H)
				if(adt->T->pixel_type->co[r][c]>=10) triplet->Tx[cnt] +=  adt->P->Kch_b * adt->P->w_dx;
			}
		}else{
			dz = adt->S->pa->co[sy][jdz][l];
			triplet->Tx[cnt] += dtheta_dpsi_from_psi(H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin)*dz/Dt;
		}
		

		//vertical hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				k = k_from_psi(jKv, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				k = k_from_psi(jKv, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}

		//flux from cell below
		if (l<Nl) {
			cnt++;
			I = adt->T->i_cont[l+1][r][c];	
			//triplet->Ti[cnt]=(UF_long)(i-1);
			//triplet->Tj[cnt]=(UF_long)(I-1);
			
			if(l==0){	//overland flow
				
				if( H->co[i]< H->co[I] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];
				
			}else{	//subsurface flow
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l+1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
				
			}
			
			triplet->Tx[cnt] = -kn/dzn;
			triplet->Tx[cntdiag] += kn/dzn;
		
		}		
			

		//flux from cell above
		if(l>0){
			cnt++;
			I = adt->T->i_cont[l-1][r][c];	
			//triplet->Ti[cnt]=(UF_long)(i-1);
			//triplet->Tj[cnt]=(UF_long)(I-1);
			
			if(l==1){	
				
				if( H->co[I] < H->co[i] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l];
				
			}else{
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l-1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l-1, r, c, adt->S), l-1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
				
			}
			
			triplet->Tx[cnt] = -kn/dzn;
			triplet->Tx[cntdiag] += kn/dzn;		
		
		}
		
		
		//lateral hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				k = k_from_psi(jKh, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				k = k_from_psi(jKh, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}		
		
		//lateral fluxes
		//4 neighbouring cells
		R = r-1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//1.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				cnt++;
				I=adt->T->i_cont[l][R][C];	
				//triplet->Ti[cnt]=(UF_long)(i-1);
				//triplet->Tj[cnt]=(UF_long)(I-1);		
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				triplet->Tx[cnt] = -(dz/ds)*kn/ds;
				triplet->Tx[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		R = r+1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//2.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				cnt++;
				I=adt->T->i_cont[l][R][C];	
				//triplet->Ti[cnt]=(UF_long)(i-1);
				//triplet->Tj[cnt]=(UF_long)(I-1);		
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				triplet->Tx[cnt] = -(dz/ds)*kn/ds;
				triplet->Tx[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}		
		
		R = r;
		C = c-1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//3.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				cnt++;
				I=adt->T->i_cont[l][R][C];	
				//triplet->Ti[cnt]=(UF_long)(i-1);
				//triplet->Tj[cnt]=(UF_long)(I-1);		
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				triplet->Tx[cnt] = -(dz/ds)*kn/ds;
				triplet->Tx[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}	
		
		
		R = r;
		C = c+1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//4.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				cnt++;
				I=adt->T->i_cont[l][R][C];	
				//triplet->Ti[cnt]=(UF_long)(i-1);
				//triplet->Tj[cnt]=(UF_long)(I-1);		
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				triplet->Tx[cnt] = -(dz/ds)*kn/ds;
				triplet->Tx[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
	}
	
	return 0;
	
}*/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//cnt is the counter of Li Lp Lx (lower diagonal)
int J_water(DOUBLEVECTOR *Lx, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i,l,r,c,I,R,C,sy,cnt=0,cntdiag;
	double h=0.0, dz=0.0, dzn=0.0, k=0.0, kn=0.0, klim=0.0, ds;
		
	for(i=1;i<=H->nh;i++){
						
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		
		cnt++;//the cell itself
		cntdiag=cnt;

		Lx->co[cntdiag]=0.0;
		
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			h = H->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
			if(h>0){
				Lx->co[cntdiag] += 1./Dt;
				//add subsuperficial flow to the channels (sink depending on H)
				if(adt->T->pixel_type->co[r][c]>=10) Lx->co[cntdiag] +=  adt->P->Kch_b * adt->P->w_dx;
			}
		}else{
			dz = adt->S->pa->co[sy][jdz][l];
			Lx->co[cntdiag] += dtheta_dpsi_from_psi(H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin)*dz/Dt;
		}
				
		
		//vertical hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				k = k_from_psi(jKv, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				k = k_from_psi(jKv, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}
		
		//flux from cell below
		if (l<Nl) {
			cnt++;
			I = i+1;
			
			if(l==0){	//overland flow
				
				if( H->co[i]< H->co[I] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];
				
			}else{	//subsurface flow
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l+1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
				
			}
			
			Lx->co[cnt] = -kn/dzn;
			Lx->co[cntdiag] += kn/dzn;

		}		
		
		//flux from cell above
		if(l>0){
			I = i-1;
			
			if(l==1){	
				
				if( H->co[I] < H->co[i] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l];
				
			}else{
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l-1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l-1, r, c, adt->S), l-1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
				
			}
			
			Lx->co[cntdiag] += kn/dzn;
			
		}
		
		
		//lateral hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				k = k_from_psi(jKh, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				k = k_from_psi(jKh, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}		
		
		//lateral fluxes
		//4 neighbouring cells
		R = r-1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//1.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
					
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
					
				if(adt->T->i_cont[l][R][C]>i){ 
					if(adt->T->i_cont[l][R][C]>i)cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}
					
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	

			}
		}
		
		R = r+1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//2.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					if(adt->T->i_cont[l][R][C]>i)cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}			
		R = r;
		C = c-1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//3.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					if(adt->T->i_cont[l][R][C]>i)cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		R = r;
		C = c+1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//4.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					if(adt->T->i_cont[l][R][C]>i)cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}		
	}
		
	return 0;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//cnt is the counter of Li Lp Lx (lower diagonal)
//cntt is the counter of Ui Up Ux (upper diagonal)
int J_water_with_transposed(DOUBLEVECTOR *Lx, DOUBLEVECTOR *Ux, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i,l,r,c,I,R,C,sy,cnt=0,cntt=0,cntdiag;
	double h=0.0, dz=0.0, dzn=0.0, kh=0.0, kv=0.0, kn=0.0, klim=0.0, ds;
	
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		
		//the cell itself
		cnt++;
		cntdiag=cnt;
		Lx->co[cntdiag]=0.0;
	
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			h = H->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
			if(h>0){
				Lx->co[cntdiag] += 1./Dt;
				//add subsuperficial flow to the channels (sink depending on H)
				if(adt->T->pixel_type->co[r][c]>=10) Lx->co[cntdiag] +=  adt->P->Kch_b * adt->P->w_dx;
			}
		}else{
			dz = adt->S->pa->co[sy][jdz][l];
			Lx->co[cntdiag] += dtheta_dpsi_from_psi(H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin)*dz/Dt;
		}
		
		
		//vertical hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				kv = k_from_psi(jKv, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				kv = k_from_psi(jKv, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}
		
		//flux from cell below
		if (l<Nl) {
			cnt++;
			I = i+1;
			
			if(l==0){	//overland flow
				
				if( H->co[i]< H->co[I] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];
				
			}else{	//subsurface flow
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l+1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, kv, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
				
			}
			
			Lx->co[cnt] = -kn/dzn;
			Lx->co[cntdiag] += kn/dzn;
			
		}		
		
		
		
		//lateral hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				kh = k_from_psi(jKh, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				kh = k_from_psi(jKh, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}		
		
		//lateral fluxes
		//4 neighbouring cells
		R = r-1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//1.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, kh, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}else{
					cntt++;
					Ux->co[cntt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		R = r+1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//2.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, kh, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}else{
					cntt++;
					Ux->co[cntt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		R = r;
		C = c-1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//3.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, kh, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}else{
					cntt++;
					Ux->co[cntt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		
		R = r;
		C = c+1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//4.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, kh, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}else{
					cntt++;
					Ux->co[cntt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		

		//flux from cell above
		if(l>0){
			cntt++;
			I = i-1;
			
			if(l==1){	
				
				if( H->co[I] < H->co[i] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l];
				
			}else{
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l-1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, kv, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l-1, r, c, adt->S), l-1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
				
			}
			
			Ux->co[cntt] = -kn/dzn;
			Lx->co[cntdiag] += kn/dzn;
			
		}
		
		
		//the cell itself
		cntt++;
		Ux->co[cntt] = Lx->co[cntdiag];
		
		
	}
			
	return 0;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//cnt is the counter of Li Lp Lx (lower diagonal without diagonal)
int find_matrix_K(DOUBLEVECTOR *Lx, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i, l, r, c, I, R, C, sy, cnt=0;
	double dz=0.0, dzn=0.0, k=0.0, kn=0.0, klim=0.0, ds;
	
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		
		
		//vertical hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				k = k_from_psi(jKv, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				k = k_from_psi(jKv, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
			dz = adt->S->pa->co[sy][jdz][l];
		}
		
				
		
		//flux from cell below
		if (l<Nl) {
			cnt++;
			I = i+1;
			
			if(l==0){	//overland flow
				
				if( H->co[i]< H->co[I] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];

			}else{	//subsurface flow
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l+1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
								
			}
						
			Lx->co[cnt] = -kn/dzn;
			
		}		
		
		
		//lateral hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				k = k_from_psi(jKh, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				k = k_from_psi(jKh, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}		
		
		
		//lateral fluxes
		//4 neighbouring cells
		R = r-1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//1.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV && adt->T->i_cont[l][R][C]>i){ 
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				cnt++;
				Lx->co[cnt] = -(dz/ds)*kn/ds;
				
			}
		}
		
		R = r+1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//2.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV && adt->T->i_cont[l][R][C]>i){ 
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				cnt++;
				Lx->co[cnt] = -(dz/ds)*kn/ds;
				
				
			}
		}
		
		R = r;
		C = c-1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//3.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV && adt->T->i_cont[l][R][C]>i){ 
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				cnt++;
				Lx->co[cnt] = -(dz/ds)*kn/ds;
				
				
			}
		}
		
		
		R = r;
		C = c+1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//4.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV && adt->T->i_cont[l][R][C]>i){ 
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, k, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				cnt++;
				Lx->co[cnt] = -(dz/ds)*kn/ds;
				
			}
		}
	}
	
	return 0;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//cnt is the counter of Li Lp Lx (lower diagonal including diagonal)
//cntt is the counter of Ui Up Ux (upper diagonal including diagonal)
int find_matrix_K_with_transposed(DOUBLEVECTOR *Lx, DOUBLEVECTOR *Ux, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i,l,r,c,I,R,C,sy,cnt=0,cntt=0,cntdiag;
	double dz=0.0, dzn=0.0, kh=0.0, kv=0.0, kn=0.0, klim=0.0, ds;
	
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		
		//the cell itself
		cnt++;
		cntdiag=cnt;
		Lx->co[cntdiag]=0.0;
				
		
		//vertical hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				kv = k_from_psi(jKv, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				kv = k_from_psi(jKv, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
			dz = adt->S->pa->co[sy][jdz][l];
		}
		
		//flux from cell below
		if (l<Nl) {
			cnt++;
			I = i+1;
			
			if(l==0){	//overland flow
				
				if( H->co[i]< H->co[I] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l+1];
				
			}else{	//subsurface flow
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l+1][r][c], l+1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l+1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, kv, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l+1, r, c, adt->S), l+1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
				
			}
			
			Lx->co[cnt] = -kn/dzn;
			Lx->co[cntdiag] += kn/dzn;
			
		}		
		
		
		
		//lateral hydraulic conductivity
		if(l>0){
			if(adt->P->UpdateK==1){
				kh = k_from_psi(jKh, H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}else{
				kh = k_from_psi(jKh, adt->S->P->co[l][r][c], l, r, c, adt->S, adt->P->imp );
			}
		}		
		
		//lateral fluxes
		//4 neighbouring cells
		R = r-1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//1.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, kh, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}else{
					cntt++;
					Ux->co[cntt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		R = r+1;
		C = c;
		ds = 1.E3*UV->U->co[2];	//[mm]
		//2.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, kh, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}else{
					cntt++;
					Ux->co[cntt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		R = r;
		C = c-1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//3.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, kh, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}else{
					cntt++;
					Ux->co[cntt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		
		R = r;
		C = c+1;
		ds = 1.E3*UV->U->co[1];	//[mm]
		//4.
		if(R>=1 && R<=Nr && C>=1 && C<=Nc && l>0 && adt->P->point_sim!=1){
			if(adt->L->LC->co[R][C]!=NoV){
				
				I=adt->T->i_cont[l][R][C];	
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKh, H->co[I] - adt->T->Z->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKh, adt->S->P->co[l][R][C], l, R, C, adt->S, adt->P->imp );
				}
				kn = Mean(adt->P->harm_or_arit_mean, ds, ds, kh, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l, R, C, adt->S), l, R, C, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				
				if(adt->T->i_cont[l][R][C]>i){ 
					cnt++;
					Lx->co[cnt] = -(dz/ds)*kn/ds;
				}else{
					cntt++;
					Ux->co[cntt] = -(dz/ds)*kn/ds;
				}
				
				Lx->co[cntdiag] += (dz/ds)*kn/ds;	
				
			}
		}
		
		
		//flux from cell above
		if(l>0){
			cntt++;
			I = i-1;
			
			if(l==1){	
				
				if( H->co[I] < H->co[i] ){	//upward flux
					kn = k_from_psi(jKv,  H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, adt->P->imp );
				}else{	//downward flow
					kn = k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp );
				}
				
				dzn = 0.5*adt->S->pa->co[sy][jdz][l];
				
			}else{
				
				if(adt->P->UpdateK==1){
					kn = k_from_psi(jKv,  H->co[I] - adt->T->Z->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
				}else{
					kn = k_from_psi(jKv,  adt->S->P->co[l-1][r][c], l-1, r, c, adt->S, adt->P->imp );
				}
				
				dzn = adt->S->pa->co[sy][jdz][l-1];
				kn = Mean(adt->P->harm_or_arit_mean, dz, dzn, kv, kn);
				klim = Harmonic_Mean( dz, dzn, k_from_psi(jKv,  psisat_from(l, r, c, adt->S), l, r, c, adt->S, adt->P->imp ),
									 k_from_psi(jKv,  psisat_from(l-1, r, c, adt->S), l-1, r, c, adt->S, adt->P->imp ) );
				if(kn > klim) kn=klim;
				dzn = 0.5*dz + 0.5*dzn;
				
			}
			
			Ux->co[cntt] = -kn/dzn;
			Lx->co[cntdiag] += kn/dzn;
			
		}
		
		
		//the cell itself
		cntt++;
		Ux->co[cntt] = Lx->co[cntdiag];
		
		
	}
	
	return 0;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_dfdH(DOUBLEVECTOR *df, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){
	
	long i,r,l,c,sy,ch;
	double h,dz;
		
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];		
		sy=adt->S->type->co[r][c];
		ch=adt->C->ch->co[r][c];
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			h = H->co[i] - 1.E3*adt->T->Z0dp->co[r][c];
			if(h>0){
				df->co[i] = 1./Dt;
			}else{
				df->co[i] = 0.;
			}
		}else{
			dz = adt->S->pa->co[sy][jdz][l];
			df->co[i] = dtheta_dpsi_from_psi(H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin)*dz/Dt;
		}
	
		
		//add subsuperficial flow to the channels (sink depending on H)
		if(ch!=0 && l==0){
			h = H->co[i] - adt->T->Z->co[l][r][c];
			if( h > adt->C->h_sup->co[ch]) df->co[i] +=  (adt->P->w_dx*adt->C->length->co[ch]/UV->U->co[1]) * adt->P->Kch_b;
		}
		
		
	}
	return 0;
}
		
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_f(DOUBLEVECTOR *f, ALLDATA *adt, DOUBLEVECTOR *H, double Dt){

	long i, l, r, c, sy, ch;
	double dz, h, V0, V1;
	
	for(i=1;i<=H->nh;i++){
		
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		ch=adt->C->ch->co[r][c];
						
		//hydraulic capacity (diagonal term)
		if(l==0){
			V1 = Fmax(0.0, H->co[i] - adt->T->Z->co[l][r][c]);
			V0 = Fmax(0.0, adt->S->P->co[l][r][c]);
		}else{
			dz = adt->S->pa->co[sy][jdz][l];		
			V1 = dz * theta_from_psi( H->co[i] - adt->T->Z->co[l][r][c], l, r, c, adt->S, PsiMin );
			V0 = dz * theta_from_psi( adt->S->P->co[l][r][c], l, r, c, adt->S, PsiMin );
		}
		f->co[i] = (V1-V0)/Dt;
			
		//subsuperficial flow to the channels
		if(ch!=0 && l==0){ 
			h = H->co[i] - adt->T->Z->co[l][r][c];
			if( h > adt->C->h_sup->co[ch]) f->co[i] += (adt->P->w_dx*adt->C->length->co[ch]/UV->U->co[1]) * adt->P->Kch_b * ( h - adt->C->h_sup->co[ch] );
		}
		
	
		//evaporation and precipitation
		if(l>0){
			f->co[i] += adt->S->ET->co[l][r][c];
		}else{
			f->co[i] -= adt->W->Pnet->co[r][c];	
		}
				
						
	}
	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_dt_max_channel(double Courant, DOUBLEVECTOR *h, TOPO *top, CHANNEL *cnet, PAR *par, double *dt){

	short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	long r, c, ch, r_down, c_down, ch_down;		
	double b[10];    //channel length
	double dn, Ks, q, Vmax, i;

	dn = par->w_dx * UV->U->co[1];		//transversal [m]
	
	b[1]=UV->U->co[1];					//distance between channel centers
	b[2]=UV->U->co[1];             
	b[3]=UV->U->co[1]*sqrt(2.); 
	b[4]=UV->U->co[1];            
	b[5]=UV->U->co[1]*sqrt(2.); 
	b[6]=UV->U->co[1];   
	b[7]=UV->U->co[1]*sqrt(2.); 
	b[8]=UV->U->co[1];             
	b[9]=UV->U->co[1]*sqrt(2.); 
	
	for(ch=1;ch<=cnet->r->nh;ch++){
		
		if(h->co[ch] > 0){
			Ks=cm_h(par->Ks_channel, h->co[ch], par->thres_hchannel);
			
			r = cnet->r->co[ch];
			c = cnet->c->co[ch];
			
			r_down = r+r_DD[top->DD->co[r][c]];
			c_down = c+c_DD[top->DD->co[r][c]];
			ch_down = cnet->ch->co[r_down][c_down];
			
			if(ch_down==0) ch_down=ch;
						
			i = Fmax( 0.0 , 1.E-3 * ( h->co[ch] - h->co[ch_down] ) / b[top->DD->co[r][c]+1] + top->i_DD->co[r][c] );	
			
			q = dn * Ks * pow( 1.E-3 * h->co[ch] , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
			Vmax = 1.E-3*h->co[ch]*dn*cnet->length->co[ch];	//m3
			
			if(q>0){
				if(Courant*Vmax/q<(*dt)) *dt=Courant*Vmax/q; 
			}
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void channel_flow(double Dt, DOUBLEVECTOR *h, DOUBLEVECTOR *dV, TOPO *top, CHANNEL *cnet, PAR *par, double *Vout)

{
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	long r,c,ch,r_down,c_down,ch_down;                                    
	double dn;                                         // the two dimensions of a channel pixel
	double Ks;											  // the Strickler's coefficent
	double b[10];                                         // distance between channels centers
	double i;											  // hydraulic gradient
	double q,tb,te,dt;
	
	if( par->point_sim!=1 && cnet->r->co[1]!=0 ){	//if it is not point simulation and there are channels
		
		dn = par->w_dx * UV->U->co[1];		//transversal length [m]
		
		b[1]=UV->U->co[1];					//distance between channel centers
		b[2]=UV->U->co[1];             
		b[3]=UV->U->co[1]*sqrt(2.); 
		b[4]=UV->U->co[1];            
		b[5]=UV->U->co[1]*sqrt(2.); 
		b[6]=UV->U->co[1];   
		b[7]=UV->U->co[1]*sqrt(2.); 
		b[8]=UV->U->co[1];             
		b[9]=UV->U->co[1]*sqrt(2.); 
		
		te=0.0;
		
		do{
			
			tb=te;
			dt=Dt;
			
			find_dt_max_channel(max_Courant_channel, h, top, cnet, par, &dt);
			
			te=tb+dt;
			if(te>Dt){
				te=Dt;
				dt=te-tb;
			}		
			
			for(ch=1;ch<=cnet->r->nh;ch++){
												
				if(h->co[ch] > 0){
					
					r = cnet->r->co[ch];
					c = cnet->c->co[ch];
					r_down = r+r_DD[top->DD->co[r][c]];
					c_down = c+c_DD[top->DD->co[r][c]];
					ch_down = cnet->ch->co[r_down][c_down];
					
					Ks=cm_h(par->Ks_channel, h->co[ch], par->thres_hchannel);
					
					if(ch_down==0){//outlet section
						
						q = Cd*(2./3.)*sqrt(2.*g*1.E-3*h->co[ch])*(1.E-3*h->co[ch])*dn;	//[m3/s]
						
					}else{

						i = Fmax( 0.0 , 1.E-3 * ( h->co[ch] - h->co[ch_down] ) / b[top->DD->co[r][c]+1] + top->i_DD->co[r][c] );	
						q = dn * Ks * pow( 1.E-3 * h->co[ch] , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
						
					}
					
					dV->co[ch] = Fmin( q*dt , 1.E-3*h->co[ch]*dn*cnet->length->co[ch] );	//m3
					
				}else{
					
					dV->co[ch] = 0.0;
				
				}
				
			}
								
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				r_down = r+r_DD[top->DD->co[r][c]];
				c_down = c+c_DD[top->DD->co[r][c]];
				ch_down = cnet->ch->co[r_down][c_down];
				if(ch_down==0) ch_down=ch;

				h->co[ch] -= 1.E3*dV->co[ch]/(dn*cnet->length->co[ch]);
												
				if (ch_down != ch){	//not the outlet
					h->co[ch_down] += 1.E3*dV->co[ch]/(dn*cnet->length->co[ch_down]);	//mm
				}else{	//outlet
					*Vout = *Vout + dV->co[ch];	//m3
				}
			}
			
			for(ch=1;ch<=cnet->r->nh;ch++){
				
				dV->co[ch] = (cnet->Qsup->co[ch]+cnet->Qsub->co[ch])*dt; //m3
				r = cnet->r->co[ch];
				c = cnet->c->co[ch];
				h->co[ch] += 1.E3*dV->co[ch]/(dn*cnet->length->co[ch]);	//mm
				
			}
						
			
		}while(te<Dt);
	}
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
