
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */


#include "water.balance.h"
#include "geotop_common.h"
#include "inputKeywords.h"

using namespace std;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//short water_balance(double Dt, double JD0, double JD1, double JD2, SOIL_STATE *L, SOIL_STATE *C, ALLDATA *adt, DOUBLEVECTOR *Vsub, DOUBLEVECTOR *Vsup,
//					double *Voutnet, double *Voutlandsub, double *Voutlandsup, double *Voutlandbottom){

short water_balance(double Dt, double JD0, double JD1, double JD2, SoilState *L, SoilState *C, AllData *adt, GeoVector<double>& Vsub, GeoVector<double>& Vsup,
					double *Voutnet, double *Voutlandsub, double *Voutlandsup, double *Voutlandbottom){

	clock_t start, end;
	FILE *flog;
	double Pnet, loss;
	long j;
	short a;
	
	// this below used only for checks..
	double mm1, mm2, mmo;
	
	double MM1, MM2, MMR=0.0, MMo=0.0, MS1, MS2;
	double m1=0., m2=0., mo=0.;
	double ds, area, dz;
	long r, c, l, sy;
	
	
	flog = fopen(geotop::common::Variables::logfile.c_str(), "a");
	
	if(adt->P->qin==1){
		time_interp_linear(JD0, JD1, JD2, adt->M->qinv, adt->M->qins, adt->M->qinsnr, 2, 0, 0, &(adt->M->qinline));		
	}
	
	if (adt->P->point_sim != 1) {//distributed simulations
		
		//surface flow: 1st half of time step
		start=clock();

		supflow(Dt/2., adt->I->time, L->P, &adt->W->h_sup[0], C->P, &adt->C->h_sup[0], adt->T, adt->L, adt->W, adt->C, adt->P, adt->M, Vsup, Voutnet, Voutlandsup, flog, &mm1, &mm2, &mmo );
		
		
		end=clock();

#ifdef VERBOSE		
		// this below for debugging: could be removed at later stage.. SC 
		MMo += mmo;
		MS1 = mm1;
		MS2 = mm2;	
		MM1 = 0.;
		MM2 = 0.;
		printf("water-balance: %e %e %e %e %e %e\n",MM1,MM2,MMR,MS1,MS2,MMo);
#endif 				
		
		geotop::common::Variables::t_sup += (end-start)/(double)CLOCKS_PER_SEC;
		
		//subsurface flow with time step Dt0 (decreasing if not converging)
		start = clock();

		
		/* ds=sqrt(UV->U[1]*UV->U[2]);
		for (j=1; j<adt->W->H1.size(); j++) {
			l=adt->T->lrc_cont[j][1];
			r=adt->T->lrc_cont[j][2];
			c=adt->T->lrc_cont[j][3];
			sy=adt->S->type[r][c];
			area=ds*ds/cos(adt->T->slope[r][c]*GTConst::Pi/180.);
			if(l==0){
				m1 += area * 1.E-3*Fmax(0.0, L->P[l][adt->T->j_cont[r][c]]) / cos(adt->T->slope[r][c]*GTConst::Pi/180.);
			}else {
				dz = adt->S->pa[sy][jdz][l];		
				m1 += area*1.E-3*dz * theta_from_psi(L->P[l][adt->T->j_cont[r][c]], 0, l, adt->S->pa, sy, GTConst::PsiMin);
				
			}
			if(l==0) mo += area * 1.E-3 * adt->W->Pnet[r][c];
		}*/
		
		
		
		
		a = Richards3D(Dt, L, C, adt, flog, &loss, Vsub, Voutlandbottom, Voutlandsub, &Pnet, adt->P->UpdateK);
		
		end=clock();
		geotop::common::Variables::t_sub += (end-start)/(double)CLOCKS_PER_SEC;
		if (a != 0){
			fclose(flog);
			return 1;
		}
		
		/*ds=sqrt(UV->U[1]*UV->U[2]);
		for (j=1; j<adt->W->H1.size(); j++) {
			l=adt->T->lrc_cont[j][1];
			r=adt->T->lrc_cont[j][2];
			c=adt->T->lrc_cont[j][3];
			sy=adt->S->type[r][c];
			area=ds*ds/cos(adt->T->slope[r][c]*GTConst::Pi/180.);
			if(l==0){
				m2 += area * 1.E-3*Fmax(0.0, L->P[l][adt->T->j_cont[r][c]]) / cos(adt->T->slope[r][c]*GTConst::Pi/180.);
			}else {
				dz = adt->S->pa[sy][jdz][l];		
				m2 += area*1.E-3*dz * theta_from_psi(L->P[l][adt->T->j_cont[r][c]], 0, l, adt->S->pa, sy, GTConst::PsiMin);
			}
		}
		
#ifdef VERBOSE		
		printf("water-balance SUB: m1:%e m2:%e mo:%e dm:%e\n",m1,m2,mo,fabs(m1+mo-m2));
#endif 		
		MM1 = m1;
		MM2 = m2;
		MMR += mo;*/
		

		
		//surface flow: 2nd half of time step
		start=clock();

		supflow(Dt/2., adt->I->time, L->P, &adt->W->h_sup[0], C->P, &adt->C->h_sup[0], adt->T, adt->L, adt->W, adt->C, adt->P, adt->M, Vsup, Voutnet, Voutlandsup, flog, &mm1, &mm2, &mmo );

		end=clock();
		
		MMo += mmo;
		MS1 = mm1;
		MS2 = mm2;
		
#ifdef VERBOSE	
		printf("water-balance after 2nd half: %e %e %e %e %e %e\n",MM1,MM2,MMR,MS1,MS2,MMo);
#endif 	
		
		geotop::common::Variables::t_sup += (end-start)/(double)CLOCKS_PER_SEC;
		
	}else {//point simulations
		
		start = clock();
		for (j=1; j<=adt->P->total_pixel; j++) {			
			a = Richards1D(j, Dt, L, adt, flog, &loss, Voutlandbottom, Voutlandsub, &Pnet, adt->P->UpdateK);
			if (a != 0){
				fclose(flog);
				return 1;
			}			
		
			if( L->P[0][j] > 0 ) L->P[0][j] = Fmin( L->P[0][j], Fmax(0.,-adt->T->BC_DepthFreeSurface[j])*cos(adt->T->slope[1][j]*GTConst::Pi/180.) );
		}
		end=clock();
		geotop::common::Variables::t_sub += (end-start)/(double)CLOCKS_PER_SEC;
		
	}
	
	geotop::common::Variables::odb[oopnet] = Pnet;
	geotop::common::Variables::odb[oomasserror] = loss;
	
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



short Richards3D(double Dt, SoilState *L, SoilState *C, AllData *adt, FILE *flog, double *loss, GeoVector<double>& Vsub, double *Vbottom, double *Vlatsub, double *Total_Pnet, short updateK){
	

	double res=0.0, res0[3], res_prev[MM], res_av, res00, lambda[3], epsilon, mu=0., hnew, hold=0.;
	double ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]), area, dz, dn, dD;
	double psi;
	
	long i, j, ch, l, r, c, m, bc, sy, cont, cont2, iter;
	long n=(geotop::common::Variables::Nl+1)*adt->P->total_pixel;
	long N=adt->W->H0.size();
	long cont_lambda_min=0;
	short out, out2;	
	int sux;
	FILE *f;
	
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
	
	for(i=1; i<N; i++){
		
		if (i<=n) {
			
			l = adt->T->lrc_cont[i][1];		
			r = adt->T->lrc_cont[i][2];	
			c = adt->T->lrc_cont[i][3];
			j = adt->T->j_cont[r][c];
			
			//precipitation

			if (l == 0 && adt->W->Pnet[r][c] > 0) {
				*Total_Pnet = *Total_Pnet + (adt->W->Pnet[r][c]/cos(adt->T->slope[r][c]*GTConst::Pi/180.)) / (double)adt->P->total_pixel;
			}
			
			//solution guess
			
			if (adt->W->Pnet[r][c] > 0 && l == 0) {
				adt->W->H1[i] = Fmax(0., L->P[l][j]) + (adt->W->Pnet[r][c]/cos(adt->T->slope[r][c]*GTConst::Pi/180.)) + adt->T->Z[l][r][c];
			}else {
				adt->W->H1[i] = L->P[l][j] + adt->T->Z[l][r][c];
			}
			
						
		}else {
			
			l = adt->C->lch[i-n][1];
			ch = adt->C->lch[i-n][2];
			r = adt->C->r[ch];
			c = adt->C->c[ch];
			
			//solution guess
			
			if (adt->W->Pnet[r][c] > 0 && l == 0) {
				adt->W->H1[i] = Fmax(0., C->P[l][ch]) + (adt->W->Pnet[r][c]/cos(adt->T->slope[r][c]*GTConst::Pi/180.)) + ( adt->T->Z[l][r][c] - adt->P->depr_channel );
			}else {
				adt->W->H1[i] = C->P[l][ch] + ( adt->T->Z[l][r][c] - adt->P->depr_channel );	
			}
			
		 						
		}
		
 //	  printf("H1 guess %ld,%f\n", i, adt->W->H1[i]);		
	}

	printf("Richard3D: Pnet: %f\n",*Total_Pnet);

	
	sux = find_matrix_K_3D(Dt, L, C, adt->W->Lx, adt->W->Klat, adt->W->Kbottom, adt->C->Kbottom, adt, adt->W->H1);
	
	find_f_3D(Dt, adt->W->f, adt, L, C, adt->W->H1, adt->W->Klat, adt->W->Kbottom, adt->C->Kbottom);
	
	product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
	
	res = norm_inf(adt->W->B, 1, N);
	printf("Richard3D: res: %f\n",res);
	
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
		
		for(i=1; i<N; i++){
			adt->W->H0[i] = adt->W->H1[i];
			adt->W->dH[i] = 0.;
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
				lambda[0] = GTConst::thmax;
				
			}else{
				lambda[2] = lambda[1];
				res0[2] = res0[1];
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1], lambda[2], res0[2]);
				
			}
			
			for(i=1; i<N; i++){
				
		
				adt->W->H1[i] = adt->W->H0[i] + lambda[0] * adt->W->dH[i];
				
				if(adt->W->H1[i] != adt->W->H1[i]) {
					f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
					fprintf(f, "Simulation Period:%ld\n",geotop::common::Variables::i_sim);
					fprintf(f, "Run Time:%ld\n",geotop::common::Variables::i_run);
					fprintf(f, "Number of days after start:%f\n",adt->I->time/86400.);					
					fprintf(f, "Error: no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);
					fclose(f);
					t_error("Fatal Error! Geotop is closed. See failing report.");	
				}
				
			}
			
			if(updateK == 1 && cont <= maxITER_rec_K) sux = find_matrix_K_3D(Dt, L, C, adt->W->Lx, adt->W->Klat, adt->W->Kbottom, adt->C->Kbottom, adt, adt->W->H1);		
			
			find_f_3D(Dt, adt->W->f, adt, L, C, adt->W->H1, adt->W->Klat, adt->W->Kbottom, adt->C->Kbottom);
			
			product_matrix_using_lower_part_by_vector_plus_vector(-1., adt->W->B, adt->W->f, adt->W->H1, adt->T->Li, adt->T->Lp, adt->W->Lx);	
			
			res = norm_inf(adt->W->B, 1, N);		
			
			out2=0;
			
			printf(" Richar3D: cnt:%ld res:%e lambda:%e Dt:%f P:%f\n",cont,res,lambda[0],Dt,*Total_Pnet);

			
			if(res <= (1.0 - ni_wat*lambda[0]*(1.-mu))*res_av) out2=1;
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
	for(i=1; i<N; i++){
		
		if (i<=n) {//land
			
			l = adt->T->lrc_cont[i][1];
			r = adt->T->lrc_cont[i][2];
			c = adt->T->lrc_cont[i][3];
			j = adt->T->j_cont[r][c];
			sy = adt->S->type[r][c];
			ch = adt->C->ch[r][c];
			bc = adt->T->BC_counter[r][c];
			
			L->P[l][j] = adt->W->H1[i] - adt->T->Z[l][r][c];
			
			//update variables
			if(l>0){
				adt->S->th[l][j] = theta_from_psi(L->P[l][j], L->thi[l][j], l, adt->S->pa, sy, GTConst::PsiMin);
				adt->S->Ptot[l][j] = psi_from_theta(adt->S->th[l][j]+L->thi[l][j], 0., l, adt->S->pa, sy, GTConst::PsiMin);
				adt->S->th[l][j] = Fmin( adt->S->th[l][j] , adt->S->pa(sy,jsat,l) - L->thi[l][j] );
			}
			
			//volume lost at the bottom
			if(l==geotop::common::Variables::Nl){
				area = ds*ds/cos(adt->T->slope[r][c]*GTConst::Pi/180.);
				if (ch>0) area -= adt->C->length[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
				*Vbottom = *Vbottom + area * adt->W->Kbottom[r][c] * 1.E-3 * Dt;
			}
			
			//lateral drainage at the border
			if (bc>0) {
				
				if (l>0) {
					
				
					if (adt->T->pixel_type[r][c] == 1) {
					//	The depth of the free surface is multiplied by cosine since Z's are the layer depths in vertical direction
					
						if ( adt->T->Z[0][r][c] - adt->T->Z[l][r][c] <= adt->T->BC_DepthFreeSurface[bc]*cos(adt->T->slope[r][c]*GTConst::Pi/180.) && adt->W->H1[i] - adt->T->Z[l][r][c] > 0 ) {
							
							dz = adt->S->pa(sy,jdz,l);//[mm]						
							
							if ((long)adt->L->LC[r+1][c]==geotop::input::gDoubleNoValue || (long)adt->L->LC[r-1][c]==geotop::input::gDoubleNoValue) {
						
								dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdN[r][c]));//mm
								dn = ds / cos(atan(adt->T->dzdE[r][c]));//m
							}else {
								dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdE[r][c]));
								dn = ds / cos(atan(adt->T->dzdN[r][c]));//m
							}
							
							*Vlatsub = *Vlatsub + Dt * (dn*dz*1.E-3) * 1.E-3*adt->W->Klat[bc][l]*adt->P->free_drainage_lateral*(adt->W->H1[i] - adt->T->Z[l][r][c]) / dD;
							
						}
						
				
					}else if (adt->T->pixel_type[r][c] == 2) {
						
						if ( adt->T->Z[0][r][c] - adt->T->Z[l][r][c] <= adt->T->BC_DepthFreeSurface[bc]*cos(adt->T->slope[r][c]*GTConst::Pi/180.) ) {
							
							dz = adt->S->pa(sy,jdz,l);//[mm]
							
						
							if ((long)adt->L->LC[r+1][c]==geotop::input::gDoubleNoValue || (long)adt->L->LC[r-1][c]==geotop::input::gDoubleNoValue) {							
								dn = ds / cos(atan(adt->T->dzdE[r][c]));//m
							}else {							
								dn = ds / cos(atan(adt->T->dzdN[r][c]));//m
							}
							
							*Vlatsub = *Vlatsub + Dt * (dn*dz*1.E-3) * 1.E-3*adt->W->Klat[bc][l]*adt->P->free_drainage_lateral;
						}
					}
				}
			}
			
			
		}else {//channel
			
		
			l = adt->C->lch[i-n][1];	
			ch = adt->C->lch[i-n][2];		
			r = adt->C->r[ch];	
			c = adt->C->c[ch];		
			sy = adt->C->soil_type[ch];
			
			if (l==0){
				//hold and hnew are normal
				hold = Fmax(0., C->P[l][ch]) / cos(adt->T->slope[r][c] * GTConst::Pi/180.);
			}
			
			//depr channel is defined vertical
			C->P[l][ch] = adt->W->H1[i] - ( adt->T->Z[l][r][c] - adt->P->depr_channel );
			
			if(l>0){
				adt->C->th[l][ch] = theta_from_psi(C->P[l][ch], C->thi[l][ch], l, adt->S->pa, sy, GTConst::PsiMin);
				adt->C->th[l][ch] = Fmin( adt->C->th[l][ch] , adt->S->pa(sy,jsat,l) - C->thi[l][ch] );
			}
			
			if (l==0){
			//	hold and hnew are normal
				hnew = Fmax(0., C->P[l][ch]) / cos(adt->T->slope[r][c] * GTConst::Pi/180.);
				Vsub[ch] += 1.E-3 * ( hnew - hold ) * adt->C->length[ch] * adt->P->w_dx * ds;
			}
			
			if(l==geotop::common::Variables::Nl){
				area = adt->C->length[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
				*Vbottom = *Vbottom + area * adt->C->Kbottom[ch] * 1.E-3 * Dt;
			}
		}
	}
	
	return 0;			
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//short Richards1D(long c, double Dt, SOIL_STATE *L, ALLDATA *adt, FILE *flog, double *loss, double *Vbottom, double *Vlat, double *Total_Pnet, short updateK){
short Richards1D(long c, double Dt, SoilState *L, AllData *adt, FILE *flog, double *loss, double *Vbottom, double *Vlat, double *Total_Pnet, short updateK){
	
	double res=0.0, res0[3], res_prev[MM], res_av, res00, lambda[3], epsilon, mu;
	double ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]), area, dz, dn, dD;
	
	long i, l, r=1, m, bc, sy, cont, cont2, iter;
	long N=adt->W->H0.size();
	long cont_lambda_min=0;
	short out, out2;	
	int sux;
	FILE *f;
	
	*Total_Pnet = 0.;
	
	for(i=1; i<=N; i++){//layers
		
		l = i-1;
		
	//	if (l == 0 && adt->W->Pnet->co[r][c] > 0) {
		if (l == 0 && adt->W->Pnet[r][c] > 0) {
		//	*Total_Pnet = *Total_Pnet + (adt->W->Pnet->co[r][c]/cos(Fmin(GTConst::max_slope,adt->T->slope->co[r][c])*GTConst::Pi/180.)) / (double)adt->P->total_pixel;
			*Total_Pnet = *Total_Pnet + (adt->W->Pnet[r][c]/cos(Fmin(GTConst::max_slope,adt->T->slope[r][c])*GTConst::Pi/180.)) / (double)adt->P->total_pixel;
		}
		
	//	solution guess
	//	adt->W->H1->co[i] = L->P->co[l][c] + adt->T->Z->co[l][r][c];
		adt->W->H1[i] = L->P[l][c] + adt->T->Z[l][r][c];
		
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
		//	adt->W->H0->co[i] = adt->W->H1->co[i];
			adt->W->H0[i] = adt->W->H1[i];
		//	adt->W->dH->co[i] = 0.;
			adt->W->dH[i] = 0.;
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
				lambda[0] = GTConst::thmax;
				
			}else{
				lambda[2] = lambda[1];
				res0[2] = res0[1];
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1], lambda[2], res0[2]);
				
			}
			
			for(i=1; i<=N; i++){
				
			//	adt->W->H1->co[i] = adt->W->H0->co[i] + lambda[0] * adt->W->dH->co[i];
				adt->W->H1[i] = adt->W->H0[i] + lambda[0] * adt->W->dH[i];
				
			//	if(adt->W->H1->co[i] != adt->W->H1->co[i]) {
				if(adt->W->H1[i] != adt->W->H1[i]) {
					f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
					fprintf(f, "Simulation Period:%ld\n",geotop::common::Variables::i_sim);
					fprintf(f, "Run Time:%ld\n",geotop::common::Variables::i_run);
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
			
			if(res <= (1.0 - ni_wat*lambda[0]*(1.-mu))*res_av) out2=1;
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
	//	sy = adt->S->type->co[r][c];
		sy = adt->S->type[r][c];
	//	bc = adt->T->BC_counter->co[r][c];
		bc = adt->T->BC_counter[r][c];
		
	//	L->P->co[l][c] = adt->W->H1->co[i] - adt->T->Z->co[l][r][c];
		L->P[l][c] = adt->W->H1[i] - adt->T->Z[l][r][c];
		
		//update variables
		if(l>0){
		//	adt->S->th->co[l][c] = theta_from_psi(L->P->co[l][c], L->thi->co[l][c], l, adt->S->pa->co[sy], GTConst::PsiMin);
			adt->S->th[l][c] = theta_from_psi(L->P[l][c], L->thi[l][c], l, adt->S->pa, sy, GTConst::PsiMin);
		//	adt->S->Ptot->co[l][c] = psi_from_theta(adt->S->th->co[l][c]+L->thi->co[l][c], 0., l, adt->S->pa->co[sy], GTConst::PsiMin);
			adt->S->Ptot[l][c] = psi_from_theta(adt->S->th[l][c]+L->thi[l][c], 0., l, adt->S->pa, sy, GTConst::PsiMin);
		//	adt->S->th->co[l][c] = Fmin( adt->S->th->co[l][c] , adt->S->pa->co[sy][jsat][l]-L->thi->co[l][c] );
			adt->S->th[l][c] = Fmin( adt->S->th[l][c] , adt->S->pa(sy,jsat,l) - L->thi[l][c] );
		}
		
		//volume lost at the bottom
		if(l==geotop::common::Variables::Nl){
			area = ds*ds;
		//	*Vbottom = *Vbottom + area * adt->W->Kbottom->co[r][c] * 1.E-3 * Dt;
			*Vbottom = *Vbottom + area * adt->W->Kbottom[r][c] * 1.E-3 * Dt;
		}
		
		//lateral drainage at the border
		if (bc>0 && l>0) {
		//	if (adt->T->pixel_type->co[r][c] == 1) {
			if (adt->T->pixel_type[r][c] == 1) {
			//	if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*GTConst::Pi/180.) && adt->W->H1->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
				if ( adt->T->Z[0][r][c] - adt->T->Z[l][r][c] <= adt->T->BC_DepthFreeSurface[bc]*cos(adt->T->slope[r][c]*GTConst::Pi/180.) && adt->W->H1[i] - adt->T->Z[l][r][c] > 0 ) {
					//dz = adt->S->pa->co[sy][jdz][l];//[mm]
					dz = adt->S->pa(sy,jdz,l);//[mm]
					dn = ds;
					dD = 0.5 * 1.E3*ds;
				//	*Vlat = *Vlat + Dt * (dn*dz*1.E-3) * 1.E-3*adt->W->Klat->co[bc][l]*adt->P->free_drainage_lateral*(adt->W->H1->co[i] - adt->T->Z->co[l][r][c]) / dD;
					*Vlat = *Vlat + Dt * (dn*dz*1.E-3) * 1.E-3*adt->W->Klat[bc][l]*adt->P->free_drainage_lateral*(adt->W->H1[i] - adt->T->Z[l][r][c]) / dD;
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



int find_matrix_K_3D(double Dt, SoilState *SL, SoilState *SC, GeoVector<double>& Lx, GeoMatrix<double>& Klat, GeoMatrix<double>& Kbottom_l, GeoVector<double>& Kbottom_ch, AllData *adt, const GeoVector<double>& H){
	
	size_t i;
	size_t l, r, c, j, I, R, C, J, sy, syn, ch, cnt=0;
	size_t n=(geotop::common::Variables::Nl+1)*adt->P->total_pixel;
	double dz=0.0, dzn=0.0, dD=0.0, k=0.0, kn=0.0, kmax=0.0, kmaxn=0.0;
	double area, ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]), dn;
	double psi, ice, a, ns, res, sat, ss, Temp;
	

	for(i=1;i<H.size();i++){
		//VERTICAL FLUXES
		if( i<=n){//land
			

			l=adt->T->lrc_cont[i][1];
			r=adt->T->lrc_cont[i][2];
			c=adt->T->lrc_cont[i][3];
			j=adt->T->j_cont[r][c];
			sy=adt->S->type[r][c];
			
			ch=adt->C->ch[r][c];
			area=ds*ds/cos(adt->T->slope[r][c]*GTConst::Pi/180.);
			if (ch>0) area-=adt->C->length[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			
		//	vertical hydraulic conductivity
			if(l>0){				
				dz = adt->S->pa(sy,jdz,l);
				if (l==geotop::common::Variables::Nl && adt->P->free_drainage_bottom>0) Kbottom_l[r][c] = k_from_psi(jKn, H[i] - adt->T->Z[l][r][c], SL->thi[l][c], SL->T[l][c], l, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat	);
			}
			
			//flux from cell below
			if (l<geotop::common::Variables::Nl) {
				
				I = i+1;
				
				if(l==0){	//overland flow
					
				
					dzn = adt->S->pa(sy,jdz,l+1);
					dD = 0.5*dzn;
					
				
					if( H[i] < H[I] ){
						//upward flux
						kn = k_from_psi(jKn, H[I] - adt->T->Z[l+1][r][c], SL->thi[l+1][j], SL->T[l+1][j], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
					}else{	
						//downward flow
						kn = k_from_psi(jKn, psisat_from(SL->thi[l+1][j], l+1, adt->S->pa, sy), SL->thi[l+1][j], SL->T[l+1][j], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
					}
					
				}else{	//subsurface flow
					
					dzn = adt->S->pa[sy][jdz][l+1];
					dD = 0.5*dz + 0.5*dzn;
					
					if( H[i] < H[I] ){	
						//upward flux
						kn = k_from_psi(jKn, H[I] - adt->T->Z[l+1][r][c], SL->thi[l+1][j], SL->T[l+1][j], l+1, adt->S->pa,sy, adt->P->imp, adt->P->k_to_ksat);	
					}else{	
						//downward flow
						kn = k_from_psi(jKn, H[i] - adt->T->Z[l][r][c], SL->thi[l][j], SL->T[l][j], l, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);	
					}
					
					/*psi = Arithmetic_Mean(dz, dzn, H[i] - adt->T->Z[l][r][c], H[I] - adt->T->Z[l+1][r][c]);
					ice = Arithmetic_Mean(dz, dzn, SL->thi[l][j], SL->thi[l+1][j]);
					Temp = Arithmetic_Mean(dz, dzn, SL->T[l][j], SL->T[l+1][j]);
					a = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][ja][l], adt->S->pa[sy][ja][l+1]);
					ns = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jns][l], adt->S->pa[sy][jns][l+1]);
					res = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jres][l], adt->S->pa[sy][jres][l+1]);
					sat = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jsat][l], adt->S->pa[sy][jsat][l+1]);
					ss = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jss][l], adt->S->pa[sy][jss][l+1]);
					k = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jKn][l], adt->S->pa[sy][jKn][l+1]);
					kn = k_hydr_soil(psi, k, adt->P->imp, ice, sat, res, a, ns, 1.-1./ns, 0.5, Temp, adt->P->k_to_ksat);*/
					
				
					kmax = k_from_psi( jKn,  psisat_from( SL->thi[l][j], l, adt->S->pa, sy), SL->thi[l][j], SL->T[l][j], l, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
					kmaxn = k_from_psi( jKn,  psisat_from( SL->thi[l+1][j], l+1, adt->S->pa, sy), SL->thi[l+1][j], SL->T[l+1][j], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
					kmaxn = Fmin(kmax, kmaxn);
					
					kn = Fmin(kn, kmaxn);
					
				}
				
				cnt++;
			
				cout << "find3D H.size()=" << H.size()-1 << "  " << "Lx.size() = " << Lx.size()-1 << "   cnt=" << cnt  << endl;
				printf("find_3D kmax,kmaxn,kn,%f,%f,%f\n",kmax,kmaxn,kn);
				Lx[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
		        cout << "find3D Lx[cnt]=" << Lx[cnt] << endl;
			}
			
		}else{//channel
			
			l=adt->C->lch[i-n][1];
			ch=adt->C->lch[i-n][2];
			r=adt->C->r[ch];
			c=adt->C->c[ch];
			sy=adt->C->soil_type[ch];
			
			area=adt->C->length[ch] * adt->P->w_dx * ds;
			
			//vertical hydraulic conductivity
			if(l>0){
				dz = adt->S->pa[sy][jdz][l];
				if (l==geotop::common::Variables::Nl && adt->P->free_drainage_bottom>0) Kbottom_ch[ch] = k_from_psi(jKn, H[i] - (adt->T->Z[l][r][c]-adt->P->depr_channel), SC->thi[l][ch], SC->T[l][ch], l, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
			}
			
			//flux from cell below
			if (l<geotop::common::Variables::Nl) {
				
				I = i+1;
				
				if(l==0){	//overland flow
					
					dzn = adt->S->pa[sy][jdz][l+1];
					dD = 0.5*dzn;
					
				//	if( H->co[i] < H->co[I] ){
					if( H[i] < H[I] ){
						//upward flux
					
						kn = k_from_psi(jKn, H[I] - (adt->T->Z[l+1][r][c]-adt->P->depr_channel), SC->thi[l+1][ch], SC->T[l+1][ch], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
						
					}else{	
						//downward flow
					
						kn = k_from_psi(jKn, psisat_from(SC->thi[l+1][ch], l+1, adt->S->pa, sy), SC->thi[l+1][ch], SC->T[l+1][ch], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
					}
					
				}else{	//subsurface flow
					
					dzn = adt->S->pa[sy][jdz][l+1];
					dD = 0.5*dz + 0.5*dzn;
					
					if( H[i] < H[I] ){	
						//upward flux
						kn = k_from_psi(jKn, H[I] - (adt->T->Z[l+1][r][c]-adt->P->depr_channel), SC->thi[l+1][ch], SC->T[l+1][ch], l+1, adt->S->pa,sy, adt->P->imp, adt->P->k_to_ksat);	
					}else{	
						//downward flow
						kn = k_from_psi(jKn, H[i] - (adt->T->Z[l][r][c]-adt->P->depr_channel), SC->thi[l][ch], SC->T[l][ch], l, adt->S->pa,sy, adt->P->imp, adt->P->k_to_ksat);	
					}
					
				//	
					
				
					kmax = k_from_psi(jKn, psisat_from(SC->thi[l][ch], l, adt->S->pa, sy), SC->thi[l][ch], SC->T[l][ch], l, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
				
					kmaxn = k_from_psi(jKn, psisat_from(SC->thi[l+1][ch], l+1, adt->S->pa, sy), SC->thi[l+1][ch], SC->T[l+1][ch], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
					kmaxn = Fmin(kmax, kmaxn);
					
					kn = Fmin(kn, kmaxn);
					
				}
				
				cnt++;
			
				Lx[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
			}
		}
		
		//LATERAL FLUXES
//		printf("i, l ,r , c %ld,%ld,%ld,%ld\n",i,);
	
		if (i<=n){
			
			//lateral hydraulic conductivity
			if(l>0){
				k = k_from_psi(jKl, H[i] - adt->T->Z[l][r][c], SL->thi[l][j], SL->T[l][j], l, adt->S->pa,sy, adt->P->imp, adt->P->k_to_ksat);	
				kmax = k_from_psi(jKl, psisat_from(SL->thi[l][j], l, adt->S->pa,sy), SL->thi[l][l], SL->T[l][j], l, adt->S->pa,sy, adt->P->imp, adt->P->k_to_ksat);	
				if(adt->T->BC_counter[r][c]>0){
					if(adt->T->pixel_type[r][c] == 1 || adt->T->pixel_type[r][c] == 2 || adt->T->pixel_type[r][c] == 11 || adt->T->pixel_type[r][c] == 12) Klat[adt->T->BC_counter[r][c]][l] = k;
				}
			}
			
			//4 neighbouring cells
			R = r-1;
			C = c;
			//1.
			if(R>=1 && R<=geotop::common::Variables::Nr && C>=1 && C<=geotop::common::Variables::Nc){
				if((long)adt->L->LC[R][C]!=geotop::input::gDoubleNoValue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0[r][c]-adt->T->Z0[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE[r][c])+0.5*atan(adt->T->dzdE[R][C]));//[m]
					
					if(l>0){
						
						//Subsurface Flow
						if (H[I] > H[i]) {
							kn = k_from_psi(jKl, H[I] - adt->T->Z[l][R][C], SL->thi[l][J], SL->T[l][J], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi[l][J], l, adt->S->pa,syn), SL->thi[l][J], SL->T[l][J], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						//Surface Flow
						//if (H[I] > H[i]) {
						//	kn = adt->L->ty[(long)adt->L->LC[R][C]][jcm];
						//	dz = Fmax(0., H[I] - adt->T->Z[l][R][C]) / cos(adt->T->slope[R][C]*GTConst::Pi/180.);
						//}else {
						//	kn = adt->L->ty[(long)adt->L->LC[R][C]][jcm];
						//	dz = Fmax(0., H[i] - adt->T->Z[l][r][c]) / cos(adt->T->slope[r][c]*GTConst::Pi/180.);
						//}
						
						//if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx[cnt] = 0. ; // -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			R = r+1;
			C = c;

			//2.
			if(R>=1 && R<=geotop::common::Variables::Nr && C>=1 && C<=geotop::common::Variables::Nc){
				if((long)adt->L->LC[R][C]!=geotop::input::gDoubleNoValue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0[r][c]-adt->T->Z0[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE[r][c])+0.5*atan(adt->T->dzdE[R][C]));//[m]
					
					if(l>0){
						if (H[I] > H[i]) {
							kn = k_from_psi(jKl, H[I] - adt->T->Z[l][R][C], SL->thi[l][J], SL->T[l][J], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi[l][J], l, adt->S->pa,syn), SL->thi[l][J], SL->T[l][J], l, adt->S->pa, syn, adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						//if (H[I] > H[i]) {
						//	kn = adt->L->ty[(long)adt->L->LC[R][C]][jcm];
						//	dz = Fmax(0., H[I] - adt->T->Z[l][R][C]) / cos(adt->T->slope[R][C]*GTConst::Pi/180.);
						//}else {
						//	kn = adt->L->ty[(long)adt->L->LC[R][C]][jcm];
						//	dz = Fmax(0., H[i] - adt->T->Z[l][r][c]) / cos(adt->T->slope[r][c]*GTConst::Pi/180.);
						//}
						
						//if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx[cnt] = 0. ;//-(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
					
				}
			}
			
			R = r;
			C = c-1;

			//3.
			if(R>=1 && R<=geotop::common::Variables::Nr && C>=1 && C<=geotop::common::Variables::Nc){
				if((long)adt->L->LC[R][C]!=geotop::input::gDoubleNoValue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0[r][c]-adt->T->Z0[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE[r][c])+0.5*atan(adt->T->dzdE[R][C]));//[m]
					
					if(l>0){
						if (H[I] > H[i]) {
							kn = k_from_psi(jKl, H[I] - adt->T->Z[l][R][C], SL->thi[l][J], SL->T[l][J], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi[l][J], l, adt->S->pa,syn), SL->thi[l][J], SL->T[l][J], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						//if (H[I] > H[i]) {
						//	kn = adt->L->ty[(long)adt->L->LC[R][C]][jcm];
						//	dz = Fmax(0., H[I] - adt->T->Z[l][R][C]) / cos(adt->T->slope[R][C]*GTConst::Pi/180.);
						//}else {
						//	kn = adt->L->ty[(long)adt->L->LC[R][C]][jcm];
						//	dz = Fmax(0., H[i] - adt->T->Z[l][r][c]) / cos(adt->T->slope[r][c]*GTConst::Pi/180.);
						//}
						
						//if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx[cnt] = 0. ; //-(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			R = r;
			C = c+1;
			//4.
			if(R>=1 && R<=geotop::common::Variables::Nr && C>=1 && C<=geotop::common::Variables::Nc){
				if((long)adt->L->LC[R][C]!=geotop::input::gDoubleNoValue && adt->T->i_cont[l][R][C]>i){ 
					
					I = adt->T->i_cont[l][R][C];	
					syn = adt->S->type[R][C];
					J = adt->T->j_cont[R][C];
					
					dD = find_3Ddistance(ds, adt->T->Z0[r][c]-adt->T->Z0[R][C]) * 1.E3;//[mm]
					dn = ds/cos(0.5*atan(adt->T->dzdE[r][c])+0.5*atan(adt->T->dzdE[R][C]));//[m]
					
					if(l>0){
						if (H[I] > H[i]) {
							kn = k_from_psi(jKl, H[I] - adt->T->Z[l][R][C], SL->thi[l][J], SL->T[l][J], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
						}else {
							kn = k;
						}
						
						kmaxn = k_from_psi(jKl, psisat_from(SL->thi[l][J], l, adt->S->pa,syn), SL->thi[l][J], SL->T[l][J], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
						kmaxn = Fmin(kmax, kmaxn);
						kn = Fmin(kn, kmaxn);
						
						cnt++;
						Lx[cnt] = -(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}else {
						
						//if (H[I] > H[i]) {
						//	kn = adt->L->ty[(long)adt->L->LC[R][C]][jcm];
						//	dz = Fmax(0., H[I] - adt->T->Z[l][R][C]) / cos(adt->T->slope[R][C]*GTConst::Pi/180.);
						//}else {
						//	kn = adt->L->ty[(long)adt->L->LC[R][C]][jcm];
						//	dz = Fmax(0., H[i] - adt->T->Z[l][r][c]) / cos(adt->T->slope[r][c]*GTConst::Pi/180.);
						//}
						//
						//if (dz < adt->P->thres_hsup_1) kn = 0.;
						
						cnt++;
						Lx[cnt] = 0.0 ;//-(dn*1.E-3*dz)*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
						
					}
				}
			}
			
			//exchange with channels
			if (l>0 && adt->C->ch[r][c] > 0) {
				
				ch = adt->C->ch[r][c];
				syn = adt->C->soil_type[ch];
				I = n + adt->C->ch3[l][ch];
				
				if (H[I] > H[i]) {
					kn = k_from_psi(jKl, H[I] - (adt->T->Z[l][r][c]-adt->P->depr_channel), SC->thi[l][ch], SC->T[l][ch], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
				}else {
					kn = k;
				}
				
				kmaxn = k_from_psi(jKl, psisat_from(SC->thi[l][ch], l, adt->S->pa,syn), SC->thi[l][ch], SC->T[l][ch], l, adt->S->pa,syn, adt->P->imp, adt->P->k_to_ksat);	
				kmaxn = Fmin(kmax, kmaxn);
				kn = Fmin(kn, kmaxn);
				
				//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]	
				//Area[m2] = 2 * channel length * layer thickness = 2 * length[m] * dz[mm] * 1.E-3[m/mm]
				//dD[mm] = 0.25 * ds * (1+w_dx)
				dD = find_3Ddistance(ds * (1.0 + adt->P->w_dx) / 4.0, 1.E-3*adt->P->depr_channel) * 1.E3;//[mm]
				
				cnt++;
				Lx[cnt] = -(2.*adt->C->length[adt->C->ch[r][c]]*1.E-3*dz)*kn/dD;						
				
			}
			
			cout << "find3D: Lx[cnt]=" << Lx[cnt] << endl;

		}
		
	}									
	
	
	return 0;
	
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int find_matrix_K_1D(long c, double Dt, SoilState *L, GeoVector<double>& Lx, GeoMatrix<double>& Klat, GeoMatrix<double>& Kbottom, AllData *adt, const GeoVector<double>& H)
	{
	long i, l, r=1, I, sy, cnt=0;
	double dz=0.0, dzn=0.0, dD=0.0, k, kn=0.0, kmax=0.0, kmaxn=0.0;
	double area, ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	double psi, ice, a, ns, res, sat, ss, Temp;
	
	for(i=1;i<H.size();i++){
		
		//VERTICAL FLUXES
		l = i-1;	
	//	sy=adt->S->type->co[r][c];
		sy=adt->S->type[r][c];
		area=ds*ds;
		
		//vertical hydraulic conductivity
		if(l>0){
			dz = adt->S->pa[sy][jdz][l];
		//	if (l==Nl && adt->P->free_drainage_bottom>0) Kbottom->co[r][c] = k_from_psi(jKn, H->co[i] - adt->T->Z->co[l][r][c], L->thi->co[l][c], L->T->co[l][c], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);
			if (l==geotop::common::Variables::Nl && adt->P->free_drainage_bottom>0) Kbottom[r][c] = k_from_psi(jKn, H[i] - adt->T->Z[l][r][c], L->thi[l][c], L->T[l][c], l, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
		}
		
		//flux from cell below
		if (l<geotop::common::Variables::Nl) {
			
			I = i+1;
			
			if(l==0){	//overland flow
				
				dzn = adt->S->pa[sy][jdz][l+1];
				dD = 0.5*dzn;
				
			//	if( H->co[i] < H->co[I] ){
				if( H[i] < H[I] ){
					//upward flux
				//	kn = k_from_psi(jKn, H->co[I] - adt->T->Z->co[l+1][r][c], L->thi->co[l+1][c], L->T->co[l+1][c], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);
					kn = k_from_psi(jKn, H[I] - adt->T->Z[l+1][r][c], L->thi[l+1][c], L->T[l+1][c], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
					
				}else{	
					//downward flow
				//	kn = k_from_psi(jKn, psisat_from(L->thi->co[l+1][c], l+1, adt->S->pa->co[sy]), L->thi->co[l+1][c], L->T->co[l+1][c], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);
					kn = k_from_psi(jKn, psisat_from(L->thi[l+1][c], l+1, adt->S->pa, sy), L->thi[l+1][c], L->T[l+1][c], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
 				 }
			}else{	//subsurface flow
				
				dzn = adt->S->pa[sy][jdz][l+1];
				dD = 0.5*dz + 0.5*dzn;
				
			//	psi = Arithmetic_Mean(dz, dzn, H->co[i] - adt->T->Z->co[l][r][c], H->co[I] - adt->T->Z->co[l+1][r][c]);
				psi = Arithmetic_Mean(dz, dzn, H[i] - adt->T->Z[l][r][c], H[I] - adt->T->Z[l+1][r][c]);
			//	ice = Arithmetic_Mean(dz, dzn, L->thi->co[l][c], L->thi->co[l+1][c]);
				ice = Arithmetic_Mean(dz, dzn, L->thi[l][c], L->thi[l+1][c]);
			//	Temp = Arithmetic_Mean(dz, dzn, L->T->co[l][c], L->T->co[l+1][c]);
				Temp = Arithmetic_Mean(dz, dzn, L->T[l][c], L->T[l+1][c]);
				a = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][ja][l], adt->S->pa[sy][ja][l+1]);
				ns = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jns][l], adt->S->pa[sy][jns][l+1]);
				res = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jres][l], adt->S->pa[sy][jres][l+1]);
				sat = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jsat][l], adt->S->pa[sy][jsat][l+1]);
				ss = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jss][l], adt->S->pa[sy][jss][l+1]);
				k = Arithmetic_Mean(dz, dzn, adt->S->pa[sy][jKn][l], adt->S->pa[sy][jKn][l+1]);
				kn = k_hydr_soil(psi, k, adt->P->imp, ice, sat, res, a, ns, 1.-1./ns, 0.5, Temp, adt->P->k_to_ksat);
				
			//	kmax = k_from_psi( jKn,  psisat_from( L->thi->co[l][c], l, adt->S->pa->co[sy]), L->thi->co[l][c], L->T->co[l][c], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);
				kmax = k_from_psi( jKn,  psisat_from( L->thi[l][c], l, adt->S->pa, sy), L->thi[l][c], L->T[l][c], l, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
			//	kmaxn = k_from_psi( jKn,  psisat_from( L->thi->co[l+1][c], l+1, adt->S->pa->co[sy]), L->thi->co[l+1][c], L->T->co[l+1][c], l+1, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);
				kmaxn = k_from_psi( jKn,  psisat_from( L->thi[l+1][c], l+1, adt->S->pa, sy), L->thi[l+1][c], L->T[l+1][c], l+1, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
				kmaxn = Fmin(kmax, kmaxn);
				
				kn = Fmin(kn, kmaxn);
				
			}
			
			cnt++;
		//	Lx->co[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
			Lx[cnt] = -area*kn/dD;	//Area[m2] * k[mm/s] * dH[mm]/dD[mm], equation written in [m2*mm/s]
		//	printf("area:%e kn:%e dD:%e psi:%e k:%e Temp:%e \n",area,kn,dD,psi,k,Temp);
		//	printf("Lx %ld %e %e %ld %ld\n",cnt,-area*kn/dD,Lx->co[cnt],r,c);
		}
		
		//LATERAL FLUXES
	//	if (adt->T->pixel_type->co[r][c] == 1){
		if (adt->T->pixel_type[r][c] == 1){
			
			//lateral hydraulic conductivity
		//	if (l>0) Klat->co[adt->T->BC_counter->co[r][c]][l] = k_from_psi(jKl, H->co[i] - adt->T->Z->co[l][r][c], L->thi->co[l][c], L->T->co[l][c], l, adt->S->pa->co[sy], adt->P->imp, adt->P->k_to_ksat);
			if (l>0) Klat[adt->T->BC_counter[r][c]][l] = k_from_psi(jKl, H[i] - adt->T->Z[l][r][c], L->thi[l][c], L->T[l][c], l, adt->S->pa, sy, adt->P->imp, adt->P->k_to_ksat);
			
		}
		
	}
	
	return 0;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//int find_dfdH_3D(double Dt, DOUBLEVECTOR *df, ALLDATA *adt, SOIL_STATE *L, SOIL_STATE *C, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat){

int find_dfdH_3D(double Dt,GeoVector<double>& df, AllData *adt, SoilState *L, SoilState *C, const GeoVector<double>& H, GeoMatrix<double>& Klat){
	
	long i, l, r, c, j, sy, ch, bc;
	long n=(geotop::common::Variables::Nl+1)*adt->P->total_pixel;
	double dz, dn, dD, psi1, ice=0.0;
	double area, ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	
	for(i=1;i<H.size();i++){
		
		df[i] = 0.;
		
		if (i<=n) {
			l=adt->T->lrc_cont[i][1];
			r=adt->T->lrc_cont[i][2];
			c=adt->T->lrc_cont[i][3];
			j=adt->T->j_cont[r][c];
			sy=adt->S->type[r][c];
			bc=adt->T->BC_counter[r][c];
			ch=adt->C->ch[r][c];
			area=ds*ds/cos(adt->T->slope[r][c]*GTConst::Pi/180.);
			if (ch>0) area-=adt->C->length[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			psi1 = H[i] - adt->T->Z[l][r][c];
			if(l>0) ice = L->thi[l][j];
			
		}else {
			l=adt->C->lch[i-n][1];
			ch=adt->C->lch[i-n][2];
			r=adt->C->r[ch];
			c=adt->C->c[ch];
			sy=adt->C->soil_type[ch];
			bc=0;
			area=adt->C->length[ch] * adt->P->w_dx * ds;
			psi1 = H[i] - (adt->T->Z[l][r][c] - adt->P->depr_channel);
			if(l>0) ice = C->thi[l][ch];
		}
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
			if(psi1>0) df[i] += ( area / cos(adt->T->slope[r][c]*GTConst::Pi/180.) ) / Dt;
			
		}else{
			dz = adt->S->pa[sy][jdz][l];
			df[i] += dtheta_dpsi_from_psi(psi1, ice, l, adt->S->pa,sy, GTConst::PsiMin) * area * dz / Dt;
			
		}
		
		//lateral drainage at the border
		if (bc>0){
			
			if (l>0) {
				if (adt->T->pixel_type[r][c] == 1 || adt->T->pixel_type[r][c] == 11) {
					if ( adt->T->Z[0][r][c] - adt->T->Z[l][r][c] <= adt->T->BC_DepthFreeSurface[bc]*cos(adt->T->slope[r][c]*GTConst::Pi/180.) && H[i] - adt->T->Z[l][r][c] > 0 ) {
						if ((long)adt->L->LC[r+1][c]==geotop::input::gDoubleNoValue || (long)adt->L->LC[r-1][c]==geotop::input::gDoubleNoValue) {
							dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdN[r][c]));//mm
							dn = ds / cos(atan(adt->T->dzdE[r][c]));//m
						}else {
							dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdE[r][c]));
							dn = ds / cos(atan(adt->T->dzdN[r][c]));//m
						}
						
						dz = adt->S->pa[sy][jdz][l];//[mm]
						df[i] += (dn*dz*1.E-3) * Klat[bc][l]*adt->P->free_drainage_lateral / dD;
						
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

//int find_dfdH_1D(long c, double Dt, SOIL_STATE *L, DOUBLEVECTOR *df, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat){
  int find_dfdH_1D(long c, double Dt, SoilState *L, GeoVector<double>& df, AllData *adt, const GeoVector<double>& H, GeoMatrix<double>& Klat){

    long r=1, l, sy, bc;
	double dz, dn, dD, psi1, ice=0.0;
//	double area, ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	double area, ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	
//	for(i=1;i<=H->nh;i++){
	for(size_t i=1;i<H.size();i++){
		
	//	df->co[i] = 0.;
		df[i] = 0.;
		
		l = i-1;	
	//	sy=adt->S->type->co[r][c];
		sy=adt->S->type[r][c];
	//	bc=adt->T->BC_counter->co[r][c];
		bc=adt->T->BC_counter[r][c];
		area=ds*ds;
	//	psi1 = H->co[i] - adt->T->Z->co[l][r][c];
		psi1 = H[i] - adt->T->Z[l][r][c];
	//	if(l>0) ice = L->thi->co[l][c];
		if(l>0) ice = L->thi[l][c];
		
		//hydraulic capacity (diagonal term) = (dV/dH)/(Ah*Dt)
		if(l==0){
		//	if(psi1>0) df->co[i] += ( area / cos( Fmin(GTConst::max_slope,adt->T->slope->co[r][c])*GTConst::Pi/180.) ) / Dt;
			if(psi1>0) df[i] += ( area / cos( Fmin(GTConst::max_slope,adt->T->slope[r][c])*GTConst::Pi/180.) ) / Dt;
		}else{
			dz = adt->S->pa[sy][jdz][l];
//			df->co[i] += dteta_dpsi(psi1, ice, adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
//									adt->S->pa->co[sy][jns][l], 1.-1./adt->S->pa->co[sy][jns][l], GTConst::PsiMin, adt->S->pa->co[sy][jss][l])*
//			area * dz / Dt;

			df[i] += dteta_dpsi(psi1, ice, adt->S->pa[sy][jsat][l], adt->S->pa[sy][jres][l], adt->S->pa[sy][ja][l],
									adt->S->pa[sy][jns][l], 1.-1./adt->S->pa[sy][jns][l], GTConst::PsiMin, adt->S->pa[sy][jss][l])* area * dz / Dt;
		}
		
		//lateral drainage at the border
		if (bc>0 && l>0) {

		//	if (adt->T->pixel_type->co[r][c] == 1) {
			if (adt->T->pixel_type[r][c] == 1) {
			//	if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*GTConst::Pi/180.) && H->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
				if ( adt->T->Z[0][r][c] - adt->T->Z[l][r][c] <= adt->T->BC_DepthFreeSurface[bc]*cos(adt->T->slope[r][c]*GTConst::Pi/180.) && H[i] - adt->T->Z[l][r][c] > 0 ) {
					dz = adt->S->pa[sy][jdz][l];//[mm]
					dn = ds;
					dD = 0.5 * 1.E3*ds;
				//	df->co[i] += (dn*dz*1.E-3) * Klat->co[bc][l]*adt->P->free_drainage_lateral / dD;
					df[i] += (dn*dz*1.E-3) * Klat[bc][l]*adt->P->free_drainage_lateral / dD;
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

//int find_f_3D(double Dt, DOUBLEVECTOR *f, ALLDATA *adt, SOIL_STATE *L, SOIL_STATE *C, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom_l, DOUBLEVECTOR *Kbottom_ch){

int find_f_3D(double Dt, GeoVector<double>& f, AllData *adt, SoilState *L, SoilState *C, const GeoVector<double>& H, GeoMatrix<double>& Klat, GeoMatrix<double>& Kbottom_l, const GeoVector<double>& Kbottom_ch){

	
	
	
    long l, r, c, j, sy, ch, bc;
	long n=(geotop::common::Variables::Nl+1)*adt->P->total_pixel;
	double dz, dn, dD, V0, V1, psi1, psi0, ice=0.0;

	double area, ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	
	for(size_t i=1;i<H.size();i++){
		
		if (i<=n) {

			l=adt->T->lrc_cont[i][1];
			r=adt->T->lrc_cont[i][2];
			c=adt->T->lrc_cont[i][3];
			j=adt->T->j_cont[r][c];
			sy=adt->S->type[r][c];
			bc=adt->T->BC_counter[r][c];
			ch=adt->C->ch[r][c];
			area=ds*ds/cos(adt->T->slope[r][c]*GTConst::Pi/180.);
			if (ch>0) area-=adt->C->length[ch] * adt->P->w_dx * ds; //area of the pixel[m2]
			psi0 = L->P[l][j];
			psi1 = H[i] - adt->T->Z[l][r][c];
			if(l>0) ice = L->thi[l][j];
			
		}else {
			l=adt->C->lch[i-n][1];
			ch=adt->C->lch[i-n][2];
			r=adt->C->r[ch];
			c=adt->C->c[ch];
			sy=adt->C->soil_type[ch];
			bc=0;
			area=adt->C->length[ch] * adt->P->w_dx * ds;
			psi0 = C->P[l][ch];
			psi1 = H[i] - (adt->T->Z[l][r][c] - adt->P->depr_channel);
			if(l>0) ice = C->thi[l][ch];
		}
		
		//hydraulic capacity (diagonal term)
		if(l==0){
			V1 = area * Fmax(0.0, psi1) / cos(adt->T->slope[r][c]*GTConst::Pi/180.);
			V0 = area * Fmax(0.0, psi0) / cos(adt->T->slope[r][c]*GTConst::Pi/180.);
		}else{
			dz = adt->S->pa[sy][jdz][l];		
			V1 = area*dz * theta_from_psi(psi1, ice, l, adt->S->pa, sy, GTConst::PsiMin);
			V0 = area*dz * theta_from_psi(psi0, ice, l, adt->S->pa, sy, GTConst::PsiMin);
		}
		
		f[i] = (V1-V0)/Dt;
	
		//drainage at the bottom
		if (l==geotop::common::Variables::Nl){
			if (i<=n) {
				f[i] += area*Kbottom_l[r][c];
			}else {
				f[i] += area*Kbottom_ch[ch];
			}
		}
		//lateral drainage at the border
		if (bc>0) {
			
			if (l>0) {
				if (adt->T->pixel_type[r][c] == 1 || adt->T->pixel_type[r][c] == 11){
					if (adt->T->Z[0][r][c] - adt->T->Z[l][r][c] <= adt->T->BC_DepthFreeSurface[bc]*cos(adt->T->slope[r][c]*GTConst::Pi/180.) && H[i] - adt->T->Z[l][r][c] > 0 ){
						if ((long)adt->L->LC[r+1][c]==geotop::input::gDoubleNoValue || (long)adt->L->LC[r-1][c]==geotop::input::gDoubleNoValue) {
						
							dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdN[r][c]));//mm
							dn = ds / cos(atan(adt->T->dzdE[r][c]));//m
						}else {
							dD = 0.5 * 1.E3*ds / cos(atan(adt->T->dzdE[r][c]));
							dn = ds / cos(atan(adt->T->dzdN[r][c]));//m
						}
						
						dz = adt->S->pa[sy][jdz][l];//[mm]
	
						f[i] += (dn*dz*1.E-3) * Klat[bc][l]*adt->P->free_drainage_lateral*(H[i] - adt->T->Z[l][r][c]) / dD;
					}
					
				}else if (adt->T->pixel_type[r][c] == 2 || adt->T->pixel_type[r][c] == 12) {
				
					if ( adt->T->Z[0][r][c] - adt->T->Z[l][r][c] <= adt->T->BC_DepthFreeSurface[bc]*cos(adt->T->slope[r][c]*GTConst::Pi/180.) ) {
			
						if ((long)adt->L->LC[r+1][c]==geotop::input::gDoubleNoValue || (long)adt->L->LC[r-1][c]==geotop::input::gDoubleNoValue) {
							dn = ds / cos(atan(adt->T->dzdE[r][c]));//m
						}else {
					
							dn = ds / cos(atan(adt->T->dzdN[r][c]));//m
						}
						
						dz = adt->S->pa[sy][jdz][l];//[mm]
						f[i] += (dn*dz*1.E-3) * Klat[bc][l]*adt->P->free_drainage_lateral;
					}
				}
				
			}
				
		}
		//evaporation and precipitation
		if(l>0){
			if(i<=n){
				f[i] += area*adt->S->ET[l][r][c]/Dt;
			}else {
				ch=adt->C->ch[r][c];
				f[i] += area*adt->C->ET[l][ch]/Dt;
			}
		}else {
					f[i] -= area*adt->W->Pnet[r][c]/Dt;
		}
 // 	printf("ff[%ld]=,%f \n",i,f[i]);	
	}
	
	return 0;
	
/// this is ok... SC 14.12.2013... 

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//int find_f_1D(long c, double Dt, SOIL_STATE *L, DOUBLEVECTOR *f, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom){
int find_f_1D(long c, double Dt, SoilState *L, GeoVector<double>& f, AllData *adt, const GeoVector<double>& H, GeoMatrix<double>& Klat, GeoMatrix<double>& Kbottom){
	
    long l, r=1, sy, bc;
	double dz, dn, dD, V0, V1, psi1, psi0, ice=0.0;
	double area, ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	
//	for(i=1;i<=H->nh;i++){
	for(size_t i=1;i<H.size();i++){
		
		l = i-1;	
	//	sy=adt->S->type->co[r][c];
		sy=adt->S->type[r][c];
	//	bc=adt->T->BC_counter->co[r][c];
		bc=adt->T->BC_counter[r][c];
		area=ds*ds;
	//	psi0 = L->P->co[l][c];
		psi0 = L->P[l][c];
	//	psi1 = H->co[i] - adt->T->Z->co[l][r][c];
		psi1 = H[i] - adt->T->Z[l][r][c];
	//	if(l>0) ice = L->thi->co[l][c];
		if(l>0) ice = L->thi[l][c];
		
	//	hydraulic capacity (diagonal term)
		if(l==0){
		//	V1 = area * Fmax(0.0, psi1) / cos( Fmin(GTConst::max_slope,adt->T->slope->co[r][c])*GTConst::Pi/180.);
			V1 = area * Fmax(0.0, psi1) / cos( Fmin(GTConst::max_slope,adt->T->slope[r][c])*GTConst::Pi/180.);
		//	V0 = area * Fmax(0.0, psi0) / cos( Fmin(GTConst::max_slope,adt->T->slope->co[r][c])*GTConst::Pi/180.);
			V0 = area * Fmax(0.0, psi0) / cos( Fmin(GTConst::max_slope,adt->T->slope[r][c])*GTConst::Pi/180.);
			
		}else{
			dz = adt->S->pa[sy][jdz][l];					
			V1 = area*dz * theta_from_psi(psi1, ice, l, adt->S->pa, sy, GTConst::PsiMin);
			V0 = area*dz * theta_from_psi(psi0, ice, l, adt->S->pa, sy, GTConst::PsiMin);
		}
		
	//	f->co[i] = (V1-V0)/Dt;
		f[i] = (V1-V0)/Dt;

		// drainage at the bottom
		if (l==geotop::common::Variables::Nl){
		//	f->co[i] += area*Kbottom->co[r][c];
			f[i] += area*Kbottom[r][c];
		}		
		
		//lateral drainage at the border
		if (bc>0 && l>0) {
		//	if (adt->T->pixel_type->co[r][c] == 1) {
			if (adt->T->pixel_type[r][c] == 1) {
			//	if ( adt->T->Z->co[0][r][c] - adt->T->Z->co[l][r][c] <= adt->T->BC_DepthFreeSurface->co[bc]*cos(adt->T->slope->co[r][c]*GTConst::Pi/180.) && H->co[i] - adt->T->Z->co[l][r][c] > 0 ) {
				if ( adt->T->Z[0][r][c] - adt->T->Z[l][r][c] <= adt->T->BC_DepthFreeSurface[bc]*cos(adt->T->slope[r][c]*GTConst::Pi/180.) && H[i] - adt->T->Z[l][r][c] > 0 ) {
					dz = adt->S->pa[sy][jdz][l];//[mm]
					dn = ds;
					dD = 0.5 * 1.E3*ds;
				//	f->co[i] += (dn*dz*1.E-3) * Klat->co[bc][l]*adt->P->free_drainage_lateral*(H->co[i] - adt->T->Z->co[l][r][c]) / dD;
					f[i] += (dn*dz*1.E-3) * Klat[bc][l]*adt->P->free_drainage_lateral*(H[i] - adt->T->Z[l][r][c]) / dD;
				}
			}
		}
		
		
		//evaporation 
		if(l>0){
		//	f->co[i] += area*adt->S->ET->co[l][r][c]/Dt;
			f[i] += area*adt->S->ET[l][r][c]/Dt;
		}else {
		//	f->co[i] -= area*adt->W->Pnet->co[r][c]/Dt;
			f[i] -= area*adt->W->Pnet[r][c]/Dt;
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

//void find_dt_max(double Courant, double *h, LAND *land, TOPO *top, CHANNEL *cnet, PAR *par, METEO *met, double t, double *dt){


  void find_dt_max(double Courant, GeoMatrix<double>& h, Land *land, Topo *top, Channel *cnet, Par *par, Meteo *met, double t, double *dt){

    // t parameter is NOT used ??? 
	double q, ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]), area, Vmax, H;
	short d;
	long r, c, j, ch;
	
	for (j=1; j<=par->total_pixel; j++) {
		
		r = top->rc_cont[j][1];
		c = top->rc_cont[j][2];
		

		H = Fmax(0.0 , h(0,j)) / cos(top->slope[r][c]*GTConst::Pi/180.); //h[i] is the pressure at the surface, H is the depth of water normal to the surface
		
		if(H > par->min_hsup_land){
			
		//  to fix the following: (SC 10.12.2013)
		//	if(DD==1){
		//		draining_land(1., j, top, land, par, cnet, h, top->Jdown->co[j], top->Qdown->co[j]);
		//	}else {
		//		draining_land(0., j, top, land, par, cnet, h, top->Jdown->co[j], top->Qdown->co[j]);
		//	}
			
			draining_land(1., j, top, land, par, h, top->Jdown, top->Qdown, j);
			
		
			area = ds*ds/cos(top->slope[r][c]*GTConst::Pi/180.);
			
			ch = cnet->ch[r][c];
		
			if (ch>0) area -= cnet->length[ch] * par->w_dx * ds; //area of the pixel[m2]
			
			Vmax = area * 1.E-3 * H;
			
			q = 0.;
			for (d=1; d<=4; d++) {
			
				q += top->Qdown[j][d];	//outgoing discharge
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

//void supflow(double Dt, double t, double *h, double *dV, double *hch, double *dhch, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet,
//			 PAR *par, METEO *met, DOUBLEVECTOR *Vsup, double *Voutnet, double *Voutland, FILE *flog)

void supflow(double Dt, double t, GeoMatrix<double>& h, double *dV, GeoMatrix<double>& hch, double *dhch, Topo *top, Land *land, Water *wat, Channel *cnet,
			 Par *par, Meteo *met, GeoVector<double>& Vsup, double *Voutnet, double *Voutland, FILE *flog, double *mm1, double *mm2, double *mmo )

{
	
	long d, r, c, j, R, C, ch;                                    
	double ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]), area, Vmax, H;
	double q, q0, tb, te=0.0, dt;
	long cnt=0,cnt2=0,cnt3=0;
	
	double m1=0.,m2=0.,mo=0.;
	

	for (j=1; j<=par->total_pixel; j++) {
		r = top->rc_cont[j][1];
		c = top->rc_cont[j][2];
		H = Fmax(0., h(0,j)) / cos(top->slope[r][c]*GTConst::Pi/180.);
		area = ds*ds;
		area /= cos(top->slope[r][c]*GTConst::Pi/180.);
		m1 += H*1.E-3*area;
	}
	
	printf("supflow: te %ld\n",te);
	
	do{
		
		tb=te;
		dt=Dt;
		
		find_dt_max(par->max_courant_land, h, land, top, cnet, par, met, t, &dt);
		cnt++;
		
		te=tb+dt;
		if(te>Dt){
			te=Dt;
			dt=te-tb;
		}		
		
		for (j=1; j<=par->total_pixel; j++) {
			r = top->rc_cont[j][1];	
			c = top->rc_cont[j][2];
			
			H = Fmax(0., h(0,j)) / cos(top->slope[r][c]*GTConst::Pi/180.);
			
			dV[j] = 0.0;
			
			if(H > 0){
				

				area = ds*ds/cos(top->slope[r][c]*GTConst::Pi/180.);
				ch = cnet->ch[r][c];
				if (ch>0) area -= cnet->length[ch] * par->w_dx * ds; //area of the pixel[m2]
				
				Vmax = area * 1.E-3 * H;
				
				q = 0.;
				for (d=1; d<=4; d++) {
					q += top->Qdown[j][d];
				}
				
				if (q*dt > Vmax){
					q0 = q;
					q = Vmax/dt;
					for (d=1; d<=4; d++) {
						top->Qdown[j][d] *= (q/q0);
					}
				}
				
				for (d=1; d<=4; d++) {
					if(top->Jdown[j][d]==0 ){
					    if (top->BC_counter[r][c] > 0 && (top->pixel_type[r][c] == 1 || top->pixel_type[r][c] == 2 || top->pixel_type[r][c] == 11 || top->pixel_type[r][c] == 12)){
					
						*Voutland = *Voutland + top->Qdown[j][d]*dt;
						mo += top->Qdown[j][d]*dt;
					    }
					}
					if (q > 0) top->Qdown[j][d] /= q;
				}
				
				dV[j] = q*dt;
				
			}
		}
		
		for (j=1; j<=par->total_pixel; j++) {	
			r = top->rc_cont[j][1];
			c = top->rc_cont[j][2];
			
			area = ds*ds/cos(top->slope[r][c]*GTConst::Pi/180.);
			ch = cnet->ch[r][c];
			if (ch>0) area -= cnet->length[ch] * par->w_dx * ds; //area of the pixel[m2]
			
		
			h(0,j) -= (1.E3*dV[j]/area) * cos(top->slope[r][c]*GTConst::Pi/180.);
			
		 
			if (top->BC_counter[r][c] > 0 && top->pixel_type[r][c] == -1) {

				if(h(0,j) > 0){
									
					h(0,j) += (1.E3*met->qinv[1]*ds*Dt/area) * cos(top->slope[r][c]*GTConst::Pi/180.);
				}else{
					if( met->qinv[1]*ds > 0) h(0,j) = (1.E3*met->qinv[1]*ds*Dt/area) * cos(top->slope[r][c]*GTConst::Pi/180.);
				}
			}
			
			for (d=1; d<=4; d++) {
		
				if(top->Jdown[j][d]>0){
					
				
					R = top->rc_cont[top->Jdown[j][d]][1];
					C = top->rc_cont[top->Jdown[j][d]][2];
								
					area = ds*ds/cos(top->slope[R][C]*GTConst::Pi/180.);
					ch = cnet->ch[R][C];
					
					if (ch>0) area -= cnet->length[ch] * par->w_dx * ds; //area of the pixel[m2]

				
					if(h(0,top->Jdown[j][d])>0){
					
						h(0,top->Jdown[j][d]) += (1.E3*dV[j]*top->Qdown[j][d]/area) * cos(top->slope[R][C]*GTConst::Pi/180.);
					}else{
						if( dV[j]*top->Qdown[j][d] > 0) h(0,top->Jdown[j][d]) = (1.E3*dV[j]*top->Qdown[j][d]/area) * cos(top->slope[R][C]*GTConst::Pi/180.);
					}
					
				}
			}
		}			
		
		supflow_chla(dt, t, h, hch, top, wat, cnet, par, Vsup, flog, &cnt2);
		channel_flow(dt, t, 1, hch, dhch, top, cnet, par, land, Voutnet, flog, &cnt3);
		
	}while(te<Dt);

	
	// this below control lines... to be removed later
	
    printf("sup-flow: %f %f %f\n",Dt/cnt,Dt/cnt2,Dt/cnt3);
	
	for (j=1; j<=par->total_pixel; j++) {
		r = top->rc_cont[j][1];
		c = top->rc_cont[j][2];
		H = Fmax(0., h(0,j)) / cos(top->slope[r][c]*GTConst::Pi/180.);
		area = ds*ds;
		area /= cos(top->slope[r][c]*GTConst::Pi/180.);
		m2 += H*1.E-3*area;
	}	
	
	*mm1 = m1;
	*mm2 = m2;
	*mmo = mo;
	
	printf("SUP: m1:%e m2:%e mo:%e dm:%e\n",m1,m2,mo,fabs(m1-mo-m2));
	
	
	
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void find_dt_max_chla(double Courant, double *h, double *hch, TOPO *top, CHANNEL *cnet, PAR *par, double t, double *dt){


  void find_dt_max_chla(double Courant, GeoMatrix<double>& h, GeoMatrix<double>& hch, Topo *top, Channel *cnet, Par *par, double t, double *dt){
	
	double q, ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]), area, areach, Vmax, H, Hch, DH;
	long r, c, ch;
	
	for (ch=1; ch<=par->total_channel; ch++) {
		
		r = cnet->r[ch];
		c = cnet->c[ch];
		H = Fmax(0.0 , h(0,ch)) / cos(top->slope[r][c]*GTConst::Pi/180.);//h[i] is the pressure at the surface, H is the depth of water normal to the surface
		
	//	area = ds*ds/cos(top->slope->co[r][c]*GTConst::Pi/180.) - cnet->length->co[ch] * par->w_dx * ds;
		area = ds*ds/cos(top->slope[r][c]*GTConst::Pi/180.) - cnet->length[ch] * par->w_dx * ds;
		
	//	Hch = Fmax(0., hch[ch] ) / cos(top->slope->co[r][c]*GTConst::Pi/180.) - par->depr_channel * cos(top->slope->co[r][c]*GTConst::Pi/180.);
		Hch = Fmax(0., hch(0,ch) ) / cos(top->slope[r][c]*GTConst::Pi/180.) - par->depr_channel * cos(top->slope[r][c]*GTConst::Pi/180.);
	//	areach = cnet->length->co[ch] * par->w_dx * ds;
		areach = cnet->length[ch] * par->w_dx * ds;
		
		if(H > par->min_hsup_land){
			
			if( Hch < -par->min_dhsup_land_channel_in ){	//free flow
				
				DH = H;
			//	q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
				q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length[ch])*1.E-3*H;//m3/s

				Vmax = Fmin( areach*1.E-3*(-Hch) , area*1.E-3*H );
				
				Vmax = Fmax(Vmax, 1.E-10);
				q = Fmax(q, 1.E-30);
				
				if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
				
			}else if( H - Hch > par->min_dhsup_land_channel_in ){//submerged flow towards channel
				
				DH = H - Hch;
			//	q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
				q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length[ch])*1.E-3*H;//m3/s
				Vmax = 1.E-3*DH / (1./area + 1./areach);
				
				Vmax = Fmax(Vmax, 1.E-10);
				q = Fmax(q, 1.E-30);
				
				if(Courant*Vmax/q < (*dt)) *dt = Courant*Vmax/q; 
				
			}
		}
		
		if ( Hch - H > par->min_dhsup_land_channel_out ) {
			
			DH = Hch - H;
		//	q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*Hch;//m3/s
			q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length[ch])*1.E-3*Hch;//m3/s
			
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

//void supflow_chla(double Dt, double t, double *h, double *hch, TOPO *top, WATER *wat, CHANNEL *cnet, PAR *par, DOUBLEVECTOR *Vsup, FILE *flog, long *cnt){
// piccola divergenza con codice c in questa routine.. 

  void supflow_chla(double Dt, double t, GeoMatrix<double>& h, GeoMatrix<double>& hch, Topo *top, Water *wat, Channel *cnet, Par *par, GeoVector<double>& Vsup, FILE *flog, long *cnt){
	
	long ch, r, c;
//	double ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	double ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	double H, Hch, DH, area, areach, q, tb, te=0., dt, Vmax;
	
	do{
		
		tb=te;
		dt=Dt;
		
		find_dt_max_chla(par->max_courant_land, h, hch, top, cnet, par, t, &dt);
		*cnt = *cnt + 1;
		
		te=tb+dt;
		if(te>Dt){
			te=Dt;
			dt=te-tb;
		}		
		
		//Superficial flow to the channels
		for(ch=1;ch<=par->total_channel;ch++){
			
		//	r = cnet->r->co[ch];
			r = cnet->r[ch];
		//	c = cnet->c->co[ch];
			c = cnet->c[ch];
			
		//	H = Fmax(0., h[top->j_cont[r][c]]) / cos(top->slope->co[r][c]*GTConst::Pi/180.);
			H = Fmax(0., h(0,top->j_cont[r][c])) / cos(top->slope[r][c]*GTConst::Pi/180.);
		//	Hch = Fmax(0., hch[ch] ) / cos(top->slope->co[r][c]*GTConst::Pi/180.) - par->depr_channel * cos(top->slope->co[r][c]*GTConst::Pi/180.);
			Hch = Fmax(0., hch(0,ch) ) / cos(top->slope[r][c]*GTConst::Pi/180.) - par->depr_channel * cos(top->slope[r][c]*GTConst::Pi/180.);
			
		//	area = ds*ds/cos(top->slope->co[r][c]*GTConst::Pi/180.) - cnet->length->co[ch] * par->w_dx * ds;
			area = ds*ds/cos(top->slope[r][c]*GTConst::Pi/180.) - cnet->length[ch] * par->w_dx * ds;
		//	areach = cnet->length->co[ch] * par->w_dx * ds;
			areach = cnet->length[ch] * par->w_dx * ds;
			
			if( Hch < 0 ){	//free flow
				
				DH = H;
			//	q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
				q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length[ch])*1.E-3*H;//m3/s
				
				Vmax = Fmin( areach*1.E-3*(-Hch) , area*1.E-3*H );
				if (q*dt > Vmax) q = Vmax/dt;
				
			//	Vsup->co[ch] += q*dt;
				Vsup[ch] += q*dt;
				
			//	h[top->j_cont[r][c]] -= (1.E3 * q*dt/area) * cos(top->slope->co[r][c]*GTConst::Pi/180.);
				h(0,top->j_cont[r][c]) -= (1.E3 * q*dt/area) * cos(top->slope[r][c]*GTConst::Pi/180.);
				
			//	if(hch[ch]>0){
				if(hch(0,ch)>0){
				//	hch[ch] += (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*GTConst::Pi/180.);	//mm;
					hch(0,ch) += (1.E3 * q*dt/areach) * cos(top->slope[r][c]*GTConst::Pi/180.);	//mm;
				}else{
				//	if( q > 0 ) hch[ch] = (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*GTConst::Pi/180.);	//mm;
					if( q > 0 ) hch(0,ch) = (1.E3 * q*dt/areach) * cos(top->slope[r][c]*GTConst::Pi/180.);	//mm;
				}				
				
			}else if( H - Hch > 0 ){//submerged flow towards channel
				
				DH = H - Hch;
			//	q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*H;//m3/s
				q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length[ch])*1.E-3*H;//m3/s
				
				Vmax = 1.E-3*DH / (1./area + 1./areach);
				if (q*dt > Vmax) q = Vmax/dt;
				
			//	Vsup->co[ch] += q*dt;
				Vsup[ch] += q*dt;
				
			//	h[top->j_cont[r][c]] -= (1.E3 * q*dt/area) * cos(top->slope->co[r][c]*GTConst::Pi/180.);
				h(0,top->j_cont[r][c]) -= (1.E3 * q*dt/area) * cos(top->slope[r][c]*GTConst::Pi/180.);
				
			//	if(hch[ch]>0){
				if(hch(0,ch)>0){
				//	hch[ch] += (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*GTConst::Pi/180.);	//mm;
					hch(0,ch) += (1.E3 * q*dt/areach) * cos(top->slope[r][c]*GTConst::Pi/180.);	//mm;
				}else{
				//	if( q > 0 ) hch[ch] = (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*GTConst::Pi/180.);	//mm;
					if( q > 0 ) hch(0,ch) = (1.E3 * q*dt/areach) * cos(top->slope[r][c]*GTConst::Pi/180.);	//mm;
				}				
				
			}else if ( Hch - H > 0 ) {
				
				DH = Hch - H;
			//	q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length->co[ch])*1.E-3*Hch;//m3/s
				q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*DH)*(2.*cnet->length[ch])*1.E-3*Hch;//m3/s
				
				Vmax = 1.E-3*DH / (1./area + 1./areach);
				if (q*dt > Vmax) q = Vmax/dt;
				
			//	Vsup->co[ch] -= q*dt;
				Vsup[ch] -= q*dt;
				
			//	hch[ch] -= (1.E3 * q*dt/areach) * cos(top->slope->co[r][c]*GTConst::Pi/180.);
				hch(0,ch) -= (1.E3 * q*dt/areach) * cos(top->slope[r][c]*GTConst::Pi/180.);
				
			//	if(h[top->j_cont[r][c]]>0){
				if(h(0,top->j_cont[r][c])>0){
				//	h[top->j_cont[r][c]] += (1.E3 * q*dt/area) * cos(top->slope->co[r][c]*GTConst::Pi/180.);	//mm;
					h(0,top->j_cont[r][c]) += (1.E3 * q*dt/area) * cos(top->slope[r][c]*GTConst::Pi/180.);	//mm;
				}else{
				//	if( q > 0 ) h[top->j_cont[r][c]] = (1.E3 * q*dt/area) * cos(top->slope->co[r][c]*GTConst::Pi/180.);	//mm;
					if( q > 0 ) h(0,top->j_cont[r][c]) = (1.E3 * q*dt/area) * cos(top->slope[r][c]*GTConst::Pi/180.);	//mm;
				}
			}
		}
		
	}while(te<Dt);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void find_dt_max_channel(short DDcomplex, double Courant, double *h, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double t, double *dt){


  void find_dt_max_channel(short DDcomplex, double Courant, GeoMatrix<double>& h, Topo *top, Channel *cnet, Par *par, Land *land, double t, double *dt){
	
	long r, c, ch, R, C;		
	double Ks, q, Vmax, i, H, dn, dD, ds;
	
	ds = sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	dn = par->w_dx*ds;
	
	for(ch=1;ch<=par->total_channel;ch++){
		
	//	if(DDcomplex!=1 && t==0) draining_channel(0., ch, top->Z0, h, cnet, &(cnet->ch_down->co[ch]));
	//	if(DDcomplex!=1 && t==0) draining_channel(0., ch, top->Z0, h, cnet, &(cnet->ch_down[ch]));

		r = cnet->r[ch];
		c = cnet->c[ch];
		H = Fmax(0., h(0,ch)) / cos(top->slope[r][c]*GTConst::Pi/180.);
		
		if(H > par->min_hsup_channel){
			
				draining_channel(1., ch, top->Z0, h, cnet, &(cnet->ch_down[ch]));


	
			Vmax = 1.E-3*H*dn*cnet->length[ch];	//m3
			
			if(top->is_on_border[r][c] == 1 && cnet->ch_down[ch]==ch){//outlet section

				q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*H)*(1.E-3*H)*dn;	//[m3/s]
				
			}else{
				
			
				R = cnet->r[cnet->ch_down[ch]];
				C = cnet->c[cnet->ch_down[ch]];
				
				if( (R-r==1 || R-r==-1) && (C-c==1 || C-c==-1) ){
					dD = find_3Ddistance(ds*sqrt(2.), top->Z0[r][c] - top->Z0[R][C]);
				}else {
					dD = find_3Ddistance(ds, top->Z0[r][c] - top->Z0[R][C]);
				}
				
				Ks = cm_h(par->Ks_channel, H, 1., par->thres_hchannel);
				
				i = ( (top->Z0[r][c] - top->Z0[R][C] ) + 1.E-3*(Fmax(0.0, h(0,ch)) - Fmax(0.0, h(0,cnet->ch_down[ch]))) ) / dD;
				
				if(i<0) i=0.;
				
				q = dn * Ks * pow( 1.E-3 * H , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
				
			}
			
			Vmax = Fmax(Vmax, 1.E-10);
			q = Fmax(q, 1.E-10);
			if(Courant*Vmax/q<(*dt)) *dt=Courant*Vmax/q; 
			if(*dt<par->dtmin_sup) *dt=par->dtmin_sup;
			
		}else{
			
			cnet->ch_down[ch] = ch;
			
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void channel_flow(double Dt, double t, short DDcomplex, double *h, double *dV, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double *Vout, FILE *f, long *cnt)



  void channel_flow(double Dt, double t, short DDcomplex, GeoMatrix<double>& h, double *dV, Topo *top, Channel *cnet, Par *par, Land *land, double *Vout, FILE *f, long *cnt){

	long r,c,ch,R,C;                                    
	double ds, dn, dD;
	double Ks;											  // the Strickler's coefficent
	double i;											  // hydraulic gradient
	double q,tb,te,dt,H;
	
//	ds = sqrt(UV->U->co[1]*UV->U->co[2]);
	ds = sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	dn = par->w_dx*ds;
	
//	if( par->point_sim!=1 && cnet->r->co[1]!=0 ){	//if it is not point simulation and there are channels
	if( par->point_sim!=1 && cnet->r[1]!=0 ){	//if it is not point simulation and there are channels
		
	//	dn = par->w_dx * UV->U->co[1];		//transversal length [m]
		dn = par->w_dx * geotop::common::Variables::UV->U[1];		//transversal length [m]
		
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
			
		//	for(ch=1;ch<=cnet->r->nh;ch++){
			for(ch=1;ch<cnet->r.size();ch++){
				
			//	r = cnet->r->co[ch];
				r = cnet->r[ch];
			//	c = cnet->c->co[ch];
				c = cnet->c[ch];
				
				dV[ch] = 0.0;
				
			//	H = Fmax(0., h[ch]) / cos(top->slope->co[r][c]*GTConst::Pi/180.);
				H = Fmax(0., h(0,ch)) / cos(top->slope[r][c]*GTConst::Pi/180.);
				
				if(H > 0){
					
				//	R = cnet->r->co[cnet->ch_down->co[ch]];
					R = cnet->r[cnet->ch_down[ch]];
				//	C = cnet->c->co[cnet->ch_down->co[ch]];
					C = cnet->c[cnet->ch_down[ch]];
					
				//	if(top->is_on_border->co[r][c] == 1 && cnet->ch_down->co[ch]==ch){//outlet section
					if(top->is_on_border[r][c] == 1 && cnet->ch_down[ch]==ch){//outlet section
						
						q = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*H)*(1.E-3*H)*dn;	//[m3/s]
						
					}else{
						
						if( (R-r==1 || R-r==-1) && (C-c==1 || C-c==-1) ){
						//	dD = find_3Ddistance(ds*sqrt(2.), top->Z0->co[r][c] - top->Z0->co[R][C]);
							dD = find_3Ddistance(ds*sqrt(2.), top->Z0[r][c] - top->Z0[R][C]);
						}else {
						//	dD = find_3Ddistance(ds, top->Z0->co[r][c] - top->Z0->co[R][C]);
							dD = find_3Ddistance(ds, top->Z0[r][c] - top->Z0[R][C]);
						}
						
						Ks = cm_h(par->Ks_channel, H, 1., par->thres_hchannel);
						
					//	i= ( (top->Z0->co[r][c] - top->Z0->co[R][C] ) + 1.E-3*(Fmax(0.0, h[ch]) - Fmax(0.0, h[cnet->ch_down->co[ch]])) ) / dD;
						i= ( (top->Z0[r][c] - top->Z0[R][C] ) + 1.E-3*(Fmax(0.0, h(0,ch)) - Fmax(0.0, h(0,cnet->ch_down[ch]))) ) / dD;
						
						if(i<0) i=0.;
						
						q = dn * Ks * pow( 1.E-3 * H , 1.0+par->gamma_m ) * sqrt(i);	//m3/s
						
					}
					
				//	dV[ch] = Fmin( q*dt , 1.E-3*H*dn*cnet->length->co[ch] );	//m3
					dV[ch] = Fmin( q*dt , 1.E-3*H*dn*cnet->length[ch] );	//m3
					
				}				
			}
			
		//	for(ch=1;ch<=cnet->r->nh;ch++){
			for(ch=1;ch<cnet->r.size();ch++){
				
			//	r = cnet->r->co[ch];
				r = cnet->r[ch];
			//	c = cnet->c->co[ch];
				c = cnet->c[ch];
				
			//	h[ch] -= (1.E3*dV[ch]/(dn*cnet->length->co[ch])) * cos(top->slope->co[r][c]*GTConst::Pi/180.);
				h(0,ch) -= (1.E3*dV[ch]/(dn*cnet->length[ch])) * cos(top->slope[r][c]*GTConst::Pi/180.);
				
			//	if(top->is_on_border->co[r][c] == 1 && cnet->ch_down->co[ch]==ch){//outlet section
				if(top->is_on_border[r][c] == 1 && cnet->ch_down[ch]==ch){//outlet section
					*Vout = *Vout + dV[ch];	//m3
				}else {
				//	R = cnet->r->co[cnet->ch_down->co[ch]];
					R = cnet->r[cnet->ch_down[ch]];
				//	C = cnet->c->co[cnet->ch_down->co[ch]];
					C = cnet->c[cnet->ch_down[ch]];

				//	if(h[cnet->ch_down->co[ch]]>0){
					if(h(0,cnet->ch_down[ch])>0){
					//	h[cnet->ch_down->co[ch]] += (1.E3*dV[ch]/(dn*cnet->length->co[cnet->ch_down->co[ch]])) * cos(top->slope->co[R][C]*GTConst::Pi/180.);	//mm;
						h(0,cnet->ch_down[ch]) += (1.E3*dV[ch]/(dn*cnet->length[cnet->ch_down[ch]])) * cos(top->slope[R][C]*GTConst::Pi/180.);	//mm;
					}else{
					//	if( dV[ch] > 0) h[cnet->ch_down->co[ch]] = (1.E3*dV[ch]/(dn*cnet->length->co[cnet->ch_down->co[ch]])) * cos(top->slope->co[R][C]*GTConst::Pi/180.);	//mm;
						if( dV[ch] > 0) h(0,cnet->ch_down[ch]) = (1.E3*dV[ch]/(dn*cnet->length[cnet->ch_down[ch]])) * cos(top->slope[R][C]*GTConst::Pi/180.);	//mm;
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

//void draining_land(double alpha, long i, TOPO *T, Land *L, PAR *P, double *h, long *I, double *Q){
  void draining_land(double alpha, long i, Topo *T, Land *L, Par *P, GeoMatrix<double>& h, GeoMatrix<long>& I, GeoMatrix<double>& Q, long row){
	double H, p, pn, dD, dn, Ks;
//	double ds=sqrt(UV->U->co[1]*UV->U->co[2]);
	double ds=sqrt(geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
	
	long d, r, c;
	long ir[5] = {0, -1, 1, 0,  0};
	long ic[5] = {0,  0, 0, -1, 1};
	
//	if (h[i] > 0) {
	if (h(0,i) > 0) {
		
	//	r = T->rc_cont->co[i][1];
		r = T->rc_cont[i][1];
	//	c = T->rc_cont->co[i][2];
		c = T->rc_cont[i][2];
	//	H = Fmax(h[i], 0.)/cos(T->slope->co[r][c]*GTConst::Pi/180.);
		H = Fmax(h(0,i), 0.)/cos(T->slope[r][c]*GTConst::Pi/180.);
	//	p = T->Z0->co[r][c] + alpha*1.E-3*Fmax(h[i], 0.);
		p = T->Z0[r][c] + alpha*1.E-3*Fmax(h(0,i), 0.);
	//	Ks = cm_h(L->ty->co[(short)L->LC->co[r][c]][jcm], H, P->thres_hsup_1,  P->thres_hsup_2);
		Ks = cm_h(L->ty[(short)L->LC[r][c]][jcm], H, P->thres_hsup_1,  P->thres_hsup_2);
		
		
		for (d=1; d<=4; d++) {
			if (r+ir[d]>=1 && r+ir[d]<=geotop::common::Variables::Nr && c+ic[d]>=1 && c+ic[d]<=geotop::common::Variables::Nc) {
				
			//	I[d] = T->j_cont[r+ir[d]][c+ic[d]];
				I(row,d) = T->j_cont[r+ir[d]][c+ic[d]];
				
			//	if (I[d]>0) {
				if (I(row,d)>0) {
					
				//	dD = find_3Ddistance(ds, T->Z0->co[r][c] - T->Z0->co[r+ir[d]][c+ic[d]]);
					dD = find_3Ddistance(ds, T->Z0[r][c] - T->Z0[r+ir[d]][c+ic[d]]);
					if(ir[d]==1 || ir[d]==-1){
					//	dn = ds/cos(0.5*atan(T->dzdE->co[r][c])+0.5*atan(T->dzdE->co[r+ir[d]][c+ic[d]]));
						dn = ds/cos(0.5*atan(T->dzdE[r][c])+0.5*atan(T->dzdE[r+ir[d]][c+ic[d]]));
					}else {
					//	dn = ds/cos(0.5*atan(T->dzdN->co[r][c])+0.5*atan(T->dzdN->co[r+ir[d]][c+ic[d]]));
						dn = ds/cos(0.5*atan(T->dzdN[r][c])+0.5*atan(T->dzdN[r+ir[d]][c+ic[d]]));
					}		
					
				//	pn = T->Z0->co[r+ir[d]][c+ic[d]] + alpha*1.E-3*Fmax(h[I[d]], 0.);
					pn = T->Z0[r+ir[d]][c+ic[d]] + alpha*1.E-3*Fmax(h(0,I(row, d)), 0.);
					
					if (pn < p) {
					//	Q[d] = Ks*dn*pow(1.E-3*H, 1.0+P->gamma_m)*sqrt((p-pn)/dD);
						Q(row, d) = Ks*dn*pow(1.E-3*H, 1.0+P->gamma_m)*sqrt((p-pn)/dD);
					}else {
					//	Q[d] = 0.;
						Q(row,d) = 0.;
					}
					
				}else {
					
				//	if(T->BC_counter->co[r][c] > 0){
					if(T->BC_counter[r][c] > 0){
					//	if (H >= -T->BC_DepthFreeSurface->co[T->BC_counter->co[r][c]] ){
						if (H >= -T->BC_DepthFreeSurface[T->BC_counter[r][c]] ){
						//	I[d] = -1;
							I(row,d) = -1;
						//	Q[d] = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*H)*(1.E-3*H)*ds;
							Q(row,d) = GTConst::Cd*(2./3.)*sqrt(2.*GTConst::GRAVITY*1.E-3*H)*(1.E-3*H)*ds;
						}else {
						//	I[d] = i;
							I(row,d) = i;
						//	Q[d] = 0.;
							Q(row, d) = 0.;
						}
					}else {
					//	I[d] = i;
						I(row, d) = i;
					//	Q[d] = 0.;
						Q(row, d) = 0.;
					}
					
				}
				
			}else {
				
			//	I[d] = i;
				I(row, d) = i;
			//	Q[d] = 0.;
				Q(row, d) = 0.;
				
			}
		}
		
	}else {
		
		for (d=1; d<=4; d++) {
		//	I[d] = i;
			I(row, d) = i;
		//	Q[d] = 0.;
			Q(row, d) = 0.;
		}
	}
}



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void draining_channel(double alpha, long ch, DOUBLEMATRIX *Z, double *h, CHANNEL *cnet, long *CH){
	
void draining_channel(double alpha, long ch, GeoMatrix<double>& Z, GeoMatrix<double>& h, Channel *cnet, long *CH){
	
	long d, r, c;
	double elev, elev1;
	long ir[9] = {0, -1, -1, 0, 1, 1, 1, 0, -1};
	long ic[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
	
	*CH = ch;
	
//	r = cnet->r->co[ch];
	r = cnet->r[ch];
//	c = cnet->c->co[ch];
	c = cnet->c[ch];
//	elev = Z->co[r][c] + alpha*1.E-3*Fmax(h[ch], 0.);
	elev = Z[r][c] + alpha*1.E-3*Fmax(h(0,ch), 0.);
	
	for (d=1; d<=8; d++) {
		if (r+ir[d]>=1 && r+ir[d]<=geotop::common::Variables::Nr && c+ic[d]>=1 && c+ic[d]<=geotop::common::Variables::Nc) {
		//	if(cnet->ch->co[r+ir[d]][c+ic[d]] > 0){
			if(cnet->ch[r+ir[d]][c+ic[d]] > 0){
			//	elev1 = Z->co[r+ir[d]][c+ic[d]] + alpha*1.E-3*Fmax(h[cnet->ch->co[r+ir[d]][c+ic[d]]], 0.);
				elev1 = Z[r+ir[d]][c+ic[d]] + alpha*1.E-3*Fmax(h(0,cnet->ch[r+ir[d]][c+ic[d]]), 0.);
				if( elev1 < elev){
					elev = elev1;
				//	*CH = cnet->ch->co[r+ir[d]][c+ic[d]];
					*CH = cnet->ch[r+ir[d]][c+ic[d]];
				}
			}
		}
	}
}

