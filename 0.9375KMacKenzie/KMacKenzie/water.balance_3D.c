
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "water.balance_3D.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "lu.h"

extern T_INIT *UV;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;
extern char *error_file_name;
#define Kred_exf 1.

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void water_balance_3D(ALLDATA *adt, double time){

	Richards_3D(adt);
	supflow2(adt->T, adt->L, adt->W, adt->C, adt->P);
	routing2(adt->C);
	if(adt->P->state_pixel==1) output_waterbalance2(adt->W, adt->S, adt->P, adt->L->LC);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Richards_3D(ALLDATA *adt){

	DOUBLEVECTOR *psi, *dpsi, *b, *psi0;
	long i, l, r, c, cont, cont2, iter,diff_bc, iter_tot=0;
	short sy=0;/* modified by Emanuele Cordano on 24/9/9 */
	double mass0, mass1, massinf, massloss0, massloss,massinf2, Infmax, Infpot, nw, resnum, k, resden, tol=adt->P->max_tol_grad_conj;
	FILE *f;
	short out=0;
	long n=Nl*adt->P->total_pixel;
	static double cum_losses;

	f=fopen(error_file_name,"a");

	if(adt->I->time==0.0) cum_losses=0.0;

	psi=new_doublevector(n);
	dpsi=new_doublevector(n);
	b=new_doublevector(n);
	psi0=new_doublevector(n);
	//ap=new_doublevector(adt->T->Ai->nh);
	//V=new_doublevector(n);

	massloss=1.E99;
	mass0=0.0;
	for(i=1;i<=n;i++){
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		psi0->co[i]=adt->S->P->co[l][r][c];
		psi->co[i]=psi0->co[i];
		mass0+=adt->S->pa->co[sy][jdz][l]*teta_psi(adt->S->P->co[l][r][c],adt->S->thice->co[l][r][c],adt->S->pa->co[sy][jsat][l],adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l],adt->S->pa->co[sy][jns][l],1-1/adt->S->pa->co[sy][jns][l],adt->P->psimin, adt->P->Esoil);
	}
	mass0/=(double)adt->P->total_pixel;	//in [mm]

	cont=0;

	do{

		cont++;

		for(i=1;i<=n;i++){
			b->co[i]=Find_b(i, psi0, adt);
		}

		/*c=1;
		for(i=1;i<=adt->T->Ai->nh;i++){

			do{
				if(i>adt->T->Ax->co[c]) c++;
			}while(i>adt->T->Ax->co[c]);

			initialize_doublevector(V, 0.0);
			V->co[c]=1.0;
			ap->co[i]=Solve_Richards_3D(adt->T->Ai->co[i], V, adt);
		}

		mat_lsolve(ap, adt->T->Ax, adt->T->Ai, b, psi);*/
		iter=tridiag_preconditioned_conjugate_gradient_search(tol, psi, b, Solve_Richards_3D_p, adt);
		iter_tot+=iter;

		for(i=1;i<=n;i++){
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			dpsi->co[i]=psi->co[i]-adt->S->P->co[l][r][c];
		}

		cont2=0.0;
		nw=adt->P->underrelax;
		massloss0=massloss;

		do{

			mass1=0.0;
			massinf=0.0;
			diff_bc=0;
			for(i=1;i<=n;i++){
				l=adt->T->lrc_cont->co[i][1];
				r=adt->T->lrc_cont->co[i][2];
				c=adt->T->lrc_cont->co[i][3];
				psi->co[i]=adt->S->P->co[l][r][c]+nw*dpsi->co[i];
				if(psi->co[i]!=psi->co[i]) printf("no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);

				if(l==1){
					sy=adt->S->type->co[r][c];
					k=K(adt->S->P->co[l][r][c], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][1l], adt->S->pa->co[sy][ja][l],
						adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][r][c]);
					Infmax=k + k*(adt->W->h_sup->co[r][c]-psi->co[i])/(0.5*adt->S->pa->co[sy][jdz][l]);
					if((adt->W->h_sup->co[r][c]-psi->co[i])<0) Infmax=k*Kred_exf + k*Kred_exf*(adt->W->h_sup->co[r][c]-psi->co[i])/(0.5*adt->S->pa->co[sy][jdz][l]);
					Infpot=adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/adt->P->Dt;
					adt->S->Jinf->co[r][c]=Fmin(Infpot, Infmax);
					massinf+=adt->S->Jinf->co[r][c]*adt->P->Dt;

					//check consistency of boundary condition
					if( (adt->S->bc->co[r][c]==0 && Infpot>Infmax) || (adt->S->bc->co[r][c]==1 && Infpot<=Infmax) ) diff_bc++;
				}

				mass1+=adt->S->pa->co[sy][jdz][l]*teta_psi(psi->co[i], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
			   		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
			}
			mass1/=(double)adt->P->total_pixel;//[mm]
			massinf/=(double)adt->P->total_pixel;

			massloss=mass0-(mass1-massinf);

			//printf("nw:%f massloss:%e massloss0:%e\n",nw,massloss,massloss0);

			cont2++;
			nw/=4.0;

		}while(fabs(massloss)>fabs(massloss0) && cont2<3);

		resnum=0.0;
		resden=0.0;
		for(i=1;i<=n;i++){
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			if(psi->co[i]<adt->P->psimin) psi->co[i]=adt->P->psimin;
			resnum+=fabs(psi->co[i]-adt->S->P->co[l][r][c]);
			resden+=fabs(adt->S->P->co[l][r][c]-psi0->co[i]);
			adt->S->P->co[l][r][c]=psi->co[i];
		}
		if(resden<1.E-20) resden=1.E-20;

		//printf("->iterCG:%ld res:%e mass0:%f mass1:%f massloss:%f massinf:%f nw:%f diff_bc:%ld tol:%e..  iterNR:%ld/%ld\n",iter,resnum/resden,mass0,mass1,massloss,massinf,nw,diff_bc,tol,cont,adt->P->MaxiterTol);

		out=0;

		if(resnum/resden<=adt->P->TolVWb) out=1;
		if(cont>=adt->P->MaxiterTol) out=1;

		if(fabs(massloss)>adt->P->MaxErrWb) out=0;

		//if(diff_bc>0) out=0;

		//if(resnum/resden<=adt->P->TolVWb && fabs(massloss)>adt->P->MaxErrWb) tol*=0.1;

		if(iter==0){
			//tol*=0.1;
			out=1;
		}

		if(cont>=adt->P->MaxiterErr) out=1;

	}while(out==0);

	cum_losses+=massloss;

	fprintf(f,"->iterCG:%ld res:%e nw:%f mass0:%f mass1:%f massloss:%f massinf:%f diff_bc:%ld tol:%e..  iterNR:%ld/%ld\n",iter_tot,resnum/resden,nw,mass0,mass1,massloss,massinf,diff_bc,tol,cont,adt->P->MaxiterTol);
	fprintf(f,"cum_loss:%f\n",cum_losses);

	massinf2=0.0;
	for(i=1;i<=n;i++){
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];

		if(l==1){
			adt->W->h_sup->co[r][c]+=(adt->W->Pn->co[r][c]-adt->S->Jinf->co[r][c])*adt->P->Dt;
			if(adt->W->h_sup->co[r][c]<0) adt->W->h_sup->co[r][c]=0.0;
		}

	}
	massinf2/=(double)adt->P->total_pixel;

	fprintf(f,"mass0:%f mass1:%f massinf:%f massloss:%f factor:%f\n",mass0,mass1,massinf2,mass0-(mass1-massinf2),(massinf-massloss)/massinf);

	adt->W->out2->co[8]+=massloss;

	free_doublevector(psi);
	free_doublevector(b);
	free_doublevector(dpsi);
	free_doublevector(psi0);
	//free_doublevector(ap);
	//free_doublevector(V);

	fclose(f);


}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Solve_Richards_3D(long i, DOUBLEVECTOR *psi, ALLDATA *adt){

	long l, r, c;
	short sy;
	double a = 0.0;
	double C, dz, k, kn, dzn;
	double Infmax;

	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];

	C=dteta_dpsi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
	dz=adt->S->pa->co[sy][jdz][l];
	a+=psi->co[i]*(C*dz/adt->P->Dt);

	k=K(adt->S->P->co[l][r][c], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(l==1){
		//Dirichlet
		if(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c]<0){
			Infmax=k*Kred_exf+Kred_exf*k*(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c])/(dz/2.);
		}else{
			Infmax=k+k*(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c])/(dz/2.);
		}
		if(adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/adt->P->Dt > Infmax){
			if(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c]<0){
				a+=(k*Kred_exf/(dz/2.))*psi->co[i];
			}else{
				a+=(k/(dz/2.))*psi->co[i];
			}
		}
	}

	if(l>1){
		kn=K(adt->S->P->co[l-1][r][c], adt->S->pa->co[sy][jKv][l-1], adt->P->imp, adt->S->thice->co[l-1][r][c], adt->S->pa->co[sy][jsat][l-1], adt->S->pa->co[sy][jres][l-1],
			adt->S->pa->co[sy][ja][l-1], adt->S->pa->co[sy][jns][l-1], 1-1/adt->S->pa->co[sy][jns][l-1], adt->S->pa->co[sy][jv][l-1], adt->S->pa->co[sy][jpsimin][l-1],
			adt->S->T->co[l-1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l-1];
		kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a+=(kn/dzn)*(psi->co[i]-psi->co[adt->T->i_cont[l-1][r][c]]);
	}

	if(l<Nl){
		kn=K(adt->S->P->co[l+1][r][c], adt->S->pa->co[sy][jKv][l+1], adt->P->imp, adt->S->thice->co[l+1][r][c], adt->S->pa->co[sy][jsat][l+1], adt->S->pa->co[sy][jres][l+1],
			adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 1-1/adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], adt->S->pa->co[sy][jpsimin][l+1],
			adt->S->T->co[l+1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l+1];
		kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a+=(kn/dzn)*(psi->co[i]-psi->co[adt->T->i_cont[l+1][r][c]]);
	}

	k=K(adt->S->P->co[l][r][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(r>1){
		if(adt->L->LC->co[r-1][c]!=NoV){
			sy=adt->S->type->co[r-1][c];
			kn=K(adt->S->P->co[l][r-1][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r-1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r-1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			a+=1.E-3*(kn/UV->U->co[2])*(psi->co[i]-psi->co[adt->T->i_cont[l][r-1][c]]);
		}
	}

	if(r<Nr){
		if(adt->L->LC->co[r+1][c]!=NoV){
			sy=adt->S->type->co[r+1][c];
			kn=K(adt->S->P->co[l][r+1][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r+1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r+1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			a+=1.E-3*(kn/UV->U->co[2])*(psi->co[i]-psi->co[adt->T->i_cont[l][r+1][c]]);
		}
	}

	if(c>1){
		if(adt->L->LC->co[r][c-1]!=NoV){
			sy=adt->S->type->co[r][c-1];
			kn=K(adt->S->P->co[l][r][c-1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c-1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c-1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			a+=1.E-3*(kn/UV->U->co[1])*(psi->co[i]-psi->co[adt->T->i_cont[l][r][c-1]]);
		}
	}

	if(c<Nc){
		if(adt->L->LC->co[r][c+1]!=NoV){
			sy=adt->S->type->co[r][c+1];
			kn=K(adt->S->P->co[l][r][c+1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c+1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c+1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			a+=1.E-3*(kn/UV->U->co[1])*(psi->co[i]-psi->co[adt->T->i_cont[l][r][c+1]]);
		}
	}

	return(a);
}

double Solve_Richards_3D_p(long i, DOUBLEVECTOR *psi, void *adt) {
	return Solve_Richards_3D(i, psi, (ALLDATA *)adt);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Find_b(long i, DOUBLEVECTOR *psi0, ALLDATA *adt){
	/* function that finds the known term dell'eq. di richards*/

	long l, r, c;
	short sy;
	double a = 0.0;
	double C, dz, k, kn, dzn, theta1;
	double Infmax;
	double theta0; /* theta at the time step n*/

	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];

	C=dteta_dpsi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
	theta0=teta_psi(psi0->co[i], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
			adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
	theta1=teta_psi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
	dz=adt->S->pa->co[sy][jdz][l];
	a+=(C*adt->S->P->co[l][r][c]+theta0-theta1)*dz/adt->P->Dt;

	k=K(adt->S->P->co[l][r][c], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(l==1){
		//Dirichlet
		if(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c]<0){
			Infmax=k*Kred_exf+Kred_exf*k*(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c])/(dz/2.);
		}else{
			Infmax=k+k*(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c])/(dz/2.);
		}
		if(adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/adt->P->Dt > Infmax){
			if(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c]<0){
				a+=(k*Kred_exf + k*Kred_exf*adt->W->h_sup->co[r][c]/(dz/2.));
			}else{
				a+=(k + k*adt->W->h_sup->co[r][c]/(dz/2.));
			}
			adt->S->bc->co[r][c]=1;
		//Neumann
		}else{
			a+=( adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/adt->P->Dt );
			adt->S->bc->co[r][c]=0;
		}
	}

	if(l==Nl) a-=k*0.0;

	if(l>1){
		kn=K(adt->S->P->co[l-1][r][c], adt->S->pa->co[sy][jKv][l-1], adt->P->imp, adt->S->thice->co[l-1][r][c], adt->S->pa->co[sy][jsat][l-1], adt->S->pa->co[sy][jres][l-1],
			adt->S->pa->co[sy][ja][l-1], adt->S->pa->co[sy][jns][l-1], 1-1/adt->S->pa->co[sy][jns][l-1], adt->S->pa->co[sy][jv][l-1], adt->S->pa->co[sy][jpsimin][l-1],
			adt->S->T->co[l-1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l-1];
		kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a+=kn;
	}

	if(l<Nl){
		kn=K(adt->S->P->co[l+1][r][c], adt->S->pa->co[sy][jKv][l+1], adt->P->imp, adt->S->thice->co[l+1][r][c], adt->S->pa->co[sy][jsat][l+1], adt->S->pa->co[sy][jres][l+1],
			adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 1-1/adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], adt->S->pa->co[sy][jpsimin][l+1],
			adt->S->T->co[l+1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l+1];
		kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a-=kn;
	}


	k=K(adt->S->P->co[l][r][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(r>1){
		if(adt->L->LC->co[r-1][c]!=NoV){
			sy=adt->S->type->co[r-1][c];
			kn=K(adt->S->P->co[l][r-1][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r-1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r-1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			a+=(kn/UV->U->co[2])*(adt->T->Z0dp->co[r-1][c] - adt->T->Z0dp->co[r][c]);
		}
	}

	if(r<Nr){
		if(adt->L->LC->co[r+1][c]!=NoV){
			sy=adt->S->type->co[r+1][c];
			kn=K(adt->S->P->co[l][r+1][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r+1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r+1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			a+=(kn/UV->U->co[2])*(adt->T->Z0dp->co[r+1][c] - adt->T->Z0dp->co[r][c]);
		}
	}

	if(c>1){
		if(adt->L->LC->co[r][c-1]!=NoV){
			sy=adt->S->type->co[r][c-1];
			kn=K(adt->S->P->co[l][r][c-1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c-1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c-1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			a+=(kn/UV->U->co[1])*(adt->T->Z0dp->co[r][c-1] - adt->T->Z0dp->co[r][c]);
		}
	}

	if(c<Nc){
		if(adt->L->LC->co[r][c+1]!=NoV){
			sy=adt->S->type->co[r][c+1];
			kn=K(adt->S->P->co[l][r][c+1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c+1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c+1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			a+=(kn/UV->U->co[1])*(adt->T->Z0dp->co[r][c+1] - adt->T->Z0dp->co[r][c]);
		}
	}

	return(a);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void supflow2(TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par)
{
	long r,c,R,C,ch,s;                                    /* counters*/
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; /* differential of number-pixel for rows and*/
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; /* columns, depending on Drainage Directions*/
	double dx,dy;                                         /* the two dimensions of a pixel*/
	double Ks;											  /* the Strickler's coefficent calculated with a standard deviation*/
	double b[10];                                         /* area perpendicular to the superficial flow divided by h_sup*/
	double i;											  /*hydraulic gradient*/
	double q,tb,te,dt;
	short lu;

if(par->point_sim==0){	//distributed simulations

	initialize_doublevector(cnet->Qsup,0.0);
	initialize_doublevector(cnet->Qsup_spread,0.0);
	initialize_doublevector(cnet->Qsub,0.0);
	initialize_doublevector(cnet->Qsub_spread,0.0);

	dx=UV->U->co[1];                                     /*cell side [m]*/
	dy=UV->U->co[2];                                     /*cell side [m]*/

	b[1]=0.0;  b[2]=dy;             b[3]=dx/2.0+dy/2.0;  b[4]=dx;            b[5]=dx/2.0+dy/2.0;
	b[6]=dy;   b[7]=dx/2.0+dy/2.0;  b[8]=dx;             b[9]=dx/2.0+dy/2.0;

	te=0.0;

	do{

		tb=te;
		dt=par->Dt;

		//printf("tb:%f dt:%f\n",tb,dt);

		//find dt min
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(wat->h_sup->co[r][c]>0 && top->DD->co[r][c]>=1 && top->DD->co[r][c]<=8 && top->pixel_type->co[r][c]==0){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
					//i=top->i_DD->co[r][c];
					if(i<0) i=0;  //correct false slopes
					q=b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;	//mm/s
					if(q>0){
						if(wat->h_sup->co[r][c]/q<dt) dt=wat->h_sup->co[r][c]/q;
					}
				}
			}
		}

		te=tb+dt;
		if(te>par->Dt){
			te=par->Dt;
			dt=te-tb;
		}

		/*wat->h_sup(=height of water over the land-surface) is updated adding the new runoff(=precipita-
		tion not infiltrated of the actual time); then it is calculated q_sub and it is checked
		that its value is not greater than the avaible water on the land-surface:*/
		/*Remember the units of measure: q_sup=[mm/s],b[m],Ks[m^(1/3)/s],wat->h_sup[mm],dx[m],dy[m]*/
		initialize_doublematrix(wat->q_sup,0.0);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(wat->h_sup->co[r][c]>0 && top->DD->co[r][c]>=1 && top->DD->co[r][c]<=8 && top->pixel_type->co[r][c]==0){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
					//i=top->i_DD->co[r][c];
					if(i<0) i=0.0;
					wat->q_sup->co[r][c]=b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;

					if(wat->q_sup->co[r][c]!=wat->q_sup->co[r][c]){
						printf("NO VALUE SUP:%ld %ld %ld %ld %f %f\n",r,c,r+r_DD[top->DD->co[r][c]],c+c_DD[top->DD->co[r][c]],wat->h_sup->co[r][c],wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]]);
						printf("i:%f Dh:%f Ks:%f pow:%f iF:%f\n",i,wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]],Ks,pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m),top->i_DD->co[r][c]);
					}
				}
			}
		}

		/*After the computation of the surface flow,these flows are moved trough D8 scheme:*/
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->pixel_type->co[r][c]==0){

					R=r+r_DD[top->DD->co[r][c]];
					C=c+c_DD[top->DD->co[r][c]];

					wat->h_sup->co[r][c]-=wat->q_sup->co[r][c]*dt;
					wat->h_sup->co[r][c]=Fmax(wat->h_sup->co[r][c], 0.0);

					//the superficial flow is added to the land pixels (code 0):
					if (top->pixel_type->co[R][C]==0){
						wat->h_sup->co[R][C]+=wat->q_sup->co[r][c]*dt;

					//the superficial flow is added to flow which flows into the channel pixels (code 10):
					}else if (top->pixel_type->co[R][C]==10){
						for(ch=1;ch<=cnet->r->nh;ch++){
							if(R==cnet->r->co[ch] && C==cnet->c->co[ch]){
								cnet->Qsup->co[ch]+=wat->q_sup->co[r][c]*dt/par->Dt;
								if(cnet->Qsup->co[ch]!=cnet->Qsup->co[ch]){
									printf("qsup no value: r:%ld c:%ld ch:%ld R:%ld C:%ld qsup:%f hsup:%f\n",r,c,ch,R,C,wat->q_sup->co[r][c],wat->h_sup->co[r][c]);
								}
							}
						}
					}

				}else if(top->pixel_type->co[r][c]==10){
					for(ch=1;ch<=cnet->r->nh;ch++){
						if(r==cnet->r->co[ch] && c==cnet->c->co[ch]){
							cnet->Qsub->co[ch]=wat->h_sup->co[r][c]/par->Dt; /*[mm/s]*/
							wat->h_sup->co[r][c]=0.0;
						}
					}
				}

			}
		}

		for(ch=1;ch<=cnet->r->nh;ch++){
			for(s=1;s<=cnet->Qsup_spread->nh;s++){
				cnet->Qsup_spread->co[s]+=(cnet->Qsup->co[ch])*(cnet->fraction_spread->co[ch][s])*0.001*dx*dy; /*in mc/s*/
				cnet->Qsub_spread->co[s]+=(cnet->Qsub->co[ch])*(cnet->fraction_spread->co[ch][s])*0.001*dx*dy; /*in mc/s*/
			}

		}


	}while(te<par->Dt);


}else{	//point simulation

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				q=0.0;
				lu=(short)land->LC->co[r][c];
				Ks=land->ty->co[lu][jcm];
				i=pow(pow(top->dz_dx->co[r][c],2.0)+pow(top->dz_dy->co[r][c],2.0),0.5);
				if(wat->h_sup->co[r][c]>0) q=Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0;	//mm/s
				wat->h_sup->co[r][c]-=q*par->Dt;
				if(wat->h_sup->co[r][c]<0) wat->h_sup->co[r][c]=0.0;
			}
		}
	}
}


}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void routing2(CHANNEL *cnet){

	long s;

	/*The just calculated Qsup_spread and Qsub_spread are added to the flow already
	present in the channel-network(Q_sup_s and Q_sub_s) and at once it is made a
	traslation of Q_sup_s e Q_sub_s towards the outlet:*/
	for(s=1;s<=cnet->Qsup_spread->nh-1;s++){
		cnet->Q_sub_s->co[s]=cnet->Q_sub_s->co[s+1]+cnet->Qsub_spread->co[s];
		cnet->Q_sup_s->co[s]=cnet->Q_sup_s->co[s+1]+cnet->Qsup_spread->co[s];
	}
	cnet->Q_sub_s->co[cnet->Qsub_spread->nh]=cnet->Qsub_spread->co[cnet->Qsub_spread->nh];
	cnet->Q_sup_s->co[cnet->Qsup_spread->nh]=cnet->Qsup_spread->co[cnet->Qsup_spread->nh];
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_waterbalance2(WATER *wat, SOIL *sl, PAR *par, DOUBLEMATRIX *Z){

	long i,l,r,c;

	for(i=1;i<=par->chkpt->nrh;i++){
		r=par->rc->co[i][1];
		c=par->rc->co[i][2];

		wat->out1->co[5][i]+=wat->Pn->co[r][c]*par->Dt;
		wat->out1->co[6][i]+=(wat->Pn->co[r][c]-sl->Jinf->co[r][c])*par->Dt;
		wat->out1->co[7][i]+=(wat->h_sup->co[r][c] - wat->out1->co[9][i] - (wat->Pn->co[r][c]-sl->Jinf->co[r][c])*par->Dt);    /*positive if entering into the pixel*/
		wat->out1->co[8][i]=0.0;
		wat->out1->co[9][i]=wat->h_sup->co[r][c];

		for(l=1;l<=Nl;l++){
			wat->out1->co[10][i]-=0.0*par->Dt;/*positive if entering into the pixel*/
		}
	}

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(Z->co[r][c]!=NoV){ /*if the pixel is not a novalue*/
				wat->out2->co[3]+=wat->Pn->co[r][c]*par->Dt;								/*[mm]*/
				wat->out2->co[4]+=sl->Jinf->co[r][c]*par->Dt;							/*[mm]*/
			}
		}
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
