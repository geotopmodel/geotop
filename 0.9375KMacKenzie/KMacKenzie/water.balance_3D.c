
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


#include "constant.h"
#include "struct.geotop.09375.h"
#include "water.balance_3D.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"

extern T_INIT *UV;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void water_balance_3D(ALLDATA *adt){

	DOUBLEMATRIX *h;
	DOUBLETENSOR *PSI;

	long n, r, c, l;
	double Dt0, Dt, tb, te, Loss;
	static double cum_losses;

	//FILE *f;

	//if(adt->I->time==0.0) cum_losses=0.0;
	cum_losses=0.0;

	h=new_doublematrix(Nr, Nc);
	PSI=new_doubletensor(Nl, Nr, Nc);

	initialize_doublevector(adt->C->Qsup_spread,0.0);
	initialize_doublevector(adt->C->Qsub_spread,0.0);

	Dt0=adt->P->Dt;
	n=1;

	te=0.0;

	do{

		tb=te;

		do{

			Dt=Dt0/(double)n;

			Richards_3D(Dt, PSI, h, &Loss, adt);

			//f=fopen(files->co[ferr]+1,"a");
			printf("n:%ld Dt:%f tb:%f %e\n\n\n",n,Dt,tb,Loss);
			//fclose(f);

			//printf("Loss:%f MaxLoss:%f\n",Loss,adt->P->MaxErrWb);

			if(fabs(Loss) > adt->P->MaxErrWb && Dt>=DtminWb) n*=nredDtWb;

		}while( fabs(Loss) > adt->P->MaxErrWb && Dt>=DtminWb);

		te=tb+Dt;

		supflow2(Dt, Dt0, h, adt->T, adt->L, adt->W, adt->C, adt->P);

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(adt->L->LC->co[r][c]!=NoV){
					adt->W->h_sup->co[r][c] = h->co[r][c];
					for(l=1;l<=Nl;l++){
						adt->S->P->co[l][r][c] = PSI->co[l][r][c];
					}
				}
			}
		}

		//printf("te/Dt0:%f\n",Dt0/te);
		//if(fabs(Dt0/te-27.)<0.001) n=27;

		adt->W->out2->co[8]+=Loss;

		cum_losses+=Loss;

	}while(te<Dt0);

	routing2(adt->C);

	free_doublematrix(h);
	free_doubletensor(PSI);

	//f=fopen(files->co[ferr]+1,"a");
	printf("Loss:%e\n\n",cum_losses);
	//fclose(f);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Richards_3D(double Dt, DOUBLETENSOR *P, DOUBLEMATRIX *h, double *loss, ALLDATA *adt){

	DOUBLEVECTOR *H00, *H0, *H1, *dH, *B;
	long i, l, r, c, cont, cont2, diff_bc, iter_tot=0, iter;
	short sy;
	double mass0=0.0, mass1, massinf, massloss0, massloss=1.E99, error=1.E99, error0, Infmax, Infpot, k, nw, resnum, resden, psi;
	short out;
	long n=Nl*adt->P->total_pixel;

	H00=new_doublevector(n);
	H0=new_doublevector(n);
	dH=new_doublevector(n);
	H1=new_doublevector(n);
	//Ad=new_doublevector(n);
	//Asup=new_doublevector(n-1);
	B=new_doublevector(n);

	for(i=1;i<=n;i++){
		l = adt->T->lrc_cont->co[i][1];
		r = adt->T->lrc_cont->co[i][2];
		c = adt->T->lrc_cont->co[i][3];
		sy = adt->S->type->co[r][c];

		H1->co[i] = adt->S->P->co[l][r][c] + adt->T->Z->co[l][r][c];
		H0->co[i] = H1->co[i];
		H00->co[i] = H0->co[i];

		mass0 += adt->S->pa->co[sy][jdz][l]*teta_psi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
	}
	mass0/=(double)adt->P->total_pixel;	//in [mm]


	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(adt->L->LC->co[r][c]!=NoV){
				sy=adt->S->type->co[r][c];
				l=1;
				i=adt->T->i_cont[l][r][c];

				k=K(adt->S->P->co[l][r][c], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l],
					adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l],
					adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][r][c]);

				Infmax=k+k*(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c])/(adt->S->pa->co[sy][jdz][l]/2.);
				Infpot=adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/Dt;

				if(Infpot > Infmax){
					adt->S->bc->co[r][c]=1;
				}else{
					adt->S->bc->co[r][c]=0;
				}
			}
		}
	}

	cont=0;

	do{

		cont++;

		for(i=1;i<=n;i++){
			//H00->co[i] = H0->co[i];
			H0->co[i] = H1->co[i];
		}

		for(i=1;i<=n;i++){
			B->co[i] = Find_b(i, H0, H00, Dt, adt);
		}

		iter=tridiag_preconditioned_conjugate_gradient_search(1.E-6, H1, H0, H00, Dt, B, Solve_Richards_3D_p, adt);
		iter_tot+=iter;

		//find_coeff_Richards_3D(Asup, Ad, B, H, Dt, adt);
		//tridiag(0, 0, 0, n, Asup, Ad, Asup, B, H1);

		for(i=1;i<=n;i++){
			dH->co[i] = H1->co[i]-H0->co[i];
		}

		cont2=0.0;
		nw=1.0;

		massloss0=massloss;
		error0=error;

		do{

			mass1=0.0;
			for(i=1;i<=n;i++){
				H1->co[i] = H0->co[i] + nw*dH->co[i];
				if(H1->co[i] != H1->co[i]) printf("no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);
				if(H1->co[i] - adt->T->Z->co[l][r][c] < adt->P->psimin) H1->co[i] = adt->P->psimin + adt->T->Z->co[l][r][c];

				l=adt->T->lrc_cont->co[i][1];
				r=adt->T->lrc_cont->co[i][2];
				c=adt->T->lrc_cont->co[i][3];
				sy=adt->S->type->co[r][c];

				psi = H1->co[i] - adt->T->Z->co[l][r][c];

				mass1 += adt->S->pa->co[sy][jdz][l]*teta_psi(psi, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
			   		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
			}
			mass1/=(double)adt->P->total_pixel;//[mm]

			massinf=0.0;
			diff_bc=0;
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(adt->L->LC->co[r][c]!=NoV){
						sy=adt->S->type->co[r][c];
						l=1;
						i=adt->T->i_cont[l][r][c];

						psi = H00->co[i] - adt->T->Z->co[l][r][c];
						k=K(psi, adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l],
							adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l],
							adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][r][c]);

						Infmax=k*(adt->W->h_sup->co[r][c]+1.E3*adt->T->Z0dp->co[r][c]-H1->co[i])/(adt->S->pa->co[sy][jdz][l]/2.);
						Infpot=adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/Dt;

						adt->S->Jinf->co[r][c]=Fmin(Infmax,Infpot);
						massinf += adt->S->Jinf->co[r][c]*Dt;

						if( (Infpot > Infmax && adt->S->bc->co[r][c]==0) || (Infpot < Infmax && adt->S->bc->co[r][c]==1) ) diff_bc++;

						if(Infpot > Infmax){
							adt->S->bc->co[r][c]=1;
						}else{
							adt->S->bc->co[r][c]=0;
						}
					}
				}
			}
			massinf/=(double)adt->P->total_pixel;
			massloss=mass0-(mass1-massinf);

			/*find_coeff_Richards_3D(Asup, Ad, B, psi, Dt, adt);
			error=0.0;
			for(i=1;i<=n;i++){
				if(i==1){
					error += pow(Ad->co[i]*psi1->co[i] + Asup->co[i]*psi1->co[i+1] - B->co[i], 2.0);
				}else if(i>1 && i<n){
					error += pow(Asup->co[i-1]*psi1->co[i-1] + Ad->co[i]*psi1->co[i] + Asup->co[i]*psi1->co[i+1] - B->co[i], 2.0);
				}else{
					error += pow(Asup->co[i-1]*psi1->co[i-1] + Ad->co[i]*psi1->co[i] - B->co[i], 2.0);
				}
			}
			error=pow(error,0.5);*/

			cont2++;
			nw/=2.0;

			printf("cont:%ld cont2:%ld error:%e error0:%e massloss:%e massloss0:%e mass0:%e mass1:%e massinf:%e diff_bc:%ld\n",
				cont,cont2,error,error0,massloss,massloss0,mass0,mass1,massinf,diff_bc);

		//}while( (fabs(massloss)>fabs(massloss0) || fabs(error)>fabs(error0)) && cont2<2 );
		}while( fabs(massloss)>fabs(massloss0) && cont2<10 );

		resnum=0.0;
		resden=0.0;
		for(i=1;i<=n;i++){
			resnum+=fabs(H1->co[i] - H0->co[i]);
			resden+=fabs(H0->co[i] - (adt->S->P->co[l][r][c]+adt->T->Z->co[l][r][c]));
		}
		if(resden<1.E-20) resden=1.E-20;

		out=0;
		if(resnum/resden<=adt->P->TolVWb) out=1;
		if(cont>=adt->P->MaxiterTol) out=1;
		if(fabs(massloss)>adt->P->MaxErrWb) out=0;
		if(diff_bc>0) out=0;
		if(cont==1) out=0;
		if(cont>=adt->P->MaxiterErr) out=1;

		printf("res:%e out:%d\n",resnum/resden,out);

	}while(out==0);

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(adt->L->LC->co[r][c]!=NoV){
				sy=adt->S->type->co[r][c];
				h->co[r][c] = adt->W->h_sup->co[r][c] + (adt->W->Pn->co[r][c]-adt->S->Jinf->co[r][c])*Dt;
				if(h->co[r][c]<0) h->co[r][c]=0.0;
			}
		}
	}

	for(i=1;i<n;i++){
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		P->co[l][r][c] = H1->co[i] - adt->T->Z->co[l][r][c];
	}

	*loss = massloss;

	free_doublevector(H00);
	free_doublevector(H0);
	free_doublevector(dH);
	free_doublevector(H1);
	//free_doublevector(Ad);
	//free_doublevector(Asup);
	free_doublevector(B);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Solve_Richards_3D(long i, DOUBLEVECTOR *H1, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, ALLDATA *adt){

	long l, r, c, I;
	short sy;
	double a = 0.0;
	double C, dz, k, kn, dzn, psi;

	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];

	psi = H0->co[i] - adt->T->Z->co[l][r][c];

	C = dteta_dpsi(psi, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
	dz = adt->S->pa->co[sy][jdz][l];
	a += H1->co[i]*(C*dz/Dt);

	psi = H00->co[i] - adt->T->Z->co[l][r][c];
	k = K(psi, adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],adt->S->T->co[l][r][c]);

	//Dirichlet
	if(l==1){ if(adt->S->bc->co[r][c]==1) a += (k/(dz/2.))*H1->co[i]; }

	if(l>1){
		I=adt->T->i_cont[l-1][r][c];

		psi = H00->co[I] - adt->T->Z->co[l-1][r][c];

		kn=K(psi, adt->S->pa->co[sy][jKv][l-1], adt->P->imp, adt->S->thice->co[l-1][r][c], adt->S->pa->co[sy][jsat][l-1], adt->S->pa->co[sy][jres][l-1],
			adt->S->pa->co[sy][ja][l-1], adt->S->pa->co[sy][jns][l-1], 1-1/adt->S->pa->co[sy][jns][l-1], adt->S->pa->co[sy][jv][l-1], adt->S->pa->co[sy][jpsimin][l-1],
			adt->S->T->co[l-1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l-1];
		kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a+=(kn/dzn)*(H1->co[i]-H1->co[I]);
	}

	if(l<Nl){
		I=adt->T->i_cont[l+1][r][c];

		psi = H00->co[I] - adt->T->Z->co[l+1][r][c];

		kn=K(psi, adt->S->pa->co[sy][jKv][l+1], adt->P->imp, adt->S->thice->co[l+1][r][c], adt->S->pa->co[sy][jsat][l+1], adt->S->pa->co[sy][jres][l+1],
			adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 1-1/adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], adt->S->pa->co[sy][jpsimin][l+1],
			adt->S->T->co[l+1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l+1];
		kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a+=(kn/dzn)*(H1->co[i]-H1->co[I]);
	}

	psi = H00->co[i] - adt->T->Z->co[l][r][c];
	k=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(r>1){
		if(adt->L->LC->co[r-1][c]!=NoV){
			I=adt->T->i_cont[l][r-1][c];

			psi = H00->co[I] - adt->T->Z->co[l][r-1][c];

			sy=adt->S->type->co[r-1][c];
			kn=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r-1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r-1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			a+=1.E-3*(kn/UV->U->co[2])*(H1->co[i]-H1->co[I]);
		}
	}

	if(r<Nr){
		if(adt->L->LC->co[r+1][c]!=NoV){
			I=adt->T->i_cont[l][r+1][c];

			psi = H00->co[I] - adt->T->Z->co[l][r+1][c];

			sy=adt->S->type->co[r+1][c];
			kn=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r+1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r+1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			a+=1.E-3*(kn/UV->U->co[2])*(H1->co[i]-H1->co[I]);
		}
	}

	if(c>1){
		if(adt->L->LC->co[r][c-1]!=NoV){
			I=adt->T->i_cont[l][r][c-1];

			psi = H00->co[I] - adt->T->Z->co[l][r][c-1];

			sy=adt->S->type->co[r][c-1];
			kn=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c-1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c-1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			a+=1.E-3*(kn/UV->U->co[1])*(H1->co[i]-H1->co[I]);
		}
	}

	if(c<Nc){
		if(adt->L->LC->co[r][c+1]!=NoV){
			I=adt->T->i_cont[l][r][c+1];

			psi = H00->co[I] - adt->T->Z->co[l][r][c+1];

			sy=adt->S->type->co[r][c+1];
			kn=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c+1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c+1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			a+=1.E-3*(kn/UV->U->co[1])*(H1->co[i]-H1->co[I]);
		}
	}

	return(a);
}

double Solve_Richards_3D_p(long i, DOUBLEVECTOR *H1, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, void *adt) {
	return Solve_Richards_3D(i, H1, H0, H00, Dt, (ALLDATA *)adt);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Find_b(long i, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, ALLDATA *adt){

	long l, r, c;
	short sy;
	double a = 0.0;
	double C, dz, k, theta1, theta0, psi;

	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];

	psi = H0->co[i] - adt->T->Z->co[l][r][c];

	C=dteta_dpsi(psi, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);

	theta1=teta_psi(psi, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->P->Esoil, adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
	theta0=teta_psi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);

	dz=adt->S->pa->co[sy][jdz][l];
	a+=(C*H0->co[i]+theta0-theta1)*dz/Dt;

	psi = H00->co[i] - adt->T->Z->co[l][r][c];

	k=K(psi, adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(l==1){
		//Dirichlet
		if(adt->S->bc->co[r][c]==1){
			a += k*( adt->W->h_sup->co[r][c] + 1.E3*adt->T->Z0dp->co[r][c] )/(dz/2.);

		//Neumann
		}else{
			a += ( adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/Dt );
		}
	}

	/*if(l==Nl) a-=k*0.0;

	if(l>1){
		I=adt->T->i_cont[l-1][r][c];

		psi = H0->co[I] - 1.E3*adt->T->Z->co[l-1][r][c];

		kn=K(psi, adt->S->pa->co[sy][jKv][l-1], adt->P->imp, adt->S->thice->co[l-1][r][c], adt->S->pa->co[sy][jsat][l-1], adt->S->pa->co[sy][jres][l-1],
			adt->S->pa->co[sy][ja][l-1], adt->S->pa->co[sy][jns][l-1], 1-1/adt->S->pa->co[sy][jns][l-1], adt->S->pa->co[sy][jv][l-1], adt->S->pa->co[sy][jpsimin][l-1],
			adt->S->T->co[l-1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l-1];
		kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a+=kn;
	}

	if(l<Nl){
		I=adt->T->i_cont[l+1][r][c];

		psi = H0->co[I] - 1.E3*adt->T->Z->co[l+1][r][c];

		kn=K(psi, adt->S->pa->co[sy][jKv][l+1], adt->P->imp, adt->S->thice->co[l+1][r][c], adt->S->pa->co[sy][jsat][l+1], adt->S->pa->co[sy][jres][l+1],
			adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 1-1/adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], adt->S->pa->co[sy][jpsimin][l+1],
			adt->S->T->co[l+1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l+1];
		kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a-=kn;
	}

	psi = H0->co[i] - 1.E3*adt->T->Z->co[l][r][c];
	k=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(r>1){
		if(adt->L->LC->co[r-1][c]!=NoV){
			I=adt->T->i_cont[l][r-1][c];
			sy=adt->S->type->co[r-1][c];

			psi = H0->co[I] - 1.E3*adt->T->Z->co[l][r-1][c];

			kn=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r-1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r-1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			a+=(kn/UV->U->co[2])*(adt->T->Z0dp->co[r-1][c] - adt->T->Z0dp->co[r][c]);
		}
	}

	if(r<Nr){
		if(adt->L->LC->co[r+1][c]!=NoV){
			I=adt->T->i_cont[l][r+1][c];
			sy=adt->S->type->co[r+1][c];

			psi = H0->co[I] - 1.E3*adt->T->Z->co[l][r+1][c];

			kn=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r+1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r+1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			a+=(kn/UV->U->co[2])*(adt->T->Z0dp->co[r+1][c] - adt->T->Z0dp->co[r][c]);
		}
	}

	if(c>1){
		if(adt->L->LC->co[r][c-1]!=NoV){
			I=adt->T->i_cont[l][r][c-1];
			sy=adt->S->type->co[r][c-1];

			psi = H0->co[I] - 1.E3*adt->T->Z->co[l][r][c-1];

			kn=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c-1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c-1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			a+=(kn/UV->U->co[1])*(adt->T->Z0dp->co[r][c-1] - adt->T->Z0dp->co[r][c]);
		}
	}

	if(c<Nc){
		if(adt->L->LC->co[r][c+1]!=NoV){
			I=adt->T->i_cont[l][r][c+1];
			sy=adt->S->type->co[r][c+1];

			psi = H0->co[I] - 1.E3*adt->T->Z->co[l][r][c+1];

			kn=K(psi, adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c+1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c+1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			a+=(kn/UV->U->co[1])*(adt->T->Z0dp->co[r][c+1] - adt->T->Z0dp->co[r][c]);
		}
	}*/

	return(a);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void supflow2(double Dt, double Dtmax, DOUBLEMATRIX *h, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par)

{
	long r,c,R,C,ch,s;
	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; // differential of number-pixel for rows and
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; // columns, depending on Drainage Directions
	double dx,dy;                                         // the two dimensions of a pixel
	double Ks;											  // the Strickler's coefficent calculated with a standard deviation
	double b[10];                                         // area perpendicular to the superficial flow divided by h_sup
	double i;											  //hydraulic gradient
	double q,tb,te,dt;
	short lu;

if(par->point_sim==0){	//distributed simulations

	initialize_doublevector(cnet->Qsup,0.0);
	initialize_doublevector(cnet->Qsub,0.0);

	dx=UV->U->co[1];
	dy=UV->U->co[2];

	b[1]=0.0;  b[2]=dy;             b[3]=dx/2.0+dy/2.0;  b[4]=dx;            b[5]=dx/2.0+dy/2.0;
	b[6]=dy;   b[7]=dx/2.0+dy/2.0;  b[8]=dx;             b[9]=dx/2.0+dy/2.0;

	te=0.0;

	do{

		tb=te;
		dt=Dt;

		//find dt min
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(h->co[r][c]>0 && top->DD->co[r][c]>=1 && top->DD->co[r][c]<=8 && top->pixel_type->co[r][c]==0){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( h->co[r][c] - h->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
					//i=top->i_DD->co[r][c];
					if(i<0) i=0;  //correct false slopes
					q=b[top->DD->co[r][c]+1]*Ks*pow(h->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;	//mm/s
					if(q>0){
						if(h->co[r][c]/q<dt) dt=h->co[r][c]/q;
					}
				}
			}
		}

		te=tb+dt;
		if(te>Dt){
			te=Dt;
			dt=te-tb;
		}

		//h(=height of water over the land-surface) is updated adding the new runoff(=precipita-
		//tion not infiltrated of the actual time); then it is calculated q_sub and it is checked
		//that its value is not greater than the avaible water on the land-surface:
		//Remember the units of measure: q_sup=[mm/s],b[m],Ks[m^(1/3)/s],h[mm],dx[m],dy[m]
		initialize_doublematrix(wat->q_sup,0.0);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(h->co[r][c]>0 && top->DD->co[r][c]>=1 && top->DD->co[r][c]<=8 && top->pixel_type->co[r][c]==0){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( h->co[r][c] - h->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
					//i=top->i_DD->co[r][c];
					if(i<0) i=0.0;
					wat->q_sup->co[r][c]=b[top->DD->co[r][c]+1]*Ks*pow(h->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;

					if(wat->q_sup->co[r][c]!=wat->q_sup->co[r][c]){
						printf("NO VALUE SUP:%ld %ld %ld %ld %f %f\n",r,c,r+r_DD[top->DD->co[r][c]],c+c_DD[top->DD->co[r][c]],h->co[r][c],h->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]]);
						printf("i:%f Dh:%f Ks:%f pow:%f iF:%f\n",i,h->co[r][c] - h->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]],Ks,pow(h->co[r][c]/1000.0,1.0+par->gamma_m),top->i_DD->co[r][c]);
					}
				}
			}
		}

		//After the computation of the surface flow,these flows are moved trough D8 scheme:
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->pixel_type->co[r][c]==0){

					R=r+r_DD[top->DD->co[r][c]];
					C=c+c_DD[top->DD->co[r][c]];

					h->co[r][c]-=wat->q_sup->co[r][c]*dt;
					h->co[r][c]=Fmax(h->co[r][c], 0.0);

					//the superficial flow is added to the land pixels (code 0):
					if (top->pixel_type->co[R][C]==0){
						h->co[R][C]+=wat->q_sup->co[r][c]*dt;

					//the superficial flow is added to flow which flows into the channel pixels (code 10):
					}else if (top->pixel_type->co[R][C]==10){
						for(ch=1;ch<=cnet->r->nh;ch++){
							if(R==cnet->r->co[ch] && C==cnet->c->co[ch]){
								cnet->Qsup->co[ch]+=wat->q_sup->co[r][c]*dt/Dtmax;
								if(cnet->Qsup->co[ch]!=cnet->Qsup->co[ch]){
									printf("qsup no value: r:%ld c:%ld ch:%ld R:%ld C:%ld qsup:%f hsup:%f\n",r,c,ch,R,C,wat->q_sup->co[r][c],h->co[r][c]);
								}
							}
						}
					}

				}else if(top->pixel_type->co[r][c]==10){
					for(ch=1;ch<=cnet->r->nh;ch++){
						if(r==cnet->r->co[ch] && c==cnet->c->co[ch]){
							cnet->Qsub->co[ch]=h->co[r][c]/Dtmax; //[mm/s]
							h->co[r][c]=0.0;
						}
					}
				}

			}
		}

		for(ch=1;ch<=cnet->r->nh;ch++){
			for(s=1;s<=cnet->Qsup_spread->nh;s++){
				cnet->Qsup_spread->co[s]+=(cnet->Qsup->co[ch])*(cnet->fraction_spread->co[ch][s])*0.001*dx*dy; //in m3/s
				cnet->Qsub_spread->co[s]+=(cnet->Qsub->co[ch])*(cnet->fraction_spread->co[ch][s])*0.001*dx*dy; //in m3/s
			}

		}


	}while(te<Dt);


}else{	//point simulation

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				/*q=0.0;
				lu=(short)land->LC->co[r][c];
				Ks=land->ty->co[lu][jcm];
				i=pow(pow(top->dz_dx->co[r][c],2.0)+pow(top->dz_dy->co[r][c],2.0),0.5);
				if(h->co[r][c]>0) q=Ks*pow(h->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0;	//mm/s
				h->co[r][c]-=q*Dt;
				if(h->co[r][c]<0) h->co[r][c]=0.0;*/
				h->co[r][c]=0.0;
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

double calc_V_water(long l, long r, long c, DOUBLEVECTOR *psi, ALLDATA *adt){

	short sy;
	double dz, k, kn, a=0.0;
	double ig, ip;

	sy=adt->S->type->co[r][c];
	dz=adt->S->pa->co[sy][jdz][l];
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
			ig=(adt->T->Z0dp->co[r-1][c] - adt->T->Z0dp->co[r][c])/UV->U->co[2];
			ip=1.E-3*(psi->co[adt->T->i_cont[l][r-1][c]]  - psi->co[adt->T->i_cont[l][r][c]])/UV->U->co[2];
			a+=kn*(ip+ig)*adt->P->Dt;
		}
	}

	if(r<Nr){
		if(adt->L->LC->co[r+1][c]!=NoV){
			sy=adt->S->type->co[r+1][c];
			kn=K(adt->S->P->co[l][r+1][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r+1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r+1][c]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);
			ig=(adt->T->Z0dp->co[r+1][c] - adt->T->Z0dp->co[r][c])/UV->U->co[2];
			ip=1.E-3*(psi->co[adt->T->i_cont[l][r+1][c]]  - psi->co[adt->T->i_cont[l][r][c]])/UV->U->co[2];
			a+=kn*(ip+ig)*adt->P->Dt;
		}
	}

	if(c>1){
		if(adt->L->LC->co[r][c-1]!=NoV){
			sy=adt->S->type->co[r][c-1];
			kn=K(adt->S->P->co[l][r][c-1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c-1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c-1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			ig=(adt->T->Z0dp->co[r][c-1] - adt->T->Z0dp->co[r][c])/UV->U->co[1];
			ip=1.E-3*(psi->co[adt->T->i_cont[l][r][c-1]]  - psi->co[adt->T->i_cont[l][r][c]])/UV->U->co[1];
			a+=kn*(ip+ig)*adt->P->Dt;
		}
	}

	if(c<Nc){
		if(adt->L->LC->co[r][c+1]!=NoV){
			sy=adt->S->type->co[r][c+1];
			kn=K(adt->S->P->co[l][r][c+1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c+1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c+1]);
			kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);
			ig=(adt->T->Z0dp->co[r][c+1] - adt->T->Z0dp->co[r][c])/UV->U->co[1];
			ip=1.E-3*(psi->co[adt->T->i_cont[l][r][c+1]]  - psi->co[adt->T->i_cont[l][r][c]])/UV->U->co[1];
			a+=kn*(ip+ig)*adt->P->Dt;
		}
	}

	return(a);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_coeff_Richards_3D(DOUBLEVECTOR *Asup, DOUBLEVECTOR *Ad, DOUBLEVECTOR *B, DOUBLEVECTOR *psi, double Dt, ALLDATA *adt){

	//psi: previous iteration
	//adt->S->P: previous time step

	long i,I,l,r,c,L,R,C;
	long n=psi->nh;
	short sy;
	double dz, HyC, theta1, theta0, k, kn, dzn;

	initialize_doublevector(Asup, 0.0);
	initialize_doublevector(Ad, 0.0);
	initialize_doublevector(B, 0.0);

	for(i=1;i<=n;i++){

		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];
		dz=adt->S->pa->co[sy][jdz][l];

		HyC=dteta_dpsi(psi->co[i], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
			adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);

		theta1=teta_psi(psi->co[i], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
			adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);
		theta0=teta_psi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
			adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin, adt->P->Esoil);

		Ad->co[i] += HyC*dz/Dt;
		B->co[i] += (HyC*psi->co[i]+theta0-theta1)*dz/Dt;

		k=K(psi->co[i], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
			adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
			adt->S->T->co[l][r][c]);

		//Boundary conditions
		if(l==1){
			//Dirichlet
			if(adt->S->bc->co[r][c]==1){
				Ad->co[i] += (k/(dz/2.));
				B->co[i] += (k + k*adt->W->h_sup->co[r][c]/(dz/2.));

			//Neumann
			}else{
				B->co[i] += ( adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/Dt );
			}
		}

		if(l==Nl) B->co[i] -= k*0.0;

		//Vertical fluxes
		if(l>1){
			L=l-1;
			I=adt->T->i_cont[L][r][c];
			kn=K(psi->co[I], adt->S->pa->co[sy][jKv][L], adt->P->imp, adt->S->thice->co[L][r][c], adt->S->pa->co[sy][jsat][L],
				adt->S->pa->co[sy][jres][L], adt->S->pa->co[sy][ja][L], adt->S->pa->co[sy][jns][L], 1-1/adt->S->pa->co[sy][jns][L],
				adt->S->pa->co[sy][jv][L], adt->S->pa->co[sy][jpsimin][L], adt->S->T->co[L][r][c]);
			dzn=adt->S->pa->co[sy][jdz][L];
			kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
			dzn=0.5*dz + 0.5*dzn;

			Ad->co[i] += (kn/dzn);
			B->co[i] += kn;
		}

		if(l<Nl){
			L=l+1;
			I=adt->T->i_cont[L][r][c];
			kn=K(psi->co[I], adt->S->pa->co[sy][jKv][L], adt->P->imp, adt->S->thice->co[L][r][c], adt->S->pa->co[sy][jsat][L],
				adt->S->pa->co[sy][jres][L], adt->S->pa->co[sy][ja][L], adt->S->pa->co[sy][jns][L], 1-1/adt->S->pa->co[sy][jns][L],
				adt->S->pa->co[sy][jv][L], adt->S->pa->co[sy][jpsimin][L], adt->S->T->co[L][r][c]);
			dzn=adt->S->pa->co[sy][jdz][L];
			kn=Mean(adt->P->harm_or_arit_mean, dz, dzn, k, kn);
			dzn=0.5*dz + 0.5*dzn;

			Ad->co[i] += (kn/dzn);
			Asup->co[i] -= (kn/dzn);
			B->co[i] -= kn;
		}

		//Lateral fluxes
		k=K(psi->co[i], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
			adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
			adt->S->T->co[l][r][c]);

		R=r-1;
		C=c;
		if(R>=1){
			if(adt->L->LC->co[R][C]!=NoV){
				sy=adt->S->type->co[R][C];
				I=adt->T->i_cont[l][R][C];
				kn=K(psi->co[I], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][R][C], adt->S->pa->co[sy][jsat][l],
					adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l],
					adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][R][C]);
				kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);

				B->co[i] += (kn/UV->U->co[2])*(1.E-3*(psi->co[I] - psi->co[i]) + (adt->T->Z0dp->co[R][C] - adt->T->Z0dp->co[r][c]));
				//B->co[i] += (kn/UV->U->co[2])*(1.E-3*(adt->S->P->co[l][R][C] - adt->S->P->co[l][r][c]) + (adt->T->Z0dp->co[R][C] - adt->T->Z0dp->co[r][c]));

			}
		}

		R=r+1;
		C=c;
		if(R<=Nr){
			if(adt->L->LC->co[R][C]!=NoV){
				sy=adt->S->type->co[R][C];
				I=adt->T->i_cont[l][R][C];
				kn=K(psi->co[I], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][R][C], adt->S->pa->co[sy][jsat][l],
					adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l],
					adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][R][C]);
				kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[2], UV->U->co[2], k, kn);

				B->co[i] += (kn/UV->U->co[2])*(1.E-3*(psi->co[I] - psi->co[i]) + (adt->T->Z0dp->co[R][C] - adt->T->Z0dp->co[r][c]));
				//B->co[i] += (kn/UV->U->co[2])*(1.E-3*(adt->S->P->co[l][R][C] - adt->S->P->co[l][r][c]) + (adt->T->Z0dp->co[R][C] - adt->T->Z0dp->co[r][c]));

			}
		}

		R=r;
		C=c-1;
		if(C>=1){
			if(adt->L->LC->co[R][C]!=NoV){
				sy=adt->S->type->co[R][C];
				I=adt->T->i_cont[l][R][C];
				kn=K(psi->co[I], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][R][C], adt->S->pa->co[sy][jsat][l],
					adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l],
					adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][R][C]);
				kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);

				B->co[i] += (kn/UV->U->co[1])*(1.E-3*(psi->co[I] - psi->co[i]) + (adt->T->Z0dp->co[R][C] - adt->T->Z0dp->co[r][c]));
				//B->co[i] += (kn/UV->U->co[1])*(1.E-3*(adt->S->P->co[l][R][C] - adt->S->P->co[l][r][c]) + (adt->T->Z0dp->co[R][C] - adt->T->Z0dp->co[r][c]));

			}
		}

		R=r;
		C=c+1;
		if(C<=Nc){
			if(adt->L->LC->co[R][C]!=NoV){
				sy=adt->S->type->co[R][C];
				I=adt->T->i_cont[l][R][C];
				kn=K(psi->co[I], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][R][C], adt->S->pa->co[sy][jsat][l],
					adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l],
					adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][R][C]);
				kn=Mean(adt->P->harm_or_arit_mean, UV->U->co[1], UV->U->co[1], k, kn);

				B->co[i] += (kn/UV->U->co[1])*(1.E-3*(psi->co[I] - psi->co[i]) + (adt->T->Z0dp->co[R][C] - adt->T->Z0dp->co[r][c]));
				//B->co[i] += (kn/UV->U->co[1])*(1.E-3*(adt->S->P->co[l][R][C] - adt->S->P->co[l][r][c]) + (adt->T->Z0dp->co[R][C] - adt->T->Z0dp->co[r][c]));

			}
		}
	}
}


