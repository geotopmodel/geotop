
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

extern T_INIT *UV;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void water_balance_3D(ALLDATA *adt, double time){

	Richards_3D(adt);
	supflow2(adt->T, adt->L, adt->W, adt->C, adt->P);
	routing2(adt->C);
	output_waterbalance2(adt->W, adt->S, adt->P, adt->L->LC);
		
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Richards_3D(ALLDATA *adt){

	DOUBLEVECTOR *psi, *dpsi, *b, *theta0;
	long i, l, r, c, cont, cont2, iter;
	short sy=0;/* modified by Emanuele Cordano on 24/9/9 */
	double mass0, mass1, massinf, massloss0, massloss, Inf, nw, res, k;
	//FILE *f;
	short out=0;
	
	psi=new_doublevector(Nl*adt->P->total_pixel);
	dpsi=new_doublevector(Nl*adt->P->total_pixel);
	b=new_doublevector(Nl*adt->P->total_pixel);
	theta0=new_doublevector(Nl*adt->P->total_pixel);

	massloss=1.E99;
	
	mass0=0.0;
	for(i=1;i<=Nl*adt->P->total_pixel;i++){
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];
		sy=adt->S->type->co[r][c];	
		mass0+=adt->S->pa->co[sy][jdz][l]*teta_psi(adt->S->P->co[l][r][c],adt->S->thice->co[l][r][c],adt->S->pa->co[sy][jsat][l],adt->S->pa->co[sy][jres][l],
		   adt->S->pa->co[sy][ja][l],adt->S->pa->co[sy][jns][l],1-1/adt->S->pa->co[sy][jns][l],adt->P->psimin2, adt->P->Esoil);
		psi->co[i]=adt->S->P->co[l][r][c];
		theta0->co[i]=teta_psi(psi->co[i], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
			adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin2, adt->P->Esoil);
		if(l==1) adt->S->Jinf->co[r][c]=adt->W->Pn->co[r][c] + adt->W->h_sup->co[r][c]/adt->P->Dt;
	}
		
	mass0/=(double)adt->P->total_pixel;	//in [mm]
	printf("mass0:%f mm\n",mass0);
	
	cont=0;
	
	do{
		
		cont++;

		/*for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(adt->L->LC->co[r][c]!=NoV){
		
					i=adt->T->i_cont[1][r][c];
					sy=adt->S->type->co[r][c];
					k=K(psi->co[i], adt->S->pa->co[sy][jKv][1], adt->P->imp, adt->S->thice->co[1][r][c], adt->S->pa->co[sy][jsat][1], 
						adt->S->pa->co[sy][jres][1], adt->S->pa->co[sy][ja][1], adt->S->pa->co[sy][jns][1], 1-1/adt->S->pa->co[sy][jns][1],
						adt->S->pa->co[sy][jv][1], adt->S->pa->co[sy][jpsimin][1], adt->S->T->co[1][r][c]);
					
					if(adt->S->Jinf->co[r][c]<=k && psi->co[i]<0 && adt->W->h_sup->co[r][c]<0.1){
						adt->W->bc->co[r][c]=0;	//Neumann
					}else{
						adt->W->bc->co[r][c]=1;	//Dirichlet
					}
				}
			}
		}*/

		for(i=1;i<=Nl*adt->P->total_pixel;i++){			
			b->co[i]=Find_b(i, theta0, adt);
		}
				
		iter=tridiag_preconditioned_conjugate_gradient_search(1.E-4, psi, b, Solve_Richards_3D_p, adt);

		/*for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(adt->L->LC->co[r][c]!=NoV){
					
					occur=0;
					
					i=adt->T->i_cont[1][r][c];
					sy=adt->S->type->co[r][c];
					k=K(psi->co[i], adt->S->pa->co[sy][jKv][1], adt->P->imp, adt->S->thice->co[1][r][c], adt->S->pa->co[sy][jsat][1], 
						adt->S->pa->co[sy][jres][1], adt->S->pa->co[sy][ja][1], adt->S->pa->co[sy][jns][1], 1-1/adt->S->pa->co[sy][jns][1],
						adt->S->pa->co[sy][jv][1], adt->S->pa->co[sy][jpsimin][1], adt->S->T->co[1][r][c]);
					
					//Neumann
					if(k + k*(adt->W->h_sup->co[r][c]-psi->co[i])/(0.5*adt->S->pa->co[sy][jdz][1]) > adt->S->Jinf->co[r][c]){
						if(adt->W->bc->co[r][c]==1){
							occur=1;
							adt->W->bc->co[r][c]=0;
						}
					//Dirichlet	
					}else{
						if(adt->W->bc->co[r][c]==0){
							occur=1;
							adt->W->bc->co[r][c]=1;
						}
					}
						
				}
			}
		}
		
		if(occur==1) iter=jacobi_preconditioned_conjugate_gradient_search(1.E-5, psi, b, Solve_Richards_3D_p, adt);*/
		
		/*for(i=1;i<=Nl*adt->P->total_pixel;i++){	
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];					
			printf("i:%ld l:%ld r:%ld c:%ld P:%f P0:%f\n",i,l,r,c,psi->co[i],adt->S->P->co[l][r][c]);
		}*/
		//stop_execution();
			
		for(i=1;i<=Nl*adt->P->total_pixel;i++){	
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];		
			dpsi->co[i]=psi->co[i]-adt->S->P->co[l][r][c];
		}
		
		cont2=0.0;
		nw=1.0;
		massloss0=massloss;
		
		do{
			
			mass1=0.0;
			massinf=0.0;
			for(i=1;i<=Nl*adt->P->total_pixel;i++){
				l=adt->T->lrc_cont->co[i][1];
				r=adt->T->lrc_cont->co[i][2];
				c=adt->T->lrc_cont->co[i][3];				
				psi->co[i]=adt->S->P->co[l][r][c]+nw*dpsi->co[i];
				if(l==1){
					sy=adt->S->type->co[r][c];
					k=K(psi->co[i], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][1l], adt->S->pa->co[sy][ja][l], 
						adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][r][c]);
					Inf=Fmin(adt->S->Jinf->co[r][c], k + k*(adt->W->h_sup->co[r][c]-psi->co[i])/(0.5*adt->S->pa->co[sy][jdz][l]));
					massinf+=Inf*adt->P->Dt;
				}
				mass1+=adt->S->pa->co[sy][jdz][l]*teta_psi(psi->co[i], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
			   		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin2, adt->P->Esoil);
			}
			mass1/=(double)adt->P->total_pixel;//[mm]
			massinf/=(double)adt->P->total_pixel;
			printf("mass1:%f mm\n",mass1);
			printf("massinf:%f mm\n",massinf);
						
			massloss=mass1-mass0-massinf;
			
			cont2++;
			printf("cont:%ld cont2:%ld massloss0:%e massloss1:%e\n",cont,cont2,massloss0,massloss);
			nw/=4.0;
			
		}while(fabs(massloss)>fabs(massloss0) && cont2<5);
		
		res=0.0;
		for(i=1;i<=Nl*adt->P->total_pixel;i++){
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];					
			if(psi->co[i]<adt->P->psimin2) psi->co[i]=adt->P->psimin2;
			res+=pow((psi->co[i]-adt->S->P->co[l][r][c]),2.0);
			adt->S->P->co[l][r][c]=psi->co[i];
		}
		res=pow(res,0.5);
		
		printf("->res:%f massloss:%f ..  %ld %ld\n",res,fabs(massloss),cont,adt->P->MaxiterVWB);
		
		if(res<=adt->P->TolVWb && fabs(massloss)<=0.02) out=1;	//go out of the "do"
																						
	}while(cont<adt->P->MaxiterVWB && out==0);
	
	for(i=1;i<=Nl*adt->P->total_pixel;i++){
		l=adt->T->lrc_cont->co[i][1];
		r=adt->T->lrc_cont->co[i][2];
		c=adt->T->lrc_cont->co[i][3];			
		
		if(l==1){			
			sy=adt->S->type->co[r][c];
			k=K(psi->co[i], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][1l], adt->S->pa->co[sy][ja][l], 
				adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l], adt->S->T->co[l][r][c]);
			Inf=Fmin(adt->S->Jinf->co[r][c], k + k*(adt->W->h_sup->co[r][c]-psi->co[i])/(0.5*adt->S->pa->co[sy][jdz][l]));

			/*if(Inf<0 || adt->W->h_sup->co[r][c]) 
				printf("a)%ld %ld Inf:%f Inf0:%f h:%f Pn:%f psi:%f k:%f\n",r,c,Inf,adt->S->Jinf->co[r][c],adt->W->h_sup->co[r][c],
					adt->W->Pn->co[r][c],psi->co[i],k);*/

			adt->S->Jinf->co[r][c]=Inf;
			adt->W->h_sup->co[r][c]+=(adt->W->Pn->co[r][c]-adt->S->Jinf->co[r][c])*adt->P->Dt;
			
			/*if(Inf<0 || adt->W->h_sup->co[r][c]) 
				printf("b)%ld %ld Inf:%f Inf0:%f h:%f Pn:%f psi:%f k:%f\n",r,c,Inf,adt->S->Jinf->co[r][c],adt->W->h_sup->co[r][c],
					adt->W->Pn->co[r][c],psi->co[i],k);*/	
			
			/*if(psi->co[i]>0 || Inf!=0){
				printf("i:%ld l:%ld r:%ld c:%ld P:%f inf:%f hsup:%f\n",i,l,r,c,psi->co[i],adt->S->Jinf->co[r][c],adt->W->h_sup->co[r][c]);
				stop_execution();
			}*/
		}
		
		if(psi->co[i]!=psi->co[i]) printf("no value psi Richards3D l:%ld r:%ld c:%ld\n",l,r,c);
	}
	
	Assign_Psi(psi, adt->S, adt->L, adt->P, adt->T);
		
	adt->W->out2->co[8]+=massloss*3600.0/adt->P->Dt;

	/*f=fopen(error_file_name,"a");
	fprintf(f,"\nERROR MASS BALANCE: %20.18fmm/h\n",massloss*3600.0/adt->P->Dt);
	fclose(f);	*/

	free_doublevector(psi);
	free_doublevector(b);
	free_doublevector(dpsi);
	free_doublevector(theta0);
		
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Assign_Psi(DOUBLEVECTOR *psi, SOIL *sl, LAND *land, PAR *par, TOPO *top){

	long r, c, l;
	double psisat, th1, theq, th0, oversat, DW, h, ice0, T0;
	short sy;

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				for(l=1;l<=Nl;l++){
					sy=sl->type->co[r][c];
					if(sl->T->co[l][r][c]<=Tfreezing && par->en_balance==1){
						sy=sl->type->co[r][c];
						psisat=psi_saturation(sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
							1-1/sl->pa->co[sy][jns][l]);

						th1=teta_psi(Fmin(psi->co[top->i_cont[l][r][c]],psisat),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
							sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil);

						theq=teta_psi(Psif(sl->T->co[l][r][c]),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],
							1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil);

						if(fabs(theq-th1)>1.E-8){

							th0=teta_psi(Fmin(sl->P->co[l][r][c],psisat),sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil);

							oversat=teta_psi(psi->co[top->i_cont[l][r][c]],sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil) - th1;

							DW=(th1-th0)*(sl->pa->co[sy][jdz][l]);
							h=internal_energy_soil(th0, sl->thice->co[l][r][c], sl->T->co[l][r][c], sl->pa->co[sy][jdz][l], sl->pa->co[sy][jct][l], sl->pa->co[sy][jsat][l]);

							ice0=sl->thice->co[l][r][c];
							T0=sl->T->co[l][r][c];

							from_internal_soil_energy(r, c, l, h+Lf*DW, &th1, &(sl->thice->co[l][r][c]), &(sl->T->co[l][r][c]), sl->pa->co[sy], par->psimin);

							sl->P->co[l][r][c]=psi_teta(th1+oversat,sl->thice->co[l][r][c],sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l],sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l],1-1/sl->pa->co[sy][jns][l],par->psimin2,par->Esoil);

							if(sl->P->co[l][r][c]!=sl->P->co[l][r][c]) printf("no value psi Assign_Psi l:%ld r:%ld c:%ld\n",l,r,c);


						}else{
							sl->P->co[l][r][c]=psi->co[top->i_cont[l][r][c]];
						}

					}else{
						sl->P->co[l][r][c]=psi->co[top->i_cont[l][r][c]];
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

double Solve_Richards_3D(long i, DOUBLEVECTOR *psi, ALLDATA *adt){

	long l, r, c;
	short sy;
	double a = 0.0;
	double C, dz, k, kn, dzn;

	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];

	C=dteta_dpsi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin2, adt->P->Esoil);
	dz=adt->S->pa->co[sy][jdz][l];
	a+=psi->co[i]*(C*dz/adt->P->Dt);

	k=K(adt->S->P->co[l][r][c], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(l==1){
		//Dirichlet
		//if(adt->W->bc->co[r][c]==1) a+=(k/(dz/2.))*psi->co[i];
		if(adt->S->Jinf->co[r][c]>k+k*(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c])/(dz/2.)) a+=(k/(dz/2.))*psi->co[i];
	}

	if(l>1){
		kn=K(adt->S->P->co[l-1][r][c], adt->S->pa->co[sy][jKv][l-1], adt->P->imp, adt->S->thice->co[l-1][r][c], adt->S->pa->co[sy][jsat][l-1], adt->S->pa->co[sy][jres][l-1],
			adt->S->pa->co[sy][ja][l-1], adt->S->pa->co[sy][jns][l-1], 1-1/adt->S->pa->co[sy][jns][l-1], adt->S->pa->co[sy][jv][l-1], adt->S->pa->co[sy][jpsimin][l-1],
			adt->S->T->co[l-1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l-1];
		kn=Harmonic_Mean(dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a+=(kn/dzn)*(psi->co[i]-psi->co[adt->T->i_cont[l-1][r][c]]);
	}

	if(l<Nl){
		kn=K(adt->S->P->co[l+1][r][c], adt->S->pa->co[sy][jKv][l+1], adt->P->imp, adt->S->thice->co[l+1][r][c], adt->S->pa->co[sy][jsat][l+1], adt->S->pa->co[sy][jres][l+1],
			adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 1-1/adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], adt->S->pa->co[sy][jpsimin][l+1],
			adt->S->T->co[l+1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l+1];
		kn=Harmonic_Mean(dz, dzn, k, kn);
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
			kn=Harmonic_Mean(UV->U->co[2], UV->U->co[2], k, kn);
			a+=1.E-3*(kn/UV->U->co[2])*(psi->co[i]-psi->co[adt->T->i_cont[l][r-1][c]]);
		}
	}

	if(r<Nr){
		if(adt->L->LC->co[r+1][c]!=NoV){
			sy=adt->S->type->co[r+1][c];
			kn=K(adt->S->P->co[l][r+1][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r+1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r+1][c]);
			kn=Harmonic_Mean(UV->U->co[2], UV->U->co[2], k, kn);
			a+=1.E-3*(kn/UV->U->co[2])*(psi->co[i]-psi->co[adt->T->i_cont[l][r+1][c]]);
		}
	}

	if(c>1){
		if(adt->L->LC->co[r][c-1]!=NoV){
			sy=adt->S->type->co[r][c-1];
			kn=K(adt->S->P->co[l][r][c-1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c-1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c-1]);
			kn=Harmonic_Mean(UV->U->co[1], UV->U->co[1], k, kn);
			a+=1.E-3*(kn/UV->U->co[1])*(psi->co[i]-psi->co[adt->T->i_cont[l][r][c-1]]);
		}
	}

	if(c<Nc){
		if(adt->L->LC->co[r][c+1]!=NoV){
			sy=adt->S->type->co[r][c+1];
			kn=K(adt->S->P->co[l][r][c+1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c+1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c+1]);
			kn=Harmonic_Mean(UV->U->co[1], UV->U->co[1], k, kn);
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

double Find_b(long i, DOUBLEVECTOR *theta0, ALLDATA *adt){

	long l, r, c;
	short sy;
	double a = 0.0;
	double C, dz, k, kn, dzn, theta1;

	l=adt->T->lrc_cont->co[i][1];
	r=adt->T->lrc_cont->co[i][2];
	c=adt->T->lrc_cont->co[i][3];
	sy=adt->S->type->co[r][c];

	C=dteta_dpsi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin2, adt->P->Esoil);
	theta1=teta_psi(adt->S->P->co[l][r][c], adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l], adt->S->pa->co[sy][ja][l],
		adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->P->psimin2, adt->P->Esoil);
	dz=adt->S->pa->co[sy][jdz][l];
	a+=(C*adt->S->P->co[l][r][c]+theta0->co[i]-theta1)*dz/adt->P->Dt;

	k=K(adt->S->P->co[l][r][c], adt->S->pa->co[sy][jKv][l], adt->P->imp, adt->S->thice->co[l][r][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
		adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
		adt->S->T->co[l][r][c]);

	if(l==1){
		//Dirichlet
		//if(adt->W->bc->co[r][c]==1){
		if(adt->S->Jinf->co[r][c]>k+k*(adt->W->h_sup->co[r][c]-adt->S->P->co[l][r][c])/(dz/2.)){
			a+=(k + k*adt->W->h_sup->co[r][c]/(dz/2.));
		//Neumann
		}else{
			a+=adt->S->Jinf->co[r][c];
		}
	}


	if(l==Nl) a-=k*0.0;

	if(l>1){
		kn=K(adt->S->P->co[l-1][r][c], adt->S->pa->co[sy][jKv][l-1], adt->P->imp, adt->S->thice->co[l-1][r][c], adt->S->pa->co[sy][jsat][l-1], adt->S->pa->co[sy][jres][l-1],
			adt->S->pa->co[sy][ja][l-1], adt->S->pa->co[sy][jns][l-1], 1-1/adt->S->pa->co[sy][jns][l-1], adt->S->pa->co[sy][jv][l-1], adt->S->pa->co[sy][jpsimin][l-1],
			adt->S->T->co[l-1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l-1];
		kn=Harmonic_Mean(dz, dzn, k, kn);
		dzn=0.5*dz + 0.5*dzn;
		a+=kn;
	}

	if(l<Nl){
		kn=K(adt->S->P->co[l+1][r][c], adt->S->pa->co[sy][jKv][l+1], adt->P->imp, adt->S->thice->co[l+1][r][c], adt->S->pa->co[sy][jsat][l+1], adt->S->pa->co[sy][jres][l+1],
			adt->S->pa->co[sy][ja][l+1], adt->S->pa->co[sy][jns][l+1], 1-1/adt->S->pa->co[sy][jns][l+1], adt->S->pa->co[sy][jv][l+1], adt->S->pa->co[sy][jpsimin][l+1],
			adt->S->T->co[l+1][r][c]);
		dzn=adt->S->pa->co[sy][jdz][l+1];
		kn=Harmonic_Mean(dz, dzn, k, kn);
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
			kn=Harmonic_Mean(UV->U->co[2], UV->U->co[2], k, kn);
			a+=1.E-3*(kn/UV->U->co[2])*(adt->T->Z0->co[r-1][c] - adt->T->Z0->co[r][c]);
		}
	}

	if(r<Nr){
		if(adt->L->LC->co[r+1][c]!=NoV){
			sy=adt->S->type->co[r+1][c];
			kn=K(adt->S->P->co[l][r+1][c], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r+1][c], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r+1][c]);
			kn=Harmonic_Mean(UV->U->co[2], UV->U->co[2], k, kn);
			a+=1.E-3*(kn/UV->U->co[2])*(adt->T->Z0->co[r+1][c] - adt->T->Z0->co[r][c]);
		}
	}

	if(c>1){
		if(adt->L->LC->co[r][c-1]!=NoV){
			sy=adt->S->type->co[r][c-1];
			kn=K(adt->S->P->co[l][r][c-1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c-1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c-1]);
			kn=Harmonic_Mean(UV->U->co[1], UV->U->co[1], k, kn);
			a+=1.E-3*(kn/UV->U->co[1])*(adt->T->Z0->co[r][c-1] - adt->T->Z0->co[r][c]);
		}
	}

	if(c<Nc){
		if(adt->L->LC->co[r][c+1]!=NoV){
			sy=adt->S->type->co[r][c+1];
			kn=K(adt->S->P->co[l][r][c+1], adt->S->pa->co[sy][jKh][l], adt->P->imp, adt->S->thice->co[l][r][c+1], adt->S->pa->co[sy][jsat][l], adt->S->pa->co[sy][jres][l],
				adt->S->pa->co[sy][ja][l], adt->S->pa->co[sy][jns][l], 1-1/adt->S->pa->co[sy][jns][l], adt->S->pa->co[sy][jv][l], adt->S->pa->co[sy][jpsimin][l],
				adt->S->T->co[l][r][c+1]);
			kn=Harmonic_Mean(UV->U->co[1], UV->U->co[1], k, kn);
			a+=1.E-3*(kn/UV->U->co[1])*(adt->T->Z0->co[r][c+1] - adt->T->Z0->co[r][c]);
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
	//double q_sup_max;                                     /* maximum superficial flow (in the first part) or the incoming flow in a channel pixel (in the second part)*/
	double b[10];                                         /* area perpendicular to the superficial flow divided by h_sup*/
	double i;											  /*hydraulic gradient*/
	double q,tb,te,dt;
	short lu;


if(par->point_sim==0){	//distributed simulations

	initialize_doublematrix(wat->q_sup,0.0);
	initialize_doublevector(cnet->Q,0.0);
	initialize_doublevector(cnet->Q_sup_spread,0.0);

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
				if(wat->h_sup->co[r][c]>0 && top->DD->co[r][c]>0){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
					//i=top->i_DD->co[r][c];
					if(i<0) i=0;  //correct false slopes
					q=b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;	//mm/s
					//q=wat->h_sup->co[r][c]/dt;
					if(q>0) if(wat->h_sup->co[r][c]/q<dt) dt=wat->h_sup->co[r][c]/q;
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
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(wat->h_sup->co[r][c]>0 && top->DD->co[r][c]>0){
					lu=(short)land->LC->co[r][c];
					Ks=land->ty->co[lu][jcm];
					i=0.001*( wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]] )/b[top->DD->co[r][c]+1] + top->i_DD->co[r][c];
					//i=top->i_DD->co[r][c];
					if(i<0) i=0;  //correct false slopes
					wat->q_sup->co[r][c]=b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(i)*1000.0/dx/dy;
					//if(wat->h_sup->co[r][c]>0)printf("%ld %ld h:%f q:%f\n",r,c,wat->h_sup->co[r][c],wat->q_sup->co[r][c]);
					//wat->q_sup->co[r][c]=wat->h_sup->co[r][c]/dt;
					if(wat->q_sup->co[r][c]!=wat->q_sup->co[r][c]){
						printf("NO VALUE SUP:%ld %ld %ld %ld %f %f\n",r,c,r+r_DD[top->DD->co[r][c]],c+c_DD[top->DD->co[r][c]],wat->h_sup->co[r][c],wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]]);
						printf("i:%f Dh:%f Ks:%f pow:%f iF:%f\n",i,wat->h_sup->co[r][c] - wat->h_sup->co[r+r_DD[top->DD->co[r][c]]][c+c_DD[top->DD->co[r][c]]],Ks,pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m),top->i_DD->co[r][c]);
					}
				}else{
					wat->q_sup->co[r][c]=wat->h_sup->co[r][c]/dt; /*[mm/s]*/
					//if(wat->h_sup->co[r][c]>0)printf("ch. %ld %ld h:%f q:%f\n",r,c,wat->h_sup->co[r][c],wat->q_sup->co[r][c]);
				}
			}
		}

		/*After the computation of the surface flow,these flows are moved trough D8 scheme:*/
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->pixel_type->co[r][c]!=9){

					R=r+r_DD[top->DD->co[r][c]];
					C=c+c_DD[top->DD->co[r][c]];

					/*the superficial flow is subtracted to the land pixels (code 0):*/
					if (top->pixel_type->co[r][c]==0){
						wat->h_sup->co[r][c]-=wat->q_sup->co[r][c]*dt;
						//possible correction
						if(wat->h_sup->co[r][c]<0){
							wat->q_sup->co[r][c]+=wat->h_sup->co[r][c]/dt;
							wat->h_sup->co[r][c]=0.0;
							if(wat->q_sup->co[r][c]!=wat->q_sup->co[r][c]) printf("A\n");
						}
						//if(wat->q_sup->co[r][c]>0)printf("corr %ld %ld qdt:%f h:%f dt:%f\n",r,c,wat->q_sup->co[r][c]*dt,wat->h_sup->co[r][c],dt);
					}

					/*the superficial flow is added to the land pixels (code 0):*/
					if (top->pixel_type->co[R][C]==0) wat->h_sup->co[R][C]+=wat->q_sup->co[r][c]*dt;

					/*the superficial flow is added to flow which flows into the channel pixels (code 10):*/
					if (top->pixel_type->co[R][C]==10){
						for(ch=1;ch<=cnet->r->nh;ch++){
							if(R==cnet->r->co[ch] && C==cnet->c->co[ch]){
								cnet->Q->co[ch]+=wat->q_sup->co[r][c]*dt/par->Dt;
								if(cnet->Q->co[ch]!=cnet->Q->co[ch]){
									printf("qsup no value: r:%ld c:%ld ch:%ld R:%ld C:%ld qsup:%f hsup:%f\n",r,c,ch,R,C,wat->q_sup->co[r][c],wat->h_sup->co[r][c]);
								}
								//if(wat->q_sup->co[r][c]>0) printf("channel r:%ld c:%ld ch:%ld q:%f %f dt:%f\n",R,C,ch,wat->q_sup->co[r][c],cnet->Q->co[ch],dt);
							}
						}
					}
				}
			}
		}

		for(ch=1;ch<=cnet->r->nh;ch++){

			r=cnet->r->co[ch];
			c=cnet->c->co[ch];

			lu=(short)land->LC->co[r][c];
			Ks=land->ty->co[lu][jcm];

			/*computation of value of flow which can go into channel:*/
			wat->h_sup->co[r][c]=Fmax(wat->h_sup->co[r][c], 0.0);
			//q_sup_max=Fmin(2.0*b[top->DD->co[r][c]+1]*Ks*pow(wat->h_sup->co[r][c]/1000.0,1.0+par->gamma_m)*sqrt(top->i_ch->co[r][c])*1000.0/dx/dy,cnet->Q->co[ch]);
			//q_sup_max=cnet->Q->co[ch];
			//if (top->DD->co[r][c]==10) q_sup_max=cnet->Q->co[ch];	//closure section
			//if(q_sup_max!=q_sup_max) printf("q_sup_max no value: qsup:%f r:%ld c:%ld ch:%ld h:%f ich:%f\n",cnet->Q->co[ch],r,c,ch,wat->h_sup->co[r][c],top->i_ch->co[r][c]);

			/*wat of channel pixel which becames lateral flow and doesn't go into the channel;
			the variable cnet->q_sup is use in this case improperly:*/
			//cnet->Q->co[ch]-=q_sup_max;

			/*the water on surface in channel pixel is updated:*/
			//wat->h_sup->co[r][c]-=wat->q_sup->co[r][c]*dt;
			//wat->h_sup->co[r][c]+=cnet->Q->co[ch]*dt;

			/*now is assingned to cnet->q_sup is right value:*/
			//cnet->Q->co[ch]=q_sup_max;

			/*The new water arrived in the channel-network is spread immediately along a
			virtual coordinate s (=distance from the outlet) so that it can travel to the
			outlet, after simple traslation, compling with the convolution of the De
			Saint Venant equations. This is made either for superficial flow or for the
			subsuperficial flow:*/
			for(s=1;s<=cnet->Q_sup_spread->nh;s++){
				cnet->Q_sup_spread->co[s]+=(cnet->Q->co[ch])*(cnet->fraction_spread->co[ch][s])*0.001*dx*dy; /*in mc/s*/
				//if(cnet->Q_sup_spread->co[s]>0)printf("s:%ld Q:%f\n",s,cnet->Q_sup_spread->co[s]);
			}

		}


	}while(te<par->Dt);

	/*for(s=1;s<=cnet->Q_sup_spread->nh;s++){
		if(cnet->Q_sup_spread->co[s]>0) printf("Spread s:%ld Q:%f\n",s,cnet->Q_sup_spread->co[s]);
	}*/


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

	/*The just calculated Q_sup_spread and Q_sub_spread are added to the flow already
	present in the channel-network(Q_sup_s and Q_sub_s) and at once it is made a
	traslation of Q_sup_s e Q_sub_s towards the outlet:*/
	for(s=1;s<=cnet->Q_sup_spread->nh-1;s++){
		//cnet->Q_sub_s->co[s]=cnet->Q_sub_s->co[s+1]+cnet->Q_sub_spread->co[s];
		cnet->Q_sup_s->co[s]=cnet->Q_sup_s->co[s+1]+cnet->Q_sup_spread->co[s];
		//if(cnet->Q_sup_s->co[s]>0)printf("->s:%ld Q:%f\n",s,cnet->Q_sup_s->co[s]);
	}
	//cnet->Q_sub_s->co[cnet->Q_sub_spread->nh]=cnet->Q_sub_spread->co[cnet->Q_sub_spread->nh];
	cnet->Q_sup_s->co[cnet->Q_sup_spread->nh]=cnet->Q_sup_spread->co[cnet->Q_sup_spread->nh];
	//if(cnet->Q_sup_s->co[cnet->Q_sup_spread->nh]>0)printf("->s:%ld Q:%f\n",cnet->Q_sup_spread->nh,cnet->Q_sup_s->co[cnet->Q_sup_spread->nh]);
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

