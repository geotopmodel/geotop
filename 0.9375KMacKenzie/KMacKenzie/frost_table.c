
#include "constant.h"
#include "struct.geotop.09375.h"
#include "frost_table.h"

extern T_INIT *UV;
		
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_activelayer(DOUBLEVECTOR *csi, double *d)
// function that calculates the depth of the level of the interface csi=thresh in the soil
// if csi==T && thresh==0 => z_posneg is the active layer depth
// if csi=Psi && thresh==0 => z_negpos is the water table depth
// Author: Matteo Dall'Amico, July 2009 in Trento
{
	double res;
	SHORTVECTOR* index;
	index=new_shortvector(3);
	initialize_shortvector(index,0);
	double thresh=0.0;
	double c_ind;// depth of the center of the layer just above the interface
	double c_buf;// depth of the layer just below the interface
	long i;
	long j=0;// number of layer below the threshold
	//long jnp=0;// negative-positive counter

	// find the indexes of the layer where there is the interface
	for (i=1; i<=(csi->nh)-1; i++){
		if(csi->co[i]>= thresh && csi->co[i+1]<thresh) {
			index->co[j+1]=i;
			j++;
		}
	}
	//printf("ind_posneg=%ld, ind_negpos=%ld", ind_posneg, ind_negpos);stop_execution();
	if (j!=0) {// there is at least one interface
		// the upper soil column is above the threshold and the lower column is below the threshold => positive-negative interface
		i=j;
		c_ind=find_c(index->co[i],d);
		c_buf=find_c(index->co[i]+1,d);
		res=c_ind-(c_ind-c_buf)/(csi->co[index->co[i]]-csi->co[index->co[i]+1])*csi->co[index->co[i]];
	}
	else{
		if(csi->co[1]<thresh){
			res=0.0;
		}else{
			res=UV->V->co[2];
		}
	}
	free_shortvector(index);
	return res;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_watertable(DOUBLEVECTOR *csi, double *d)
// function that calculates the depth of the level of the interface csi=thresh in the soil
// if csi==T && thresh==0 => z_posneg is the active layer depth
// if csi=Psi && thresh==0 => z_negpos is the water table depth
// Author: Matteo Dall'Amico, July 2009 in Trento
{
	double res;
	SHORTVECTOR* index;
	index=new_shortvector(3);
	initialize_shortvector(index,0);
	double thresh=0.0;
	double c_ind;// depth of the center of the layer just above the interface
	double c_buf;// depth of the layer just below the interface
	long i;
	long j=0;// number of layer below the threshold
	//long jnp=0;// negative-positive counter

	// find the indexes of the layer where there is the interface
	for (i=1; i<=(csi->nh)-1; i++){
		if(csi->co[i]<= thresh && csi->co[i+1]>thresh) {
			index->co[j+1]=i;
			j++;
		}
	}
	//printf("ind_posneg=%ld, ind_negpos=%ld", ind_posneg, ind_negpos);stop_execution();
	if (j!=0) {// there is at least one interface
		// the upper soil column is below the threshold and the lower column is above the threshold => negative-positive interface
		i=j;
		c_ind=find_c(index->co[i],d);
		c_buf=find_c(index->co[i]+1,d);
		res=c_ind-(c_ind-c_buf)/(csi->co[index->co[i]]-csi->co[index->co[i]+1])*csi->co[index->co[i]];
	}
	else{
		if(csi->co[1]>thresh){
			res=0.0;
		}else{
			res=UV->V->co[2];
		}
	}
	free_shortvector(index);
	return res;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_c(long indx, double *d)
// function that calculates the depth of the center of a layer with index=indx
// Author: Matteo Dall'Amico, July 2009 in Trento
{
	long i; double res=0.0;
	if (indx==1){
		res=d[1]*0.5;
	}else{
		for (i=1; i<=indx-1; i++){
			res+=d[i];
		}
		res+=d[indx]*0.5;
	}
	return res;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
