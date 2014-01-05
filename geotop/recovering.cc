
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2011 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "recovering.h"
#include "geotop_common.h"

using namespace std;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


//void assign_recovered_map(long n, char *name, DOUBLEMATRIX *assign, PAR *par, DOUBLEMATRIX *Zdistr, DOUBLEMATRIX *Zpoint){
//  void assign_recovered_map(long n, char *name, DOUBLEMATRIX *assign, Par *par, DOUBLEMATRIX *Zdistr, DOUBLEMATRIX *Zpoint){
//
//
//	long r, c;
//	//long i;
////	DOUBLEMATRIX *M;
//	GeoMatrix<double> M;
//	char *temp;
//
//	temp = namefile_i_we2(name, n);
//
//	//if(par->point_sim == 0){
//		//M = read_map(1, temp, Zdistr, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
//		meteoio_readMap(string(temp), M);
//
//	//	for (r=1; r<=M->nrh; r++) {
//		for (r=1; r< M.getRows(); r++) {
//		//	for (c=1; c<=M->nch; c++) {
//			for (c=1; c< M.getCols(); c++) {
//			//	assign->co[r][c] = M->co[r][c];
//				assign->co[r][c] = M[r][c];
//			}
//		}
//
//	/*}else{
//		M = read_map(1, temp, Zpoint, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
//		for(i=1; i<=par->r_points->nh; i++){
//			r = par->r_points->co[i];
//			c = par->c_points->co[i];
//			assign->co[1][i] = M->co[r][c];
//		}
//	}*/
//
////	free_doublematrix(M);
//	free(temp);
//}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

  void assign_recovered_map_vector(long n, std::string name, GeoVector<double>& assign, GeoMatrix<long>& rc, Par *par, GeoMatrix<double>& Zdistr, GeoMatrix<double>& Zpoint){
	long i, r, c;
//	DOUBLEMATRIX *M;
	GeoMatrix<double> M;
    std::string temp;
	
	temp = namefile_i_we2(name, n);
	
	//if(par->point_sim == 0){
	//	M = read_map(1, temp, Zdistr, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
		meteoio_readMap(string(temp), M);
	//	for (i=1; i<=rc->nrh; i++) {
		for (i=1; i<rc.getRows(); i++) {
		//	r = rc->co[i][1];
			r = rc[i][1];
		//	c = rc->co[i][2];
			c = rc[i][2];
		//	assign->co[i] = M->co[r][c];
			assign[i] = M[r][c];
		}
		
	/*}else{
		M = read_map(1, temp, Zpoint, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
		for(i=1; i<=par->r_points->nh; i++){
			r = par->r_points->co[i];
			c = par->c_points->co[i];
			assign->co[i] = M->co[r][c];
		}
	}*/
	
//	free_doublematrix(M);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

  void assign_recovered_map_long(long n, std::string name, GeoMatrix<long>& assign, Par *par, GeoMatrix<double>& Zdistr, GeoMatrix<double>& Zpoint){

 	long r, c;
 	//long i;
 //	DOUBLEMATRIX *M;
 	GeoMatrix<double> M;
      std::string temp;

 	temp = namefile_i_we2(name, n);

 	//	if(par->point_sim == 0){
 	//	M = read_map(1, temp, Zdistr, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
 		meteoio_readMap(string(temp), M);
 	//	for (r=1; r<=M->nrh; r++) {
 		for (r=1; r< M.getRows(); r++) {
 		//	for (c=1; c<=M->nch; c++) {
 			for (c=1; c< M.getCols(); c++) {
 			//	assign->co[r][c] = (long)M->co[r][c];
 				assign[r][c] = (long)M[r][c];
 			}
 		}

 	/*}else{
 		M = read_map(1, temp, Zpoint, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
 		for(i=1; i<=par->r_points->nh; i++){
 			r = par->r_points->co[i];
 			c = par->c_points->co[i];
 			assign->co[1][i] = (long)M->co[r][c];
 		}
 	}*/

 //	free_doublematrix(M);
 }


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
  //void assign_recovered_tensor(long n, char *name, DOUBLETENSOR *assign, PAR *par, DOUBLEMATRIX *Zdistr, DOUBLEMATRIX *Zpoint){
    void assign_recovered_tensor(long lbegin, long n, std::string name, GeoTensor<double>& assign, Par *par, GeoMatrix<double>& Zdistr, GeoMatrix<double>& Zpoint){

//lbegin can be either 0 or 1 depending if the tensor is defined for l==0
		
  	long r, c, l;
//  DOUBLEMATRIX *M;
  	GeoMatrix<double> M;

        std::string temp1, temp2;

// 	for (l=assign->ndl; l<=assign->ndh; l++) {
 	for (l=lbegin; l<assign.getDh(); l++) {  //TODO: assign->ndl  is replaced by 1, need to check!

  		temp1 = namefile_i_we2(name, n);
  		temp2 = namefile_i_we(temp1, l);

  		//	if(par->point_sim == 0){
  		//	M = read_map(1, temp2, Zdistr, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
  			meteoio_readMap(string(temp2), M);
  		//	for (r=1; r<=M->nrh; r++) {
  			for (r=1; r< M.getRows(); r++) {
  			//	for (c=1; c<=M->nch; c++) {
  				for (c=1; c< M.getCols(); c++) {
  				//	assign->co[l][r][c] = M->co[r][c];
  					assign[l][r][c] = M[r][c];
  				}
  			}

  		/*}else{
  			M = read_map(1, temp2, Zpoint, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
  			for(i=1; i<=par->r_points->nh; i++){
  				r = par->r_points->co[i];
  				c = par->c_points->co[i];
  				assign->co[l][1][i] = M->co[r][c];
  			}
  		}*/

  	//	free_doublematrix(M);
  	}
  }


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void assign_recovered_tensor_vector(long n, char *name, DOUBLEMATRIX *assign, LONGMATRIX *rc, PAR *par, DOUBLEMATRIX *Zdistr, DOUBLEMATRIX *Zpoint){
void assign_recovered_tensor_vector(long lbegin, long n, std::string name, GeoMatrix<double>& assign, GeoMatrix<long>& rc, Par *par, GeoMatrix<double>& Zdistr, GeoMatrix<double>& Zpoint){

    	long r, c, i, l;
    //	DOUBLEMATRIX *M;
    	GeoMatrix<double> M;
    std::string temp1, temp2;

    //	for (l=assign->nrl; l<=assign->nrh; l++) {
    	for (l=lbegin; l<assign.getRows(); l++) {

    		temp1 = namefile_i_we2(name, n);
    		temp2 = namefile_i_we(temp1, l);

    		//	if(par->point_sim == 0){
    		//	M = read_map(1, temp2, Zdistr, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
    			meteoio_readMap(temp2, M);
    		//	for (i=1; i<=rc->nrh; i++) {
    			for (i=1; i<rc.getRows(); i++) {
    			//	r = rc->co[i][1];
    				r = rc[i][1];
    			//	c = rc->co[i][2];
    				c = rc[i][2];
    			//	assign->co[l][i] = M->co[r][c];
    				assign[l][i] = M[r][c];
    			}

    		/*}else{
    			M = read_map(1, temp2, Zpoint, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
    			for(i=1; i<=par->r_points->nh; i++){
    				r = par->r_points->co[i];
    				c = par->c_points->co[i];
    				assign->co[l][i] = M->co[r][c];
    			}
    		}*/

    	//	free_doublematrix(M);
    	}
    }


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


   // void assign_recovered_tensor_channel(long n, char *name, DOUBLEMATRIX *assign, LONGVECTOR *r, LONGVECTOR *c, DOUBLEMATRIX *Zdistr){
   void assign_recovered_tensor_channel(long lbegin, long n, std::string name, GeoMatrix<double>& assign,const GeoVector<long> r, const GeoVector<long> c, GeoMatrix<double>& Zdistr){

   	long ch, l;
//  DOUBLEMATRIX *M;
   	GeoMatrix<double> M;
       std::string temp1, temp2;

// 	for (l=assign->nrl; l<=assign->nrh; l++) {
 	for (l=lbegin; l<assign.getRows(); l++) {

   		temp1 = namefile_i_we2(name, n);
   		temp2 = namefile_i_we(temp1, l);

   		//M = read_map(1, temp2, Zdistr, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
   		meteoio_readMap(temp2, M);

   	//	for (ch=1; ch<=r->nh; ch++) {
   		for (ch=1; ch<r.size(); ch++) {
   		//	if(r->co[ch] > 0) assign->co[l][ch] = M->co[r->co[ch]][c->co[ch]];
   			if(r[ch] > 0) assign[l][ch] = M[r[ch]][c[ch]];
   		}

   	//	free_doublematrix(M);
   	}
   }

