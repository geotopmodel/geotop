/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file g_raster2plvector.h

Copyright, 2009 Emanuele Cordano and Riccardo Rigon

This file is part of BGEOMETRY.
 BGEOMETRY is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BGEOMETRY is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "turtle.h"
#include "t_alloc.h"


/* These functions allow to use basin with a generic shape.
 * WARNING: The mask is the raster which contains all the basin, the boundary of the basin must contain no-value elements.
 *
 *
 * right vertical line
 */


#include "geometry.h"
#include "geometry_freememory.h"
#include "geometry_io.h"

#include "rw_maps.h"
//#include "gridded.element.input.geotop.h"

#include "g_raster2plvector.h"

#define NO_ELEVATION 0.0
#define CENTER -0.5
#define CENTER -0.5
#define TOP -1.0
#define BOTTOM 0.0
#define LEFT -1.0
#define RIGHT 0.0

#define IVERTEX 10000

#define INIT_VALUE -999

#define FLOATING_POINT_TYPE 0
#define MAP_FORMAT 2



//LONGBIN *addresses(LONGMATRIX *lmask,DOUBLEVECTOR *V, int (* check_novalue2)(double x,double y,DOUBLEVECTOR *V),int reply,long row_shift,long col_shift){
/** \param V (DOUBLEVECTOR *) - vector regarding the query appied (searching novalue)
	 * \param (* check_novalue)(double x,DOUBLEVECTOR *V) (int) - query expressed as a function returning a long integer
	 * \param reply - expected reply from the query to save the element addresses
	 */


//int create_raster_index_matrices(LONGMATRIX *lmask,long novalue, LONGMATRIX *i_pixels,LONGMATRIX *i_horizontal_lines, LONGMATRIX *i_vertical_lines, LONGMATRIX *i_vertex){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date January 2008
	 *
*/

LONGBIN *addresses(LONGMATRIX *lmask,long row_shift,long col_shift){
	/*
	 *
	 *
	 * \param lmask (DOUBLEMATRIX *) - matrix wth 0 or 1 values (if element is 1 the element is within the basin, 0 otherwise)
	 * \param row_shitf (long) - integer parameters to be set for vertical line searching (-1: left vertical line searching, 0: pixel searching, 1: right vertical line searching)
	 * \param col_shift (long) - integer parameters to be set for horizonttal line searching (-1: top horizontal line searching, 0: pixel searching, 1: bottom vertical line searching)
	 *
	 *\return a LONGBIN containing for each rows the number of columns at wihch the element retrun the query equal to the variable reply
	 * \author Emanuele Cordano
	 * \date January 2009
	 *
	 *
	 */
	LONGVECTOR *ndata;

	LONGBIN *LB;
	long r,c,j,data_counter;

	printf("ENTRANCE:  %ld %ld %ld %ld ",row_shift,col_shift,lmask->nrh,lmask->nch);

	if ((row_shift<-1)|| (row_shift>1) || (col_shift<-1)|| (col_shift>1) ) {
		printf ("Warning in addresses function source file g_raster2plvector_novalue_manegement: parameters row_shift or col_shift exceed 1 (in absolute value) \n");
		printf ("This could cause bus errors!!\n ");
	}


	for (r=lmask->nrl;r<=lmask->nrh;r++) {
		if ((lmask->element[r][lmask->ncl]==1) || (lmask->element[r][lmask->nch]==1)) printf ("Warning in addresses function source file g_raster2plvector_novalue_manegement: at row %ld the mask of the basin reaches the border !\n",r);
	}
	for (c=lmask->ncl;c<=lmask->nch;c++) {
		if ((lmask->element[lmask->nrl][c]==1) || (lmask->element[lmask->nrh][c]==1)) printf ("Warning in addresses function source file g_raster2plvector_novalue_manegement: at column %ld the mask of the basin reaches the border!\n",c);
	}
	ndata=new_longvector(lmask->nrh);
	ndata->element[ndata->nl]=0;
//	ndata->element[ndata->nh]=0;
	for (r=lmask->nrl+1;r<=lmask->nrh;r++){
		data_counter=0;
		for (c=lmask->ncl+1;c<=lmask->nch;c++){
			if ((lmask->element[r][c]==1) || (lmask->element[r][c+col_shift]==1) || (lmask->element[r+row_shift][c]==1) || (lmask->element[r+row_shift][c+col_shift]==1))
				data_counter++;

		}
		ndata->element[r]=data_counter;
	}
//	print_longvector_elements(ndata,PRINT);
//	stop_execution();
	LB=new_longbin(ndata);
//	printf("OKK pl");
//	stop_execution();
	for (r=lmask->nrl+1;r<=lmask->nrh;r++){
		j=0;
		for (c=lmask->ncl+1;c<=lmask->nch;c++){
			if ((lmask->element[r][c]==1) || (lmask->element[r][c+col_shift]==1) || (lmask->element[r+row_shift][c]==1) || (lmask->element[r+row_shift][c+col_shift]==1)) {
				j++;
				LB->element[r][j]=c;
			}
		}

	}
	return LB;

}
//COME FARE CON I VERTICI???



long index_pixel_from_a_bin(long r, long c, LONGVECTOR *s_index){
	/*
	 *\param lmask (DOUBLEMATRIX *) - matrix with 0 or 1 values (if element is 1 the element is within the basin, 0 otherwise)
	 *\param row_shitf (long) - integer parameters to be set for vertical line searching (-1: left vertical line searching, 0: pixel searching, 1: right vertical line searching)
	 *\param col_shift (long) - integer parameters to be set for horizonttal line searching (-1: top horizontal line searching, 0: pixel searching, 1: bottom vertical line searching)

	 * \param - (long) r row
	 * \param - (long) c column
	 * \param - (LONGVECTOR *) sum of the index of the bin
	 *
	 * \auhor Emanuele Cordano
	 * \date January 2009
	 *
	 *\return index_pixel an index related
	 *
	 */
	long l;

	if (r%2==0) {
		l=s_index->element[r]-c+1;
	}else {
		if (r>1) {
			l=s_index->element[r-1]+c;
		} else {
			l=c;
		}


	}
	return l;

}

LONGMATRIX *m_indices_with_novalues(LONGBIN *laddresses, long nch, long novalue, long IBASE, long (*index_pixel_from_a_bin)(long r, long c, LONGVECTOR *s_index)){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date January 2009
	 *
	 * \param (LONGBIN *) - the bin of the addresses of pixels inside the domain
	 * \param (long)      - number of columnus for the requested longmatrix
	 * \param (long)       - integer values for NULL
	 * \param (long)      - IBASE base number
	 * \param (long)       -  (*index_pixel_from_a_bin)(long r, long c, LONGVECTOR *s_index) - function which numbers the pixel inside the domain
	 *
	 * \brief It creates a longmatrix with cell index (for pixels inside the domain)  and novalue (for pixels outside the domain)
	 */
	LONGMATRIX *m_ind;
	LONGVECTOR *s_index;
	long r,c,j,l;

	s_index=new_longvector(laddresses->index->nh);
	s_index->element[s_index->nl]=laddresses->index->element[s_index->nl];
	for (l=s_index->nl+1;l<=s_index->nh;l++){
		if (laddresses->index->element[l]>nch) printf("Error: in m_indices_with_novalues at row %ld , %ld exceeds number of columns %ld \n",l,laddresses->index->element[l],nch);
		s_index->element[l]=s_index->element[l-1]+laddresses->index->element[l];
	}
//	stop_execution();
//	printf("22novalue %ld :\n",novalue);
//	stop_execution();
	m_ind=new_longmatrix(laddresses->index->nh,nch);
	for(r=m_ind->nrl;r<=m_ind->nrh;r++) {
		for(c=m_ind->ncl;c<=m_ind->nch;c++) {
			m_ind->element[r][c]=novalue;
	//		printf("m_IND[%ld][%ld]: %ld\n",r,c,m_ind->element[r][c]);

		}
		for(j=m_ind->ncl;j<=laddresses->index->element[r];j++){
			m_ind->element[r][laddresses->element[r][j]]=(*index_pixel_from_a_bin)(r,j,s_index)+IBASE;
//			printf("m_IND[%ld][%ld]: %ld\n",r,laddresses->element[r][j],m_ind->element[r][laddresses->element[r][j]]);
		}


	}


	return m_ind;




}


/*
 * LONGMATRIX *m_indices_with_novalues(LONGBIN *addresses, long nch, long novalue, long IBASE, long (*index_pixel_from_a_bin)(long r, long c, LONGVECTOR *s_index))
 * LONGBIN *addresses(LONGMATRIX *lmask,long row_shift,long col_shift)
 */

LONGMATRIX *m_indices_from_mask(DOUBLEMATRIX *mask,long row_shift,long col_shift,long novalue,long IBASE,long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index),DOUBLEVECTOR *V, int (*check_novalues)(double x, DOUBLEVECTOR *V)){
	/*
	 *
	 * \param mask (DOUBLEMATRIX *) - mask of the basin
	 * \param row_shitf (long) - integer parameters to be set for vertical line searching (-1: left vertical line searching, 0: pixel searching, 1: right vertical line searching)
	 * \param col_shift (long) - integer parameters to be set for horizonttal line searching (-1: top horizontal line searching, 0: pixel searching, 1: bottom vertical line searching)
	 * \param (long)       - integer values for NULL
	 * \param (long)      - IBASE base number
	 * \param (long)       -  (*index_pixel_from_a_bin)(long r, long c, LONGVECTOR *s_index) - function which numbers the pixel inside the domain
	 * \param V (DOUBLEVECTOR *)- vector containing novalue information
	 * \param int (*check_novalues)(double x, DOUBLEVECTOR *V)) -
	 *
	 * \author Emanuele Cordano
	 * \date January 2008
	 *
	 *
	 */
	LONGBIN *LB;
	LONGMATRIX *m_ind;
	LONGMATRIX *lmask;
	long r,c;

	lmask=new_longmatrix(mask->nrh,mask->nch);
	printf("1: lmask %ld %ld:\n",lmask->nrh,lmask->nch);
//	stop_execution();
	for (r=mask->nrl;r<=mask->nrh;r++) {
		for (c=mask->ncl;c<=mask->nch;c++) {
	//		printf("value:%lf at  %ld,%ld \n",mask->element[r][c],r,c);
			if ((*check_novalues)(mask->element[r][c],V)==0) {
	//			printf("NONNULL value:%lf at  %ld,%ld \n",mask->element[r][c],r,c);
				lmask->element[r][c]=1;
	//			stop_execution();
			}else if ((*check_novalues)(mask->element[r][c],V)==1) {
				lmask->element[r][c]=0;
		//		printf("NULL value:%lf at  %ld,%ld \n",mask->element[r][c],r,c);
			}
		}
	}
	printf("lmask %ld %ld :\n",lmask->nrh,lmask->nch);


	//printf("novalue %ld :\n",novalue);
	//stop_execution();
	LB=addresses(lmask,row_shift,col_shift);
	m_ind=m_indices_with_novalues(LB,lmask->nch,novalue,IBASE,index_pixel_from_a_bin);
	printf("STOP!\n");
	free_longbin(LB);
	free_longmatrix(lmask);

	//stop_execution();
	return m_ind;

}
//LONGMATRIX *indices_conditioned(long nrh, long nch,long IBASE,long (*t_index)(long r, long c,long nrh, long nch)



// POINT *new_point_from_raster(long r,long c, long nrh, long nch, double lx, double ly, double nsres, double ewres, double blc_x, double blc_y, long NBASE,long (*t_index)(long r, long c,long nrh, long nch)) {

POINT *new_point_from_raster(long r,long c, long nrh, long nch, double lx, double ly, double nsres, double ewres, double blc_x, double blc_y, long index) {
	/*
	 *
	 * \param r - (long) r row
	 * \param c - (long) c column
	 * \param nrh - (long) number of rows
	 * \param nch - (long) number of columns
	 * \param lx - (double) lx dimensionless x displacement for point location
	 * \param ly - (double)  ly dimensionless y displacement for point location
	 * \param nsres - (double) north-south resolution
	 * \param ewres - (double) east-west resolution
	 * \param blc_x - (double) x coordinate for the bottom left corner
	 * \param blc_y - (double) y coordinate for the bottom left corner
	 * \param index - (long) index of the point.
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 * \return a point struct located in a specific point of the pixel r,c
	 *
	 */

	POINT *P;
	double x,y;

//	index=(*t_index)(r,c,nrh,nch)+NBASE;
	x=((double)c+lx)*ewres+blc_x;
	y=((double)(nrh-r+1)+ly)*nsres+blc_x;

	P=new_point(index,x,y,NO_ELEVATION);

	return P;

}

LINE *new_horizontal_line_from_raster(long r,long c, long nrh, long nch, double nsres, double ewres, double blc_x, double blc_y, long iline, long ivertex1, long ivertex2){


	/*
		 *
		 * \param r - (long) r row
		 * \param c - (long) c column
		 * \param nrh - number of rows
		 * \param nch - number of columns
		 * \param lx - (double) lx dimensionless x displacement for point location
		 * \param ly - (double)  ly dimensionless y displacement for point location
		 * \param nsres - (double) north-south resolution
		 * \param ewres - (double) east-west resolution
		 * \param blc_x - (double) x coordinate for the bottom left corner
		 * \param blc_y - (double) y coordinate for the bottom left corner
		 * \param iline - (long) index number for line;
		 * \param ivertex1 - (long) index number for the vertex 1
		 * \param ivertex2 - (long) index number for the vertex 2
		 *
		 *
		 * \author Emanuele Cordano
		 * \date November 2008
		 *
		 * \return a line struct for an horizontal line located in a specic point of the pixel r,c
		 *
		 */

	LINE *line;
	POINT *P1, *P2;


//	index=(*t_index)(r,c,nrh,nch)+NBASE;
	P1=new_point_from_raster(r,c,nrh,nch,LEFT,TOP,nsres,ewres,blc_x,blc_y,ivertex1);
	P2=new_point_from_raster(r,c,nrh,nch,RIGHT,TOP,nsres,ewres,blc_x,blc_y,ivertex2);

	line=new_line_from_points(iline,P1,P2);


	return line;

}


LINE *new_vertical_line_from_raster(long r,long c, long nrh, long nch, double nsres, double ewres, double blc_x, double blc_y, long iline, long ivertex1, long ivertex2){


	/*
		 *
		 * \param r - (long) r row
		 * \param c - (long) c column
		 * \param nrh - number of rows
		 * \param nch - number of columns
		 * \param lx - (double) lx dimensionless x displacement for point location
		 * \param ly - (double)  ly dimensionless y displacement for point location
		 * \param nsres - (double) north-south resolution
		 * \param ewres - (double) east-west resolution
		 * \param blc_x - (double) x coordinate for the bottom left corner
		 * \param blc_y - (double) y coordinate for the bottom left corner
		 * \param iline - (long) index number for line;
		 * \param ivertex1 - (long) index number for the vertex 1
		 * \param ivertex2 - (long) index number for the vertex 2
		 *
		 *
		 * \author Emanuele Cordano
		 * \date November 2008
		 *
		 * \return a line struct for an horizontal line located in a specic point of the pixel r,c
		 *
		 */

	LINE *line;
	POINT *P1, *P2;
	long index;

//	index=(*t_index)(r,c,nrh,nch)+NBASE;
	P1=new_point_from_raster(r,c,nrh,nch,LEFT,BOTTOM,nsres,ewres,blc_x,blc_y,ivertex1);
	P2=new_point_from_raster(r,c,nrh,nch,LEFT,TOP,nsres,ewres,blc_x,blc_y,ivertex2);
//	printf("index=%ld r=%ld c=%ld  \n",iline,r,c);
//stop_execution();
	line=new_line_from_points(iline,P1,P2);

	free_point(P1);
	free_point(P2);

	return line;

}


LINEVECTOR *get_linevector_from_raster_grid(LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,LONGMATRIX *i_vertex, double nsres, double ewres, double blc_x, double blc_y, long novalue){

/*
 * \author Emanuele Cordano
 * \date November 2008
 *
 * \param i_horizontal  - (LONGMATRIX *) matrix for horizontal lines
 * \param i_horizontal  - (LONGMATRIX *) matrix for vertical lines
 * \param i_vertex      -
 * \param lx - (double) lx dimensionless x displacement for point location
 * \param ly - (double)  ly dimensionless y displacement for point location
 * \param nsres - (double) north-south resolution
 * \param ewres - (double) east-west resolution
 * \param blc_x - (double) x coordinate for the bottom left corner
 * \param blc_y - (double) y coordinate for the bottom left corner
 * \param novalue - (long) null value
 * \return a linevector with vertical and horizontal line ordering assigned
 *
 *\\OCCORRE TENERE CONTO DEL NOVALUE!!!
 *
 */

LINEVECTOR *lines;
long r,c,i,count;
/* Counting lines */
count=0;
/* VERTICA LINES */
for(r=i_vertical->nrl;r<=i_vertical->nrh;r++){
	for (c=i_vertical->ncl;c<=i_vertical->nch;c++){
		/* vertical lines at the left edge of the grid element*/

		if (i_vertical->element[r][c]!=novalue) count++;
	}
}
//stop_execution();
/* HORIZONTAL LINES */
/* VERTICA LINES */
for(r=i_horizontal->nrl;r<=i_horizontal->nrh;r++){
	for (c=i_horizontal->ncl;c<=i_horizontal->nch;c++){
		/* horizontal lines at the upper limit of the grid element */
		if (i_horizontal->element[r][c]!=novalue) count++ ;
	}
}
/* verify of lines element */
printf("OK HERE IN LINEVECTOR");
stop_execution();
lines=new_linevector(count);
count=0;

/* VERTICA LINES */
for(r=i_vertical->nrl;r<=i_vertical->nrh;r++){
	for (c=i_vertical->ncl;c<=i_vertical->nch;c++){
		/* vertical lines at the left edge of the grid element*/
//		printf("indexxx=%ld r=%ld c=%ld  \n",i_vertical->element[r][c],r,c);
		if (i_vertical->element[r][c]!=novalue) lines->element[i_vertical->element[r][c]]=new_vertical_line_from_raster(r,c,i_vertical->nrh,i_vertical->nch,nsres,ewres,blc_x,blc_y,i_vertical->element[r][c],i_vertex->element[r+1][c],i_vertex->element[r][c]);
	}
}
//stop_execution();
/* HORIZONTAL LINES */
/* VERTICA LINES */
for(r=i_horizontal->nrl;r<=i_horizontal->nrh;r++){
	for (c=i_horizontal->ncl;c<=i_horizontal->nch;c++){
		/* horizontal lines at the upper limit of the grid element */
		if (i_horizontal->element[r][c]!=novalue) lines->element[i_horizontal->element[r][c]]=new_horizontal_line_from_raster(r,c,i_horizontal->nrh,i_horizontal->nch,nsres,ewres,blc_x,blc_y,i_horizontal->element[r][c],i_vertex->element[r][c],i_vertex->element[r][c+1]);
	}
}
/* verify of lines element */

for (i=lines->nl;i<=lines->nh;i++){
	if (!lines->element[i]) printf("Warning:: element %ld corresponding of lines was not assigned !!\n",i);
}

return lines;
}

POLYGON *new_pixel_from_raster(long index,long r, long c ,LINEVECTOR *lines, LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,double nsres, double ewres, double blc_x, double blc_y,short print) {
	/*
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 *\param index (long) - index of the element
	 *\param r     (long) - row
	 *\param c     (long) - colums
	 *\param lines (LINEVECTOR *) - array of lines
	 *\param i_orizontal (LONGMATRIX *)  - matrix of indices for horizontal lines
	 *\param i_vertical  (LONGMATRIX *)  - matrix of indices for vertical lines
	 *\param nsres - (double) north-south resolution
	 *\param ewres - (double) east-west resolution
	 *\param blc_x - (double) x coordinate for the bottom left corner
	 *\param blc_y - (double) y coordinate for the bottom left corner
	 *\param print - (short) it is a verbose modality if print==1
	 *
	 *
	 *
	 * \return a polygon (Attributes are not allocated!!)
	 *
	 */

	LINEVECTOR *edges;
	POLYGON *PO;
	LONGVECTOR *ledges;
	POINT *centroid;
	double x,y,z;
	long nrh,nch;

	nrh=i_vertical->nrh;
	nch=i_horizontal->nch;


	centroid=new_point_from_raster(r,c,nrh,nch,CENTER,CENTER,nsres,ewres,blc_x,blc_y,index);
    /* indices of the edge of the pixel */

	ledges=new_longvector(4);


	if ((!i_vertical->element[r][c]) || (i_vertical->element[r][c]<=0)) {
		printf ("Warning: left vertical line missing  (polygon %ld) at r=%ld c=%ld \n ",index,r,c);
		ledges->element[1]=i_vertical->element[r][c];
	} else {
		ledges->element[1]=i_vertical->element[r][c];
//		printf ("DEBUG ON VERTICAL LINE %ld  (polygon %ld) at r=%ld c=%ld \n ",i_vertical->element[r][c],index,r,c);

	}
	if ((!i_vertical->element[r][c+1]) || (i_vertical->element[r][c+1]<=0)) {
		printf ("Warning: right vertical line missing (polygon %ld) at r=%ld c=%ld \n ",index,r,c);
		ledges->element[2]=i_vertical->element[r][c+1];
	} else {
		ledges->element[2]=i_vertical->element[r][c+1];
//		printf ("DEBUG ON VERTICAL LINE %ld  (polygon %ld) at r=%ld c=%ld \n ",i_vertical->element[r][c+1],index,r,c+1);
	}

	if ((!i_horizontal->element[r][c]) || (i_horizontal->element[r][c]<=0)) {
		printf ("Warning: top horizontal line missing (polygon %ld) at r=%ld c=%ld \n",index,r,c);
		ledges->element[3]=i_horizontal->element[r][c];
	} else {
		ledges->element[3]=i_horizontal->element[r][c];

	}


	if ((!i_horizontal->element[r+1][c]) || (i_horizontal->element[r+1][c]<=0)) {
		printf ("Warning: bottom horizontal line missing (polygon %ld) at r=%ld c=%ld \n",index,r,c);
		ledges->element[4]=i_horizontal->element[r+1][c];
	} else {
		ledges->element[4]=i_horizontal->element[r+1][c];

	}


	if (print==1) printf(" Polygon : %ld [r= %ld , c = %ld ] lines (edges): %ld, %ld, %ld, %ld  \n",index,r,c,ledges->element[1],ledges->element[2],ledges->element[3],ledges->element[4]);

	edges=extract_linvector_from_linevector(ledges,lines);
	PO=new_polygon_from_a_linevector(edges,centroid);


	if (print==1) printf(" Polygon : %ld [r= %ld , c = %ld ] lines (edges): %ld, %ld, %ld, %ld was allocated \n",index,r,c,ledges->element[1],ledges->element[2],ledges->element[3],ledges->element[4]);
	free_point(centroid);
	free_longvector(ledges);
	free_linevector(edges);

	return PO;


}

POLYGONVECTOR *get_polygonvector_from_raster(LINEVECTOR *lines,LONGMATRIX *i_pixels,LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,double nsres, double ewres, double blc_x, double blc_y,long novalue, short print) {
	/*
	 * \author Emanuele Cordano
	 * \date November 2008
	 *

	 *\param lines (LINEVECTOR *) - array of lines
	 *\param i_pixels   (LONGMATRIX *) - matrix for pixel indices
	 *\param i_orizontal (LONGMATRIX *)  - matrix of indices for horizontal lines
	 *\param i_vertical  (LONGMATRIX *)  - matrix of indices for vertical lines
	 *\param nsres - (double) north-south resolution
	 *\param ewres - (double) east-west resolution
	 *\param blc_x - (double) x coordinate for the bottom left corner
	 *\param blc_y - (double) y coordinate for the bottom left corner
	 *\param novalue (long) NULL value
	 *\param print - (short) it is a verbose modality if print==1
	 *
	 *
	 */
	POLYGONVECTOR *pv;
	long i,r,c,count;
	/* counter for the pixels within te basin    */
	count=0;

	for(r=i_pixels->nrl;r<=i_pixels->nrh;r++){
		for(c=i_pixels->ncl;c<=i_pixels->nch;c++){
			if (i_pixels->element[r][c]!=novalue) count++;

		}
	}

	pv=new_polygonvector(count);
	count=0;

	for(r=i_pixels->nrl;r<=i_pixels->nrh;r++){
		for(c=i_pixels->ncl;c<=i_pixels->nch;c++){
			if (i_pixels->element[r][c]!=novalue) {
				i=i_pixels->element[r][c];
				pv->element[i]=new_pixel_from_raster(i,r,c,lines,i_horizontal,i_vertical,nsres, ewres,blc_x,blc_y,print);
			}
		}
	}




	/* verify of polygonvector element */

	for (i=pv->nl;i<=pv->nh;i++){
		if (!pv->element[i]) printf("Warning:: element %ld corresponding of polygons was not assigned !!\n",i);
	}

	return pv;

}

/* Functions which converts DOUBLEVECTOR and DOUBLEMATRIX according to a given function */

DOUBLEVECTOR *get_doublevector_from_doublematrix(LONGMATRIX *indices,DOUBLEMATRIX *M, double novalue){
	/*
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 * \param  indices (LONGMATRIX *) - matrix of the indices;
	 * \param  M (DOUBLEMATRIX *) - matrix to be transformed into a vector;
	 * \param  novalue (double) - null value
	 * \return a vector
	 *
	 */
	DOUBLEVECTOR *v;
	long r,c,i,cnt;

	if ((indices->nrh!=M->nrh) || (indices->nch!=M->nch)) printf("Error:: in get_doublevector_from_doublematrix indices [%ld,%ld] and M [%ld,%ld] has different sizes! \n",indices->nrh,indices->nch,M->nrh,M->nch);
	cnt=0;
	for (r=M->nrl;r<=M->nrh;r++){
		for (c=M->ncl;c<=M->nch;c++){
					if (M->element[r][c]!=novalue) cnt++;
		}
	}

	v=new_doublevector(cnt);
	for (i=v->nl;i<=v->nh;i++){
		v->element[i]=INIT_VALUE;
	}
	if (cnt>M->nrh*M->nch) printf ("Error:: Error:: in get_doublevector_from_doublematrix vector elements [%lf] execeed doublematrix elements!!\n",v->nh,M->nrh*M->nch);
	for (r=M->nrl;r<=M->nrh;r++){
		for (c=M->ncl;c<=M->nch;c++){
			i=indices->element[r][c];
			if ((i<v->nl) && (i>v->nh))  {
//				printf ("Error:: in get_doublevector_from_doublevector index %ld exceeds size of matrix [%ld,%ld] at %ld,%ld",i,M->nrh,M->nch,r,c);
			}else{
				v->element[i]=M->element[r][c];
			}
		}
	}

	for (i=v->nl;i<=v->nh;i++){
		if (v->element[i]==INIT_VALUE) printf("Error:: in get_doublevector_from_doublevector index %ld  was not assigned (%lf) !!\n",i,v->element[i]);
	}

	return v;
}

//*read_map(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref)

DOUBLEMATRIX *get_doublematrix_from_doublevector(DOUBLEVECTOR *v, LONGMATRIX *indices, DOUBLEMATRIX *Mref, double novalue) {
/*
 *
 *
 * \author Emanuele Cordano
 * \date November 2008
 *
 * \param v - (DOUBLEVECTOR *) - vector which will be trasformed into a raster
 * \param indices  (LONGVECTOR *) - matrix of indices
 * \param Mref - (DOUBLEMATRIX *) - reference raster
 * \param novalue - (double) - no value
 *
 *\return
 *
 */

DOUBLEMATRIX *M;
long i,r,c;

M=new_doublematrix(Mref->nrh,Mref->nch);

if ((indices->nrh!=M->nrh) || (indices->nch!=M->nch)) printf("Error:: in get_doublematrix_from_doublevector indices [%ld,%ld] and M [%ld,%ld] has different sizes! \n",indices->nrh,indices->nch,M->nrh,M->nch);
//if (v->nh!=nrh*nch) {
//	printf("Error::in get_doublematrix_from_doublevector number of elemets in vector %lf does not correspond to the number of rows and column %lf and %lf. The matrix was not created!!\n",v->nh,nrh,nch);
//	return M;
//}

for (r=M->nrl;r<=M->nrh;r++){
		for (c=M->ncl;c<=M->nch;c++){
			i=indices->element[r][c];
			M->element[r][c]=INIT_VALUE;
			if ((i>v->nh) || (i<v->nl))  {
				if (Mref->element[r][c]!=novalue) printf ("Error:: in get_doublematrix_from_doublevector index %ld exceeds size of matrix [%ld,%ld] at %ld,%ld \n",i,M->nrh,M->nch,r,c);
			}else{
				M->element[r][c]=v->element[i];
			}
			if ((M->element[r][c]==INIT_VALUE) && (Mref->element[r][c]==novalue)) M->element[r][c]=novalue;
	}
}

for (r=M->nrl;r<=M->nrh;r++){
		for (c=M->ncl;c<=M->nch;c++){
			if (M->element[r][c]==INIT_VALUE) printf("Error:: in get_doublematrix_from_doublevector matrix element %ld, %ld  [%ld,%ld] was not assigned \n",r,c,M->nrh,M->nch);
		}
	}


return M;

}

DOUBLEVECTOR *read_doublevector_from_raster(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref,LONGMATRIX *indices){
////*read_map(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref)
	/*
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 * \param a (short) - flag about the use of the check with the reference map see read_map comment
	 * \param filename (char *) name of the file
	 * \param Mref - (DOUBLEMATRIX *) - reference map
	 * \param UVref - (T_INIT *) T_IIT struct containg ewres and nwres information
	 * \param  (*t_index)(long r, long c,long nrh, long nch) -(long) reordering  function
	 *
	 * \brief it reads a doublematrix with "read_map" functions and it gets an re-ordered doublevector
	 *
	 * \details see paramentes for "read_map"
	 *
	 */
	 DOUBLEMATRIX *M;
	 DOUBLEVECTOR *v;

	 M=read_map(a,filename,Mref,UVref);
	 v=get_doublevector_from_doublematrix(indices,M,UVref->V->element[2]);
	 free_doublematrix(M);

	 return v;

}

int write_raster_from_doublevector(char *filename, DOUBLEVECTOR *v, T_INIT *UVref, LONGMATRIX *indices, DOUBLEMATRIX *Mref){
	/*
	 * \author Emanuele Cordano
	 *
	 * \date November 2008
	 *

	 * \param filename (char *) name of the file
	 * \param v - (DOUBLECTOR *) - vecto to be mapped
	 * \param UVref - (T_INIT *) T_IIT struct containg ewres and nwres information
	 * \param indices  (LONGVECTOR *) - matrix of indices
     * \param Mref - (DOUBLEMATRIX *) - reference raster

	 *

	 * \brief it writes a map from a vector using get_doublematrix_from_doublevecto and write_map (This function does not take into account the NOVALUE!!!)
	 */

	DOUBLEMATRIX *M;

	M=get_doublematrix_from_doublevector(v,indices,Mref,UVref->V->element[2]);
	write_map(filename, FLOATING_POINT_TYPE,MAP_FORMAT,M,UVref);

	free_doublematrix(M);

	return 0;

}

DOUBLEMATRIX *get_doublematrix_from_mapseries(LONGMATRIX *indices,DOUBLETENSOR *mapseries, T_INIT *UVref){
	/*
	 *\author Emanuele Cordano
	 *\date December 2008
	 *
	 *
	 */
	long l,c;
	DOUBLEMATRIX *map, *mv;
	DOUBLEVECTOR *v;

	for (l=mapseries->ndl;l<=mapseries->ndh;l++){
		map=extract_a_new_map(mapseries,l,UVref);
		v=get_doublevector_from_doublematrix(indices,map,UVref->V->element[2]);
/*		printf ("loop %ld of %ld \n",l,mapseries->ndh);
		printf ("Values: %lf %lf \n",map->element[1][1],v->element[1]);
		stop_execution();*/
		if (l==mapseries->ndl) mv=new_doublematrix(mapseries->ndh,v->nh);
		for (c=mv->ncl;c<=mv->nch;c++){
			mv->element[l][c]=v->element[c];
		}
		free_doublematrix(map);
		free_doublevector(v);
	}

	return mv;

}


