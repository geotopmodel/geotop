/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
BGEOMETRY Version 0.9375 KMackenzie

file bigcells.h

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
#include "rw_maps.h"
#include "geometry.h"
#include "geometry_utilities.h"
#include "geometry_attribute.h"
#include "geometry_freememory.h"
#include "geometry_io.h"
#include "g_raster2plvector.h"
#include "sorting.h"
#include "linear_span.h"
#include "bigcells2.h"
#define BOUNDARY -10
#define EWRES 1
#define NSRES 2
#define X_BOTTOMLEFTCORNER 4
#define Y_BOTTOMLEFTCORNER 3
#define XY_ORIGIN 0.0
#define A_FLAG 2
#define INTEGER_NULL -99
#define INIT_VALUE2 -990
#define MAX_COL 100
#define REFERENCE_INDEX 1
SQUARE_GRID *get_square_grid(DOUBLEMATRIX *DTM, T_INIT *UV,char *file_resume_lines, char *file_resume_polygons, char *file_resume_connections,long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index), int (*check_novalues)(double x, DOUBLEVECTOR *V),long displacement,short print) {
	/*
	 *232: LONGMATRIX *m_indices_from_mask(DOUBLEMATRIX *mask,long row_shift,long col_shift,long novalue,long IBASE,long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index),DOUBLEVECTOR *V, int (*check_novalues)(double x, DOUBLEVECTOR *V)){
	 * \author Emanuele Cordano
	 * \date April 2008
	 *
	 *\param - DTM (DOUBLEMATRIX *) - reference map
	 *\param - UV (T_INIT *) - T_INIT information
	 *\param - (char *) - name of file_resume_lines
	 *\param - (char *) - name of file_resume_polygons
	 *\param - (char *) - name of file_resume_connections;
	 *\param - long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index)
	 *\param - int (*check_novalues)(double x, DOUBLEVECTOR *V)
	 *\param   displacement_index (long), displacement amang which connected polygon are searched
	 *\param - print (short)
	 *
	 *
	 *\return a square grid
	 */

	SQUARE_GRID *sq;
	long count,r,c,vl;

	sq=(SQUARE_GRID *)malloc(sizeof(SQUARE_GRID));
	if (!sq) t_error("Square Grid sq was not allocated");
	sq->indices_pixel=m_indices_from_mask(DTM,0,0,INTEGER_NULL,0,index_pixel_from_a_bin,UV->V,(*check_novalues));
	sq->indices_vertex=m_indices_from_mask(DTM,-1,-1,INTEGER_NULL,0,(*index_pixel_from_a_bin),UV->V,(*check_novalues));
	sq->indices_horizontal_lines=m_indices_from_mask(DTM,-1,0,INTEGER_NULL,0,(*index_pixel_from_a_bin),UV->V,(*check_novalues));
	count=0;
	for(r=sq->indices_horizontal_lines->nrl;r<=sq->indices_horizontal_lines->nrh;r++){
		for (c=sq->indices_horizontal_lines->ncl;c<=sq->indices_horizontal_lines->nch;c++){
	    	if (sq->indices_horizontal_lines->element[r][c]!=INTEGER_NULL) count++;
		}
	}
	sq->nhorizontal_lines=count;
	sq->novalue=INTEGER_NULL;
	sq->indices_vertical_lines=m_indices_from_mask(DTM,0,-1,INTEGER_NULL,count,(*index_pixel_from_a_bin),UV->V,(*check_novalues));

	sq->grid=(GRID *)malloc(sizeof(GRID));
	if (!sq->grid) t_error("Grid of a Square Grid was not allocated");



	 /* build lines and polygons */
	sq->grid->lines=get_linevector_from_raster_grid(sq->indices_horizontal_lines,sq->indices_vertical_lines,sq->indices_vertex,UV->U->element[NSRES],UV->U->element[EWRES],UV->U->element[X_BOTTOMLEFTCORNER],UV->U->element[Y_BOTTOMLEFTCORNER],INTEGER_NULL);
	if (print==1) printf("number of lines: %ld \n",sq->grid->lines->nh);
	sq->grid->polygons=get_polygonvector_from_raster(sq->grid->lines,sq->indices_pixel,sq->indices_horizontal_lines,sq->indices_vertical_lines,UV->U->element[NSRES],UV->U->element[EWRES],UV->U->element[X_BOTTOMLEFTCORNER],UV->U->element[Y_BOTTOMLEFTCORNER],INTEGER_NULL,print);
	if (print==1) printf("number of polygons: %ld \n",sq->grid->polygons->nh);
	// displacement=sq->indices_vertex->nch*2
	sq->grid->links=get_connection_array(sq->grid->polygons,BOUNDARY,displacement,print);
	sq->grid->boundary_indicator=BOUNDARY;
	/* chech the simmetry of the connections */
	vl=connections_symmetry(sq->grid->links,sq->grid->boundary_indicator);
	sq->grid->file_resume_lines=copy_string(file_resume_lines);
	sq->grid->file_resume_polygons=copy_string(file_resume_polygons);
	sq->grid->file_resume_connections=copy_string(file_resume_connections);
	write_grid(sq->grid);

	return sq;

}

void write_grid(GRID *grid){
	/*
	 *\author Emanuele Cordano
	 *\date April 2009
	 *
	 */
	int l,s,r;

	//if (strcmp(filenames->element[O_INFO_LINES]+1,MISSING_FILE))
	l=fprint_linevector(grid->file_resume_lines,grid->lines);
	//if (strcmp(filenames->element[O_INFO_POLYGONS]+1,MISSING_FILE))
	s=fprint_polygonvector(grid->file_resume_polygons,grid->polygons);
	//if (strcmp(filenames->element[O_INFO_CONNECTIONS]+1,MISSING_FILE))
	r=fprint_polygonconnectionattributearray(grid->file_resume_connections,grid->links);

	printf("write_grid function: fprint_linevector %s (exit statues %d)\n",grid->file_resume_lines,l);
	printf("write_grid function: fprint_polygonvector %s (exit statues %d)\n",grid->file_resume_polygons,s);
	printf("write_grid function: fprint_polygonconnectionattributearray %s (exit statues %d)\n",grid->file_resume_connections,r);
	//r=fprint_polygonconnectionattributearray(grid->file_resume_connections,grid->links);

}

RASTER_MAP *new_raster_map(long nh) {
	/*
	 *
	 *\author Emanuele Cordano
	 *\date April 2009
	 *
	 *\long nh - number of layers;
	 *
	 *\creates and allocates a RASTER_MAP class
	 *
	 *
	 */
	RASTER_MAP *raster;

	raster=(RASTER_MAP *)malloc(sizeof(RASTER_MAP));
	if (!raster) t_error("Raster map was not allocated");


	raster->nl=NL;
	raster->nh=nh;
	raster->check_novalues=no_value_function;
	raster->reference_index_map=INTEGER_NULL;
	raster->UV=(T_INIT *)malloc(sizeof(T_INIT));
	if (!raster->UV) t_error("Raster map UV was not allocated");

	raster->layer=(DOUBLEMATRIX **)malloc((size_t)((nh-NL+1+NR_END)*sizeof(DOUBLEMATRIX *)));
	if (!raster->layer) t_error("Raster layer  was not allocated");

	return raster;

}

int add_map_to_raster(long index,char *filename, RASTER_MAP *raster, short a,short print) {
/*
 *
 *
 * \author Emanuele Cordano
 * \date April 2008
 *
 *\param index (param) - index at which he map is added in the raster map class;
 *\param filename (char *) - filename of the map in asccii format
 *\param raster (RASTER_MAP *) - raster map class
 *\param a - (short) parameter to check boundary and novalue (0: nocheck, 1: check boundary, 2: check boundary and novaule (if raster->reference_index_map is set to NULL, a=0 by default)
 *\param print (short) - print parameter

 *
 *\bief add a map from asccii format to a raster_map class
 */

DOUBLEMATRIX *M;
if(print==1) printf("\n raster [%ld] reference_map %ld: \n",index,raster->reference_index_map);

if (raster->reference_index_map==INTEGER_NULL) {
	M=new_doublematrix(2,2);
	raster->layer[index]=read_map(0,filename,M,raster->UV);
	free_doublematrix(M);
	raster->reference_index_map=index;
	if (print==1) printf ("Add_map_to_raster map: map %s was successfully read (index  %ld of %ld) (reference index: %ld) (map %ld,%ld)\n",filename,index,raster->nh,raster->reference_index_map,raster->layer[index]->nrh,raster->layer[index]->nch);
	return 0;
} else if (index==raster->reference_index_map) {
	printf ("Error in add_map_to_raster map  index of %s corresponds to the reference map index %ld    \n",filename,index);
} else {
	raster->layer[index]=read_map(2,filename,raster->layer[raster->reference_index_map],raster->UV);
	if (print==1) printf ("Add_map_to_raster map: map %s was successfully read (index  %ld of %ld) (reference index: %ld) \n",filename,index,raster->nh,raster->reference_index_map);
	return 0;
}
printf ("Error in add_map_to_raster map: map (%s) was not read (index  %ld)   \n",filename,index);

return -1;

}

LONGMATRIX_VECTOR *new_longmatrix_vector(long nh){
/*
 * \author Emanuele Cordano
 * \date April 2009
 *
 * \param nh (long) - number of elements;
 *
 *
 *
 */

	LONGMATRIX_VECTOR *vector;

		vector=(LONGMATRIX_VECTOR *)malloc(sizeof(LONGMATRIX_VECTOR));
		if (!vector) t_error("Longmatrix_vector map was not allocated");


		vector->nl=NL;
		vector->nh=nh;

		vector->element=(LONGMATRIX **)malloc((size_t)((nh-NL+1+NR_END)*sizeof(LONGMATRIX *)));
		if (!vector->element) t_error("Longmatrix vector  was not allocated");
		printf("mew_longmatrix_vector\n\n ");
		return vector;


}

LONGMATRIX_VECTOR *get_fine_indices(LONGMATRIX *fine,LONGMATRIX *coarse,long i_firstcell,long i_lastcell,long novalue,short print) {
	/*
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 *\param (LONGMATRIX*) fine;
	 *\param (LONGMATRIX*) coarse;
	 *\param (long) i_firstcell;
	 *\param (long) I_lastcells;
	 *\param (long) novalue;
	 *\param (short) print;
	 *
	 *\return a LONGMATRIX_VECTOR which has the addresses of the
	 *
	 */


	LONGMATRIX_VECTOR *lco;
	long l,r,c,rs,cs,rp,cp;
	long nrsh,ncsh;
	long i_novalue=novalue-5;
	long ncells=i_lastcell-i_firstcell+1;
	lco=new_longmatrix_vector(ncells);

	nrsh=fine->nrh/coarse->nrh;
	ncsh=fine->nch/coarse->nch;
	if (print==1) printf ("get_fine_indices: generation of a LONGMATRIX_VECTOR (ncells=%ld(%ld,%ld)) (small square dimension %ld %ld) (coarse gird %ld %ld)(fine grid %ld %ld)\n",lco->nh,i_firstcell,i_lastcell,nrsh,ncsh,coarse->nrh,coarse->nch,fine->nrh,fine->nch);
	for(l=lco->nl;l<=lco->nh;l++) {
		lco->element[l]=new_longmatrix(nrsh,ncsh);
		for (r=lco->element[l]->nrl;r<=lco->element[l]->nrh;r++){
			for (c=lco->element[l]->ncl;c<=lco->element[l]->nch;c++){
				lco->element[l]->element[r][c]=i_novalue;
			}
		}
	//	if (print==1) printf("get_fine_indices: Element %ld of %ld was initialized \n",l,lco->nh);
	}

	for(r=fine->nrl;r<=fine->nrh;r++){
		for(c=fine->ncl;c<=fine->nch;c++){
			rs=(r-1)/nrsh+1;
			cs=(c-1)/ncsh+1;
			rp=(r-1)%nrsh+1;
			cp=(c-1)%ncsh+1;
		//	if (print==1) printf("get_fine_indices: r=%ld of %ld c=%ld of %ld  rs=%ld cs=%ld rp=%ld cp=%ld element=%ld novalue=%ld \n",r,fine->nrh,c,fine->nch,rs,cs,rp,cp,coarse->element[rs][cs],novalue);
			l=coarse->element[rs][cs];
			if (l!=novalue) l=l-i_firstcell+1;
			if ((l!=novalue) && ((l<=lco->nh) && (l>=lco->nl))) lco->element[l]->element[rp][cp]=fine->element[r][c];
			if ((l!=novalue) && ((l>lco->nh) || (l<lco->nl))) printf("Error in get_fine_indices: uncorrect value %ld of %ld  for coearse->element[%ld][%ld] r=%ld c=%ld \n",l,lco->nh,rs,cs,rp,cp);


		}
	}

	for(l=lco->nl;l<=lco->nh;l++) {
		for (r=lco->element[l]->nrl;r<=lco->element[l]->nrh;r++){
			for (c=lco->element[l]->ncl;c<=lco->element[l]->nch;c++){
				if (lco->element[l]->element[r][c]==i_novalue) {
					printf("Error in get_fine_elements element %ld (row=%ld col=%ld) was not assigned (number of elements: %ld)",l,r,c,lco->nh);
					return lco;
				}
			}
		}
	}

	if (print==1) printf("Function get_fine_elements  was successfully assigned (number of elements: %ld)",lco->nh);

	return lco;

}

LONGBIN *get_fine_line_indices (LONGMATRIX *h_fine,LONGMATRIX *h_coarse,LONGMATRIX *v_fine,LONGMATRIX *v_coarse,long n_horizontal,long n_lines, long novalue,short print){
/*
 *
 *\author Emanuele Cordano
 *\date April 2009
 *
 *
 */

LONGMATRIX_VECTOR *horizontal, *vertical;

//long n_lines;
LONGVECTOR *index;
LONGBIN *lb;

long i,l;
/// vedre qui!!!
if (print==1) printf("get_fine_line_indices n_horizonal=%ld \n",n_horizontal);
horizontal=get_fine_indices(h_fine,h_coarse,1,n_horizontal,novalue,print);
//n_horizontal=horizontal->nh;
if (print==1) printf("get_fine_line_indices (vertical) n_horizonal=%ld \n",n_horizontal);
vertical=get_fine_indices(v_fine,v_coarse,n_horizontal+1,n_lines,novalue,print);
//n_lines=vertical->nh+n_horizontal;

index=new_longvector(n_lines);
for (l=index->nl;l<=n_horizontal;l++){
	index->element[l]=horizontal->element[l]->nch;
}
for (l=n_horizontal+1;l<=index->nh;l++){
	index->element[l]=vertical->element[l-n_horizontal]->nrh;
}

lb=new_longbin(index);

for(l=lb->index->nl;l<=n_horizontal;l++){
	for (i=1;i<=lb->index->element[l];i++){
		lb->element[l][i]=horizontal->element[l]->element[horizontal->element[l]->nrl][i];
		if (print==1) printf("get_fine_line_indices: horizontal line element [%ld][%ld] = %ld was assigned!\n",l,i,lb->element[l][i]);
	}

}

for(l=n_horizontal+1;l<=lb->index->nh;l++) {
	for (i=1;i<=lb->index->element[l];i++){
		lb->element[l][i]=vertical->element[l-n_horizontal]->element[i][vertical->element[l-n_horizontal]->ncl];
		if (print==1) printf("get_fine_line_indices: vertical line element [%ld][%ld] = %ld was assigned!\n",l,i,lb->element[l][i]);
	}
}



free_longvector(index);

free_longmatrix_vector(vertical);
free_longmatrix_vector(horizontal);

if (print==1) printf("Function get_fine_elements  was successfully assigned (number of elements: %ld)",lb->index->nh);

return lb;




}


DOUBLERASTER_MAP *get_doubleraster_map(long n_coarsemaps, long n_finemaps,char *coarse_mapname, char *fine_mapname,short print){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 *
	 *\param - (long) number of coarse-grid maps
	 *\param - (long) number of fine-grid maps
	 *\param - (char *) - name of the file containing the coarse grid map (reference map)
	 *\param - (char *) - name of the file containing the fine grid map (reference map);
	 *\param (short) print;
	 *
	 * \return a DOUBLERASTER_MAP
	 */

	DOUBLERASTER_MAP *draster;
	char *aa;

	draster=(DOUBLERASTER_MAP *)malloc(sizeof(DOUBLERASTER_MAP));
	if (!draster) t_error("DOUBLERASTER_MAP was not allocated");

	draster->coarse=new_raster_map(n_coarsemaps);
	draster->fine=new_raster_map(n_finemaps);


	add_map_to_raster(draster->coarse->nl,coarse_mapname,draster->coarse,0,print);
	add_map_to_raster(draster->fine->nl,fine_mapname,draster->fine,0,print);

//	if (print==1) printf("\n DOUBLERASTER_MAP was successfully read %ld \n",strlen(coarse_mapname)+1);

//	strcpy(aa,coarse_mapname);
//	draster->mapname_coarse=(char *)malloc((strlen(coarse_mapname)+1)*sizeof(char));
//	memcpy(draster->mapname_coarse,coarse_mapname,strlen(coarse_mapname)+1);
//	strcpy(draster->mapname_coarse,coarse_mapname);
//	draster->mapname_fine=(char *)malloc((strlen(fine_mapname)+1)*sizeof(char));
//	strcpy(draster->mapname_fine,fine_mapname);

	draster->mapname_coarse=copy_string(coarse_mapname);
	draster->mapname_fine=copy_string(fine_mapname);

	if (print==1) printf("\n Doubleraster_map %s (coarse_map) and %s (fine_map) was successfully read \n",draster->mapname_coarse,draster->mapname_fine); // ,coarse_mapname,fine_mapname);
//	if (print==1) printf("\n DOUBLERASTER_MAP was successfully read: %s  \n",aa);

	return draster;

}


DOUBLESQUARE_GRID *get_doublesquare_grid(DOUBLERASTER_MAP *draster, char *resume_filenames,long (*index_pixel_from_a_bin_coarse)(long r, long c,LONGVECTOR *s_index),long (*index_pixel_from_a_bin_fine)(long r, long c,LONGVECTOR *s_index),long d_coarse,long d_fine,short print){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 *\param - DOUBLERASTER_MAP *draster
	 *\param - (char *) root neme of textfiles containing the struct information.
	 *\param long (*index_pixel_from_a_bin_coarse)(long r, long c,LONGVECTOR *s_index) - equation of the filling curve for the pixel of a coarse grid
	 *\param long (*index_pixel_from_a_bin_fine)(long r, long c,LONGVECTOR *s_index),long d_coarse,long d_fine,short print) - equation of the filling curve for the pixel of a fine grid
	 *\param (long) d_coarse - displacement for coarse grid
	 *\param (long) d_fine - displacement for fine grid
	 *\param (short) print
	 *
	 */

	DOUBLESQUARE_GRID *dsq;
	if (print==1) printf("Function get_doublesquare_grid file_resume: %s \n",resume_filenames);
	char *file_resume_lines_coarse=join_strings(resume_filenames,"_lines_coarse.txt");
	char *file_resume_polygons_coarse=join_strings(resume_filenames,"_polygons_coarse.txt");
	char *file_resume_connections_coarse=join_strings(resume_filenames,"_connections_coarse.txt");

	char *file_resume_lines_fine=join_strings(resume_filenames,"_lines_fine.txt");
	char *file_resume_polygons_fine=join_strings(resume_filenames,"_polygons_fine.txt");
	char *file_resume_connections_fine=join_strings(resume_filenames,"_connections_fine.txt");
	char *file_resume_c_polygon=join_strings(resume_filenames,"_c_polygon.txt");
	char *file_resume_c_line=join_strings(resume_filenames,"_c_line.txt");

	dsq=(DOUBLESQUARE_GRID *)malloc(sizeof(DOUBLESQUARE_GRID));
	if (!dsq) t_error("DOUBLESQUARE_GRID was not allocated");


	dsq->big=get_square_grid(draster->coarse->layer[draster->coarse->reference_index_map],draster->coarse->UV,file_resume_lines_coarse,file_resume_polygons_coarse,file_resume_connections_coarse,(*index_pixel_from_a_bin_coarse),draster->coarse->check_novalues,d_coarse,print);
	dsq->fine=get_square_grid(draster->fine->layer[draster->fine->reference_index_map],draster->fine->UV,file_resume_lines_fine,file_resume_polygons_fine,file_resume_connections_fine,(*index_pixel_from_a_bin_fine),draster->fine->check_novalues,d_fine,print);

	dsq->small_content_polygon=get_fine_indices(dsq->fine->indices_pixel,dsq->big->indices_pixel,dsq->big->grid->polygons->nl,dsq->big->grid->polygons->nh,dsq->big->novalue,print);
	dsq->small_content_line=get_fine_line_indices(dsq->fine->indices_horizontal_lines,dsq->big->indices_horizontal_lines,dsq->fine->indices_vertical_lines,dsq->big->indices_vertical_lines,dsq->big->nhorizontal_lines,dsq->big->grid->lines->nh,dsq->big->novalue,print);

//	strcpy(dsq->file_resume_c_polygon,file_resume_c_polygon);
//	strcpy(dsq->file_resume_c_line,file_resume_c_line);

	dsq->file_resume_c_line=copy_string(file_resume_c_line);
	dsq->file_resume_c_polygon=copy_string(file_resume_c_polygon);

	if (print==1) printf("Function get_doublesquare_grid was executed successfully! \n");



	return dsq;
}


/* write functions */

int write_longmatrix_vector(LONGMATRIX_VECTOR *lmv,char *filename) {
	/*
	 *\param  (LONGMATRIX_VECTOR *) lmv
	 *\param  (char *)   filename
	 *
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 * \brief
	 */
	FILE *f;
	long i;

	f=t_fopen(filename,"w");
	fprintf(f,"index{%ld} \n\n",lmv->nh);
	fprintf(f," /** MATRICES OF ''SMALL'' cells   */     \n ");
	for (i=lmv->nl;i<=lmv->nh;i++) {
		fprintf(f,"\n");
		fprintf(f,"%ld: long matrix small_cell {%ld,%ld} \n",i,lmv->element[i]->nrh,lmv->element[i]->nch);
		write_longmatrix_elements(f,lmv->element[i],MAX_COL);

	}


	t_fclose(f);
	return 0;

}

int write_small_lines_indices(LONGBIN *small_lines,char *filename) {
	/*
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 *\param   (LONBIN *) small_lines
	 *\param  (char *)   filename
	 *
	 *\brief write the indices of small lines in a textfles with write_longbin_element function
	 */

	FILE *f;

	f=t_fopen(filename,"w");
	write_longbin_elements(f,small_lines,MAX_COL);
	t_fclose(f);

	return 0;

}

int write_doublesquare_grid(DOUBLESQUARE_GRID *dsq){
	/*
	 * \author Emanuele Cordano
	 * \date April 2009
	 */

	write_grid(dsq->big->grid);
	write_grid(dsq->fine->grid);
	write_longmatrix_vector(dsq->small_content_polygon,dsq->file_resume_c_polygon);
	write_small_lines_indices(dsq->small_content_line,dsq->file_resume_c_line);
	return 0;
}

/* free function */

void free_longmatrix_vector (LONGMATRIX_VECTOR *lmv) {
/*
 *
 * \author Emanuele Cordano
 * \date April 2009
 *
 */
	long i;

	for (i=lmv->nl;i<=lmv->nh;i++){
		free_longmatrix(lmv->element[i]);
	}
	free(lmv->element);
	free(lmv);


}

void free_grid(GRID *grid) {
	/*
	 *\author Emanuele Cordano
	 *\date April 2009
	 *
	 */
	//da fare...
	free_polygon_connection_attribute_array(grid->links);
	free_polygonvector(grid->polygons);
	free_linevector(grid->lines);

	free(grid);

}

void free_square_grid(SQUARE_GRID *sq) {
	/*
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 */

	free_grid(sq->grid);
	free_longmatrix(sq->indices_horizontal_lines);
	free_longmatrix(sq->indices_vertical_lines);
	free_longmatrix(sq->indices_vertex);
	free_longmatrix(sq->indices_pixel);


}

void free_raster_map(RASTER_MAP *raster) {
	/*
	 * \author Emanuele Cordano
	 * \date April 2008
	 *
	 *
	 */
	long i;
	free_doublevector(raster->UV->U);
	free_doublevector(raster->UV->V);
	free(raster->UV);
	for(i=raster->nl;i<=raster->nh;i++){
		if (!(!raster->layer[i]->element)) free_doublematrix(raster->layer[i]);
	}

	free(raster->layer);
	free(raster);



}

void free_doubleraster_map(DOUBLERASTER_MAP *draster){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 *
	 */
	 free_raster_map(draster->coarse);
	 free_raster_map(draster->fine);

	 free(draster);
}

void free_doublesquare_grid(DOUBLESQUARE_GRID *dsq) {
	/*
	 *
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 */

	free_square_grid(dsq->big);
	free_square_grid(dsq->fine);
	free_longbin(dsq->small_content_line);
	free_longmatrix_vector(dsq->small_content_polygon);

	free(dsq);

}

/* proof utilities */

char *strcpy2(char *dest, const char *src)
{
   const char *p;
   char *q;

   for(p = src, q = dest; *p != '\0'; p++, q++)
     *q = *p;

   *q = '\0';

   return dest;
}

char *copy_string(char *origin){
/*
 * \author Emanuele Cordano
 * \date April 2009
 *
 *
 */
char *dest;

dest=(char *)malloc((strlen(origin)+1)*sizeof(char));
strcpy(dest,origin);
return dest;
}

