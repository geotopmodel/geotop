
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
BGEOMETRY Version 0.9375 KMackenzie

file geometry_attribute.h

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




/* function header */

polygon_connection_attributes *get_connection(POLYGON *polygon,POLYGONVECTOR *polygons, long boundary,long displacement,  short print);

polygon_connection_attribute_array *new_connection_attributes(long nh);


polygon_connection_attribute_array *get_connection_array(POLYGONVECTOR *polygons, long boundary,long displacement,short print);

int fprint_polygonconnectionattributearray(char *filename,polygon_connection_attribute_array *pca);

int connections_symmetry(polygon_connection_attribute_array* pca, long boundary);
