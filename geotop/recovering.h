
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
#ifndef RECOVERING_H
#define RECOVERING_H
#include "struct.geotop.h"
//#include "../meteoio_plugin/meteoioplugin.h"
#include <string>
#include "../libraries/ascii/rw_maps.h"

void assign_recovered_map(long n, std::string name, DOUBLEMATRIX *assign, Par *par, DOUBLEMATRIX *Zdistr, DOUBLEMATRIX *Zpoint);

void assign_recovered_map_vector(long n, std::string name, GeoVector<double>& assign, GeoMatrix<long>& rc, Par *par, GeoMatrix<double>& Zdistr, GeoMatrix<double>& Zpoint);

void assign_recovered_map_long(long n, std::string name, GeoMatrix<long>& assign, Par *par, GeoMatrix<double>&  Zdistr, GeoMatrix<double>& Zpoint);


void assign_recovered_tensor(long lbegin, long n, std::string name, GeoTensor<double>& assign, Par *par, GeoMatrix<double>& Zdistr, GeoMatrix<double>& Zpoint);

void assign_recovered_tensor_vector(long lbegin, long n, std::string name, GeoMatrix<double>& assign, GeoMatrix<long>& rc, Par *par, GeoMatrix<double>& Zdistr, GeoMatrix<double>& Zpoint);

void assign_recovered_tensor_channel(long lbegin, long n, std::string name, GeoMatrix<double>& assign,const GeoVector<long> r, const GeoVector<long> c, GeoMatrix<double>& Zdistr);

#endif
