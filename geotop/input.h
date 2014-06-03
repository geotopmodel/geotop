/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.0.0 - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi
 
 This file is part of GEOtop 2.0.0 
 
 GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
    
#ifndef INPUT_H
#define INPUT_H
#include "struct.geotop.h"
//#include "../libraries/geomorphology/geomorphology.0875.h"
//#include "../libraries/geomorphology/geomorphology.h"
#include "pedo.funct.h"
//#include "../libraries/geomorphology/networks.h"
#include "constants.h"
//#include "../libraries/geomorphology/dtm_resolution.h"
#include "../libraries/ascii/rw_maps.h"
#include "../libraries/ascii/tabs.h"
#include "snow.h"
#include "vegetation.h"
#include "output.h"
#include "times.h"
#include "clouds.h"
#include "meteo.h"
#include "meteodata.h"
#include "channels.h"
#include "indices.h"
#include "recovering.h"
#ifdef USE_NETCDF
//#include "../gt_utilities/gt_utilities.h"
#include "../netCDF/read_command_line.h"
#endif
#include <iostream>
#include <string>
#include <vector>
#define __MATHOPTIM_H__
#include <meteoio/MeteoIO.h>

template<typename T> class DBGGeoVector : public GeoVector<T>
{
public:
    DBGGeoVector(const size_t& asize=0) {
        
    } ;
    
    /**
     * A constructor that creates a vector filled with constant values
     * @param asize size of the new array
     * @param init initial value to fill the vector with
     */
    DBGGeoVector(const size_t& asize, const T& init) {
        
    } ;
    
    ~DBGGeoVector() {
        
    }
    

    T& at(const size_t& index) {
        return GeoVector<T>::at(index) ;
    } ;

    const T at(const size_t& index) const {
        return GeoVector<T>::at(index) ;
    } ;
    
    size_t size() const {
        return GeoVector<T>::size() ;
    } ;

    void resize(const size_t& asize) {
        GeoVector<T>::resize(asize) ;
    } ;
    
    void resize(const size_t& asize, const T& init) {
        GeoVector<T>::resize(asize, init) ;
    }
    
    void clear() {
        GeoVector<T>::clear() ;
    } ;
    
    /**
     * Calling this void procedure sets all elements in GeoVector to val
     * (the size is not affected)
     * @param val Value of type T
     */
    void reset(const T& val) {
        GeoVector<T>::reset(val) ;
    };
    
    /**
     * Calling this void procedure sets all elements in GeoVector to val
     * except the ones that already have value val_to_omit (the size is not affected)
     * @param val Value of type T that will be copied to all elements, except if
     * @param val_to_omit Value of type T that will be preserved
     */
    void reset(const T& val, const T& val_to_omit){
        GeoVector<T>::reset(val, val_to_omit) ;
    } ;
    
    void setMetaData(const std::string& name="unknown", const std::string& unit="unknown", const size_t& precision=3){
        GeoVector<T>::setMetaData(name, unit, precision) ;
    };
} ;

typedef struct __INIT_TOOLS__
{
	double swe0;
	double Tsnow0;
	double agesnow0;
    double rhosnow0;
	double rhoglac0;
	double Dglac0;
	double Tglac0;
	std::vector<std::string> met_col_names;
	std::vector<std::string> soil_col_names;
	std::vector<std::string> horizon_col_names;
	std::vector<std::string> point_col_names;
    std::vector<std::string> lapserates_col_names;
	std::vector<std::string> meteostations_col_names;
//TODO: remove LU
	GeoMatrix<double> LU;

	GeoMatrix<double> bed;
	GeoTensor<double> pa_bed ;
	DBGGeoVector<double> init_water_table_depth;

} InitTools;


  void get_all_input(long argc, char *argv[], Topo *top, Soil *sl, Land *land, Meteo *met, Water *wat, Channel *cnet,
		  Par *par, Energy *egy, Snow *snow, Glacier *glac, Times *times, mio::IOManager& iomanager);

  void read_inputmaps(Topo *top, Land *land, Soil *sl, Par *par, FILE *flog, mio::IOManager& iomanager);

  void read_optionsfile_point(Par *par, Topo *top, Land *land, Soil *sl, Times *times, InitTools *IT, FILE *flog);

  void set_bedrock(InitTools *IT, Soil *sl, Channel *cnet, Par *par, Topo *top, GeoMatrix<double>& LC, FILE *flog);

  GeoTensor<double> find_Z_of_any_layer(GeoMatrix<double>& Zsurface, GeoMatrix<double>& slope, GeoMatrix<double>& LC, Soil *sl, short point);

#ifdef WITH_LOGGER
short file_exists(short key);
#else
short file_exists(short key, FILE *flog);
#endif

double peat_thickness(double dist_from_channel);

void initialize_soil_state(SoilState *S, long n, long nl);

void copy_soil_state(SoilState *from, SoilState *to);

void initialize_veg_state(StateVeg *V, long n);

void copy_veg_state(StateVeg *from, StateVeg *to);


#endif
