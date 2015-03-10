/*
 * @brief Energy Data implementation
 */

#include "meteo_class.h"

Meteo::Meteo(double novalue, size_t Z0nrows, size_t Z0ncols)
{
   tau_cl_map = GeoMatrix<double>(Z0nrows, Z0ncols, novalue);
   tau_cl_av_map = GeoMatrix<double>(Z0nrows, Z0ncols, novalue);
   tau_cl_map_yes = GeoMatrix<short>(Z0nrows, Z0ncols, (short)novalue);
   tau_cl_av_map_yes = GeoMatrix<short>(Z0nrows, Z0ncols, (short)novalue);
}

Meteo::Meteo(size_t nrows, size_t ncols, double Vmin, double Pa)
{   
   Tgrid = GeoMatrix<double>(nrows + 1, ncols + 1, 5.);
   Pgrid = GeoMatrix<double>(nrows + 1, ncols + 1, Pa);
   Vgrid = GeoMatrix<double>(nrows + 1, ncols + 1, Vmin);
   Vdir = GeoMatrix<double>(nrows + 1, ncols + 1, 0.);
   RHgrid = GeoMatrix<double>(nrows + 1, ncols + 1, 0.7);
   ILWRgrid = GeoMatrix<double>(nrows + 1, ncols + 1, 0);
   // GeoVector<long> imeteo_stations; // SERAFIN: I don't know what to do, because this variable is allocated in parameters.cc
}

//FIXME: Horrible hack needed to cope with legacy code structure
void Meteo::allocate_data(double novalue, size_t nrows, size_t ncols, size_t Z0nrows, size_t Z0ncols, double Vmin, double Pa)
{
   tau_cl_map.resize(Z0nrows, Z0ncols, novalue);
   tau_cl_av_map.resize(Z0nrows, Z0ncols, novalue);
   tau_cl_map_yes.resize(Z0nrows, Z0ncols, (short)novalue);
   tau_cl_av_map_yes.resize(Z0nrows, Z0ncols, (short)novalue);
   
   Tgrid.resize(nrows + 1, ncols + 1, 5.);
   Pgrid.resize(nrows + 1, ncols + 1, Pa);
   Vgrid.resize(nrows + 1, ncols + 1, Vmin);
   Vdir.resize(nrows + 1, ncols +1, 0.);
   RHgrid.resize(nrows + 1, ncols + 1, 0.7);
   ILWRgrid.resize(nrows + 1, ncols + 1, novalue);
}

