/*
 * @brief Water Data implementation
 */

#include "water_class.h"

Water::Water(double novalue, size_t nrows, size_t ncols, size_t total_pixel)
{
    // output variables
    PrecTot = GeoMatrix<double>(nrows + 1, ncols + 1, 0.);
    Pnet = GeoMatrix<double>(nrows + 1, ncols + 1, 0.);
    HN = GeoMatrix<double>(nrows + 1, ncols + 1, 0.);
    PrTOT_mean = GeoVector<double>(total_pixel + 1, 0.);
    PrSNW_mean = GeoVector<double>(total_pixel + 1, 0.);
    Pt = GeoVector<double>(total_pixel + 1, novalue);
    Ps = GeoVector<double>(total_pixel + 1, novalue);
    // computational variables
    h_sup = GeoVector<double>(total_pixel + 1, 0.);
    // error = GeoMatrix<double>(nlayers, total_pixel + 1, novalue);
    // Lx = GeoVector<double>(total_pixel + 1, novalue);
    // H0 = GeoVector<double>(total_pixel + 1, novalue);
    // H1 = GeoVector<double>(total_pixel + 1, novalue);
    // dH = GeoVector<double>(total_pixel + 1, novalue);
    // B = GeoVector<double>(total_pixel + 1, novalue);
    // f = GeoVector<double>(total_pixel + 1, novalue);
    // df = GeoVector<double>(total_pixel + 1, novalue);
    // Klat = GeoMatrix<double>(nlayers, total_pixel + 1, 0.);
    // Kbottom = GeoMatrix<double>(nlayers, total_pixel + 1, 0.);
}

//FIXME: Horrible hack needed to cope with legacy code structure
void Water::allocate_data(double novalue, size_t nrows, size_t ncols, size_t total_pixel)
{
    // output variables
    PrecTot.resize(nrows + 1, ncols + 1, 0.);
    Pnet.resize(nrows + 1, ncols + 1, 0.);
    HN.resize(nrows + 1, ncols + 1, 0.);
    PrTOT_mean.resize(total_pixel + 1, 0.);
    PrSNW_mean.resize(total_pixel + 1, 0.);
    Pt.resize(total_pixel + 1, novalue);
    Ps.resize(total_pixel + 1, novalue);
    // computational variables
    h_sup.resize(total_pixel + 1, 0.);
    // error.resize(nlayers, total_pixel + 1, novalue);
    // Lx.resize(total_pixel + 1, novalue);
    // H0.resize(total_pixel + 1, novalue);
    // H1.resize(total_pixel + 1, novalue);
    // dH.resize(total_pixel + 1, novalue);
    // B.resize(total_pixel + 1, novalue);
    // f.resize(total_pixel + 1, novalue);
    // df.resize(total_pixel + 1, novalue);



	
    // Klat.resize(nlayers, total_pixel + 1, 0.);
    // Kbottom.resize(nlayers, total_pixel + 1, 0.);

    Voutlandsub = 0.;
    Voutlandsup = 0.;
    Voutbottom = 0.;
}
