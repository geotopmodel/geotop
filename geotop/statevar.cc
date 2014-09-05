/**
 * @brief Snow and Glacier State Variables implementation
 * @date September 2014
 */

#include "statevar.h"

Statevar3D::Statevar3D(double novalue, size_t layers, size_t rows, size_t columns)
{
	type = GeoMatrix<short>(rows+1, columns+1, 2);
	lnum = GeoMatrix<long>(rows+1, columns+1,0);

    //FIXME: GeoTensor's constructor assumes that the order of indices is:
    //Number of Rows, Number of Columns and Number of Layers BUT the rest of the code
    //assumes that the order is layers, rows and then columns.
	Dzl = GeoTensor<double>(layers+1, rows+1, columns+1 , 0.);
	w_liq = GeoTensor<double>(layers+1, rows+1, columns+1, 0.);
	w_ice = GeoTensor<double>(layers+1, rows+1, columns+1, 0.);
	T = GeoTensor<double>(layers+1, rows+1, columns+1, novalue);
}
