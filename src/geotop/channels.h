#ifndef _GEOTOP_CHANNELS_H
#define _GEOTOP_CHANNELS_H


/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */


void enumerate_channels(CHANNEL *cnet, Matrix<double> *LC, Matrix<short> *pixel_type,
                        Matrix<double> *Z, Matrix<double> *slope, long novalue);

void next_down_channel_pixel(long r, long c, Matrix<double> *Z, Matrix<double> *LC, Matrix<short> *pixel_type,
                             Matrix<long> *CH, long novalue, long *R, long *C);

void find_max_constraint(Matrix<double> *Z, Matrix<double> *LC, Matrix<short> *pixel_type, Matrix<long> *CH,
                         long novalue, long *R, long *C);

short neighboring_down_channel_pixel(long r, long c, long ir, long ic, Matrix<double> *Z, Matrix<double> *LC,
                                     Matrix<short> *pixel_type, Matrix<long> *CH, long novalue);


#endif
