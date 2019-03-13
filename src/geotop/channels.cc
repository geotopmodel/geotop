
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

#include "struct.geotop.h"
#include "channels.h"
#include "constants.h"

extern T_INIT *UV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void enumerate_channels(CHANNEL *cnet, Matrix<double> *LC,
                        Matrix<short> *pixel_type, Matrix<double> *Z, Matrix<double> *slope, long novalue)
{
    long r, c, rnext, cnext, i=0;

    do
    {
        find_max_constraint(Z, LC, pixel_type, cnet->ch.get(), novalue, &r, &c);
        if (r>0)
        {
            i++;
            (*cnet->r)(i)=r;
            (*cnet->c)(i)=c;
            (*cnet->ch)(r,c)=i;

            do
            {
                next_down_channel_pixel(r, c, Z, LC, pixel_type, cnet->ch.get(), novalue, &rnext, &cnext);
                if (rnext>0)
                {
                    i++;
                    if (fabs(rnext-r)==1 && fabs(cnext-c)==1)
                    {
                        (*cnet->length)(i-1)+=0.5*sqrt(2.);
                        (*cnet->length)(i)+=0.5*sqrt(2.);
                    }
                    else
                    {
                        (*cnet->length)(i-1)+=0.5;
                        (*cnet->length)(i)+=0.5;
                    }

                    (*cnet->r)(i)=rnext;
                    (*cnet->c)(i)=cnext;
                    (*cnet->ch)(rnext,cnext)=i;
                    r=rnext;
                    c=cnext;
                }
            }
            while (rnext>0);

            if (rnext<0)
            {
                if (fabs(-rnext-r)==1 && fabs(-cnext-c)==1)
                {
                    (*cnet->length)(i)+=0.5*sqrt(2.);
                }
                else
                {
                    (*cnet->length)(i)+=0.5;
                }
            }
            else
            {
                (*cnet->length)(i)+=0.5;
            }
        }
    }
    while (r>0);

    for (i=1; i<=cnet->r->nh; i++)
    {
        r = (*cnet->r)(i);
        c = (*cnet->c)(i);
        (*cnet->length)(i) *= ((*UV->U)(1)/cos((*slope)(r,c)*GTConst::Pi/180.));
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void next_down_channel_pixel(long r, long c, Matrix<double> *Z, Matrix<double> *LC,
                             Matrix<short> *pixel_type, Matrix<long> *CH, long novalue, long *R, long *C)
{
    *R = 0;
    *C = 0;

    /** (1a) (ir, ic) = (1, 1) */
    if (neighboring_down_channel_pixel(r, c, 1, 1, Z, LC, pixel_type, CH, novalue) == -1)
    {
        *R=r+1;
        *C=c+1;
        *R=(*R)*(-1.);
        *C=(*C)*(-1.);
    }

    /** (2a) (ir, ic) = (1, -1) */
    if (neighboring_down_channel_pixel(r, c, 1, -1, Z, LC, pixel_type, CH, novalue) == -1)
    {
        *R=r+1;
        *C=c-1;
        *R=(*R)*(-1.);
        *C=(*C)*(-1.);
    }

    /** (3a) (ir, ic) = (-1, 1) */
    if (neighboring_down_channel_pixel(r, c, -1, 1, Z, LC, pixel_type, CH, novalue) == -1)
    {
        *R=r-1;
        *C=c+1;
        *R=(*R)*(-1.);
        *C=(*C)*(-1.);
    }

    /** (4a) (ir, ic) = (-1, -1) */
    if (neighboring_down_channel_pixel(r, c, -1, -1, Z, LC, pixel_type, CH, novalue) == -1)
    {
        *R=r-1;
        *C=c-1;
        *R=(*R)*(-1.);
        *C=(*C)*(-1.);
    }

    /** (5a) (ir, ic) = (1, 0) */
    if (neighboring_down_channel_pixel(r, c, 1, 0, Z, LC, pixel_type, CH, novalue) == -1)
    {
        *R=r+1;
        *C=c;
        *R=(*R)*(-1.);
        *C=(*C)*(-1.);
    }

    /** (6a) (ir, ic) = (0, 1) */
    if (neighboring_down_channel_pixel(r, c, 0, 1, Z, LC, pixel_type, CH, novalue) == -1)
    {
        *R=r;
        *C=c+1;
        *R=(*R)*(-1.);
        *C=(*C)*(-1.);
    }

    /** (7a) (ir, ic) = (-1, 0) */
    if (neighboring_down_channel_pixel(r, c, -1, 0, Z, LC, pixel_type, CH, novalue) == -1)
    {
        *R=r-1;
        *C=c;
        *R=(*R)*(-1.);
        *C=(*C)*(-1.);
    }

    /** (8a) (ir, ic) = (0, -1) */
    if (neighboring_down_channel_pixel(r, c, 0, -1, Z, LC, pixel_type, CH, novalue) == -1)
    {
        *R=r;
        *C=c-1;
        *R=(*R)*(-1.);
        *C=(*C)*(-1.);
    }

    /** (5b) (ir, ic) = (1, 0) */
    if (neighboring_down_channel_pixel(r, c, 1, 0, Z, LC, pixel_type, CH, novalue) == 1)
    {
        *R=r+1;
        *C=c;
    }

    /** (6b) (ir, ic) = (0, 1) */
    if (neighboring_down_channel_pixel(r, c, 0, 1, Z, LC, pixel_type, CH, novalue) == 1)
    {
        *R=r;
        *C=c+1;
    }

    /** (7b) (ir, ic) = (-1, 0) */
    if (neighboring_down_channel_pixel(r, c, -1, 0, Z, LC, pixel_type, CH, novalue) == 1)
    {
        *R=r-1;
        *C=c;
    }

    /** (8b) (ir, ic) = (0, -1) */
    if (neighboring_down_channel_pixel(r, c, 0, -1, Z, LC, pixel_type, CH, novalue) == 1)
    {
        *R=r;
        *C=c-1;
    }

    /** (1b) (ir, ic) = (1, 1) */
    if (neighboring_down_channel_pixel(r, c, 1, 1, Z, LC, pixel_type, CH, novalue) == 1)
    {
        *R=r+1;
        *C=c+1;
    }

    /** (2b) (ir, ic) = (1, -1) */
    if (neighboring_down_channel_pixel(r, c, 1, -1, Z, LC, pixel_type, CH, novalue) == 1)
    {
        *R=r+1;
        *C=c-1;
    }

    /** (3b) (ir, ic) = (-1, 1) */
    if (neighboring_down_channel_pixel(r, c, -1, 1, Z, LC, pixel_type, CH, novalue) == 1)
    {
        *R=r-1;
        *C=c+1;
    }

    /** (4b) (ir, ic) = (-1, -1) */
    if (neighboring_down_channel_pixel(r, c, -1, -1, Z, LC, pixel_type, CH, novalue) == 1)
    {
        *R=r-1;
        *C=c-1;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void find_max_constraint(Matrix<double> *Z, Matrix<double> *LC,
                         Matrix<short> *pixel_type, Matrix<long> *CH, long novalue, long *R, long *C)
{
/** find highest channel pixel that has not been enumerated yet */
    long r, c;
    double z = -1.E99;

    *R = 0;
    *C = 0;

    for (r=1; r<=CH->nrh; r++)
    {
        for (c=1; c<=CH->nch; c++)
        {
            if ( (long)(*LC)(r,c) != novalue /** land cover has a value */
                 && (*pixel_type)(r,c) >= 10 /** channel pixel */
                 && (*CH)(r,c) == 0 /** channel pixel not enumerated */
                 && (*Z)(r,c) > z )
            {
                z = (*Z)(r,c);
                *R = r;
                *C = c;
            }
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short neighboring_down_channel_pixel(long r, long c, long ir, long ic, Matrix<double> *Z, Matrix<double> *LC,
                                     Matrix<short> *pixel_type, Matrix<long> *CH, long novalue)
{
    short yes=0;
    long R=r+ir, C=c+ic;

    if (R>=1 && R<=CH->nrh && C>=1 && C<=CH->nch)
    {
        if ((long)(*LC)(R,C)!=novalue)
        {
            /** neighboring pixel is a channel */
            if ((*Z)(R,C)<=(*Z)(r,c) && (*pixel_type)(R,C)>=10)
                yes = -1;
            /** neighboring pixel is a channels that has not been enumerated yet */
            if ((*Z)(R,C)<=(*Z)(r,c) && (*pixel_type)(R,C)>=10 && (*CH)(R,C)==0)
                yes = 1;
        }
    }

    return (yes);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
