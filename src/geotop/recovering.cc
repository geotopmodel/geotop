
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.0.0 - 9 Mar 2012
 
 Copyright (c), 2011 - Stefano Endrizzi 
 
 This file is part of GEOtop 2.0.0 
 
 GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "recovering.h"
#include "geotop_common.h"

using namespace std;

void assign_recovered_map_vector   (long n,
                                    std::string name,
                                    GeoVector<double> &assign,
                                    GeoMatrix<long> &rc)
{
    long i, r, c;
    GeoMatrix<double> M;
    std::string temp;

    temp = namefile_i_we2(name, n);

    meteoio_readMap(string(temp), M);

    for (i = 1; i < rc.getRows(); i++)
    {
        r = rc[i][1];
        c = rc[i][2];
        assign[i] = M[r][c];
    }

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map_long (long n,
                                std::string name,
                                GeoMatrix<long> &assign)
{

    long r, c;
    GeoMatrix<double> M;
    std::string temp;

    temp = namefile_i_we2(name, n);

    meteoio_readMap(string(temp), M);
    for (r = 1; r < M.getRows (); r++)
    {
        for (c = 1; c < M.getCols(); c++)
        {
            assign[r][c] = (long) M[r][c];
        }
    }

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void assign_recovered_tensor   (long lbegin,
                                long n,
                                std::string name,
                                GeoTensor<double> &assign)
{

    long r, c, l;
    GeoMatrix < double >M;

    std::string temp1, temp2;

    //TODO: assign->ndl is replaced by 1, need to check!
    for (l = lbegin; l < assign.getDh(); l++) {

          temp1 = namefile_i_we2 (name, n);
          temp2 = namefile_i_we (temp1, l);

          meteoio_readMap (string (temp2), M);

          for (r = 1; r < M.getRows(); r++)
          {
              for (c = 1; c < M.getCols(); c++)
              {
                  assign[l][r][c] = M[r][c];
              }
          }
    }
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void assign_recovered_tensor_vector    (long lbegin,
                                        long n,
                                        std::string name,
                                        GeoMatrix<double> &assign,
                                        GeoMatrix<long> &rc)
{

// lbegin 1 or 0 depending of the size of the data-structure considered (psi is 0)
// this should disappear once we have all data structure correclty allocated (SC.       
    long r, c, i, l;
    GeoMatrix<double> M;
    std::string temp1, temp2;

    for (l = lbegin; l < assign.getRows(); l++)
    {

        temp1 = namefile_i_we2(name, n);
        temp2 = namefile_i_we(temp1, l);

        meteoio_readMap (temp2, M);
        for (i = 1; i < rc.getRows(); i++)
        {
            r = rc[i][1];
            c = rc[i][2];
            assign[l][i] = M[r][c];
        }

    }
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void assign_recovered_tensor_channel   (long lbegin,
                                        long n,
                                        std::string name,
                                        GeoMatrix<double> &assign,
                                        const GeoVector<long> r,
                                        const GeoVector<long> c)
{

    long l;
    size_t ch;
    GeoMatrix<double> M;
    std::string temp1, temp2;

    for (l = lbegin; l < assign.getRows(); l++)
    {

        temp1 = namefile_i_we2(name, n);
        temp2 = namefile_i_we(temp1, l);

        meteoio_readMap(temp2, M);

        for (ch = 1; ch < r.size(); ch++)
        {
            if (r[ch] > 0)
            assign[l][ch] = M[r[ch]][c[ch]];
        }
    }
}

