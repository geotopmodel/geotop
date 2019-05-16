#include "meteoio_plugin.h"
#include "turtle.h"
#include "vector.h"
#include "matrix.h"

extern T_INIT *UV;
extern long number_novalue;
// ----------------------------------------------------------------------------------------------------------------
void meteoio_initUV(mio::DEMObject& dem, Matrix<double>* matrix)
{
    /**
     * copy DEM map from MeteoIO to GEOtop
     */
    UV->V.reset(new Vector<double>{2});
    (*UV->V)(1) = -1.0;
    (*UV->V)(2) = number_novalue;  // GEOtop nodata -9999.0

    UV->U.reset(new Vector<double>{4});
    (*UV->U)(1) = dem.cellsize;
    (*UV->U)(2) = dem.cellsize;
    (*UV->U)(3) = dem.llcorner.getNorthing();
    (*UV->U)(4) = dem.llcorner.getEasting();
}
// ----------------------------------------------------------------------------------------------------------------
void copyGridToMatrix(mio::Grid2DObject& gridObject, Matrix<double>* mymatrix)
{
    /**
     * copy map from MeteoIO to GEOtop
     */
    for (std::size_t i=0; i<gridObject.getNy(); i++)
    {
        for (std::size_t j=0; j<gridObject.getNx(); j++)
        {
            if (gridObject.grid2D(j, gridObject.getNy()-1-i) == mio::IOUtils::nodata)
            {
                (*mymatrix)(i+1, j+1) = -9999;
            } else
            {
                (*mymatrix)(i+1, j+1) = gridObject(j, gridObject.getNy()-1-i);
            }
        }
    }
}

// ----------------------------------------------------------------------------------------------------------------