#ifndef STANDARD_NAMES_H
#define STANDARD_NAMES_H

#include <vector>
#include <string>

namespace geotop
{
    namespace input
    {
        /**
         * @brief Output variables known
         */
        enum Variable {
            SOIL_TEMP = 0, //Soil temperature
            SOIL_WATER_CONTENT, //Water content in soil [mm]
            SOIL_ICE_CONTENT, //Ice content in soil depth
            SOIL_WATER_PRESSURE, //Liquid water pressure in soil depth
            SOIL_TOTAL_PRESSURE, //Total water and ice pressure in soil depth
            // SNOW MAP VARIABLES
            SNOW_AGE, //Age of the snow
            SNOW_DEPTH, //Depth of the snow
            SNOW_HN, //The height of the new snow fallend in the time interval
            SNOW_MELTED, //Snowmelt
            SNOW_SUBL, //The sublimation of the snow
            SNOW_DURATION, //The duration of the snow
            // SNOW_CA, //Snow covered area
            // WATER MAPS VARIABLES
            PREC_TOTAL, //Total precipitation (snow + rain)
            PREC_LIQ, //Liquid precipitation which reaches the soil surface in the time interval
            PREC_SNOW, //Snowy part of the total precipitation
            VECTOR_TEST,
            UNKNOWN_VAR
        };


        typedef struct {

            std::string longstring;
            std::string shortstring;
            std::string unitstring;

        } stringpair_s;


        void fillStandardNames();

        std::string getLongString(Variable v);
        std::string getShortString(Variable v);
        std::string getUnitString(Variable v);

    }
}

#endif
