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
            SNOW_AGE, // Age of the snow
            SNOW_HNcum, // The height of the snow
            SNOW_MELTED, // Snowmelt
            SNOW_SUBL, // The sublimation of the snow
            SNOW_DURATION, // The duration of the snow
            // SNOW_CA, // Snow covered area
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
