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
            SOIL_TEMP = 0,       // Soil temperature
            SOIL_WATER_CONTENT,  // Water content in soil [mm]
            SOIL_ICE_CONTENT,    // Ice content in soil depth
            SOIL_WATER_PRESSURE, // Liquid water pressure in soil depth
            SOIL_TOTAL_PRESSURE, // Total water and ice pressure in soil depth
            // SOIL_ET,             // Evapotranspiration from soil
            SOIL_CAN_RAIN,       // Canopy intercepted rain
            SOIL_CAN_SNOW,       // Canopy intercepted snow
            // SNOW MAP VARIABLES
            SNOW_AGE,            // Age of the snow
            // SNOW_DEPTH,          //Depth of the snow
            SNOW_HN,             // The height of the new snow fallend in the time interval
            SNOW_MELTED,         // Snow melting
            SNOW_SUBL,           // The sublimation of the snow
            SNOW_DURATION,       // The duration of the snow
            // SNOW_CA,          // Snow covered area
            // WATER MAPS VARIABLES
            PREC_TOTAL,          // Total precipitation (snow + rain)
            PREC_LIQ,            // Liquid precipitation in the time interval
            PREC_SNOW,           // Snowy part of the total precipitation
            // ENERGY MAPS VARIABLES
            ENER_LWin,           // Incoming longwave radiation
            ENER_SW,             // Shortwave radiation
            ENER_LW,             // Longwave radiation
            ENER_LE,             // Surface latent heat flux
            ENER_H,              // Surface sensible heat flux
            ENER_G,              // Surface heat flux
            ENER_Ts,             // Surface temperature
            ENER_SWin,           // Incoming shortwave radiation
            ENER_SWinb,          // Direct incoming shortwave radiation
            // GLACIER MAPS VARIABLES
            GLAC_MELT,           // Glacier melting
            GLAC_SUBL,           // Sublimation of the glacier
            VECTOR_TEST,
            // METEO MAPS VARIABLES
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
