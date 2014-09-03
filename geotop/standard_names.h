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
            UNKNOWN_VAR
        };


			typedef struct {

					std::string longstring;
					std::string shortstring;
					std::string unitstring;
					
			} stringpair_s;


			std::vector<stringpair_s> StandardNames = std::vector<stringpair_s>((size_t)UNKNOWN_VAR);

			static void fillStandardName(size_t index, const std::string& longstring, const std::string& shortstring, const std::string& unitstring)
			{
					StandardNames[index].longstring = longstring;
					StandardNames[index].shortstring = shortstring;
					StandardNames[index].unitstring = unitstring;
			}

			void fillStandardNames()
			{
					static bool init = false;

					//Avoid multiple initializations
					if (not init)
					{
							fillStandardName(SOIL_TEMP, "SoilTemperature", "SoilT", "C");
							fillStandardName(SOIL_WATER_CONTENT, "SoilWaterContent", "Theta", "mm");
							init = true;
					}
			}

	}
}

#endif
