#include "standard_names.h"


namespace geotop
{
    namespace input
    {

        static std::vector<stringpair_s> StandardNames = std::vector<stringpair_s>((size_t)UNKNOWN_VAR);

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
                fillStandardName(SOIL_WATER_CONTENT, "SoilWaterContent", "SoilTh", "mm");
                fillStandardName(SOIL_ICE_CONTENT, "SoilIceContent", "SoilThI", "mm");
                fillStandardName(SOIL_WATER_PRESSURE, "SoilWaterPressure", "SoilP", "m");
                fillStandardName(SOIL_TOTAL_PRESSURE, "SoilTotalPressure", "SoilTotP", "m");
                init = true;
            }
        }

        std::string getLongString(Variable v)
        {
            return StandardNames[v].longstring;
        }
        
        std::string getShortString(Variable v)
        {
            return StandardNames[v].shortstring;
        }
        
        std::string getUnitString(Variable v)
        {
            return StandardNames[v].unitstring;
        }
    }
}

