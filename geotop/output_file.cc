/**
 * @file output_file.cc
 * @Author  Gianfranco Gallizia (skyglobe83@gmail.com)
 * @copyright (C) 2014 eXact lab srl
 * @date June, 2014
 * @brief Output file definition class
 * 
 * Defines a container for the output definition: output variable, dimension (1D, 2D, 3D),
 * type and period of integration.
 */

#include "output_file.h"
#include <vector>
#include <stdio.h>
#include <string.h>
#include "times.h"
#include "global_logger.h"

namespace geotop
{
    namespace input
    {

        /*=====================================================================
         * Constants
         =====================================================================*/
        const double minPeriod = 60;

        /*=====================================================================
         * Static functions
         =====================================================================*/

        /**
         * @internal
         * @brief Splits an extended key
         * @param[in] key the key to split
         * @return a std::vector with the subkeys
         */
        static std::vector<std::string> split_ext_key(std::string key)
        {
            std::vector<std::string> output;
            std::string tmp;

            size_t i = 0;

            while (i < key.size())
            {
                if (key.at(i) != ':')
                    tmp.push_back(key.at(i));
                else
                {
                    output.push_back(tmp);
                    tmp="";
                    i++; //Skip the separator character
                }

                i++;
            }

            output.push_back(tmp);

            return output;
        }

        /**
         * @internal
         * @brief Converts a double to a long with rounding
         * @param[in] d the double to convert
         * @return truncate (d + 0.5)
         */
        static long rounding(double d)
        {
            d += 0.5;
            return (long)d ;
        }

        /**
         * @internal
         * @brief Converts a string to lowercase
         * @param[in] s the string to convert
         * @return the lowercase version of s
         */
        static std::string toLower(std::string s)
        {
            std::string tmp;
            size_t i;

            for (i = 0; i < s.size(); i++)
            {
                char c = s.at(i);
                if (c >= 'A' && c <= 'Z')
                    c += ('a' - 'A');

                tmp.push_back(c);
            }

            return tmp;
        }

        /*=====================================================================
         * TemporaryValues class members
         =====================================================================*/

        /**
         * @brief Default constructor: sets eveything as invalid
         */
        TemporaryValues::TemporaryValues()
        {
            mWhatIsValid = -1;
            mDValue = geotop::input::gDoubleNoValue;
            mVValue = NULL;
            mMValue = NULL;
            mTValue = NULL;
        }

        TemporaryValues::TemporaryValues(double init)
        {
            mWhatIsValid = 0;
            mDValue = init;
            mVValue = NULL;
            mMValue = NULL;
            mTValue = NULL;
        }

        TemporaryValues::TemporaryValues(GeoVector<double>* init)
        {
            mWhatIsValid = 1;
            mDValue = geotop::input::gDoubleNoValue;
            mVValue = init;
            mMValue = NULL;
            mTValue = NULL;
        }

        TemporaryValues::TemporaryValues(GeoMatrix<double>* init)
        {
            mWhatIsValid = 2;
            mDValue = geotop::input::gDoubleNoValue;
            mVValue = NULL;
            mMValue = init;
            mTValue = NULL;
        }

        TemporaryValues::TemporaryValues(GeoTensor<double>* init)
        {
            mWhatIsValid = 3;
            mDValue = geotop::input::gDoubleNoValue;
            mVValue = NULL;
            mMValue = NULL;
            mTValue = init;
        }

        double TemporaryValues::getValueD()
        {
            assert(mWhatIsValid == 0);
            return mDValue;
        }

        GeoVector<double>* TemporaryValues::getValuesV()
        {
            assert(mWhatIsValid == 1);
            return mVValue;
        }

        GeoMatrix<double>* TemporaryValues::getValuesM()
        {
            assert(mWhatIsValid == 2);
            return mMValue;
        }

        GeoTensor<double>* TemporaryValues::getValuesT()
        {
            assert(mWhatIsValid == 3);
            return mTValue;
        }

        /*=====================================================================
         * OutputFile class members
         =====================================================================*/

        OutputFile::OutputFile(std::string extended_key, double period, long layer, std::string prefix)
        {
            mVariable = geotop::input::UNKNOWN_VAR;
            mDimension = geotop::input::UNKNOWN_DIM;
            mType = geotop::input::UNKNOWN_INTEG;
            mLayerIndex = layer;
            mPrefix = prefix;

            size_t l = mPrefix.length();

			//Construct known variables names array
			fillStandardNames();

            if (l > 0)
            {
                if (mPrefix.at(l - 1) != '/')
                    mPrefix.push_back('/');
            }

            if (period >= minPeriod)
            {
               
                mPeriod = rounding(period);

                std::vector<std::string> values = split_ext_key(extended_key);

                if (values.size() == 3)
                {
                    //Variable
                    std::string tmp = values.at(0);
                    mVariable = str2var(tmp);

                    if (mVariable == geotop::input::UNKNOWN_VAR)
                    {
                        geotop::logger::GlobalLogger* lg =
                            geotop::logger::GlobalLogger::getInstance();

                        lg->logsf(geotop::logger::CRITICAL,
                                  "Unknown output variable: '%s'. Modify the geotop.inpts. Aborting.",
                                  tmp.c_str());

						exit(1);

                    }

                    //Dimension
                    tmp = values.at(1);
                    if (tmp.compare("1Dp") == 0) mDimension = D1Dp;
                    if (tmp.compare("1Ds") == 0) mDimension = D1Ds;
                    if (tmp.compare("2D") == 0) mDimension = D2D;
                    if (tmp.compare("3D") == 0) mDimension = D3D;

                    //Integration Type
                    tmp = values.at(2);
                    if (tmp.compare("AVG") == 0) mType = AVG;
                    if (tmp.compare("CUM") == 0) mType = CUM;
                    if (tmp.compare("INS") == 0) mType = INS;
                    
                    if (isValidDimension() == false)
                    {
                        mVariable = geotop::input::UNKNOWN_VAR;
                        mDimension = geotop::input::UNKNOWN_DIM;
                        mType = geotop::input::UNKNOWN_INTEG;

                        geotop::logger::GlobalLogger* lg =
                            geotop::logger::GlobalLogger::getInstance();

                        lg->logsf(geotop::logger::ERROR,
                                  "Invalid output file definition: %s",
                                  extended_key.c_str());

                    }

                }

            }
            else
            {
                geotop::logger::GlobalLogger* lg =
                    geotop::logger::GlobalLogger::getInstance();

                lg->logsf(geotop::logger::ERROR,
                          "Invalid integration period for key '%s': %f",
                          extended_key.c_str(), period);
            }

        }

        OutputFile::~OutputFile()
        {
        }

        std::string OutputFile::getFileName(double dateeur12, long layer)
        {
            char buffer[13] = {'\0'};
            std::string output;
            long day = 0L, month = 0L, year = 0L, hour = 0L, min = 0L;

            if (mDimension == D1Dp || mDimension == D1Ds) 
            {
                output = std::string(""); //Do not prepend the date for tab files
            }
            else
            {
                convert_dateeur12_daymonthyearhourmin(dateeur12, &day, &month, &year, &hour, &min);
                sprintf(buffer, "%.4ld%.2ld%.2ld%.2ld%.2ld", year, month, day, hour, min);

                output = std::string(buffer);
                output.append("_");
            }

            output.append(var2str(mVariable));

            output.append("_");

            switch (mDimension)
            {
                case D1Dp:
                    output.append("1Dp");
                    break;
                case D1Ds:
                    output.append("1Ds");
                    break;
                case D2D:
                    output.append("2D");
                    break;
                case D3D:
                    output.append("3D");
                    break;
                default:
                    output.append("UNKNOWN");
                    break;
            }

            output.append("_");

            switch (mType)
            {
                case AVG:
                    output.append("AVG");
                    break;
                case CUM:
                    output.append("CUM");
                    break;
                case INS:
                    output.append("INS");
                    break;
                default:
                    output.append("UNKNOWN");
                    break;
            }

            output.append("_");

            //Period
            memset(buffer, 0, 13);

            if (mPeriod > 0 && mPeriod < 60)
                sprintf(buffer,"%.2lds", mPeriod);

            if (mPeriod >= 60 && mPeriod < 3600)
                sprintf(buffer,"%.2ldm", mPeriod / 60);

            if (mPeriod >= 3600 && mPeriod < 86400)
                sprintf(buffer,"%.2ldh", mPeriod / 3600);

            if (mPeriod >= 86400)
                sprintf(buffer,"%ldd", mPeriod / 86400);

            output.append(buffer);

            //Layer
            if (layer != -1 && mDimension != D1Dp && mDimension != D1Ds)
            {
                output.append("_L");
                memset(buffer, 0, 13);
                sprintf(buffer, "%.4ld", layer);
                output.append(buffer);
            }

            //Extension
            output.append(".asc");
            
            return output;
        }

        std::string OutputFile::getFilePath(double dateeur12, long layer)
        {
            std::string output = std::string(mPrefix);
            output.append(getFileName(dateeur12, layer));
            return output;
        }

       
        std::string OutputFile::toString()
        {
            char buffer[256] = {'\0'};

            sprintf(buffer,"%ld", mPeriod);

            std::string output(var2str(mVariable));

            output.append("::");

            switch (mDimension)
            {
                case D1Dp:
                    output.append("1Dp");
                    break;
                case D1Ds:
                    output.append("1Ds");
                    break;
                case D2D:
                    output.append("2D");
                    break;
                case D3D:
                    output.append("3D");
                    break;
                default:
                    output.append("UNKNOWN");
                    break;
            }

            output.append("::");

            switch (mType)
            {
                case AVG:
                    output.append("AVG");
                    break;
                case CUM:
                    output.append("CUM");
                    break;
                case INS:
                    output.append("INS");
                    break;
                default:
                    output.append("UNKNOWN");
                    break;
            }

            output.append("=");

            output.append(buffer);

            return output;
        }

        std::string OutputFile::var2str(Variable v)
        {
            std::string output = "";

            switch (v)
            {
                case SOIL_TEMP:
                case SOIL_WATER_CONTENT:
                case SOIL_ICE_CONTENT:
                case SOIL_WATER_PRESSURE:
                case SOIL_TOTAL_PRESSURE:
                // case SOIL_ET:
                case SOIL_CAN_RAIN:
                case SOIL_CAN_SNOW:
                // SNOW MAP VARIABLES
                case SNOW_AGE:
                // case SNOW_DEPTH:
                case SNOW_HN:
                case SNOW_MELTED:
                case SNOW_SUBL:
                case SNOW_DURATION:
                // case SNOW_CA:
                // WATER MAP VARIABLES
                case PREC_TOTAL:
                case PREC_LIQ:
                case PREC_SNOW:
                // ENERGY MAP VARIABLES
                case ENER_LWin:
                case ENER_SW:
                case ENER_LW:
                case ENER_LE:
                case ENER_H:
                case ENER_G:
                case ENER_Ts:
                case ENER_SWin:
                case ENER_SWinb:
                case GLAC_MELT:
                case GLAC_SUBL:
                case VECTOR_TEST:
                    output.append(getLongString(v));
                    break;
                default:
                    output.append("UNKNOWN");
                    break;
            }

            return output;
        }

        Variable OutputFile::str2var(std::string v)
        {
            std::string tmp = toLower(v);

            if (tmp.compare(toLower(getLongString(SOIL_TEMP))) == 0) return SOIL_TEMP;
            if (tmp.compare(toLower(getLongString(SOIL_WATER_CONTENT))) == 0) return SOIL_WATER_CONTENT;
            if (tmp.compare(toLower(getLongString(SOIL_ICE_CONTENT))) == 0) return SOIL_ICE_CONTENT;
            if (tmp.compare(toLower(getLongString(SOIL_WATER_PRESSURE))) == 0) return SOIL_WATER_PRESSURE;
            if (tmp.compare(toLower(getLongString(SOIL_TOTAL_PRESSURE))) == 0) return SOIL_TOTAL_PRESSURE;
            // if (tmp.compare(toLower(getLongString(SOIL_ET))) == 0) return SOIL_ET;
            if (tmp.compare(toLower(getLongString(SOIL_CAN_RAIN))) == 0) return SOIL_CAN_RAIN;
            if (tmp.compare(toLower(getLongString(SOIL_CAN_SNOW))) == 0) return SOIL_CAN_SNOW;
            // SNOW MAP VARIABLES
            if (tmp.compare(toLower(getLongString(SNOW_AGE))) == 0) return SNOW_AGE;
            // if (tmp.compare(toLower(getLongString(SNOW_DEPTH))) == 0) return SNOW_DEPTH;
            if (tmp.compare(toLower(getLongString(SNOW_HN))) == 0) return SNOW_HN;
            if (tmp.compare(toLower(getLongString(SNOW_MELTED))) == 0) return SNOW_MELTED;
            if (tmp.compare(toLower(getLongString(SNOW_SUBL))) == 0) return SNOW_SUBL;
            if (tmp.compare(toLower(getLongString(SNOW_DURATION))) == 0) return SNOW_DURATION;
            // if (tmp.compare(toLower(getLongString(SNOW_CA))) == 0) return SNOW_CA;
            // WATER MAP VARIABLES
            if (tmp.compare(toLower(getLongString(PREC_TOTAL))) == 0) return PREC_TOTAL;
            if (tmp.compare(toLower(getLongString(PREC_LIQ))) == 0) return PREC_LIQ;
            if (tmp.compare(toLower(getLongString(PREC_SNOW))) == 0) return PREC_SNOW;
            if (tmp.compare(toLower(getLongString(ENER_LWin))) == 0) return ENER_LWin;
            if (tmp.compare(toLower(getLongString(ENER_SW))) == 0) return ENER_SW;
            if (tmp.compare(toLower(getLongString(ENER_LW))) == 0) return ENER_LW;
            if (tmp.compare(toLower(getLongString(ENER_LE))) == 0) return ENER_LE;
            if (tmp.compare(toLower(getLongString(ENER_H))) == 0) return ENER_H;
            if (tmp.compare(toLower(getLongString(ENER_G))) == 0) return ENER_G;
            if (tmp.compare(toLower(getLongString(ENER_Ts))) == 0) return ENER_Ts;
            if (tmp.compare(toLower(getLongString(ENER_SWin))) == 0) return ENER_SWin;
            if (tmp.compare(toLower(getLongString(ENER_SWinb))) == 0) return ENER_SWinb;
            if (tmp.compare(toLower(getLongString(GLAC_MELT))) == 0) return GLAC_MELT;
            if (tmp.compare(toLower(getLongString(GLAC_SUBL))) == 0) return GLAC_SUBL;
            if (tmp.compare(toLower(getLongString(VECTOR_TEST))) == 0) return VECTOR_TEST;

			
            return UNKNOWN_VAR;
        }

        bool OutputFile::isValidDimension()
        {
            bool output = false;

            if (mDimension == UNKNOWN_DIM)
                return false;

            switch (mVariable)
            {
                case SOIL_TEMP:
                    output = true;
                    break;
                case SOIL_WATER_CONTENT:
                    output = true;
                    break;
                case SOIL_ICE_CONTENT:
                    output = true;
                    break;
                case SOIL_WATER_PRESSURE:
                    output = true;
                    break;
                case SOIL_TOTAL_PRESSURE:
                    output = true;
                    break;
                // case SOIL_ET:
                //     output = true;
                //     break;
                case SOIL_CAN_RAIN:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case SOIL_CAN_SNOW:
                    output = (mDimension == D3D) ? false : true;
                    break;
                // SNOW MAP VARIABLES
                //If mDimension is D3D then false else true
                case SNOW_AGE:
                    output = (mDimension == D3D) ? false : true;
                    break;
                // case SNOW_DEPTH:// Dzl is a GeoTensor that has to be cumulate in space
                //     output = true;
                //     break;
                case SNOW_HN:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case SNOW_MELTED:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case SNOW_SUBL:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case SNOW_DURATION:
                    output = (mDimension == D3D) ? false : true;
                    break;
                // case SNOW_CA:
                //     output = true;
                //     break;
                case PREC_TOTAL:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case VECTOR_TEST:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case PREC_LIQ:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case PREC_SNOW:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_LWin:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_SW:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_LW:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_LE:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_H:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_G:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_Ts:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_SWin:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case ENER_SWinb:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case GLAC_MELT:
                    output = (mDimension == D3D) ? false : true;
                    break;
                case GLAC_SUBL:
                    output = (mDimension == D3D) ? false : true;
                    break;
                default:
                    output = false;
                    break;
            }

            return output;
        }
    }
}

