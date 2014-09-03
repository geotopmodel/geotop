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
                default:
                    output = false;
                    break;
            }

            return output;
        }
    }
}

