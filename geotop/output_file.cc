/* STATEMENT:

   GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
   GEOtop 2.0.0 - 9 Mar 2012

   Copyright (c), 2012 - Stefano Endrizzi 

   This file is part of GEOtop 2.0.0 

   GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

   GEOtop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
   If you just use the code, please give feedback to the authors and the community.
   Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.

   If you have satisfactorily used the code, please acknowledge the authors.

*/

/**
 * @file output_file.cc
 * @Author  Gianfranco Gallizia (skyglobe83@gmail.com)
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
            mCount = 0;
            mWhatIsValid = -1;
            mDValue = geotop::input::gDoubleNoValue;
            mVValue = NULL;
            mMValue = NULL;
            mTValue = NULL;
        }

        TemporaryValues::TemporaryValues(double init)
        {
            mCount = 1;
            mWhatIsValid = 0;
            mDValue = init;
            mVValue = NULL;
            mMValue = NULL;
            mTValue = NULL;
        }

        TemporaryValues::TemporaryValues(GeoVector<double>* init)
        {
            mCount = 1;
            mWhatIsValid = 1;
            mDValue = NULL;
            mVValue = init;
            mMValue = NULL;
            mTValue = NULL;
        }

        TemporaryValues::TemporaryValues(GeoMatrix<double>* init)
        {
            mCount = 1;
            mWhatIsValid = 2;
            mDValue = geotop::input::gDoubleNoValue;
            mVValue = NULL;
            mMValue = init;
            mTValue = NULL;
        }

        TemporaryValues::TemporaryValues(GeoTensor<double>* init)
        {
            mCount = 1;
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

                        lg->logsf(geotop::logger::WARNING,
                                  "Unknown output variable: '%s'.",
                                  tmp.c_str());
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

            convert_dateeur12_daymonthyearhourmin(dateeur12, &day, &month, &year, &hour, &min);

            sprintf(buffer, "%.4ld%.2ld%.2ld%.2ld%.2ld", year, month, day, hour, min);

            if (mDimension == D1Dp || mDimension == D1Ds) 
            {
                output = std::string(""); //Do not prepend the date for tab files
            }
            else
            {
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
                    output.append("SoilTemperature");
                    break;
                case SOIL_WATER_CONTENT:
                    output.append("SoilWaterContent");
                    break;
                default:
                    output.append("UNKNOWN");
                    break;
            }

            return output;
        }

        Variable OutputFile::str2var(std::string v)
        {
            Variable lVar = UNKNOWN_VAR;
            std::string tmp = toLower(v);

            if (tmp.compare("soiltemperature") == 0) lVar = SOIL_TEMP;
            if (tmp.compare("soilwatercontent") == 0) lVar = SOIL_WATER_CONTENT;

            return lVar;
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
                default:
                    output = false;
                    break;
            }

            return output;
        }
    }
}

