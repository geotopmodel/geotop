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

namespace geotop
{
    namespace input
    {
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

        static long rounding(double d)
        {
            d += 0.5;
            return (long)d ;
        }

        static long floor(double d)
        {
            return (long)d;
        }

        static void convert_dateeur12_daymonthyearhourmin(double date, long *day, long *month, long *year, long *hour, long *min)
        {

            *day = floor(date/1.E10);
            *month = floor(date/1.E8 - (*day)*1.E2);
            *year = floor(date/1.E4 - (*day)*1.E6 - (*month)*1.E4);
            *hour = floor(date/1.E2 - (*day)*1.E8 - (*month)*1.E6 - (*year)*1.E2 );
            *min = floor(date/1.E0 - (*day)*1.E10 - (*month)*1.E8 - (*year)*1.E4 - (*hour)*1.E2);

            if (*day<1 || *day>31 || *month<1 || *month>12 || *year<1700 || *year>2900 || *hour<0 || *hour>23 || *min<0 || *min>59) {
                geotop::input::DateEur12Exception e(date);
                throw e;
            }
        }

        const char* DateEur12Exception::what()
        {
            char buffer[512] = {'\0'};

            sprintf(buffer, "Invalid date: %f", mDate);

            return (const char*)strdup(buffer);
        }

        OutputFile::OutputFile(std::string extended_key, double period)
        {
            mVariable = geotop::input::UNKNOWN_VAR;
            mDimension = geotop::input::UNKNOWN_DIM;
            mType = geotop::input::UNKNOWN_INTEG;
            mPeriod = rounding(period);
            std::vector<std::string> values = split_ext_key(extended_key);

            if (values.size() == 3)
            {
                //Dimension
                std::string tmp = values.at(1);
                if (tmp.compare("1Dp") == 0) mDimension = D1Dp;
                if (tmp.compare("1Ds") == 0) mDimension = D1Ds;
                if (tmp.compare("2D") == 0) mDimension = D2D;
                if (tmp.compare("3D") == 0) mDimension = D3D;

                //Integration Type
                tmp = values.at(2);
                if (tmp.compare("AVG") == 0) mType = AVG;
                if (tmp.compare("CUM") == 0) mType = CUM;
                if (tmp.compare("INS") == 0) mType = INS;
            }

        }

        OutputFile::~OutputFile()
        {
        }

        std::string OutputFile::getFileName(double dateeur12)
        {
            char buffer[13] = {'\0'};
            long day = 0L, month = 0L, year = 0L, hour = 0L, min = 0L;

            convert_dateeur12_daymonthyearhourmin(dateeur12, &day, &month, &year, &hour, &min);

            sprintf(buffer, "%.4ld%.2ld%.2ld%.2ld%.2ld", year, month, day, hour, min);

            std::string output(buffer);
            output.append("_");

            switch (mVariable)
            {
                default:
                    output.append("UNKNOWN");
                    break;
            }

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

            output.append(".asc");
            
            return output;
        }
    }
}

