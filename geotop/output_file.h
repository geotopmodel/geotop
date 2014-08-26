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
 * @file output_file.h
 * @Author  Gianfranco Gallizia (skyglobe83@gmail.com)
 * @date June, 2014
 * @brief Output file definition class
 * 
 * Defines a container for the output definition: output variable, dimension (1D, 2D, 3D),
 * type and period of integration.
 */

#ifndef OUTPUT_FILE_H
#define OUTPUT_FILE_H

#include <string>
#include <exception>
#include "struct.geotop.h"

namespace geotop
{
    namespace input
    {
        /**
         * @brief Output Dimensions
         *
         * - D1Dp Single point table
         * - D1Ds Whole domain mean table
         * - D2D  Single layer map
         * - D3D  Multiple maps (one per layer) 
         */
        enum Dimension {
            D1Dp,
            D1Ds,
            D2D,
            D3D,
            UNKNOWN_DIM
        };

        /**
         * @brief Type of integration
         *
         * - AVG Time Average (value += current_value/period)
         * - CUM Cumulate (value += current_value)
         * - INS Instant (value = current_value)
         */
        enum IntegrationType {
            AVG,
            CUM,
            INS,
            UNKNOWN_INTEG
        };

        /**
         * @brief Output variables known
         */
        enum Variable {
            SOIL_TEMP, //Soil temperature
            SOIL_WATER_CONTENT, //Water content in soil [mm]
            UNKNOWN_VAR
        };

        /**
         * @brief Holds the temporary values for cumulates and averages
         */
        class TemporaryValues
        {
            public:
                TemporaryValues();
                TemporaryValues(double init);
                TemporaryValues(GeoMatrix<double>* init);
                TemporaryValues(GeoTensor<double>* init);
                int whatIsValid() { return mWhatIsValid; }
                double getValueD();
                GeoMatrix<double>* getValuesM();
                GeoTensor<double>* getValuesT();
            private:
                int mWhatIsValid;
                double mDValue;
                GeoMatrix<double>* mMValue;
                GeoTensor<double>* mTValue;

        } ;

        class OutputFile
        {
        public:
            /*=================================================================
             * Constructor and Destructor
             =================================================================*/
            /**
             * @brief Constructor
             * @param[in] extended_key the key in VARIABLE::DIMENSION::INTEGRATION format
             * @param[in] period time of integration
             * @param[in] layer optional layer index
             */
            OutputFile(std::string extended_key, double period, long layer = 0L, std::string prefix = std::string(""));
            virtual ~OutputFile();

            /*=================================================================
             * Methods
             =================================================================*/
            /**
             * @brief retrieves the output file name
             * @param[in] dateeur12 the date in dateeur12 format (see times.cc)
             * @param[in] layer optional layer index (set to -1 to omit)
             * @return a string with the file's name
             */
            std::string getFileName(double dateeur12, long layer = -1L);
            
            /**
             * @brief retrieves the output file path (based on prefix)
             * @param[in] dateeur12 the date in dateeur12 format (see times.cc)
             * @param[in] layer optional layer index (set to -1 to omit)
             * @return a string with the file's name
             */
            std::string getFilePath(double dateeur12, long layer = -1L);

            /*=================================================================
             * Read-only Properties
             =================================================================*/
            geotop::input::Variable getVariable() { return mVariable; }
            geotop::input::Dimension getDimension() { return mDimension; }
            geotop::input::IntegrationType getIntegrationType() { return mType; }
            long getPeriod() { return mPeriod; }
            long getLayer() { return mLayerIndex; }
            std::string getPrefix() { return std::string(mPrefix); }

            /*=================================================================
             * Public Fields
             =================================================================*/
            geotop::input::TemporaryValues values;

            /*=================================================================
             * Static methods
             =================================================================*/
            /**
             * @brief Converts a string to a Variable
             */
            static Variable str2var(std::string v);
            /**
             * @brief Converts a Variable to a string
             */
            static std::string var2str(Variable v);
        private:
            bool isValidDimension();
            std::string mPrefix;
            geotop::input::Variable mVariable;
            geotop::input::Dimension mDimension;
            geotop::input::IntegrationType mType;
            long mPeriod;
            long mLayerIndex;
        };
    }
}
#endif

