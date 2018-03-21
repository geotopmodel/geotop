/**
 * @file output_file.h
 * @Author  Gianfranco Gallizia (skyglobe83@gmail.com)
 * @date June, 2014
 * @copyright (C) 2014 eXact lab srl
 * @brief Output file definition class
 *
 * Defines a container for the output definition: output variable, dimension
 * (1D, 2D, 3D), type and period of integration.
 */

#ifndef OUTPUT_FILE_H
#define OUTPUT_FILE_H

#include <string>
#include <exception>
#include "struct.geotop.h"
#include "standard_names.h"

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
    enum Dimension { D1Dp, D1Ds, D2D, D3D, UNKNOWN_DIM };

    /**
     * @brief Type of integration
     *
     * - AVG Time Average (value += current_value/period)
     * - CUM Cumulate (value += current_value)
     * - INS Instant (value = current_value)
     */
    enum IntegrationType { AVG, CUM, INS, UNKNOWN_INTEG };

    /**
     * @brief Holds the temporary values for cumulates and averages
     */
    class TemporaryValues
    {
    public:
      TemporaryValues();
      TemporaryValues(double init);
      TemporaryValues(GeoVector<double> *init);
      TemporaryValues(GeoMatrix<double> *init);
      TemporaryValues(GeoTensor<double> *init);
      int whatIsValid() { return mWhatIsValid; }
      double getValueD();
      GeoVector<double> *getValuesV();
      GeoMatrix<double> *getValuesM();
      GeoTensor<double> *getValuesT();

    private:
      int mWhatIsValid;
      double mDValue;
      GeoVector<double> *mVValue;
      GeoMatrix<double> *mMValue;
      GeoTensor<double> *mTValue;
    };

    class OutputFile
    {
    public:
      /*=================================================================
       * Constructor and Destructor
       =================================================================*/
      /**
       * @brief Constructor
       * @param[in] extended_key the key in VARIABLE::DIMENSION::INTEGRATION
       * format
       * @param[in] period time of integration
       * @param[in] layer optional layer index
       */
      OutputFile(std::string extended_key,
                 double period,
                 long layer = 0L,
                 std::string prefix = std::string(""));
      virtual ~OutputFile();

      /*=================================================================
       * Methods
       =================================================================*/
      /**
       * @brief retrieves the output file name
       * @param[in] dateeur12 the date in dateeur12 format (see times.cc)
       * defaults to 1/1/1900 00:00
       * @param[in] layer optional layer index (set to -1 to omit)
       * @return a string with the file's name
       */
      std::string getFileName(double dateeur12 = 10119000000.0,
                              long layer = -1L);

      /**
       * @brief retrieves the output file path (based on prefix)
       * @param[in] dateeur12 the date in dateeur12 format (see times.cc)
       * @param[in] layer optional layer index (set to -1 to omit)
       * @return a string with the file's name
       */
      std::string getFilePath(double dateeur12 = 10119000000.0,
                              long layer = -1L);

      /**
       * @brief converts the output file to a std::string
       * @return a string that contains the output file specification
       */
      std::string toString();

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
  }  // namespace input
}  // namespace geotop
#endif
