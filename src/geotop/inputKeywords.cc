/**
 * @file   inputKeywords.cc
 * @Author Angelo Leto (angleto@gmail.com)
 * @date   November, 2013
 * @brief  generic configuration store class
 *
 * Parse configuration file and store parameters on single container
 */

#include "inputKeywords.h"
#include <boost/spirit/include/classic.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/period.hpp>

#include <boost/date_time/local_time/local_time.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/date.hpp>
#include <boost/date_time/date_facet.hpp>
#include <boost/date_time/time_parsing.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "output_file.h"

#include "logger.h"
#include "global_logger.h"
using namespace geotop::logger;

//Uncomment the following line to enable the trace log.
//#define TRACELOG

#ifdef TRACELOG
static std::ofstream tlf("inputKeywords_trace.log");
static Logger trace_log((std::ostream*)(&tlf), TRACE);
#endif


using namespace boost::assign;

/** @internal
 * @brief convert a string to a double
 */
bool stringToDouble(std::string const& pString, double &pValue)
{
    std::istringstream lStringStream(pString);
    if (!(lStringStream >> pValue))
        return false ;
    return true ;
}

enum Sections
{
    GENERAL_SEC,
    OUTPUT_SEC
};

/** @internal
 * @brief configuration file parsing class, this class
 *  must not be used, in the future if more configuration format will be
 *  available, this class should be handled using the decorator pattern.
 */
class ConfGrammar : public boost::spirit::classic::grammar<ConfGrammar>
{
public:
    
    ConfGrammar (boost::shared_ptr< std::map<std::string, boost::any> > pMap){
        mCurrent_section = boost::shared_ptr<Sections>(new Sections);
        *mCurrent_section = GENERAL_SEC;
        mMap = pMap ;
        mUnsupportedKeys.push_back("NumSimulationTimes") ;
        mKey = boost::shared_ptr< std::string >( new std::string() ) ;
        mCurrentPrefix  = boost::shared_ptr< std::string >( new std::string("") ) ;
        mCurrentLayerIndex = boost::shared_ptr<double>(new double);
        *mCurrentLayerIndex  = 0.;
    };
    
    ~ConfGrammar (){
        
    };

    void actionSection(char const * const pBegin, char const * const pEnd) const
    {
        std::string lKey = std::string(pBegin, pEnd);
        std::string lLowercaseKey ( lKey );
        boost::algorithm::to_lower(lLowercaseKey);
        *mCurrent_section = GENERAL_SEC; //By default set the section to General

#ifdef TRACELOG
        trace_log.logsf(TRACE, "actionSection: %s", lKey.c_str());
#endif

        if (lLowercaseKey.compare("output") == 0)
        {
            *mCurrent_section = OUTPUT_SEC;
            if (mMap->count("OUTPUT_SEC") == 0)
            {
                (*mMap)["OUTPUT_SEC"] = boost::shared_ptr< std::vector<geotop::input::OutputFile> >(new std::vector<geotop::input::OutputFile>());
#ifdef TRACELOG
        trace_log.log("actionSection: created OUTPUT_SEC key", TRACE);
#endif
            }
        }
    }

    void actionOutputKey(char const * const pBegin, char const * const pEnd) const
    {
        std::string lKey = std::string(pBegin, pEnd);

#ifdef TRACELOG
        trace_log.logsf(TRACE, "actionOutputKey: %s", lKey.c_str());
#endif

        if (*mCurrent_section == OUTPUT_SEC)
            mKey->assign(lKey);
        else
            mKey->assign("");
    }

    void actionKey(char const * const pBegin, char const * const pEnd) const
    {
        std::string lKey = std::string(pBegin, pEnd) ;
        std::string lLowercaseKey ( lKey );
        boost::algorithm::to_lower(lLowercaseKey);

#ifdef TRACELOG
        trace_log.logsf(TRACE, "actionKey: %s", lKey.c_str());
#endif
        for (size_t it = 0; it < mUnsupportedKeys.size(); ++it)
        {
            std::string ukl(mUnsupportedKeys.at(it));
            boost::algorithm::to_lower(ukl);
            if (lKey.compare(ukl) == 0)
            {
                GlobalLogger* lg = GlobalLogger::getInstance();
                lg->logsf(WARNING, "The %s parameter is no longer supported. Please update your configuration file.",
                          lKey.c_str());
                return;
            }
        }

        mKey->assign( lLowercaseKey ) ;
    };
    
    void actionValueString(char const * const pBegin, char const * const pEnd) const
    {
        std::string lValue = std::string(pBegin, pEnd) ;
        
        //std::cout << "StringValue: " << *mKey << ":" << lValue << std::endl;
#ifdef TRACELOG
        trace_log.logsf(TRACE, "StringValue: %s:%s", ((std::string)*mKey).c_str(), lValue.c_str());
#endif
        switch (*mCurrent_section)
        {
            case GENERAL_SEC:
                if( mKey->compare ("") != 0 )
                {
                    (*mMap)[*mKey] = lValue;
                } else {
                    GlobalLogger* lg = GlobalLogger::getInstance();
                    lg->log("actionValueString : no key was pushed for the value, value will be discarded", ERROR);
                }
                break;
            case OUTPUT_SEC:
                if(mKey->compare("pathprefix") == 0)
                {
                    mCurrentPrefix->assign(lValue);
                }
                if(mKey->compare("") == 0)
                {
                    GlobalLogger* lg = GlobalLogger::getInstance();
                    lg->log("actionValueString : no key was pushed for the value, value will be discarded", ERROR);
                }
                break;
        }

        mKey->assign( "" );

    };
    
    void actionValueDate(char const * const pBegin, char const * const pEnd) const
    {
        std::string lValue = std::string(pBegin, pEnd) ;
     
        GlobalLogger* lg = GlobalLogger::getInstance();
#ifdef TRACELOG
        trace_log.logsf(TRACE, "DateValue: %s:%s", ((std::string)*mKey).c_str(), lValue.c_str());
#endif

        if (*mCurrent_section == OUTPUT_SEC)
        {
            lg->log("Invalid input file format: date value not allowed for output file period", CRITICAL);
            exit(1);
        }

        std::stringstream lStringStream;
        lStringStream.imbue(std::locale(lStringStream.getloc(), new boost::posix_time::time_input_facet("%d/%m/%Y %H:%M")));
        lStringStream.exceptions(std::ios_base::failbit);
        lStringStream << lValue;
        boost::posix_time::ptime lDate(boost::date_time::not_a_date_time) ;
        lStringStream >> lDate ;

        std::stringstream lOutStringStream;
        lOutStringStream.imbue(std::locale(lStringStream.getloc(), new boost::posix_time::time_facet("%d%m%Y%H%M")));
        lOutStringStream << lDate ;
        
        double lDoubleEncodedDate = geotop::input::gDoubleNoValue ;
        if ( not stringToDouble(lOutStringStream.str(), lDoubleEncodedDate) )
        {
            lg->logsf(ERROR, "Unable to convert string to double: %s", lOutStringStream.str().c_str());
            return ;
        }

        if( mKey->compare ("") != 0 )
        {
            boost::any lAny = (*mMap)[*mKey] ;
            if(lAny.type() == typeid(std::vector<double>))
            {
                std::vector<double> lValue ;
                lValue.push_back(lDoubleEncodedDate) ;
                (*mMap)[*mKey] = lValue;
            } else {
                (*mMap)[*mKey] = lDoubleEncodedDate;
            }
        } else {
            lg->log("actionValueDate : no key was pushed for the value, the value will be discarded", ERROR);
        }
        mKey->assign( "" ) ;
    };
    
    void actionValueDouble(const double pValue) const
    {
        GlobalLogger* lg = GlobalLogger::getInstance();
#ifdef TRACELOG
        trace_log.logsf(TRACE, "actionValueDouble: %f", pValue);
#endif
        switch(*mCurrent_section)
        {
            case GENERAL_SEC:
                if( mKey->compare ("") != 0 )
                {
                    boost::any lAny = (*mMap)[*mKey] ;
                    if(lAny.type() == typeid(std::vector<double>))
                    {
                        std::vector<double> lValue ;
                        lValue.push_back(pValue) ;
                        (*mMap)[*mKey] = lValue;
                    } else {
                        (*mMap)[*mKey] = pValue;
                    }
                } else {
                    lg->log("actionValueDouble : no key was pushed for the value, the value will be discarded", ERROR);
                }
                break;
            case OUTPUT_SEC:
                if (mKey->compare("") != 0)
                {
                    //Fetch the key and convert it to lower case
                    std::string lcKey = std::string(*mKey);
                    boost::algorithm::to_lower(lcKey);

                    //If it's not a LayerIndex key create a new OutputFile
                    if (lcKey.compare("layerindex") != 0)
                    {
                        boost::shared_ptr< std::vector<geotop::input::OutputFile> > files = 
                            boost::any_cast< boost::shared_ptr< std::vector<geotop::input::OutputFile> > >(mMap->at("OUTPUT_SEC"));
                        files->push_back(geotop::input::OutputFile(*mKey, pValue, (long)(*mCurrentLayerIndex), (std::string)(*mCurrentPrefix)));
#ifdef TRACELOG
                        trace_log.logsf(TRACE, "actionValueDouble: created new output file: %s %f", mKey->c_str(), pValue);
                        trace_log.logsf(TRACE, "actionValueDouble: files' size: %u", files->size());
#endif
                    }
                    else
                    {
                        //update  the current layer index
                        *mCurrentLayerIndex = pValue;
                    }
                }
                else
                {
                    lg->log("actionValueDouble : no key was pushed for the value, the value will be discarded", ERROR);
                }
                break;
        }
        
        mKey->assign( "" ) ;
    };
    
    void actionValueDoubleArray(char const * const pBegin, char const * const pEnd) const
    {
        GlobalLogger* lg = GlobalLogger::getInstance();
        std::string lString = std::string(pBegin, pEnd) ;
        boost::algorithm::erase_all(lString, " ");

        if (*mCurrent_section == OUTPUT_SEC)
        {
            lg->log("Invalid input file format: too many values for output file period", CRITICAL);
            exit(1);
        }
        
        std::vector<std::string> lTokenArray;
        boost::split(lTokenArray, lString, boost::algorithm::is_any_of(","));
        
        std::vector<double> lDoubleArray ;
        BOOST_FOREACH(std::string lToken, lTokenArray)
        {
            lDoubleArray.push_back( boost::lexical_cast<double>(lToken) );
        }
        
        if( mKey->compare ("") != 0 )
        {
            (*mMap)[*mKey] = lDoubleArray;
        } else {
            lg->log("actionValueDoubleArray : no key was pushed for the value, the value will be discarded", ERROR);
        }
        mKey->assign( "" ) ;
        
        //std::cout << "actionValueDoubleArray Size: " << lDoubleArray.size() << std::endl;
#ifdef TRACELOG
        trace_log.logsf(TRACE, "actionValueDoubleArray Size: %u", lDoubleArray.size());
#endif
    };

    void actionTrace(char const * const pBegin, char const * const pEnd) const
    {
        std::string lValue = std::string(pBegin, pEnd);
        GlobalLogger* lg = GlobalLogger::getInstance();
        lg->logsf(TRACE, "actionTrace: %s", lValue.c_str());
#ifdef TRACELOG
        trace_log.logsf(TRACE, "actionTrace: %s", lValue.c_str());
#endif
    };
    
    /** @internal
     * @brief specification of the Extended Bakus Naur form grammar
     */
    template <typename Scanner>
    struct definition
    {
        definition(const ConfGrammar &self)
        {
            configfile = +(row);
            row = blanks >> (section | comment | parameter ) >> blanks | +(boost::spirit::classic::space_p);
            blanks = *(boost::spirit::classic::space_p) ;
            section = "[" >> section_ident[boost::bind(&ConfGrammar::actionSection, self, _1, _2)] >> "]";
            section_ident = (+(boost::spirit::classic::alpha_p) >> *(boost::spirit::classic::alpha_p | boost::spirit::classic::ch_p('_')));
            parameter = (output_key[boost::bind(&ConfGrammar::actionOutputKey, self, _1, _2)] | standard_key[boost::bind(&ConfGrammar::actionKey, self, _1, _2)]) >> blanks >> "=" >>
                blanks >> value ;
            standard_key = +(boost::spirit::classic::alpha_p|boost::spirit::classic::alnum_p) ;
            output_key = output_var >> "::" >> output_dim >> "::" >> output_type ;
            output_var = boost::spirit::classic::alpha_p >> *(boost::spirit::classic::alnum_p);
            output_dim = (boost::spirit::classic::str_p("1Dp") | boost::spirit::classic::str_p("1Ds") | boost::spirit::classic::str_p("2D") | boost::spirit::classic::str_p("3D"));
            output_type = (boost::spirit::classic::str_p("AVG") | boost::spirit::classic::str_p("CUM") | boost::spirit::classic::str_p("INS"));
            value = date | array | num | string ;
            num = (boost::spirit::classic::real_p)[boost::bind(&ConfGrammar::actionValueDouble, self, _1)] ;
            date_array = (date >> +(blanks >> "," >>
                                                     blanks >> date))[boost::bind(&ConfGrammar::actionValueDate, self, _1, _2)] ;
            array = (boost::spirit::classic::real_p >> +(blanks >> "," >>
                                                blanks >> boost::spirit::classic::real_p))[boost::bind(&ConfGrammar::actionValueDoubleArray, self, _1, _2)] ;
            string = "\"" >> (*(~(boost::spirit::classic::ch_p("\""))))[boost::bind(&ConfGrammar::actionValueString, self, _1, _2)] >> "\"" ;
            comment = boost::spirit::classic::comment_p("!") | boost::spirit::classic::comment_p("#");
            date = (DD >> "/" >> MM >> "/" >> YYYY >>
                    blanks >> hh >> ":" >> mm)[boost::bind(&ConfGrammar::actionValueDate, self, _1, _2)] ;
            DD = boost::spirit::classic::digit_p >> boost::spirit::classic::digit_p ;
            MM = boost::spirit::classic::digit_p >> boost::spirit::classic::digit_p ;
            YYYY = boost::spirit::classic::digit_p >> boost::spirit::classic::digit_p >>
            boost::spirit::classic::digit_p >> boost::spirit::classic::digit_p ;
            hh = boost::spirit::classic::digit_p >> boost::spirit::classic::digit_p ;
            mm = boost::spirit::classic::digit_p >> boost::spirit::classic::digit_p ;
        }
        
        const boost::spirit::classic::rule<Scanner> &start()
        {
            return configfile;
        }
        
    private:
        
        boost::spirit::classic::rule<Scanner> configfile,
        row,
        section,
        section_ident,
        parameter,
        standard_key,
        output_key,
        output_var,
        output_dim,
        output_type,
        value,
        date_array,
        array,
        string,
        num,
        blanks,
        date, DD, MM, YYYY, hh, mm,
        comment ;
    };
    
private:
    
    std::vector<std::string> mUnsupportedKeys; //No longer supported keys
    boost::shared_ptr<std::string> mKey ;
    boost::shared_ptr<std::string> mCurrentPrefix ;
    boost::shared_ptr< std::map<std::string, boost::any> > mMap ;
    boost::shared_ptr<Sections> mCurrent_section;
    boost::shared_ptr<double> mCurrentLayerIndex;
};

geotop::input::ConfigStore::ConfigStore()
{
    mValueMap = boost::shared_ptr<std::map<std::string, boost::any> > (new std::map<std::string, boost::any> () ) ;
}

bool geotop::input::ConfigStore::parse(const std::string pFileName) {

    GlobalLogger* lg = GlobalLogger::getInstance();
    std::ifstream fs(pFileName.c_str());
    std::ostringstream ss;
    ss << fs.rdbuf();
    
    init() ;
    
    ConfGrammar lConfGrammar (mValueMap) ;
    boost::spirit::classic::parse_info<> lParserInfo = boost::spirit::classic::parse(ss.str().c_str(), lConfGrammar);
    if (lParserInfo.hit)
    {
        lg->logsf(NOTICE, "Configuration file parsing completed. %u characters parsed.", lParserInfo.length);
    }
    else
    {
        lg->log("Parsing failed.", ERROR);
        return false ;
    }

    if (lParserInfo.full)
    {
        lg->log("Full match!", NOTICE);
    }
    else
    {
        char tmpbuf[21] = {0};
        size_t i;
        for (i=0; i<21; i++) {
            if (lParserInfo.stop && lParserInfo.stop[i] != '\n')
                tmpbuf[i] = lParserInfo.stop[i];
            else
                break;
        }
        lg->logsf(WARNING, "Partial parsing, matched length(%u) stopped at '%.20s'", lParserInfo.length, tmpbuf);
        return false ;
    }

    return true ;
}

/** @internal
 * @brief initialization of the configuration parameters
 */
void geotop::input::ConfigStore::init()
{
    mValueMap->clear() ;
    
    //BEGIN INITIALIZATION OF STRING PARAMETERS
    initValue("SpecificPlotSurfaceSensibleHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverGlacierLayerThick", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderMeteoStationLongitude", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("NetLongwaveRadiationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("GlacierProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLEBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderMeteoStationElevation", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSnowOnCanopy", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilAveragedTempProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWbeamPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("Headerz0vegPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWglobal", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWinPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWvBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWatContentSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointLandCoverType", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLEvPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLEvBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SubfolderRecoveryFiles", std::string("rec")) ;
    
    initValue("RecoverRunSoilMaximumTemperatureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilAveragedLiqContentProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SuccessfulRunFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverVegTemp", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLEPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderCoordinatePointY", std::string("Ywgs")) ;
    
    initValue("HeaderCoordinatePointX", std::string("Xwgs")) ;
    
    initValue("HeaderDateGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPeriodGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotVegTempMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderN", std::string("n")) ;
    
    initValue("HeaderV", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("NetRadiationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("PointOutputFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderCanopyFractionPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderRunSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLapseRateTemp", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHorizonAngle", std::string("azimuth")) ;
    
    initValue("HeaderLWNetPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderMeteoStationLatitude", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderIceContentSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSoilDz", std::string("Dz")) ;
    
    initValue("EvapotranspirationFromSoilMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWdirect", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("AspectMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("GlacierMeltedMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTimeFromStartSoil", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowCoveredAreaFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWESublBlownPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotSurfaceLatentHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilTempProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSnowTempPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSoilTemp", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("NetPrecipitationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HNMapFile", std::string(geotop::input::gStringNoValue)) ;//TODO mattiu

    initValue("SurfaceHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPeriodPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilAveragedIceContentTensorFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDateDDMMYYYYhhmmLapseRates", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("AirTempMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointLatitude", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTSurfBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("MeteoStationsListFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotSurfaceWaterContentMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("CanopyInterceptedWaterMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLapseRateDewTemp", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSurfaceTemperature", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilLiqWaterPressProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderMeteoStationCoordinateY", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverRunSoilAveragedSnowWaterEquivalentFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSnowLiqMass", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLEgUnvegPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilAveragedIceContentProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWupPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("InitWaterTableDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTvegBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("LandCoverMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("LapseRateFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverRunSoilMinimumTotalSoilMoistureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("WindSpeedMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderJulianDayFromYear0Glac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSnowTemp", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("InLongwaveRadiationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("GlacierProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderAirPressPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLapseRatePrec", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointCurvatureWestEastDirection", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWinBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SkyViewFactorMapFile", std::string("input_maps/sky")) ;
    
    initValue("HeaderIDPointPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotIncomingShortwaveRadMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowLiqContentProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderEvapSurfaceBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLEupPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTraspCanopyPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SurfaceSensibleHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDepthSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWvBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWvPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderGlacDensityPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPeriodBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPRainBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWindVelocity", std::string("WindSp")) ;
    
    initValue("LandSurfaceWaterDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("PointOutputFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSurfaceEBPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilSaturationRatioProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDecayKCanopyPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilTempTensorFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTCanopyAirPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointSoilType", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTempGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotVegSensibleHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotVegLatentHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverGlacierTemp", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWupPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLateralHydrConductivity", std::string("Kh")) ;
    
    initValue("HeaderSoilInitTemp", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("DemFile", std::string("input_maps/pit")) ;
    
    initValue("RecoverRunSoilMinimumTemperatureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotCanopyAirTempMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderCloudSWTransmissivity",std::string(geotop::input::gStringNoValue) ) ;
    
    initValue("SoilAveragedLiqContentTensorFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilTotWaterPressProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWNetBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderNormalHydrConductivity", std::string("Kv")) ;
    
    initValue("HeaderWindDirection", std::string("WindDir")) ;
    
    initValue("HeaderDateSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderIPrec", std::string("Iprec")) ;
    
    initValue("InitSnowDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLObukhovPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("InitGlacierDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderMeanTimeStep", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilLiqWaterPressProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotTotalSensibleHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderJulianDayfrom0Meteo", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderCloudFactor", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointMaxSWE", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("WaterTableDepthFromAboveMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWEPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPSnowBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotWindDirMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverRunSoilAveragedTemperatureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointLongitude", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDewTemp", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("BedrockDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RunSoilMinimumTotalSoilMoistureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPeriodSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSnowMeltedPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSnowIceMass", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderJulianDayFromYear0Snow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("InitSnowAgeMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWindSpeedTopCanopyPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilLiqWaterPressTensorFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverGlacierLayerNumber", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderThetaRes", std::string("res")) ;
    
    initValue("HorizonMeteoStationFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWinBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverRunSoilAveragedInternalEnergyFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderEstoredCanopyPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWEBlownPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSnowSublPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("InitSWEMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RunSoilMaximumTemperatureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPrainPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderGlacDepthPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("DaysDelayMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("TimeDependentIncomingDischargeFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWdiffuse", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("FirstSoilLayerIceContentMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("FailedRunFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilIceContentProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotWindSpeedMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWinMaxPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderAirTempBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("InShortwaveRadiationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHighestWaterTableDepthPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTimeFromStartSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("ChannelSurfaceWaterDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPsnowPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderCthSoilSolids", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHighestThawedSoilDepthPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowIceContentProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPSnowNetBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSoilInitPres", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderRunPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SWEMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHupPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("WaterTableDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderGlacMeltedPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointSlope", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTempSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("WindDirMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTimeStepAverage", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("PrecipitationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLEgVegPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWiltingPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSoilIceContChannel", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWinPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDateDDMMYYYYhhmmMeteo", std::string("Date")) ;
    
    initValue("DirectInShortwaveRadiationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWindSpeedPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("PointFile", std::string("ListPoints")) ;
    
    initValue("RecoverTime", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotSnowDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RunSoilAveragedTotalSoilMoistureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowDurationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SurfaceTempMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLSAIPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RunSoilAveragedInternalEnergyFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderQSurfPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWindDirPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPrainOnSnowPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilIceContentProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilTotWaterPressTensorFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWaterOnCanopyPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilLiqContentProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderAirTempPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotNetVegShortwaveRadMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilTempProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotNetSurfaceShortwaveRadMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSoilWatPres", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotSurfaceHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWNetBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointID", std::string("ID")) ;
    
    initValue("SnowTempProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotSurfaceTempMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWatContentGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RunSoilMinimumTemperatureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("ShadowFractionTimeMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("GlacierDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("BasinOutputFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSoilTempChannel", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("FirstSoilLayerTempMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSnowOnCanopyPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSoilWatPresChannel", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderIDMeteoStation", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTDewPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSnowLayerNumber", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SurfaceLatentHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderGWEPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointCurvatureNorthwestSoutheastDirection", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWnet", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotNetVegLongwaveRadMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("ThawedSoilDepthFromAboveMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowLiqContentProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHgUnvegPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowDepthLayersFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotIncomingLongwaveRadMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderRH", std::string("RH")) ;
    
    initValue("SoilAveragedTempTensorFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilTotWaterPressProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("GlacierWaterEqMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("FirstSoilLayerLiqContentMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointHorizon", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverGlacierLiqMass", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SlopeMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderMeteoStationSkyViewFactor", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RelHumMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("GlacierSublimatedMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotAboveVegAirTempMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("TimeDependentVegetationParameterFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RunSoilMaximumTotalSoilMoistureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderMeteoStationStandardTime", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPrainNetPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLowestThawedSoilDepthPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTimeFromStartBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverRunSoilAveragedTotalSoilMoistureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSnowDensityPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPrec", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderQVegPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderIDPointGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTsurfPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("TimeStepsFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("ThawedSoilDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderJulianDayFromYear0Soil", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLObukhovCanopyPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderJulianDayFromYear0Basin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTimeFromStartPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSnowLayerThick", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderIDPointSnow", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointCurvatureNorthSouthDirection", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HorizonPointFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTimeFromStartGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SuccessfulRecoveryFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("FirstSoilLayerAveragedTempMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderRunGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowSublMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowDepthMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotRelHumMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointSkyViewFactor", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("MeteoFile", std::string("meteo/meteo")) ;
    
    initValue("HeaderPointElevation", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("BasinOutputFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilAveragedLiqContentProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTraspCanopyBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RunSoilAveragedSnowWaterEquivalentFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWinMinPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHvPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWindX", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderQCanopyAirPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDepthGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointBedrockDepth", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderFieldCapacity", std::string("fc")) ;
    
    initValue("HeaderQAirPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowIceContentProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderRunBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("NetShortwaveRadiationMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointCurvatureNortheastSouthwestDirection", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSoilHeatFluxPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotNetSurfaceLongwaveRadMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSpecificStorativity", std::string("SS")) ;
    
    initValue("HeaderAlpha", std::string("a")) ;
    
    initValue("HeaderMeteoStationCoordinateX", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointDepthFreeSurface", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSnowDepthPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderIDPointSoil", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("CurvaturesMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderJulianDayFromYear0Point", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPRainNetBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverRunSoilMaximumTotalSoilMoistureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPsnowNetPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWNetPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderIceContentGlac", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderTvegPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowMeltedMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderGlacTempPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderAirTemp", std::string("AirT")) ;
    
    initValue("RecoverLiqWaterOnCanopy", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("DischargeFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverSoilIceCont", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHgVegPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPeriodSoil", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilLiqContentProfileFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderRunSoil", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDateSoil", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHvBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("Headerd0vegPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SpecificPlotTotalLatentHeatFluxMapFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilLiqContentTensorFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverGlacierIceMass", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderThetaSat", std::string("sat")) ;
    
    initValue("SoilAveragedIceContentProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderLWvPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderKthSoilSolids", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RecoverNonDimensionalSnowAge", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RunSoilAveragedTemperatureFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderSWdiffPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("RiverNetwork", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderWindY", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDateBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPNetBasin", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderHorizonHeight", std::string("horizon_ele")) ;
    
    initValue("HeaderLowestWaterTableDepthPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilIceContentTensorFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderPointAspect", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilAveragedTempProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderGlacSublPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderDatePoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderEvapSurfacePoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowTempProfileFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("HeaderRHPoint", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SoilParFile", std::string(geotop::input::gStringNoValue)) ;
    
    initValue("SnowDepthLayersFileWriteEnd", std::string(geotop::input::gStringNoValue)) ;
    //END INITIALIZATION OF STRING PARAMETERS
    
    //BEGIN INITIALIZATION OF NUMERIC PARAMETERS
    initValue("SWNetBasin", double(-1)) ;
    
    initValue("BasinAll", double(0)) ;
    
    initValue("TraspCanopyPoint", double(-1)) ;
    
    initValue("CanopyMaxIter", double(3)) ;
    
    std::vector<double> lAlphaVanGenuchten ;
    lAlphaVanGenuchten += 0.004,0.004,0.004,0.004,0.004 ;
    initValue("AlphaVanGenuchten", lAlphaVanGenuchten) ;
    
    initValue("SnowRoughness", double(0.1)) ;
    
    initValue("MaxPrecDecreaseFactorWithElev", double(4.4)) ;
    
    initValue("DecayKCanopyPoint", double(-1)) ;
    
    initValue("EvapSurfacePoint", double(-1)) ;
    
    std::vector<double> lLeafAngles ;
    lLeafAngles += 0 ;
    initValue("LeafAngles", lLeafAngles) ;
    
    initValue("OutputDepthsVertical", double(0)) ;
    
    std::vector<double> lOutputVegetationMaps ;
    lOutputVegetationMaps += 0 ;
    initValue("OutputVegetationMaps", lOutputVegetationMaps) ;
    
    initValue("ExitMinLambdaEnergy", double(0)) ;
    
    std::vector<double> lSpinUpLayerBottom ;
    lSpinUpLayerBottom += 10000 ;
    initValue("SpinUpLayerBottom", lSpinUpLayerBottom) ;
    
    initValue("FreeDrainageAtLateralBorder", double(1)) ;
    
    initValue("SWEPoint", double(-1)) ;
    
    initValue("DepthFreeSurfaceAtTheBoundary", double(0)) ;
    
    initValue("MaxCourantSupFlowChannel", double(0.1)) ;
    
    initValue("DateBasin", double(-1)) ;
    
    initValue("GlacDensityPoint", double(-1)) ;
    
    initValue("ExitMinLambdaWater", double(1)) ;
    
    initValue("InitSnowTemp", double(0)) ;
    
    std::vector<double> lNVanGenuchtenBedrock ;
    lNVanGenuchtenBedrock += geotop::input::gDoubleNoValue,
        geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,
        geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("NVanGenuchtenBedrock", lNVanGenuchtenBedrock) ;

    initValue("DateSoil", double(-1)) ;

    initValue("JulianDayFromYear0Soil", double(-1)) ;

    initValue("TimeFromStartSoil", double(-1)) ;

    initValue("PeriodSoil", double(-1)) ;

    initValue("RunSoil", double(-1)) ;
    
    initValue("IDPointSoil", double(-1)) ;
    
    initValue("SoilHeatFluxPoint", double(-1)) ;
    
    initValue("WindAsWindXAndWindY", double(0)) ;
    
    initValue("z0vegPoint", double(-1)) ;
    
    initValue("PSnowBasin", double(-1)) ;
    
    std::vector<double> lSoilEmissiv ;
    lSoilEmissiv += 0.96;
    initValue("SoilEmissiv", lSoilEmissiv) ;
    
    initValue("DDChannel", double(1)) ;
    
    initValue("HeatEqMaxIter", double(700)) ;
    
    initValue("MaxGlacLayersMiddle", double(0)) ;
    
    initValue("NumberDayIntervalsToCalculateCloudiness", double(3)) ;
    
    initValue("TimeStepBlowingSnow", double(3600)) ;
    
    std::vector<double> lNormalHydrConductivityBedrock ;
    lNormalHydrConductivityBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("NormalHydrConductivityBedrock", lNormalHydrConductivityBedrock) ;
    
    initValue("NumLowPassFilterOnDemForCurv", double(0)) ;
    
    std::vector<double> lCanopyFraction ;
    lCanopyFraction += 0;
    initValue("CanopyFraction", lCanopyFraction) ;
    
    std::vector<double> lVegTransNIR ;
    lVegTransNIR += 0.2 ;
    initValue("VegTransNIR", lVegTransNIR) ;
    
    initValue("MaxCourantSupFlowChannelLand", double(0.1)) ;
    
    initValue("JulianDayFromYear0Glac", double(-1)) ;
    
    std::vector<double> lOutputGlacierMaps ;
    lOutputGlacierMaps += 0 ;
    initValue("OutputGlacierMaps", lOutputGlacierMaps) ;

    std::vector<double> lIDPointGlac ;
    lIDPointGlac += -1 ;
    initValue("IDPointGlac", lIDPointGlac) ;

    std::vector<double> lTempGlac ;
    lTempGlac += -1 ;
    initValue("TempGlac", lTempGlac) ;

    std::vector<double> lIceContentGlac ;
    lIceContentGlac += -1 ;
    initValue("IceContentGlac", lIceContentGlac) ;

    std::vector<double> lWatContentGlac ;
    lWatContentGlac += -1 ;
    initValue("WatContentGlac", lWatContentGlac) ;

    std::vector<double> lThetaSat ;
    lThetaSat += 0.5;
    initValue("ThetaSat", lThetaSat) ;
    
    std::vector<double> lVegSnowBurying ;
    lVegSnowBurying += 1;
    initValue("VegSnowBurying", lVegSnowBurying) ;
    
    std::vector<double> lGlacPlotDepths ;
    lGlacPlotDepths += geotop::input::gDoubleNoValue ;
    initValue("GlacPlotDepths", lGlacPlotDepths) ;
    
    initValue("NumLandCoverTypes", double(1)) ;
    
    initValue("SnowDensityPoint", double(-1)) ;
    
    std::vector<double> lLinearInterpolation ;
    lLinearInterpolation += 0 ;
    initValue("LinearInterpolation", lLinearInterpolation) ;
    
    initValue("PointSim", double(0)) ;
    
    initValue("LEvBasin", double(-1)) ;
    
    initValue("SWbeamPoint", double(-1)) ;
    
    initValue("SurfaceEnergyFlux", double(geotop::input::gDoubleNoValue)) ;
    
    initValue("DrySnowDefRate", double(1)) ;
    
    std::vector<double> lOutputSoilMaps ;
    lOutputSoilMaps += 1 ;
    initValue("OutputSoilMaps", lOutputSoilMaps) ;

    std::vector<double> lInitWaterTableDepth ;
    lInitWaterTableDepth += 5000.0 ;
    initValue("InitWaterTableDepth", lInitWaterTableDepth) ;
    
    initValue("Lozone", double(0.3)) ;
    
    initValue("RicalculateCloudiness", double(0)) ;
    
    std::vector<double> lSoilPlotDepths ;
    lSoilPlotDepths += geotop::input::gDoubleNoValue ;
    initValue("SoilPlotDepths", lSoilPlotDepths) ;
    
    initValue("BaseIPrec", double(0)) ;
    
    initValue("DDLand", double(1)) ;
    
    initValue("MinTimeStep", double(10)) ;
    
    std::vector<double> lCanDensSurface ;
    lCanDensSurface += 2 ;
    initValue("CanDensSurface", lCanDensSurface) ;
    
    initValue("RunBasin", double(-1)) ;
    
    std::vector<double> lNormalHydrConductivity ;
    lNormalHydrConductivity += 0.0001;
    initValue("NormalHydrConductivity", lNormalHydrConductivity) ;
    
    initValue("MinIceContentForBlowingSnow", double(8)) ;
    
    initValue("SnowThermalConductivityPar", double(1)) ;
    
    initValue("HPoint", double(-1)) ;
    
    std::vector<double> lSoilLayerThicknesses ;
    lSoilLayerThicknesses += geotop::input::gDoubleNoValue ;
    initValue("SoilLayerThicknesses", lSoilLayerThicknesses) ;
    
    initValue("SnowCorrFactor", double(1.3)) ;
    
    initValue("ThresWaterDepthLandInf", double(0)) ;
    
    initValue("DateSnow", double(-1)) ;

    initValue("JulianDayFromYear0Snow", double(-1)) ;
    
    initValue("TimeFromStartSnow", double(-1)) ;
    
    initValue("PeriodSnow", double(-1)) ;
    
    initValue("RunSnow", double(-1)) ;
    
    initValue("IDPointSnow", double(-1)) ;
    
    initValue("SlopeWeight", double(0)) ;
    
    initValue("WindDirPoint", double(-1)) ;
    
    initValue("FreshSnowReflNIR", double(0.65)) ;
    
    initValue("HighestThawedSoilDepthPoint", double(-1)) ;
    
    initValue("PRainNetBasin", double(-1)) ;
    
    initValue("CanopyFractionPoint", double(-1)) ;
    
    initValue("EnergyBalance", double(1)) ;
    
    initValue("LEgUnvegPoint", double(-1)) ;
    
    initValue("HvBasin", double(-1)) ;
    
    initValue("ZeroTempAmplitDepth", double(1e+20)) ;
    
    initValue("MaxSnowLayersMiddle", double(10)) ;
    
    initValue("AirTempBasin", double(-1)) ;
    
    initValue("CanopyStabCorrection", double(1)) ;
    
    initValue("LWinMaxPoint", double(-1)) ;
    
    initValue("TDewPoint", double(-1)) ;
    
    initValue("SnowMeltedPoint", double(-1)) ;
    
    initValue("SWinPoint", double(-1)) ;
    
    initValue("SnowDepthPoint", double(-1)) ;

    std::vector<double> lMeteoStationsID ;
    lMeteoStationsID += double(geotop::input::gDoubleNoValue) ;
    initValue("MeteoStationsID", lMeteoStationsID) ;

    initValue("QSurfPoint", double(-1)) ;
    
    std::vector<double> lMeteoStationWindVelocitySensorHeight ;
    lMeteoStationWindVelocitySensorHeight += 10 ;
    initValue("MeteoStationWindVelocitySensorHeight", lMeteoStationWindVelocitySensorHeight) ;
    
    initValue("ThresTempRain", double(3)) ;
    
    initValue("HighestWaterTableDepthPoint", double(-1)) ;
    
    initValue("EvapSurfaceBasin", double(-1)) ;
    
    initValue("GlacAll", double(0)) ;
    
    std::vector<double> lSnowPlotDepths ;
    lSnowPlotDepths += geotop::input::gDoubleNoValue ;
    initValue("SnowPlotDepths", lSnowPlotDepths) ;
    
    initValue("RichardTol", double(1e-06)) ;
    
    initValue("WindCompaction1D", double(0)) ;
    
    initValue("InitGlacierDepth", double(0)) ;

    std::vector<double> lOutputSnowMaps ;
    lOutputSnowMaps += 1 ;
    initValue("OutputSnowMaps",lOutputSnowMaps) ;
    
    initValue("SoilLayerNumber", double(5)) ;
    
    initValue("InitInNewPeriods", double(0)) ;
    
    initValue("ThresWaterDepthChannel", double(50)) ;
    
    initValue("MoninObukhov", double(1)) ;
    
    std::vector<double> lSpecificStorativityBedrock ;
    lSpecificStorativityBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("SpecificStorativityBedrock", lSpecificStorativityBedrock) ;
    
    std::vector<double> lDtPlotDischarge ;
    lDtPlotDischarge += 0 ;
    initValue("DtPlotDischarge", lDtPlotDischarge) ;
    
    std::vector<double> lVMualem ;
    lVMualem += 0.5;
    initValue("VMualem", lVMualem) ;
    
    std::vector<double> lWiltingPointBedrock ;
    lWiltingPointBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("WiltingPointBedrock", lWiltingPointBedrock) ;
    
    initValue("Longitude", double(11.7)) ;
    
    initValue("LWinBasin", double(-1)) ;
    
    initValue("EstoredCanopyPoint", double(-1)) ;
    
    initValue("IDPointPoint", double(-1)) ;
    
    std::vector<double> lThetaRes ;
    lThetaRes += 0.05;
    initValue("ThetaRes", lThetaRes) ;
    
    initValue("AlbExtParSnow", double(10)) ;

    std::vector<double> lThresSnowVegDown ;
    lThresSnowVegDown += 1000 ;
    initValue("ThresSnowVegDown", lThresSnowVegDown) ;
    
    std::vector<double> lSurFlowResLand ;
    lSurFlowResLand += 0.5;
    initValue("SurFlowResLand", lSurFlowResLand) ;
    
    initValue("SnowTempPoint", double(-1)) ;
    
    initValue("SWEBlownPoint", double(-1)) ;
    
    initValue("MaxTimesMinLambdaWater", double(0)) ;
    
    initValue("LWinPoint", double(-1)) ;
    
    initValue("LObukhovCanopyPoint", double(-1)) ;
    
    std::vector<double> lLSAI ;
    lLSAI += 1 ;
    initValue("LSAI", lLSAI) ;
    
    initValue("SurFlowResChannel", double(20)) ;
    
    initValue("PeriodPoint", double(-1)) ;
    
    initValue("AirPressPoint", double(-1)) ;
    
    initValue("DEMRotationAngle", double(0)) ;
    
    std::vector<double> lMeteoStationStandardTime ;
    lMeteoStationStandardTime += 0 ;
    initValue("MeteoStationStandardTime", lMeteoStationStandardTime) ;
    
    std::vector<double> lSavingPoints ;
    lSavingPoints += 0;
    initValue("SavingPoints", lSavingPoints) ;
    
    initValue("Iobsint", double(1)) ;
    
    initValue("SnowEmissiv", double(0.98)) ;
    
    initValue("MaxSnowPorosity", double(0.7)) ;
    
    initValue("KonzelmannB", double(8)) ;
    
    initValue("FlagSkyViewFactor", double(0)) ;
    
    initValue("KonzelmannA", double(0.484)) ;
    
    initValue("LWinMinPoint", double(-1)) ;
    
    initValue("HupPoint", double(-1)) ;
    
    initValue("ThresWaterDepthLandSup", double(0)) ;
    
    std::vector<double> lMeteoStationCoordinateY ;
    lMeteoStationCoordinateY += geotop::input::gDoubleNoValue ;
    initValue("MeteoStationCoordinateY", lMeteoStationCoordinateY) ;
    
    std::vector<double> lMeteoStationCoordinateX ;
    lMeteoStationCoordinateX += geotop::input::gDoubleNoValue ;
    initValue("MeteoStationCoordinateX", lMeteoStationCoordinateX) ;
    
    initValue("MeanTimeStep", double(-1)) ;
    
    std::vector<double> lTimeStepEnergyAndWater ;
    lTimeStepEnergyAndWater += 3600 ;
    initValue("TimeStepEnergyAndWater", lTimeStepEnergyAndWater) ;
    
    std::vector<double> lSoilRoughness ;
    lSoilRoughness += 10;
    initValue("SoilRoughness", lSoilRoughness) ;
    
    initValue("WetSnowDefRate", double(1.5)) ;
    
    std::vector<double> lMeteoStationLatitude ;
    lMeteoStationLatitude += 45 ;
    initValue("MeteoStationLatitude", lMeteoStationLatitude) ;
    
    initValue("MinSupWaterDepthLand", double(1)) ;
    
    std::vector<double> lFieldCapacity ;
    lFieldCapacity += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("FieldCapacity", lFieldCapacity) ;
    
    initValue("StandardTimeSimulation", double(0)) ;
    
    std::vector<double> lThresSnowVegUp ;
    lThresSnowVegUp += 0,200,0,200,0,1900,1900,800 ;
    initValue("ThresSnowVegUp", lThresSnowVegUp) ;
    
    initValue("BottomBoundaryHeatFlux", double(0)) ;
    
    initValue("GlacMeltedPoint", double(-1)) ;
    
    initValue("SurfaceTemperature", double(geotop::input::gDoubleNoValue)) ;
    
    initValue("WaterOnCanopyPoint", double(-1)) ;
    
    initValue("PrainNetPoint", double(-1)) ;
    
    initValue("DefaultSoilTypeLand", double(1)) ;
    
    initValue("PrecAsIntensity", double(0)) ;
    
    initValue("QVegPoint", double(-1)) ;
    
    initValue("TimeFromStartBasin", double(-1)) ;
    
    std::vector<double> lVegTransVis ;
    lVegTransVis += 0.2;
    initValue("VegTransVis", lVegTransVis) ;
    
    initValue("TsMaxIter", double(2)) ;
    
    initValue("RHmin", double(10)) ;
    
    initValue("WindSpeedTopCanopyPoint", double(-1)) ;
    
    initValue("PrainPoint", double(-1)) ;
    
    initValue("LWNetPoint", double(-1)) ;
    
    initValue("WaterBalance", double(0)) ;
    
    initValue("ZeroTempAmplitTemp", double(20)) ;
    
    initValue("BaseWindDirection", double(0)) ;
    
    initValue("LWNetBasin", double(-1)) ;
    
    initValue("RunIfAnOldRunIsPresent", double(1)) ;
    
    initValue("Surroundings", double(0)) ;
    
    initValue("Dn", double(1)) ;
    
    initValue("RichardMaxIter", double(100)) ;
    
    initValue("GlacDepthPoint", double(-1)) ;
    
    initValue("SWvBasin", double(-1)) ;
    
    initValue("LWvPoint", double(-1)) ;
    
    initValue("BlowingSnow", double(0)) ;
    
    std::vector<double> lThermalConductivitySoilSolids ;
    lThermalConductivitySoilSolids += 2.5;
    initValue("ThermalConductivitySoilSolids", lThermalConductivitySoilSolids) ;
    
    std::vector<double> lMeteoStationElevation ;
    lMeteoStationElevation += 0 ;
    initValue("MeteoStationElevation", lMeteoStationElevation) ;
    
    std::vector<double> lSoilAlbNIRDry ;
    lSoilAlbNIRDry += 0.33;
    initValue("SoilAlbNIRDry", lSoilAlbNIRDry) ;
    
    initValue("BaseRelativeHumidity", double(70)) ;

    std::vector<double> lInitDateDDMMYYYYhhmm ;
    lInitDateDDMMYYYYhhmm += double(010119000000.) ;
    initValue("InitDateDDMMYYYYhhmm", lInitDateDDMMYYYYhhmm) ;

    std::vector<double> lEndDateDDMMYYYYhhmm ;
    lEndDateDDMMYYYYhhmm += double(010119000000.) ;
    initValue("EndDateDDMMYYYYhhmm", lEndDateDDMMYYYYhhmm) ;

    std::vector<double> lSpecialPlotBegin ;
    lSpecialPlotBegin += double(0.) ;
    initValue("SpecialPlotBegin", lSpecialPlotBegin) ;
    
    std::vector<double> lSpecialPlotEnd ;
    lSpecialPlotEnd += double(0.) ;
    initValue("SpecialPlotEnd", lSpecialPlotEnd) ;

    initValue("CurvatureWeightI", double(0)) ;
    
    initValue("HBasin", double(-1)) ;
    
    std::vector<double> lVMualemBedrock ;
    lVMualemBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("VMualemBedrock", lVMualemBedrock) ;
    
    initValue("QAirPoint", double(-1)) ;
    
    initValue("PsnowPoint", double(-1)) ;
    
    initValue("CurvatureWeightD", double(0)) ;
    
    initValue("NumberOfMeteoStations", double(1)) ;
    
    initValue("TvegBasin", double(-1)) ;
    
    initValue("GlacTempPoint", double(-1)) ;

    std::vector<double> lDtPlotBasin ;
    lDtPlotBasin += 0;
    initValue("DtPlotBasin", lDtPlotBasin) ;
    
    initValue("UpwindBorderBlowingSnow", double(0)) ;
    
    std::vector<double> lThetaResBedrock ;
    lThetaResBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("ThetaResBedrock", lThetaResBedrock) ;
    
    std::vector<double> lAlphaVanGenuchtenBedrock ;
    lAlphaVanGenuchtenBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("AlphaVanGenuchtenBedrock", lAlphaVanGenuchtenBedrock) ;
    
    std::vector<double> lCoordinatePointX ;
    lCoordinatePointX += geotop::input::gDoubleNoValue ;
    initValue("CoordinatePointX", lCoordinatePointX) ;

    std::vector<double> lCoordinatePointY ;
    lCoordinatePointY += geotop::input::gDoubleNoValue ;
    initValue("CoordinatePointY", lCoordinatePointY) ;

    std::vector<double> lPointElevation ;
    lPointElevation += geotop::input::gDoubleNoValue ;
    initValue("PointElevation", lPointElevation) ;

    std::vector<double> lPointLandCoverType ;
    lPointLandCoverType += geotop::input::gDoubleNoValue ;
    initValue("PointLandCoverType", lPointLandCoverType) ;

    std::vector<double> lPointSoilType ;
    lPointSoilType += geotop::input::gDoubleNoValue ;
    initValue("PointSoilType", lPointSoilType) ;

    std::vector<double> lPointSlope ;
    lPointSlope += geotop::input::gDoubleNoValue ;
    initValue("PointSlope", lPointSlope) ;

    std::vector<double> lPointAspect ;
    lPointAspect += geotop::input::gDoubleNoValue ;
    initValue("PointAspect", lPointAspect) ;

    std::vector<double> lPointSkyViewFactor ;
    lPointSkyViewFactor += geotop::input::gDoubleNoValue ;
    initValue("PointSkyViewFactor", lPointSkyViewFactor) ;

    std::vector<double> lPointCurvatureNorthSouthDirection ;
    lPointCurvatureNorthSouthDirection += geotop::input::gDoubleNoValue ;
    initValue("PointCurvatureNorthSouthDirection", lPointCurvatureNorthSouthDirection) ;

    std::vector<double> lPointCurvatureWestEastDirection ;
    lPointCurvatureWestEastDirection += geotop::input::gDoubleNoValue ;
    initValue("PointCurvatureWestEastDirection", lPointCurvatureWestEastDirection) ;

    std::vector<double> lPointCurvatureNorthwestSoutheastDirection ;
    lPointCurvatureNorthwestSoutheastDirection += geotop::input::gDoubleNoValue ;
    initValue("PointCurvatureNorthwestSoutheastDirection", lPointCurvatureNorthwestSoutheastDirection) ;

    std::vector<double> lPointCurvatureNortheastSouthwestDirection ;
    lPointCurvatureNortheastSouthwestDirection += geotop::input::gDoubleNoValue ;
    initValue("PointCurvatureNortheastSouthwestDirection", lPointCurvatureNortheastSouthwestDirection) ;

    std::vector<double> lPointDepthFreeSurface ;
    lPointDepthFreeSurface += geotop::input::gDoubleNoValue ;
    initValue("PointDepthFreeSurface", lPointDepthFreeSurface) ;

    std::vector<double> lPointHorizon ;
    lPointHorizon += geotop::input::gDoubleNoValue ;
    initValue("PointHorizon", lPointHorizon) ;

    std::vector<double> lPointMaxSWE ;
    lPointMaxSWE += geotop::input::gDoubleNoValue ;
    initValue("PointMaxSWE", lPointMaxSWE) ;

    std::vector<double> lPointLatitude ;
    lPointLatitude += geotop::input::gDoubleNoValue ;
    initValue("PointLatitude", lPointLatitude) ;

    std::vector<double> lPointLongitude ;
    lPointLongitude += geotop::input::gDoubleNoValue ;
    initValue("PointLongitude", lPointLongitude) ;
    
    std::vector<double> lPointBedrock ;
    lPointBedrock += geotop::input::gDoubleNoValue ;
    initValue("PointBedrock", lPointBedrock) ;

    std::vector<double> lVegReflectVis ;
    lVegReflectVis += 0.2 ;
    initValue("VegReflectVis", lVegReflectVis) ;
    
    initValue("IrriducibleWatSatGlacier", double(0.02)) ;
    
    initValue("GWEPoint", double(-1)) ;
    
    initValue("FetchUp", double(1000)) ;
    
    initValue("BlowingSnowSoftLayerIceContent", double(0)) ;
    
    initValue("TvegPoint", double(-1)) ;
    
    initValue("ThresTempSnow", double(0)) ;
    
    initValue("GWEtop", double(0)) ;
    
    initValue("HgUnvegPoint", double(-1)) ;
    
    initValue("TimeFromStartPoint", double(-1)) ;
    
    initValue("CurvatureWeight", double(0)) ;

    initValue("MaxWaterEqSnowLayerContent", double(5)) ;
    
    initValue("SoilAll", double(0)) ;
    
    initValue("LocMaxIter", double(3)) ;
    
    initValue("SWNetPoint", double(-1)) ;
    
    initValue("LWupPoint", double(-1)) ;
    
    initValue("TimeFromStartGlac", double(-1)) ;
    
    initValue("FrozenSoilHydrCondReduction", double(7)) ;
    
    initValue("ConsiderAlbedoInSWin", double(0)) ;
    
    initValue("LEBasin", double(-1)) ;
    
    initValue("DewTemperatureAsRH", double(1)) ;
    
    initValue("GlacSublPoint", double(-1)) ;
    
    initValue("SnowSublPoint", double(-1)) ;
    
    initValue("LEPoint", double(-1)) ;
    
    initValue("MassErrorBasin", double(-1)) ;
    
    initValue("FetchDown", double(100)) ;
    
    initValue("JulianDayFromYear0Basin", double(-1)) ;
    
    initValue("MinPrecToRestoreFreshSnowAlbedo", double(10)) ;
    
    std::vector<double> lOutputMeteoMaps ;
    lOutputMeteoMaps += 0 ;
    initValue("OutputMeteoMaps", lOutputMeteoMaps) ;
    
    initValue("Vmin", double(0.5)) ;
    
    initValue("SWESublBlownPoint", double(-1)) ;
    
    initValue("InitGlacierDensity", double(800)) ;
    
    initValue("MinDiffSupWaterDepthChannelLand", double(1)) ;
    
    initValue("CalculateCastShadow", double(1)) ;
    
    initValue("GWEbottom", double(0)) ;
    
    initValue("MinRatioKactualToKSat", double(0)) ;
    
    initValue("InitSnowAge", double(0)) ;
    
    std::vector<double> lMeteoStationLongitude ;
    lMeteoStationLongitude += 0 ;
    initValue("MeteoStationLongitude", lMeteoStationLongitude) ;

    
    initValue("AirTempPoint", double(-1)) ;
    
    std::vector<double> lWiltingPoint ;
    lWiltingPoint += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("WiltingPoint", lWiltingPoint) ;
    
    initValue("AngstromBeta", double(0.1)) ;
    
    initValue("SnowAgingCoeffVis", double(0.2)) ;
    
    initValue("IrriducibleWatSatSnow", double(0.02)) ;
    
    initValue("DatePoint", double(-1)) ;
    
    std::vector<double> lThetaSatBedrock ;
    lThetaSatBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("ThetaSatBedrock", lThetaSatBedrock) ;
    
    initValue("BaseWindSpeed", double(0.5)) ;
    
    initValue("InitSWE", double(0)) ;
    
    initValue("SoilLayerTypes", double(1)) ;
    
    initValue("SlopeWeightD", double(0)) ;
    
    initValue("SlopeWeightI", double(0)) ;
    
    std::vector<double> lInitSoilPressureBedrock ;
    lInitSoilPressureBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("InitSoilPressureBedrock", lInitSoilPressureBedrock) ;
    
    std::vector<double> lRoughElemXUnitArea ;
    lRoughElemXUnitArea += 0 ;
    initValue("RoughElemXUnitArea", lRoughElemXUnitArea) ;
    
    std::vector<double> lThermalCapacitySoilSolidsBedrock ;
    lThermalCapacitySoilSolidsBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("ThermalCapacitySoilSolidsBedrock", lThermalCapacitySoilSolidsBedrock) ;
    
    initValue("LWvBasin", double(-1)) ;
    
    initValue("MinLambdaWater", double(1e-07)) ;
    
    initValue("ConvectiveHeatTransferCoefficient", double(geotop::input::gDoubleNoValue)) ;
    
    initValue("BusingerMaxIter", double(5)) ;
    
    initValue("SnowViscosity", double(1e+06)) ;
    
    initValue("SWdiffPoint", double(-1)) ;
    
    initValue("SimulationHours", double(1)) ;
    
    initValue("MinDiffSupWaterDepthLandChannel", double(1)) ;
    
    initValue("HvPoint", double(-1)) ;
    
    initValue("SnowSMIN", double(30)) ;
    
    initValue("PRainBasin", double(-1)) ;
    
    initValue("RichardInitForc", double(0.01)) ;
    
    initValue("ActualOrProjectedArea", double(0)) ;
    
    std::vector<double> lMinStomatalRes ;
    lMinStomatalRes += 60;
    initValue("MinStomatalRes", lMinStomatalRes) ;
    
    std::vector<double> lLateralHydrConductivity ;
    lLateralHydrConductivity += 0.0001;
    initValue("LateralHydrConductivity", lLateralHydrConductivity) ;
    
    initValue("RunGlac", double(-1)) ;
    
    initValue("SnowAgingCoeffNIR", double(0.5)) ;
    
    initValue("PsnowNetPoint", double(-1)) ;
    
    std::vector<double> lPointID ;
    lPointID += geotop::input::gDoubleNoValue;
    initValue("PointID", lPointID) ;
    
    initValue("SWupPoint", double(-1)) ;
    
    std::vector<double> lDecayCoeffCanopy ;
    lDecayCoeffCanopy += 2.5;
    initValue("DecayCoeffCanopy", lDecayCoeffCanopy) ;
    
    initValue("HighestNodeCorrespondsToLayer", double(0)) ;
    
    initValue("TsurfPoint", double(-1)) ;
    
    std::vector<double> lNVanGenuchten ;
    lNVanGenuchten += 1.3 ;
    initValue("NVanGenuchten", lNVanGenuchten) ;
    
    initValue("PNetBasin", double(-1)) ;
    
    initValue("AlphaSnow", double(100000)) ;
    
    initValue("Latitude", double(46.3)) ;
    
    std::vector<double> lRoughElemDiam ;
    lRoughElemDiam += 50 ;
    initValue("RoughElemDiam", lRoughElemDiam) ;
    
    initValue("LWinParameterization", double(9)) ;
    
    initValue("LEgVegPoint", double(-1)) ;
    
    initValue("SnowDensityCutoff", double(100)) ;
    
    initValue("DefaultSoilTypeBedrock", double(1)) ;
    
    initValue("PrainOnSnowPoint", double(-1)) ;
    
    initValue("ConsiderMicrometeorology", double(1)) ;
    
    std::vector<double> lLateralHydrConductivityBedrock ;
    lLateralHydrConductivityBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("LateralHydrConductivityBedrock", lLateralHydrConductivityBedrock) ;
    
    std::vector<double> lLapseRateDewTemp ;
    lLapseRateDewTemp += geotop::input::gDoubleNoValue ;
    initValue("LapseRateDewTemp", lLapseRateDewTemp) ;
    
    initValue("SWEbottom", double(20)) ;
    
    std::vector<double> lSoilAlbNIRWet ;
    lSoilAlbNIRWet += 0.16;
    initValue("SoilAlbNIRWet", lSoilAlbNIRWet) ;
    
    initValue("JulianDayFromYear0Point", double(-1)) ;
    
    initValue("FreeDrainageAtBottom", double(0)) ;
    
    initValue("MinTimeStepSupFlow", double(0.01)) ;
    
    std::vector<double> lMeteoStationSkyViewFactor ;
    lMeteoStationSkyViewFactor += 1 ;
    initValue("MeteoStationSkyViewFactor", lMeteoStationSkyViewFactor) ;
    
    initValue("PSnowNetBasin", double(-1)) ;
    
    initValue("MinLambdaEnergy", double(1e-05)) ;
    
    initValue("SnowOnCanopyPoint", double(-1)) ;
    
    initValue("MaxTimesMinLambdaEnergy", double(0)) ;
    
    std::vector<double> lLapseRateTemp ;
    lLapseRateTemp += 6.5 ;
    initValue("LapseRateTemp", lLapseRateTemp) ;
    
    initValue("AngstromAlpha", double(1.3)) ;
    
    initValue("FormatOutputMaps", double(3)) ;
    
    initValue("HgVegPoint", double(-1)) ;
    
    initValue("SWvPoint", double(-1)) ;
    
    std::vector<double> lInitSoilTempBedrock ;
    lInitSoilTempBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ;
    initValue("InitSoilTempBedrock", lInitSoilTempBedrock) ;
    
    std::vector<double> lSoilAlbVisDry ;
    lSoilAlbVisDry += 0.16;
    initValue("SoilAlbVisDry", lSoilAlbVisDry) ;
    
    initValue("LEupPoint", double(-1)) ;
    
    std::vector<double> lRootDepth ;
    lRootDepth += 300 ;
    initValue("RootDepth", lRootDepth) ;
    
    initValue("RecoverSim", double(0)) ;
    
    initValue("SurfaceEBPoint", double(-1)) ;
    
    initValue("PeriodGlac", double(-1)) ;
    
    initValue("LowestThawedSoilDepthPoint", double(-1)) ;
    
    initValue("RHAsDewTemperature", double(0)) ;
    
    initValue("BaseAirTemperature", double(5)) ;
    
    std::vector<double> lThermalCapacitySoilSolids ;
    lThermalCapacitySoilSolids += 2.3e+06,2.3e+06,2.3e+06,2.3e+06,2.3e+06 ;
    initValue("ThermalCapacitySoilSolids", lThermalCapacitySoilSolids) ;
    
    initValue("LEvPoint", double(-1)) ;
    
    std::vector<double> lVegHeight ;
    lVegHeight += 1000;
    initValue("VegHeight", lVegHeight) ;
    
    initValue("ChannelDepression", double(500)) ;
    
    initValue("RainCorrFactor", double(1)) ;
    
    initValue("SnowAll", double(1)) ;
    
    std::vector<double> lLapseRatePrec ;
    lLapseRatePrec += -0.2 ;
    initValue("LapseRatePrec", lLapseRatePrec) ;
    
    initValue("InitGlacierTemp", double(-3)) ;
    
    std::vector<double> lMeteoStationTemperatureSensorHeight ;
    lMeteoStationTemperatureSensorHeight += 2 ;
    initValue("MeteoStationTemperatureSensorHeight", lMeteoStationTemperatureSensorHeight) ;
    
    std::vector<double> lVegReflNIR ;
    lVegReflNIR += 0.2 ;
    initValue("VegReflNIR", lVegReflNIR) ; 
    
    initValue("ContinuousRecovery", double(0)) ;
    
    initValue("d0vegPoint", double(-1)) ;
    
    initValue("LObukhovPoint", double(-1)) ;
    
    initValue("WindAsSpeedAndDirection", double(1)) ;
    
    initValue("NumLowPassFilterOnDemForAll", double(0)) ;
    
    initValue("TSurfBasin", double(-1)) ;
    
    initValue("MinSupWaterDepthChannel", double(1)) ;
    
    initValue("WindSpeedPoint", double(-1)) ;
    
    initValue("RHPoint", double(-1)) ;
    
    initValue("QCanopyAirPoint", double(-1)) ;
    
    initValue("RatioChannelWidthPixelWidth", double(0.1)) ;
    
    std::vector<double> lThermalConductivitySoilSolidsBedrock ;
    lThermalConductivitySoilSolidsBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ; 
    initValue("ThermalConductivitySoilSolidsBedrock", lThermalConductivitySoilSolidsBedrock) ; 
    
    initValue("LSAIPoint", double(-1)) ;
    
    initValue("FreshSnowReflVis", double(0.95)) ;
    
    std::vector<double> lOutputSurfEBALMaps ;
    lOutputSurfEBALMaps += 0 ;
    initValue("OutputSurfEBALMaps", lOutputSurfEBALMaps) ;
    
    initValue("PointAll", double(0)) ;
    
    initValue("DewTempOrNormTemp", double(0)) ;
    
    initValue("PeriodBasin", double(-1)) ;
    
    std::vector<double> lFieldCapacityBedrock ;
    lFieldCapacityBedrock += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ; 
    initValue("FieldCapacityBedrock", lFieldCapacityBedrock) ; 
    
    initValue("RunPoint", double(-1)) ;
    
    initValue("MaxCourantSupFlowLand", double(0.1)) ;
    
    initValue("LowestWaterTableDepthPoint", double(-1)) ;

    std::vector<double> lDtPlotPoint ;
    lDtPlotPoint += 1 ;
    initValue("DtPlotPoint", lDtPlotPoint) ;
    
    std::vector<double> lSpecificStorativity ;
    lSpecificStorativity += 1e-07 ;
    initValue("SpecificStorativity", lSpecificStorativity) ; 
    
    initValue("HeatEqTol", double(1e-06)) ;
    
    initValue("DateGlac", double(-1)) ;
    
    std::vector<double> lSoilAlbVisWet ;
    lSoilAlbVisWet += 0.16 ;
    initValue("SoilAlbVisWet", lSoilAlbVisWet) ; 
    
    initValue("SnowSMAX", double(80)) ;
    
    initValue("SWinBasin", double(-1)) ;
    
    std::vector<double> lInitSoilTemp ;
    lInitSoilTemp += 5 ;
    initValue("InitSoilTemp", lInitSoilTemp) ; 
    
    initValue("TCanopyAirPoint", double(-1)) ;
    
    initValue("DefaultSoilTypeChannel", double(1)) ;
    
    initValue("InitSnowDensity", double(250)) ;
    
    initValue("SnowCURV", double(-200)) ;
    
    initValue("TraspCanopyBasin", double(-1)) ;
    
    initValue("SWEtop", double(50)) ;
    
    initValue("UpdateHydraulicConductivity", double(0)) ;
    
    std::vector<double> lThresSnowSoilRough ;
    lThresSnowSoilRough += 10 ;
    initValue("ThresSnowSoilRough", lThresSnowSoilRough) ; 
    
    initValue("MaxWaterEqGlacLayerContent", double(5)) ;
    
    std::vector<double> lInitSoilPressure ;
    lInitSoilPressure += geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue,geotop::input::gDoubleNoValue ; 
    initValue("InitSoilPressure", lInitSoilPressure) ; 
    
    initValue("SurFlowResExp", double(0.666666666667)) ;
    
    initValue("MinPrecIncreaseFactorWithElev", double(0.1)) ;
    //END INITIALIZATION OF NUMERIC PARAMETERS
}

geotop::input::ConfigStore::~ConfigStore()
{
    
}

boost::shared_ptr<geotop::input::ConfigStore> geotop::input::ConfigStoreSingletonFactory::getInstance() {
    if ( ! mInstance ) {
        mMutex.lock ();
        if ( ! mInstance ) {
            boost::shared_ptr<geotop::input::ConfigStore> lTemp ( new geotop::input::ConfigStore() );
            mInstance = lTemp ;
        }
        mMutex.unlock ();
    }
    return mInstance;
}

boost::shared_ptr<geotop::input::ConfigStore> geotop::input::ConfigStoreSingletonFactory::mInstance ;
boost::signals2::mutex geotop::input::ConfigStoreSingletonFactory::mMutex ;
