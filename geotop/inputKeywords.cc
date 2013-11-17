/**
 * @file   inputKeywords.cc
 * @Author Angelo Leto (angleto@gmail.com)
 * @date   November, 2013
 * @brief  generic configuration store class
 *
 * Parse configuration file and store parameters on single container
 */

#include <inputKeywords.h>
#include <boost/spirit.hpp>
#include <boost/spirit/include/qi_real.hpp>
#include <boost/spirit/include/qi_char.hpp>
#include <boost/spirit/include/qi_eol.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/assign/std/vector.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace boost::assign;

/** @internal
 * @brief configuration file parsing class, this class
 *  must not be used, in the future if more configuration format will be
 *  available, this class should be handled using the decorator pattern.
 */
class ConfGrammar
    : public boost::spirit::grammar<ConfGrammar>
{
public:
    
    ConfGrammar (boost::shared_ptr< std::map<std::string, boost::any> > pMap){
        mMap = pMap ;
        mKey = boost::shared_ptr< std::string >( new std::string() ) ;
    };
    
    ~ConfGrammar (){
        
    };

    void actionKey(char const * const pBegin, char const * const pEnd) const
    {
        std::string lKey = std::string(pBegin, pEnd) ;
        std::cout << "actionKey: " << lKey << std::endl ;
        mKey->assign( lKey ) ;
    };
    
    void actionValueString(char const * const pBegin, char const * const pEnd) const
    {
        std::string lValue = std::string(pBegin, pEnd) ;
        
        std::cout << "actionValueString: " << lValue << std::endl;
        
        if( mKey->compare ("") != 0 )
        {
            (*mMap)[*mKey] = lValue;
        } else {
            std::cerr << "Error: actionValueString : no key was pushed for the value, value will be discarded" << std::endl ;
        }
        mKey->assign( "" ) ;
    };
    
    void actionValueDate(char const * const pBegin, char const * const pEnd) const
    {
        std::string lValue = std::string(pBegin, pEnd) ;
     
        std::cout << "actionValueDate: " << lValue << std::endl;
        
        if( mKey->compare ("") != 0 )
        {
            (*mMap)[*mKey] = lValue;
        } else {
            std::cerr << "Error: actionValueDate : no key was pushed for the value, value will be discarded" << std::endl ;
        }
        mKey->assign( "" ) ;
    };
    
    void actionValueDouble(const double pValue) const
    {
        std::cout << "actionValueDouble: " << pValue << std::endl;
        
        if( mKey->compare ("") != 0 )
        {
            (*mMap)[*mKey] = pValue;
        } else {
            std::cerr << "Error: actionValueDouble : no key was pushed for the value, value will be discarded" << std::endl ;
        }
        mKey->assign( "" ) ;
    };
    
    void actionValueDoubleArray(char const * const pBegin, char const * const pEnd) const
    {
        std::string lString = std::string(pBegin, pEnd) ;
        boost::algorithm::erase_all(lString, " ");
        
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
            std::cerr << "Error: actionValueDoubleArray : no key was pushed for the value, value will be discarded" << std::endl ;
        }
        mKey->assign( "" ) ;
        
        std::cout << "actionValueDoubleArray Size: " << lDoubleArray.size() << std::endl;
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
            row = blanks >> (comment | parameter) ;
            blanks = *(boost::spirit::blank_p) ;
            parameter = key[boost::bind(&ConfGrammar::actionKey, self, _1, _2)] >> blanks >> "=" >>
            blanks >> value ;
            key = +(boost::spirit::alpha_p|boost::spirit::alnum_p) ;
            value = date | array | num | string ;
            num = (boost::spirit::real_p)[boost::bind(&ConfGrammar::actionValueDouble, self, _1)] ;
            array = (boost::spirit::real_p >> +(blanks >> "," >>
                                                blanks >> boost::spirit::real_p))[boost::bind(&ConfGrammar::actionValueDoubleArray, self, _1, _2)] ;
            string = "\"" >> (*(~(boost::spirit::ch_p("\""))))[boost::bind(&ConfGrammar::actionValueString, self, _1, _2)] >> "\"" ;
            comment = boost::spirit::comment_p("!");
            date = (DD >> "/" >> MM >> "/" >> YYYY >>
                    blanks >> hh >> ":" >> mm)[boost::bind(&ConfGrammar::actionValueDate, self, _1, _2)] ;
            DD = boost::spirit::digit_p >> boost::spirit::digit_p ;
            MM = boost::spirit::digit_p >> boost::spirit::digit_p ;
            YYYY = boost::spirit::digit_p >> boost::spirit::digit_p >>
            boost::spirit::digit_p >> boost::spirit::digit_p ;
            hh = boost::spirit::digit_p >> boost::spirit::digit_p ;
            mm = boost::spirit::digit_p >> boost::spirit::digit_p ;
        }
        
        const boost::spirit::rule<Scanner> &start()
        {
            return configfile;
        }
        
    private:
        
        boost::spirit::rule<Scanner> configfile,
        row,
        parameter,
        key,
        value,
        array,
        string,
        num,
        blanks,
        date, DD, MM, YYYY, hh, mm,
        comment ;
    };
    
private:
    
    boost::shared_ptr<std::string> mKey ;
    boost::shared_ptr< std::map<std::string, boost::any> > mMap ;
};

geotop::input::ConfigStore::ConfigStore()
{
    mValueMap = boost::shared_ptr<std::map<std::string, boost::any> > (new std::map<std::string, boost::any> () ) ;
}

bool geotop::input::ConfigStore::parse(const std::string pFileName) {
    std::ifstream fs(pFileName);
    std::ostringstream ss;
    ss << fs.rdbuf();
    
    init() ;
    
    ConfGrammar lConfGrammar (mValueMap) ;
    boost::spirit::parse_info<> lParserInfo = boost::spirit::parse(ss.str().c_str(), lConfGrammar, boost::spirit::space_p);
    if (lParserInfo.hit)
    {
        std::cout << "Info: configuration file parsing completed" << std::endl;
        std::cout << lParserInfo.length << " characters parsed" << std::endl;
    }
    else
    {
        std::cout << "parsing failed: stopped at '" << lParserInfo.stop << "'" << std::endl;
    }
}

/** @internal
 * @brief initialization of the configuration parameters
 */
void geotop::input::ConfigStore::init()
{
    mValueMap->clear() ;
    
    //BEGIN INITIALIZATION OF STRING PARAMETERS
    initValue("SpecificPlotSurfaceSensibleHeatFluxMapFile", std::string("none")) ;
    
    initValue("RecoverGlacierLayerThick", std::string("none")) ;
    
    initValue("HeaderMeteoStationLongitude", std::string("none")) ;
    
    initValue("NetLongwaveRadiationMapFile", std::string("none")) ;
    
    initValue("GlacierProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderLEBasin", std::string("none")) ;
    
    initValue("HeaderMeteoStationElevation", std::string("none")) ;
    
    initValue("RecoverSnowOnCanopy", std::string("none")) ;
    
    initValue("SoilAveragedTempProfileFile", std::string("none")) ;
    
    initValue("HeaderSWbeamPoint", std::string("none")) ;
    
    initValue("Headerz0vegPoint", std::string("none")) ;
    
    initValue("HeaderSWglobal", std::string("Swglob")) ;
    
    initValue("HeaderLWinPoint", std::string("none")) ;
    
    initValue("HeaderLWvBasin", std::string("none")) ;
    
    initValue("HeaderWatContentSnow", std::string("none")) ;
    
    initValue("HeaderPointLandCoverType", std::string("none")) ;
    
    initValue("HeaderLEvPoint", std::string("none")) ;
    
    initValue("HeaderLEvBasin", std::string("none")) ;
    
    initValue("SubfolderRecoveryFiles", std::string("rec")) ;
    
    initValue("RecoverRunSoilMaximumTemperatureFile", std::string("none")) ;
    
    initValue("SoilAveragedLiqContentProfileFile", std::string("none")) ;
    
    initValue("SuccessfulRunFile", std::string("none")) ;
    
    initValue("RecoverVegTemp", std::string("none")) ;
    
    initValue("HeaderLEPoint", std::string("none")) ;
    
    initValue("HeaderCoordinatePointY", std::string("Ywgs")) ;
    
    initValue("HeaderCoordinatePointX", std::string("Xwgs")) ;
    
    initValue("HeaderDateGlac", std::string("none")) ;
    
    initValue("HeaderPeriodGlac", std::string("none")) ;
    
    initValue("SpecificPlotVegTempMapFile", std::string("none")) ;
    
    initValue("HeaderN", std::string("n")) ;
    
    initValue("HeaderV", std::string("none")) ;
    
    initValue("NetRadiationMapFile", std::string("output_maps/RadNet")) ;
    
    initValue("PointOutputFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderCanopyFractionPoint", std::string("none")) ;
    
    initValue("HeaderRunSnow", std::string("none")) ;
    
    initValue("HeaderLapseRateTemp", std::string("none")) ;
    
    initValue("HeaderHorizonAngle", std::string("azimuth")) ;
    
    initValue("HeaderLWNetPoint", std::string("none")) ;
    
    initValue("HeaderMeteoStationLatitude", std::string("none")) ;
    
    initValue("HeaderIceContentSnow", std::string("none")) ;
    
    initValue("HeaderSoilDz", std::string("Dz")) ;
    
    initValue("EvapotranspirationFromSoilMapFile", std::string("none")) ;
    
    initValue("HeaderSWdirect", std::string("none")) ;
    
    initValue("AspectMapFile", std::string("input_maps/aspect")) ;
    
    initValue("GlacierMeltedMapFile", std::string("none")) ;
    
    initValue("HeaderTimeFromStartSoil", std::string("none")) ;
    
    initValue("SnowCoveredAreaFile", std::string("none")) ;
    
    initValue("HeaderSWESublBlownPoint", std::string("none")) ;
    
    initValue("SpecificPlotSurfaceLatentHeatFluxMapFile", std::string("none")) ;
    
    initValue("SoilTempProfileFile", std::string("output_tabs/soiltemp")) ;
    
    initValue("HeaderSnowTempPoint", std::string("none")) ;
    
    initValue("RecoverSoilTemp", std::string("none")) ;
    
    initValue("NetPrecipitationMapFile", std::string("none")) ;
    
    initValue("SurfaceHeatFluxMapFile", std::string("none")) ;
    
    initValue("HeaderHPoint", std::string("none")) ;
    
    initValue("HeaderPeriodPoint", std::string("none")) ;
    
    initValue("SoilAveragedIceContentTensorFile", std::string("none")) ;
    
    initValue("HeaderDateDDMMYYYYhhmmLapseRates", std::string("none")) ;
    
    initValue("AirTempMapFile", std::string("none")) ;
    
    initValue("HeaderPointLatitude", std::string("none")) ;
    
    initValue("HeaderTSurfBasin", std::string("none")) ;
    
    initValue("MeteoStationsListFile", std::string("none")) ;
    
    initValue("SpecificPlotSurfaceWaterContentMapFile", std::string("none")) ;
    
    initValue("CanopyInterceptedWaterMapFile", std::string("none")) ;
    
    initValue("HeaderLapseRateDewTemp", std::string("none")) ;
    
    initValue("HeaderSurfaceTemperature", std::string("none")) ;
    
    initValue("SoilLiqWaterPressProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderMeteoStationCoordinateY", std::string("none")) ;
    
    initValue("RecoverRunSoilAveragedSnowWaterEquivalentFile", std::string("none")) ;
    
    initValue("RecoverSnowLiqMass", std::string("none")) ;
    
    initValue("HeaderLEgUnvegPoint", std::string("none")) ;
    
    initValue("SoilAveragedIceContentProfileFile", std::string("none")) ;
    
    initValue("HeaderLWupPoint", std::string("none")) ;
    
    initValue("InitWaterTableDepthMapFile", std::string("none")) ;
    
    initValue("HeaderTvegBasin", std::string("none")) ;
    
    initValue("LandCoverMapFile", std::string("input_maps/landcover")) ;
    
    initValue("LapseRateFile", std::string("none")) ;
    
    initValue("RecoverRunSoilMinimumTotalSoilMoistureFile", std::string("none")) ;
    
    initValue("WindSpeedMapFile", std::string("none")) ;
    
    initValue("HeaderJulianDayFromYear0Glac", std::string("none")) ;
    
    initValue("RecoverSnowTemp", std::string("none")) ;
    
    initValue("InLongwaveRadiationMapFile", std::string("output_maps/LWin")) ;
    
    initValue("GlacierProfileFile", std::string("none")) ;
    
    initValue("HeaderAirPressPoint", std::string("none")) ;
    
    initValue("HeaderLapseRatePrec", std::string("none")) ;
    
    initValue("HeaderPointCurvatureWestEastDirection", std::string("none")) ;
    
    initValue("HeaderSWinBasin", std::string("none")) ;
    
    initValue("SkyViewFactorMapFile", std::string("input_maps/sky")) ;
    
    initValue("HeaderIDPointPoint", std::string("none")) ;
    
    initValue("SpecificPlotIncomingShortwaveRadMapFile", std::string("none")) ;
    
    initValue("SnowLiqContentProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderEvapSurfaceBasin", std::string("none")) ;
    
    initValue("HeaderLEupPoint", std::string("none")) ;
    
    initValue("HeaderTraspCanopyPoint", std::string("none")) ;
    
    initValue("SurfaceSensibleHeatFluxMapFile", std::string("none")) ;
    
    initValue("HeaderDepthSnow", std::string("none")) ;
    
    initValue("HeaderSWvBasin", std::string("none")) ;
    
    initValue("HeaderSWvPoint", std::string("none")) ;
    
    initValue("HeaderGlacDensityPoint", std::string("none")) ;
    
    initValue("HeaderPeriodBasin", std::string("none")) ;
    
    initValue("HeaderPRainBasin", std::string("none")) ;
    
    initValue("HeaderWindVelocity", std::string("WindSp")) ;
    
    initValue("LandSurfaceWaterDepthMapFile", std::string("none")) ;
    
    initValue("PointOutputFile", std::string("output_tabs/point")) ;
    
    initValue("HeaderSurfaceEBPoint", std::string("none")) ;
    
    initValue("SoilSaturationRatioProfileFile", std::string("none")) ;
    
    initValue("HeaderDecayKCanopyPoint", std::string("none")) ;
    
    initValue("SoilTempTensorFile", std::string("none")) ;
    
    initValue("HeaderTCanopyAirPoint", std::string("none")) ;
    
    initValue("HeaderPointSoilType", std::string("none")) ;
    
    initValue("HeaderTempGlac", std::string("none")) ;
    
    initValue("SpecificPlotVegSensibleHeatFluxMapFile", std::string("none")) ;
    
    initValue("SpecificPlotVegLatentHeatFluxMapFile", std::string("none")) ;
    
    initValue("RecoverGlacierTemp", std::string("none")) ;
    
    initValue("HeaderSWupPoint", std::string("none")) ;
    
    initValue("HeaderLateralHydrConductivity", std::string("Kh")) ;
    
    initValue("HeaderSoilInitTemp", std::string("none")) ;
    
    initValue("DemFile", std::string("input_maps/pit")) ;
    
    initValue("RecoverRunSoilMinimumTemperatureFile", std::string("none")) ;
    
    initValue("SpecificPlotCanopyAirTempMapFile", std::string("none")) ;
    
    initValue("HeaderCloudSWTransmissivity", std::string("CloudTrans")) ;
    
    initValue("SoilAveragedLiqContentTensorFile", std::string("none")) ;
    
    initValue("SoilTotWaterPressProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderSWNetBasin", std::string("none")) ;
    
    initValue("HeaderNormalHydrConductivity", std::string("Kv")) ;
    
    initValue("HeaderWindDirection", std::string("WindDir")) ;
    
    initValue("HeaderDateSnow", std::string("none")) ;
    
    initValue("HeaderIPrec", std::string("Iprec")) ;
    
    initValue("InitSnowDepthMapFile", std::string("none")) ;
    
    initValue("HeaderLObukhovPoint", std::string("none")) ;
    
    initValue("InitGlacierDepthMapFile", std::string("none")) ;
    
    initValue("HeaderMeanTimeStep", std::string("none")) ;
    
    initValue("SoilLiqWaterPressProfileFile", std::string("none")) ;
    
    initValue("SpecificPlotTotalSensibleHeatFluxMapFile", std::string("none")) ;
    
    initValue("HeaderJulianDayfrom0Meteo", std::string("none")) ;
    
    initValue("HeaderCloudFactor", std::string("none")) ;
    
    initValue("HeaderPointMaxSWE", std::string("none")) ;
    
    initValue("WaterTableDepthFromAboveMapFile", std::string("none")) ;
    
    initValue("HeaderSWEPoint", std::string("none")) ;
    
    initValue("HeaderPSnowBasin", std::string("none")) ;
    
    initValue("SpecificPlotWindDirMapFile", std::string("none")) ;
    
    initValue("RecoverRunSoilAveragedTemperatureFile", std::string("none")) ;
    
    initValue("HeaderPointLongitude", std::string("none")) ;
    
    initValue("HeaderDewTemp", std::string("none")) ;
    
    initValue("BedrockDepthMapFile", std::string("none")) ;
    
    initValue("RunSoilMinimumTotalSoilMoistureFile", std::string("none")) ;
    
    initValue("HeaderPeriodSnow", std::string("none")) ;
    
    initValue("HeaderSnowMeltedPoint", std::string("none")) ;
    
    initValue("RecoverSnowIceMass", std::string("none")) ;
    
    initValue("HeaderJulianDayFromYear0Snow", std::string("none")) ;
    
    initValue("InitSnowAgeMapFile", std::string("none")) ;
    
    initValue("HeaderWindSpeedTopCanopyPoint", std::string("none")) ;
    
    initValue("SoilLiqWaterPressTensorFile", std::string("output_maps/psiz")) ;
    
    initValue("RecoverGlacierLayerNumber", std::string("none")) ;
    
    initValue("HeaderThetaRes", std::string("res")) ;
    
    initValue("HorizonMeteoStationFile", std::string("hor_meteo/horizon")) ;
    
    initValue("HeaderLWinBasin", std::string("none")) ;
    
    initValue("RecoverRunSoilAveragedInternalEnergyFile", std::string("none")) ;
    
    initValue("HeaderEstoredCanopyPoint", std::string("none")) ;
    
    initValue("HeaderSWEBlownPoint", std::string("none")) ;
    
    initValue("HeaderSnowSublPoint", std::string("none")) ;
    
    initValue("InitSWEMapFile", std::string("none")) ;
    
    initValue("RunSoilMaximumTemperatureFile", std::string("none")) ;
    
    initValue("HeaderPrainPoint", std::string("none")) ;
    
    initValue("HeaderGlacDepthPoint", std::string("none")) ;
    
    initValue("DaysDelayMapFile", std::string("none")) ;
    
    initValue("TimeDependentIncomingDischargeFile", std::string("none")) ;
    
    initValue("HeaderSWdiffuse", std::string("none")) ;
    
    initValue("FirstSoilLayerIceContentMapFile", std::string("none")) ;
    
    initValue("FailedRunFile", std::string("none")) ;
    
    initValue("SoilIceContentProfileFileWriteEnd", std::string("none")) ;
    
    initValue("SpecificPlotWindSpeedMapFile", std::string("none")) ;
    
    initValue("HeaderLWinMaxPoint", std::string("none")) ;
    
    initValue("HeaderAirTempBasin", std::string("none")) ;
    
    initValue("InShortwaveRadiationMapFile", std::string("output_maps/SWin")) ;
    
    initValue("HeaderHighestWaterTableDepthPoint", std::string("none")) ;
    
    initValue("HeaderTimeFromStartSnow", std::string("none")) ;
    
    initValue("ChannelSurfaceWaterDepthMapFile", std::string("none")) ;
    
    initValue("HeaderPsnowPoint", std::string("none")) ;
    
    initValue("HeaderCthSoilSolids", std::string("none")) ;
    
    initValue("HeaderHighestThawedSoilDepthPoint", std::string("none")) ;
    
    initValue("SnowIceContentProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderPSnowNetBasin", std::string("none")) ;
    
    initValue("HeaderSoilInitPres", std::string("none")) ;
    
    initValue("HeaderRunPoint", std::string("none")) ;
    
    initValue("SWEMapFile", std::string("output_maps/SWE")) ;
    
    initValue("HeaderHupPoint", std::string("none")) ;
    
    initValue("WaterTableDepthMapFile", std::string("none")) ;
    
    initValue("HeaderGlacMeltedPoint", std::string("none")) ;
    
    initValue("HeaderPointSlope", std::string("none")) ;
    
    initValue("HeaderTempSnow", std::string("none")) ;
    
    initValue("WindDirMapFile", std::string("none")) ;
    
    initValue("HeaderTimeStepAverage", std::string("none")) ;
    
    initValue("PrecipitationMapFile", std::string("none")) ;
    
    initValue("HeaderLEgVegPoint", std::string("none")) ;
    
    initValue("HeaderWiltingPoint", std::string("none")) ;
    
    initValue("RecoverSoilIceContChannel", std::string("none")) ;
    
    initValue("HeaderSWinPoint", std::string("none")) ;
    
    initValue("HeaderDateDDMMYYYYhhmmMeteo", std::string("Date")) ;
    
    initValue("DirectInShortwaveRadiationMapFile", std::string("none")) ;
    
    initValue("HeaderWindSpeedPoint", std::string("none")) ;
    
    initValue("PointFile", std::string("ListPoints")) ;
    
    initValue("RecoverTime", std::string("none")) ;
    
    initValue("SpecificPlotSnowDepthMapFile", std::string("none")) ;
    
    initValue("RunSoilAveragedTotalSoilMoistureFile", std::string("none")) ;
    
    initValue("SnowDurationMapFile", std::string("none")) ;
    
    initValue("SurfaceTempMapFile", std::string("none")) ;
    
    initValue("HeaderLSAIPoint", std::string("none")) ;
    
    initValue("RunSoilAveragedInternalEnergyFile", std::string("none")) ;
    
    initValue("HeaderQSurfPoint", std::string("none")) ;
    
    initValue("HeaderWindDirPoint", std::string("none")) ;
    
    initValue("HeaderPrainOnSnowPoint", std::string("none")) ;
    
    initValue("SoilIceContentProfileFile", std::string("none")) ;
    
    initValue("SoilTotWaterPressTensorFile", std::string("none")) ;
    
    initValue("HeaderWaterOnCanopyPoint", std::string("none")) ;
    
    initValue("SoilLiqContentProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderAirTempPoint", std::string("none")) ;
    
    initValue("SpecificPlotNetVegShortwaveRadMapFile", std::string("none")) ;
    
    initValue("SoilTempProfileFileWriteEnd", std::string("none")) ;
    
    initValue("SpecificPlotNetSurfaceShortwaveRadMapFile", std::string("none")) ;
    
    initValue("RecoverSoilWatPres", std::string("none")) ;
    
    initValue("SpecificPlotSurfaceHeatFluxMapFile", std::string("none")) ;
    
    initValue("HeaderLWNetBasin", std::string("none")) ;
    
    initValue("HeaderPointID", std::string("ID")) ;
    
    initValue("SnowTempProfileFile", std::string("output_tabs/snowtemp")) ;
    
    initValue("SpecificPlotSurfaceTempMapFile", std::string("none")) ;
    
    initValue("HeaderWatContentGlac", std::string("none")) ;
    
    initValue("RunSoilMinimumTemperatureFile", std::string("none")) ;
    
    initValue("SoilMapFile", std::string("input_maps/soiltype")) ;
    
    initValue("ShadowFractionTimeMapFile", std::string("none")) ;
    
    initValue("GlacierDepthMapFile", std::string("none")) ;
    
    initValue("BasinOutputFile", std::string("none")) ;
    
    initValue("RecoverSoilTempChannel", std::string("none")) ;
    
    initValue("FirstSoilLayerTempMapFile", std::string("none")) ;
    
    initValue("HeaderSnowOnCanopyPoint", std::string("none")) ;
    
    initValue("RecoverSoilWatPresChannel", std::string("none")) ;
    
    initValue("HeaderIDMeteoStation", std::string("none")) ;
    
    initValue("HeaderTDewPoint", std::string("none")) ;
    
    initValue("RecoverSnowLayerNumber", std::string("none")) ;
    
    initValue("SurfaceLatentHeatFluxMapFile", std::string("none")) ;
    
    initValue("HeaderGWEPoint", std::string("none")) ;
    
    initValue("HeaderPointCurvatureNorthwestSoutheastDirection", std::string("none")) ;
    
    initValue("HeaderSWnet", std::string("none")) ;
    
    initValue("SpecificPlotNetVegLongwaveRadMapFile", std::string("none")) ;
    
    initValue("ThawedSoilDepthFromAboveMapFile", std::string("none")) ;
    
    initValue("SnowLiqContentProfileFile", std::string("output_tabs/snowliq")) ;
    
    initValue("HeaderHgUnvegPoint", std::string("none")) ;
    
    initValue("SnowDepthLayersFile", std::string("output_tabs/snowly")) ;
    
    initValue("SpecificPlotIncomingLongwaveRadMapFile", std::string("none")) ;
    
    initValue("HeaderHBasin", std::string("none")) ;
    
    initValue("HeaderRH", std::string("RH")) ;
    
    initValue("SoilAveragedTempTensorFile", std::string("output_maps/T")) ;
    
    initValue("SoilTotWaterPressProfileFile", std::string("none")) ;
    
    initValue("GlacierWaterEqMapFile", std::string("none")) ;
    
    initValue("FirstSoilLayerLiqContentMapFile", std::string("none")) ;
    
    initValue("HeaderPointHorizon", std::string("none")) ;
    
    initValue("RecoverGlacierLiqMass", std::string("none")) ;
    
    initValue("SlopeMapFile", std::string("input_maps/slope")) ;
    
    initValue("HeaderMeteoStationSkyViewFactor", std::string("none")) ;
    
    initValue("RelHumMapFile", std::string("none")) ;
    
    initValue("GlacierSublimatedMapFile", std::string("none")) ;
    
    initValue("SpecificPlotAboveVegAirTempMapFile", std::string("none")) ;
    
    initValue("TimeDependentVegetationParameterFile", std::string("none")) ;
    
    initValue("RunSoilMaximumTotalSoilMoistureFile", std::string("none")) ;
    
    initValue("HeaderMeteoStationStandardTime", std::string("none")) ;
    
    initValue("HeaderPrainNetPoint", std::string("none")) ;
    
    initValue("HeaderLowestThawedSoilDepthPoint", std::string("none")) ;
    
    initValue("HeaderTimeFromStartBasin", std::string("none")) ;
    
    initValue("RecoverRunSoilAveragedTotalSoilMoistureFile", std::string("none")) ;
    
    initValue("HeaderSnowDensityPoint", std::string("none")) ;
    
    initValue("HeaderPrec", std::string("none")) ;
    
    initValue("HeaderQVegPoint", std::string("none")) ;
    
    initValue("HeaderIDPointGlac", std::string("none")) ;
    
    initValue("HeaderTsurfPoint", std::string("none")) ;
    
    initValue("TimeStepsFile", std::string("none")) ;
    
    initValue("ThawedSoilDepthMapFile", std::string("none")) ;
    
    initValue("HeaderJulianDayFromYear0Soil", std::string("none")) ;
    
    initValue("HeaderLObukhovCanopyPoint", std::string("none")) ;
    
    initValue("HeaderJulianDayFromYear0Basin", std::string("none")) ;
    
    initValue("HeaderTimeFromStartPoint", std::string("none")) ;
    
    initValue("RecoverSnowLayerThick", std::string("none")) ;
    
    initValue("HeaderIDPointSnow", std::string("none")) ;
    
    initValue("HeaderPointCurvatureNorthSouthDirection", std::string("none")) ;
    
    initValue("HorizonPointFile", std::string("hor_points/horizon")) ;
    
    initValue("HeaderTimeFromStartGlac", std::string("none")) ;
    
    initValue("SuccessfulRecoveryFile", std::string("none")) ;
    
    initValue("FirstSoilLayerAveragedTempMapFile", std::string("output_maps/MMGST")) ;
    
    initValue("HeaderRunGlac", std::string("none")) ;
    
    initValue("SnowSublMapFile", std::string("none")) ;
    
    initValue("SnowDepthMapFile", std::string("output_maps/snowdepth")) ;
    
    initValue("SpecificPlotRelHumMapFile", std::string("none")) ;
    
    initValue("HeaderPointSkyViewFactor", std::string("none")) ;
    
    initValue("MeteoFile", std::string("meteo/meteo")) ;
    
    initValue("HeaderPointElevation", std::string("none")) ;
    
    initValue("BasinOutputFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderLWin", std::string("none")) ;
    
    initValue("SoilAveragedLiqContentProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderTraspCanopyBasin", std::string("none")) ;
    
    initValue("RunSoilAveragedSnowWaterEquivalentFile", std::string("none")) ;
    
    initValue("HeaderLWinMinPoint", std::string("none")) ;
    
    initValue("HeaderHvPoint", std::string("none")) ;
    
    initValue("HeaderWindX", std::string("none")) ;
    
    initValue("HeaderQCanopyAirPoint", std::string("none")) ;
    
    initValue("HeaderDepthGlac", std::string("none")) ;
    
    initValue("HeaderPointBedrockDepth", std::string("none")) ;
    
    initValue("HeaderFieldCapacity", std::string("fc")) ;
    
    initValue("HeaderQAirPoint", std::string("none")) ;
    
    initValue("SnowIceContentProfileFile", std::string("output_tabs/snowice")) ;
    
    initValue("HeaderRunBasin", std::string("none")) ;
    
    initValue("NetShortwaveRadiationMapFile", std::string("none")) ;
    
    initValue("HeaderPointCurvatureNortheastSouthwestDirection", std::string("none")) ;
    
    initValue("HeaderSoilHeatFluxPoint", std::string("none")) ;
    
    initValue("SpecificPlotNetSurfaceLongwaveRadMapFile", std::string("none")) ;
    
    initValue("HeaderSpecificStorativity", std::string("SS")) ;
    
    initValue("HeaderAlpha", std::string("a")) ;
    
    initValue("HeaderMeteoStationCoordinateX", std::string("none")) ;
    
    initValue("HeaderPointDepthFreeSurface", std::string("none")) ;
    
    initValue("HeaderSnowDepthPoint", std::string("none")) ;
    
    initValue("HeaderIDPointSoil", std::string("none")) ;
    
    initValue("CurvaturesMapFile", std::string("none")) ;
    
    initValue("HeaderJulianDayFromYear0Point", std::string("none")) ;
    
    initValue("HeaderPRainNetBasin", std::string("none")) ;
    
    initValue("RecoverRunSoilMaximumTotalSoilMoistureFile", std::string("none")) ;
    
    initValue("HeaderPsnowNetPoint", std::string("none")) ;
    
    initValue("HeaderSWNetPoint", std::string("none")) ;
    
    initValue("HeaderIceContentGlac", std::string("none")) ;
    
    initValue("HeaderTvegPoint", std::string("none")) ;
    
    initValue("SnowMeltedMapFile", std::string("none")) ;
    
    initValue("HeaderGlacTempPoint", std::string("none")) ;
    
    initValue("HeaderAirTemp", std::string("AirT")) ;
    
    initValue("RecoverLiqWaterOnCanopy", std::string("none")) ;
    
    initValue("DischargeFile", std::string("none")) ;
    
    initValue("RecoverSoilIceCont", std::string("none")) ;
    
    initValue("HeaderHgVegPoint", std::string("none")) ;
    
    initValue("HeaderPeriodSoil", std::string("none")) ;
    
    initValue("SoilLiqContentProfileFile", std::string("none")) ;
    
    initValue("HeaderRunSoil", std::string("none")) ;
    
    initValue("HeaderDateSoil", std::string("none")) ;
    
    initValue("HeaderHvBasin", std::string("none")) ;
    
    initValue("Headerd0vegPoint", std::string("none")) ;
    
    initValue("SpecificPlotTotalLatentHeatFluxMapFile", std::string("none")) ;
    
    initValue("SoilLiqContentTensorFile", std::string("output_maps/thetaliq")) ;
    
    initValue("RecoverGlacierIceMass", std::string("none")) ;
    
    initValue("HeaderThetaSat", std::string("sat")) ;
    
    initValue("SoilAveragedIceContentProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderLWvPoint", std::string("none")) ;
    
    initValue("HeaderKthSoilSolids", std::string("none")) ;
    
    initValue("RecoverNonDimensionalSnowAge", std::string("none")) ;
    
    initValue("RunSoilAveragedTemperatureFile", std::string("none")) ;
    
    initValue("HeaderSWdiffPoint", std::string("none")) ;
    
    initValue("RiverNetwork", std::string("none")) ;
    
    initValue("HeaderWindY", std::string("none")) ;
    
    initValue("HeaderDateBasin", std::string("none")) ;
    
    initValue("HeaderPNetBasin", std::string("none")) ;
    
    initValue("HeaderHorizonHeight", std::string("horizon_ele")) ;
    
    initValue("HeaderLowestWaterTableDepthPoint", std::string("none")) ;
    
    initValue("SoilIceContentTensorFile", std::string("none")) ;
    
    initValue("HeaderPointAspect", std::string("none")) ;
    
    initValue("SoilAveragedTempProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderGlacSublPoint", std::string("none")) ;
    
    initValue("HeaderDatePoint", std::string("none")) ;
    
    initValue("HeaderEvapSurfacePoint", std::string("none")) ;
    
    initValue("SnowTempProfileFileWriteEnd", std::string("none")) ;
    
    initValue("HeaderRHPoint", std::string("none")) ;
    
    initValue("SoilParFile", std::string("soil/soil")) ;
    
    initValue("SnowDepthLayersFileWriteEnd", std::string("none")) ;
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
    lLeafAngles += 0,0.3,0,0.3,0,0.1,0.01,0.01 ;
    initValue("LeafAngles", lLeafAngles) ;
    
    initValue("OutputDepthsVertical", double(0)) ;
    
    initValue("OutputVegetationMaps", double(0)) ;
    
    initValue("ExitMinLambdaEnergy", double(0)) ;
    
    initValue("SpinUpLayerBottom", double(10000)) ;
    
    initValue("FreeDrainageAtLateralBorder", double(1)) ;
    
    initValue("SWEPoint", double(4)) ;
    
    initValue("DepthFreeSurfaceAtTheBoundary", double(0)) ;
    
    initValue("MaxCourantSupFlowChannel", double(0.1)) ;
    
    initValue("DateBasin", double(-1)) ;
    
    initValue("GlacDensityPoint", double(-1)) ;
    
    initValue("ExitMinLambdaWater", double(1)) ;
    
    initValue("InitSnowTemp", double(0)) ;
    
    std::vector<double> lNVanGenuchtenBedrock ;
    lNVanGenuchtenBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("NVanGenuchtenBedrock", lNVanGenuchtenBedrock) ;
    
    initValue("SoilHeatFluxPoint", double(8)) ;
    
    initValue("WindAsWindXAndWindY", double(0)) ;
    
    initValue("z0vegPoint", double(-1)) ;
    
    initValue("PSnowBasin", double(-1)) ;
    
    std::vector<double> lSoilEmissiv ;
    lSoilEmissiv += 0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99 ;
    initValue("SoilEmissiv", lSoilEmissiv) ;
    
    initValue("DDChannel", double(1)) ;
    
    initValue("HeatEqMaxIter", double(700)) ;
    
    initValue("MaxGlacLayersMiddle", double(0)) ;
    
    initValue("NumberDayIntervalsToCalculateCloudiness", double(3)) ;
    
    initValue("TimeStepBlowingSnow", double(3600)) ;
    
    std::vector<double> lNormalHydrConductivityBedrock ;
    lNormalHydrConductivityBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("NormalHydrConductivityBedrock", lNormalHydrConductivityBedrock) ;
    
    initValue("NumLowPassFilterOnDemForCurv", double(0)) ;
    
    std::vector<double> lCanopyFraction ;
    lCanopyFraction += 0,0,0,0,0,0,0,0 ;
    initValue("CanopyFraction", lCanopyFraction) ;
    
    std::vector<double> lVegTransNIR ;
    lVegTransNIR += 0,0.32,0,0.32,0,0.22,0.09,0.32 ;
    initValue("VegTransNIR", lVegTransNIR) ;
    
    initValue("MaxCourantSupFlowChannelLand", double(0.1)) ;
    
    initValue("JulianDayFromYear0Glac", double(-1)) ;
    
    initValue("OutputGlacierMaps", double(0)) ;
    
    std::vector<double> lThetaSat ;
    lThetaSat += 0.5,0.5,0.5,0.5,0.5 ;
    initValue("ThetaSat", lThetaSat) ;
    
    std::vector<double> lVegSnowBurying ;
    lVegSnowBurying += 1,1,1,1,1,1,1,1 ;
    initValue("VegSnowBurying", lVegSnowBurying) ;
    
    std::vector<double> lGlacPlotDepths ;
    lGlacPlotDepths += -9999,-9999 ;
    initValue("GlacPlotDepths", lGlacPlotDepths) ;
    
    initValue("NumLandCoverTypes", double(8)) ;
    
    initValue("SnowDensityPoint", double(-1)) ;
    
    std::vector<double> lLinearInterpolation ;
    lLinearInterpolation += 0,0,0 ;
    initValue("LinearInterpolation", lLinearInterpolation) ;
    
    initValue("PointSim", double(0)) ;
    
    initValue("LEvBasin", double(-1)) ;
    
    initValue("SWbeamPoint", double(10)) ;
    
    initValue("SurfaceEnergyFlux", double(-9999)) ;
    
    initValue("DrySnowDefRate", double(1)) ;
    
    initValue("OutputSoilMaps", double(1)) ;
    
    initValue("InitWaterTableDepth", double(1000)) ;
    
    initValue("Lozone", double(0.3)) ;
    
    initValue("RicalculateCloudiness", double(0)) ;
    
    std::vector<double> lSoilPlotDepths ;
    lSoilPlotDepths += -9999,-9999 ;
    initValue("SoilPlotDepths", lSoilPlotDepths) ;
    
    initValue("BaseIPrec", double(0)) ;
    
    initValue("DDLand", double(1)) ;
    
    initValue("MinTimeStep", double(10)) ;
    
    std::vector<double> lCanDensSurface ;
    lCanDensSurface += 0,1,0,1,0,20,20,5 ;
    initValue("CanDensSurface", lCanDensSurface) ;
    
    initValue("RunBasin", double(-1)) ;
    
    std::vector<double> lNormalHydrConductivity ;
    lNormalHydrConductivity += 0.0001,0.0001,0.0001,0.0001,0.0001 ;
    initValue("NormalHydrConductivity", lNormalHydrConductivity) ;
    
    initValue("MinIceContentForBlowingSnow", double(8)) ;
    
    initValue("SnowThermalConductivityPar", double(1)) ;
    
    initValue("HPoint", double(-1)) ;
    
    initValue("SoilLayerThicknesses", double(-9999)) ;
    
    initValue("SnowCorrFactor", double(1.3)) ;
    
    initValue("ThresWaterDepthLandInf", double(0)) ;
    
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
    
    std::vector<double> lMaxSnowLayersMiddle ;
    lMaxSnowLayersMiddle += 5,5 ;
    initValue("MaxSnowLayersMiddle", lMaxSnowLayersMiddle) ;
    
    initValue("AirTempBasin", double(-1)) ;
    
    initValue("CanopyStabCorrection", double(1)) ;
    
    initValue("LWinMaxPoint", double(-1)) ;
    
    initValue("TDewPoint", double(14)) ;
    
    initValue("SnowMeltedPoint", double(-1)) ;
    
    initValue("SWinPoint", double(9)) ;
    
    initValue("SnowDepthPoint", double(3)) ;
    
    initValue("MeteoStationsID", double(-9999)) ;
    
    initValue("QSurfPoint", double(-1)) ;
    
    initValue("MeteoStationWindVelocitySensorHeight", double(2)) ;
    
    initValue("ThresTempRain", double(3)) ;
    
    initValue("HighestWaterTableDepthPoint", double(-1)) ;
    
    initValue("EvapSurfaceBasin", double(-1)) ;
    
    initValue("GlacAll", double(0)) ;
    
    std::vector<double> lSnowPlotDepths ;
    lSnowPlotDepths += -9999,-9999 ;
    initValue("SnowPlotDepths", lSnowPlotDepths) ;
    
    initValue("RichardTol", double(1e-06)) ;
    
    initValue("WindCompaction1D", double(0)) ;
    
    initValue("InitGlacierDepth", double(0)) ;
    
    initValue("OutputSnowMaps", double(1)) ;
    
    initValue("SoilLayerNumber", double(5)) ;
    
    initValue("InitInNewPeriods", double(0)) ;
    
    initValue("ThresWaterDepthChannel", double(50)) ;
    
    initValue("MoninObukhov", double(1)) ;
    
    std::vector<double> lSpecificStorativityBedrock ;
    lSpecificStorativityBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("SpecificStorativityBedrock", lSpecificStorativityBedrock) ;
    
    initValue("DtPlotDischarge", double(0)) ;
    
    std::vector<double> lVMualem ;
    lVMualem += 0.5,0.5,0.5,0.5,0.5 ;
    initValue("VMualem", lVMualem) ;
    
    std::vector<double> lWiltingPointBedrock ;
    lWiltingPointBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("WiltingPointBedrock", lWiltingPointBedrock) ;
    
    initValue("Longitude", double(11.7)) ;
    
    initValue("LWinBasin", double(-1)) ;
    
    initValue("EstoredCanopyPoint", double(-1)) ;
    
    initValue("IDPointPoint", double(-1)) ;
    
    std::vector<double> lThetaRes ;
    lThetaRes += 0.05,0.05,0.05,0.05,0.05 ;
    initValue("ThetaRes", lThetaRes) ;
    
    initValue("AlbExtParSnow", double(10)) ;
    
    initValue("NumSimulationTimes", double(1)) ;
    
    std::vector<double> lThresSnowVegDown ;
    lThresSnowVegDown += 0,200,0,200,0,1900,1900,800 ;
    initValue("ThresSnowVegDown", lThresSnowVegDown) ;
    
    std::vector<double> lSurFlowResLand ;
    lSurFlowResLand += 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 ;
    initValue("SurFlowResLand", lSurFlowResLand) ;
    
    initValue("SnowTempPoint", double(-1)) ;
    
    initValue("SWEBlownPoint", double(-1)) ;
    
    initValue("MaxTimesMinLambdaWater", double(0)) ;
    
    initValue("LWinPoint", double(12)) ;
    
    initValue("LObukhovCanopyPoint", double(-1)) ;
    
    std::vector<double> lLSAI ;
    lLSAI += 0,2,0,2,0,4,4,2 ;
    initValue("LSAI", lLSAI) ;
    
    initValue("SurFlowResChannel", double(20)) ;
    
    initValue("PeriodPoint", double(-1)) ;
    
    initValue("AirPressPoint", double(-1)) ;
    
    initValue("DEMRotationAngle", double(0)) ;
    
    std::vector<double> lMeteoStationStandardTime ;
    lMeteoStationStandardTime += 1,1,1 ;
    initValue("MeteoStationStandardTime", lMeteoStationStandardTime) ;
    
    std::vector<double> lSavingPoints ;
    lSavingPoints += 10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,165 ;
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
    lMeteoStationCoordinateY += 5.00857e+06,5.01001e+06,5.0065e+06 ;
    initValue("MeteoStationCoordinateY", lMeteoStationCoordinateY) ;
    
    std::vector<double> lMeteoStationCoordinateX ;
    lMeteoStationCoordinateX += 641960,643962,644544 ;
    initValue("MeteoStationCoordinateX", lMeteoStationCoordinateX) ;
    
    initValue("MeanTimeStep", double(-1)) ;
    
    initValue("TimeStepEnergyAndWater", double(3600)) ;
    
    std::vector<double> lSoilRoughness ;
    lSoilRoughness += 10,10,10,10,10,10,10,10 ;
    initValue("SoilRoughness", lSoilRoughness) ;
    
    initValue("WetSnowDefRate", double(1.5)) ;
    
    std::vector<double> lMeteoStationLatitude ;
    lMeteoStationLatitude += 46.5563,46.461,46.45 ;
    initValue("MeteoStationLatitude", lMeteoStationLatitude) ;
    
    initValue("MinSupWaterDepthLand", double(1)) ;
    
    std::vector<double> lFieldCapacity ;
    lFieldCapacity += -9999,-9999,-9999,-9999,-9999 ;
    initValue("FieldCapacity", lFieldCapacity) ;
    
    initValue("StandardTimeSimulation", double(0)) ;
    
    std::vector<double> lThresSnowVegUp ;
    lThresSnowVegUp += 0,200,0,200,0,1900,1900,800 ;
    initValue("ThresSnowVegUp", lThresSnowVegUp) ;
    
    initValue("BottomBoundaryHeatFlux", double(0)) ;
    
    initValue("GlacMeltedPoint", double(-1)) ;
    
    initValue("SurfaceTemperature", double(-9999)) ;
    
    initValue("WaterOnCanopyPoint", double(-1)) ;
    
    initValue("PrainNetPoint", double(-1)) ;
    
    initValue("DefaultSoilTypeLand", double(1)) ;
    
    initValue("PrecAsIntensity", double(0)) ;
    
    initValue("QVegPoint", double(-1)) ;
    
    initValue("TimeFromStartBasin", double(-1)) ;
    
    std::vector<double> lVegTransVis ;
    lVegTransVis += 0,0.07,0,0.07,0,0.04,0.04,0.07 ;
    initValue("VegTransVis", lVegTransVis) ;
    
    initValue("TsMaxIter", double(2)) ;
    
    initValue("RHmin", double(10)) ;
    
    initValue("WindSpeedTopCanopyPoint", double(-1)) ;
    
    initValue("PrainPoint", double(-1)) ;
    
    initValue("LWNetPoint", double(16)) ;
    
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
    lThermalConductivitySoilSolids += 2.5,2.5,2.5,2.5,2.5 ;
    initValue("ThermalConductivitySoilSolids", lThermalConductivitySoilSolids) ;
    
    std::vector<double> lMeteoStationElevation ;
    lMeteoStationElevation += 1200,1450,990 ;
    initValue("MeteoStationElevation", lMeteoStationElevation) ;
    
    std::vector<double> lSoilAlbNIRDry ;
    lSoilAlbNIRDry += 0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33 ;
    initValue("SoilAlbNIRDry", lSoilAlbNIRDry) ;
    
    initValue("BaseRelativeHumidity", double(70)) ;
    
    initValue("InitDateDDMMYYYYhhmm", double(2.0102e+11)) ;
    initValue("CurvatureWeightI", double(0)) ;
    
    initValue("HBasin", double(-1)) ;
    
    std::vector<double> lVMualemBedrock ;
    lVMualemBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("VMualemBedrock", lVMualemBedrock) ;
    
    initValue("EndDateDDMMYYYYhhmm", double(2.9102e+11)) ;
    initValue("QAirPoint", double(-1)) ;
    
    initValue("PsnowPoint", double(-1)) ;
    
    initValue("CurvatureWeightD", double(0)) ;
    
    initValue("NumberOfMeteoStations", double(3)) ;
    
    initValue("TvegBasin", double(-1)) ;
    
    initValue("GlacTempPoint", double(-1)) ;
    
    initValue("DtPlotBasin", double(0)) ;
    
    initValue("UpwindBorderBlowingSnow", double(0)) ;
    
    std::vector<double> lThetaResBedrock ;
    lThetaResBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("ThetaResBedrock", lThetaResBedrock) ;
    
    std::vector<double> lAlphaVanGenuchtenBedrock ;
    lAlphaVanGenuchtenBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("AlphaVanGenuchtenBedrock", lAlphaVanGenuchtenBedrock) ;
    
    initValue("CoordinatePointX", double(-9999)) ;
    
    initValue("CoordinatePointY", double(-9999)) ;
    
    std::vector<double> lVegReflectVis ;
    lVegReflectVis += 0,0.15,0,0.15,0,0.12,0.09,0.15 ;
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
    
    std::vector<double> lMaxWaterEqSnowLayerContent ;
    lMaxWaterEqSnowLayerContent += 10,10 ;
    initValue("MaxWaterEqSnowLayerContent", lMaxWaterEqSnowLayerContent) ;
    
    initValue("SoilAll", double(1)) ;
    
    initValue("LocMaxIter", double(3)) ;
    
    initValue("SWNetPoint", double(17)) ;
    
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
    
    initValue("OutputMeteoMaps", double(0)) ;
    
    initValue("Vmin", double(0.5)) ;
    
    initValue("SWESublBlownPoint", double(-1)) ;
    
    initValue("InitGlacierDensity", double(800)) ;
    
    initValue("MinDiffSupWaterDepthChannelLand", double(1)) ;
    
    initValue("CalculateCastShadow", double(1)) ;
    
    initValue("GWEbottom", double(0)) ;
    
    initValue("MinRatioKactualToKSat", double(0)) ;
    
    initValue("InitSnowAge", double(0)) ;
    
    std::vector<double> lMeteoStationLongitude ;
    lMeteoStationLongitude += 12.4259,12.4103,12.417 ;
    initValue("MeteoStationLongitude", lMeteoStationLongitude) ;
    
    initValue("AirTempPoint", double(5)) ;
    
    std::vector<double> lWiltingPoint ;
    lWiltingPoint += -9999,-9999,-9999,-9999,-9999 ;
    initValue("WiltingPoint", lWiltingPoint) ;
    
    initValue("AngstromBeta", double(0.1)) ;
    
    initValue("SnowAgingCoeffVis", double(0.2)) ;
    
    initValue("IrriducibleWatSatSnow", double(0.02)) ;
    
    initValue("DatePoint", double(1)) ;
    
    std::vector<double> lThetaSatBedrock ;
    lThetaSatBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("ThetaSatBedrock", lThetaSatBedrock) ;
    
    initValue("BaseWindSpeed", double(0.5)) ;
    
    initValue("InitSWE", double(0)) ;
    
    initValue("SoilLayerTypes", double(1)) ;
    
    initValue("SlopeWeightD", double(0)) ;
    
    initValue("SlopeWeightI", double(0)) ;
    
    std::vector<double> lInitSoilPressureBedrock ;
    lInitSoilPressureBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("InitSoilPressureBedrock", lInitSoilPressureBedrock) ;
    
    std::vector<double> lRoughElemXUnitArea ;
    lRoughElemXUnitArea += 0,0,0,0,0,0,0,0 ;
    initValue("RoughElemXUnitArea", lRoughElemXUnitArea) ;
    
    std::vector<double> lThermalCapacitySoilSolidsBedrock ;
    lThermalCapacitySoilSolidsBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("ThermalCapacitySoilSolidsBedrock", lThermalCapacitySoilSolidsBedrock) ;
    
    initValue("LWvBasin", double(-1)) ;
    
    initValue("MinLambdaWater", double(1e-07)) ;
    
    initValue("ConvectiveHeatTransferCoefficient", double(-9999)) ;
    
    initValue("BusingerMaxIter", double(5)) ;
    
    initValue("SnowViscosity", double(1e+06)) ;
    
    initValue("SWdiffPoint", double(11)) ;
    
    initValue("SimulationHours", double(1)) ;
    
    initValue("MinDiffSupWaterDepthLandChannel", double(1)) ;
    
    initValue("HvPoint", double(-1)) ;
    
    initValue("SnowSMIN", double(30)) ;
    
    initValue("PRainBasin", double(-1)) ;
    
    initValue("RichardInitForc", double(0.01)) ;
    
    initValue("ActualOrProjectedArea", double(0)) ;
    
    std::vector<double> lMinStomatalRes ;
    lMinStomatalRes += 0,60,0,60,0,60,60,60 ;
    initValue("MinStomatalRes", lMinStomatalRes) ;
    
    std::vector<double> lLateralHydrConductivity ;
    lLateralHydrConductivity += 0.0001,0.0001,0.0001,0.0001,0.0001 ;
    initValue("LateralHydrConductivity", lLateralHydrConductivity) ;
    
    initValue("RunGlac", double(-1)) ;
    
    initValue("SnowAgingCoeffNIR", double(0.5)) ;
    
    initValue("PsnowNetPoint", double(-1)) ;
    
    initValue("PointID", double(-9999)) ;
    
    initValue("SWupPoint", double(-1)) ;
    
    std::vector<double> lDecayCoeffCanopy ;
    lDecayCoeffCanopy += 0,2.5,0,2.5,0,4,4,2.5 ;
    initValue("DecayCoeffCanopy", lDecayCoeffCanopy) ;
    
    initValue("HighestNodeCorrespondsToLayer", double(0)) ;
    
    initValue("TsurfPoint", double(6)) ;
    
    std::vector<double> lNVanGenuchten ;
    lNVanGenuchten += 1.3,1.3,1.3,1.3,1.3 ;
    initValue("NVanGenuchten", lNVanGenuchten) ;
    
    initValue("PNetBasin", double(-1)) ;
    
    initValue("AlphaSnow", double(100000)) ;
    
    initValue("Latitude", double(46.3)) ;
    
    std::vector<double> lRoughElemDiam ;
    lRoughElemDiam += 50,50,50,50,50,50,50,50 ;
    initValue("RoughElemDiam", lRoughElemDiam) ;
    
    initValue("LWinParameterization", double(9)) ;
    
    initValue("LEgVegPoint", double(-1)) ;
    
    initValue("SnowDensityCutoff", double(100)) ;
    
    initValue("DefaultSoilTypeBedrock", double(1)) ;
    
    initValue("PrainOnSnowPoint", double(-1)) ;
    
    initValue("ConsiderMicrometeorology", double(1)) ;
    
    std::vector<double> lLateralHydrConductivityBedrock ;
    lLateralHydrConductivityBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("LateralHydrConductivityBedrock", lLateralHydrConductivityBedrock) ;
    
    initValue("LapseRateDewTemp", double(-9999)) ;
    
    initValue("SWEbottom", double(20)) ;
    
    std::vector<double> lSoilAlbNIRWet ;
    lSoilAlbNIRWet += 0.16,0.16,0.16,0.16,0.16,0.16,0.16,0.16 ;
    initValue("SoilAlbNIRWet", lSoilAlbNIRWet) ;
    
    initValue("JulianDayFromYear0Point", double(2)) ;
    
    initValue("FreeDrainageAtBottom", double(0)) ;
    
    initValue("MinTimeStepSupFlow", double(0.01)) ;
    
    std::vector<double> lMeteoStationSkyViewFactor ;
    lMeteoStationSkyViewFactor += 1,1,1 ;
    initValue("MeteoStationSkyViewFactor", lMeteoStationSkyViewFactor) ;
    
    initValue("PSnowNetBasin", double(-1)) ;
    
    initValue("MinLambdaEnergy", double(1e-05)) ;
    
    initValue("SnowOnCanopyPoint", double(-1)) ;
    
    initValue("MaxTimesMinLambdaEnergy", double(0)) ;
    
    initValue("LapseRateTemp", double(6.5)) ;
    
    initValue("AngstromAlpha", double(1.3)) ;
    
    initValue("FormatOutputMaps", double(3)) ;
    
    initValue("HgVegPoint", double(-1)) ;
    
    initValue("SWvPoint", double(-1)) ;
    
    std::vector<double> lInitSoilTempBedrock ;
    lInitSoilTempBedrock += -9999,-9999,-9999,-9999,-9999 ;
    initValue("InitSoilTempBedrock", lInitSoilTempBedrock) ;
    
    std::vector<double> lSoilAlbVisDry ;
    lSoilAlbVisDry += 0.16,0.16,0.16,0.16,0.16,0.16,0.16,0.16 ;
    initValue("SoilAlbVisDry", lSoilAlbVisDry) ;
    
    initValue("LEupPoint", double(-1)) ;
    
    std::vector<double> lRootDepth ;
    lRootDepth += 0,30,0,30,0,2000,2000,300 ;
    initValue("RootDepth", lRootDepth) ;
    
    initValue("RecoverSim", double(0)) ;
    
    initValue("SurfaceEBPoint", double(7)) ;
    
    initValue("PeriodGlac", double(-1)) ;
    
    initValue("LowestThawedSoilDepthPoint", double(-1)) ;
    
    initValue("RHAsDewTemperature", double(0)) ;
    
    initValue("BaseAirTemperature", double(5)) ;
    
    std::vector<double> lThermalCapacitySoilSolids ;
    lThermalCapacitySoilSolids += 2.3e+06,2.3e+06,2.3e+06,2.3e+06,2.3e+06 ;
    initValue("ThermalCapacitySoilSolids", lThermalCapacitySoilSolids) ;
    
    initValue("LEvPoint", double(-1)) ;
    
    std::vector<double> lVegHeight ;
    lVegHeight += 0,200,0,200,0,1900,1900,800 ;
    initValue("VegHeight", lVegHeight) ;
    
    initValue("ChannelDepression", double(500)) ;
    
    initValue("RainCorrFactor", double(1)) ;
    
    std::vector<double> lSnowAll ;
    lSnowAll += 1,1 ;
    initValue("SnowAll", lSnowAll) ; 
    
    initValue("LapseRatePrec", double(-0.2)) ;
    
    initValue("InitGlacierTemp", double(-3)) ;
    
    initValue("MeteoStationTemperatureSensorHeight", double(5)) ;
    
    std::vector<double> lVegReflNIR ;
    lVegReflNIR += 0,0.4,0,0.4,0,0.43,0.36,0.4 ; 
    initValue("VegReflNIR", lVegReflNIR) ; 
    
    initValue("ContinuousRecovery", double(0)) ;
    
    initValue("d0vegPoint", double(-1)) ;
    
    initValue("LObukhovPoint", double(-1)) ;
    
    initValue("WindAsSpeedAndDirection", double(1)) ;
    
    initValue("NumLowPassFilterOnDemForAll", double(0)) ;
    
    initValue("TSurfBasin", double(-1)) ;
    
    initValue("MinSupWaterDepthChannel", double(1)) ;
    
    initValue("WindSpeedPoint", double(15)) ;
    
    initValue("RHPoint", double(13)) ;
    
    initValue("QCanopyAirPoint", double(-1)) ;
    
    initValue("RatioChannelWidthPixelWidth", double(0.1)) ;
    
    std::vector<double> lThermalConductivitySoilSolidsBedrock ;
    lThermalConductivitySoilSolidsBedrock += -9999,-9999,-9999,-9999,-9999 ; 
    initValue("ThermalConductivitySoilSolidsBedrock", lThermalConductivitySoilSolidsBedrock) ; 
    
    initValue("LSAIPoint", double(-1)) ;
    
    initValue("FreshSnowReflVis", double(0.95)) ;
    
    initValue("OutputSurfEBALMaps", double(0)) ;
    
    initValue("PointAll", double(0)) ;
    
    initValue("DewTempOrNormTemp", double(0)) ;
    
    initValue("PeriodBasin", double(-1)) ;
    
    std::vector<double> lFieldCapacityBedrock ;
    lFieldCapacityBedrock += -9999,-9999,-9999,-9999,-9999 ; 
    initValue("FieldCapacityBedrock", lFieldCapacityBedrock) ; 
    
    initValue("RunPoint", double(-1)) ;
    
    initValue("MaxCourantSupFlowLand", double(0.1)) ;
    
    initValue("LowestWaterTableDepthPoint", double(-1)) ;
    
    initValue("DtPlotPoint", double(1)) ;
    
    std::vector<double> lSpecificStorativity ;
    lSpecificStorativity += 1e-07,1e-07,1e-07,1e-07,1e-07 ; 
    initValue("SpecificStorativity", lSpecificStorativity) ; 
    
    initValue("HeatEqTol", double(1e-06)) ;
    
    initValue("DateGlac", double(-1)) ;
    
    std::vector<double> lSoilAlbVisWet ;
    lSoilAlbVisWet += 0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08 ; 
    initValue("SoilAlbVisWet", lSoilAlbVisWet) ; 
    
    initValue("SnowSMAX", double(80)) ;
    
    initValue("SWinBasin", double(-1)) ;
    
    std::vector<double> lInitSoilTemp ;
    lInitSoilTemp += 5,5,5,5,5 ; 
    initValue("InitSoilTemp", lInitSoilTemp) ; 
    
    initValue("TCanopyAirPoint", double(-1)) ;
    
    initValue("DefaultSoilTypeChannel", double(1)) ;
    
    initValue("InitSnowDensity", double(250)) ;
    
    initValue("SnowCURV", double(-200)) ;
    
    initValue("TraspCanopyBasin", double(-1)) ;
    
    initValue("SWEtop", double(50)) ;
    
    initValue("UpdateHydraulicConductivity", double(0)) ;
    
    std::vector<double> lThresSnowSoilRough ;
    lThresSnowSoilRough += 10,10,10,10,10,10,10,10 ; 
    initValue("ThresSnowSoilRough", lThresSnowSoilRough) ; 
    
    initValue("MaxWaterEqGlacLayerContent", double(5)) ;
    
    std::vector<double> lInitSoilPressure ;
    lInitSoilPressure += -9999,-9999,-9999,-9999,-9999 ; 
    initValue("InitSoilPressure", lInitSoilPressure) ; 
    
    initValue("SurFlowResExp", double(0.666667)) ;
    
    initValue("MinPrecIncreaseFactorWithElev", double(0.1)) ;
    //END INITIALIZATION OF NUMERIC PARAMETERS
}

geotop::input::ConfigStore::~ConfigStore()
{
    
}