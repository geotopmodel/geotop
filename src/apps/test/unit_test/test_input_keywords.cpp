#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE Main
#endif
#include <boost/test/unit_test.hpp>

#include "../../../geotop/inputKeywords.h"

using namespace boost::assign;

struct ToolsFixture
{

    geotop::input::ConfigStore mConfigStore ;

    ToolsFixture()
    {
        BOOST_TEST_MESSAGE("setup tools");
        
        const std::string lFilePath = "test_data/test_input_keywords/geotop.inpts.real_case" ;
        bool lResult = mConfigStore.parse(lFilePath) ;
        BOOST_CHECK(lResult) ;
    }

    ~ToolsFixture()
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(test_queue, ToolsFixture)

BOOST_AUTO_TEST_CASE(parse_files)
{
    // parse file test
    {
        geotop::input::ConfigStore lConfigStore ;
        const std::string lFilePath = "test_data/test_input_keywords/geotop.inpts.0" ;
        bool lResult = lConfigStore.parse(lFilePath) ;
        
        BOOST_CHECK( lResult );
    }

    // parse file test
    {
        geotop::input::ConfigStore lConfigStore ;
        const std::string lFilePath = "test_data/test_input_keywords/geotop.inpts.1" ;
        bool lResult = lConfigStore.parse(lFilePath) ;
        
        BOOST_CHECK( lResult );
    }

    // parse file test
    {
        geotop::input::ConfigStore lConfigStore ;
        const std::string lFilePath = "test_data/test_input_keywords/geotop.inpts.2" ;
        bool lResult = lConfigStore.parse(lFilePath) ;
        
        BOOST_CHECK( lResult );
    }

    // parse file test
    {
        geotop::input::ConfigStore lConfigStore ;
        const std::string lFilePath = "test_data/test_input_keywords/geotop.inpts.3" ;
        bool lResult = lConfigStore.parse(lFilePath) ;
        
        BOOST_CHECK( lResult );
    }

    // parse file test
    {
        geotop::input::ConfigStore lConfigStore ;
        const std::string lFilePath = "test_data/test_input_keywords/geotop.inpts.4" ;
        bool lResult = lConfigStore.parse(lFilePath) ;
        
        BOOST_CHECK( lResult );
    }

}

BOOST_AUTO_TEST_CASE(set_get_doublearray_value)
{
    {
        // set double array value test
        std::vector<double> lDecayCoeffCanopy ;
        lDecayCoeffCanopy += 0,2.5,0,2.8,0,4,4,2.5 ;
        BOOST_TEST_MESSAGE("setup tools");
        
        bool lResult = mConfigStore.set("DecayCoeffCanopy", lDecayCoeffCanopy) ;
        
        BOOST_CHECK( lResult );
    }
    
    {
        // get double array value test
        std::vector<double> lDecayCoeffCanopyTest ;
        lDecayCoeffCanopyTest += 0,2.5,0,2.8,0,4,4,2.5 ;
        std::vector<double> lDecayCoeffCanopy ;
        bool lGetResult = mConfigStore.get("DecayCoeffCanopy", lDecayCoeffCanopy) ;
        BOOST_CHECK( lGetResult );
        
        bool lDataCheck = std::equal ( lDecayCoeffCanopyTest.begin(), lDecayCoeffCanopyTest.end(), lDecayCoeffCanopy.begin() ) ;
        BOOST_CHECK( lDataCheck );
    }
}


BOOST_AUTO_TEST_SUITE_END()
