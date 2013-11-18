/**
 * @file   testInputKeywords.cc
 * @Author Angelo Leto (angleto@gmail.com)
 * @date   November, 2013
 * @brief  test program for generic configuration store class
 *
 * test program for generic configuration store class
 */

#include <inputKeywords.h>
#include <boost/filesystem.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/version.hpp>
#include <string>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace boost::assign;

int main(int argc, char *argv[])
{

    boost::program_options::options_description usage("Usage"); // name of help function
    usage.add_options() //detailed specification of command line interface
    ("help,h", "this help") // this option return boolean variable and is used to print command line help
    ("file,f", boost::program_options::value<std::string>(), "input file path");

    // map for options/value
    boost::program_options::variables_map vm;

    // option parsing statements.
    try{
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(usage).run(), vm);
        boost::program_options::notify(vm);
    }
    catch ( boost::program_options::unknown_option &u ){
        std::cerr << "Option parsing error: " << u.what() << ": please use only valid flags"<< std::endl;
        return 10 ;
    }
    catch ( boost::program_options::invalid_command_line_syntax &u ){
        std::cerr << "Option parsing error: " << u.what() << ": please specificy a valid value for this option "<< std::endl;
        return 11 ;
    }
    catch ( boost::program_options::error &u ){
        std::cout << "Option parsing error: " << u.what() << std::endl;
        return 12 ;
    }
    
    // if no options or --help or -h print help and exit.
    if(argc <= 1 || vm.count("help"))
    {
        std::cout << usage << std::endl;
        return 13 ;
    }

    std::string lInputFile ;
    if (vm.count("file"))
    {
        lInputFile = vm["file"].as<std::string>();
    } else {
        std::cerr << "The input file mut be specified, see --input option" << std::endl ;
        return 20  ;
    }

    boost::filesystem::path lInputFilePath (lInputFile.c_str());
    if(not boost::filesystem::exists(lInputFile))
    {
        std::cerr << "Input file not found: " << lInputFilePath.string() << std::endl ;
        return 21 ;
    }

    geotop::input::ConfigStore lConfigStore ;
    lConfigStore.parse(lInputFilePath.string()) ;

    std::vector<double> lDecayCoeffCanopyOrig ;
    lDecayCoeffCanopyOrig += 0,2.5,0,2,5,0,4,4,2.5 ;
    std::vector<double> lDecayCoeffCanopy ;
    lConfigStore.get("DecayCoeffCanopy", lDecayCoeffCanopy) ;
    if(not std::equal ( lDecayCoeffCanopyOrig.begin(), lDecayCoeffCanopyOrig.end(), lDecayCoeffCanopy.begin() )) {
        std::cerr << "Error: DecayCoeffCanopy : FAILED" << std::endl ;
    } else {
        std::cout << "Info: DecayCoeffCanopy : OK" << std::endl ;
    }

    return 0 ;
}
