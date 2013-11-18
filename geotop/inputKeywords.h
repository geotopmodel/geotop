/**
 * @file   inputKeywords.h
 * @Author Angelo Leto (angleto@gmail.com)
 * @date   November, 2013
 * @brief  generic configuration store class
 *
 * Parse configuration file and store parameters on single container
 */

#ifndef __INPUT_KEYWORDS__
#define __INPUT_KEYWORDS__

#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/shared_ptr.hpp>
#include <map>
#include <algorithm>
#include <iterator>
#include <string>
#include <iostream>

namespace geotop
{
    namespace input
    {
        /** @brief Configuration file parser and storage class
         */
        class ConfigStore
        {
        public:
            
            /** @brief ConfigStore constructor
             */
            ConfigStore();
            
            /** @brief ConfigStore constructor
             */
            virtual ~ConfigStore();
            
            /** @brief configuration parameters getter accessor, you must parse a
             *      file before to use this method otherwise you will obtain no values.
             *      The method infer the right type by the parameter passed by reference
             *      and handle the cast consequently
             *  @param[in] pName the name of the parameter to get, note that this name is case insensitive
             *  @param[out] rValue the return value of the requested parameter
             *  @return true if the value stored in rValue is valid, otherwise return false.
             *      False will be returned also if type of the return parameters passed is not
             *      compatible with the type of the requested
             *      parameter.
             */
            template <typename T> bool get(const std::string pName, T &rValue) const {
                
                std::string lName ( pName );
                boost::algorithm::to_lower(lName);

                if( not mValueMap->count(lName) )
                {
                    std::cerr << "Error: configuration item not found: " << pName << std::endl ;
                    return false ;
                }

                try
                {
                    rValue = boost::any_cast<T>(mValueMap->at(lName));
                    return true;
                }
                catch(const boost::bad_any_cast &)
                {
                    std::cerr << "Error: runtime typechecking failed for item: " << pName << std::endl ;
                    return false;
                }
            } ;

            /** @brief configuration parameters setter accessor, you must parse a
             *      file before to use this method otherwise you will no set any.
             *      The method infer the right type by the parameter passed by reference
             *      and handle the cast consequently.
             *  @param[in] pName the name of the parameter to modify, note that this name is case insensitive
             *  @param[in] pValue to set
             *  @return true if the value vas successfully modified, otherwise return false.
             *      If the keyword does not exixts the function will return false.
             *      False will be returned also if type of the return parameters passed is not
             *      compatible with the type of the parameter to be modifyed.
             */
            template <typename T> bool set(const std::string pName, const T &pValue) {
                
                std::string lName ( pName );
                boost::algorithm::to_lower(lName);
                
                if( not mValueMap->count(lName) )
                {
                    std::cerr << "Error: configuration item not present on configuration dictionary: " << pName << std::endl ;
                    return false ;
                }
                
                try
                {
                    boost::any_cast<T>(mValueMap->at(lName));
                }
                catch(const boost::bad_any_cast &)
                {
                    std::cerr << "Error: runtime typechecking failed for item: " << pName << std::endl ;
                    return false;
                }
                
                (*mValueMap)[lName] = pValue ;

                return true ;
            } ;

            /** @brief parse a configuration file
             * @param[in] pFileName the configuration file path
             * @return true if the configuration file was successfully parsed, otherwise return false
             */
            bool parse(const std::string pFileName) ;

            
        private:

            /** @internal
             * @brief this function is a relaxed version of the set
             * @param[in] pName the name of the parameter to modify, note that this name is case insensitive
             * @param[in] pValue to set
             * @return true if the value vas successfully initialized, otherwise return false.
             */
            template <typename T> bool initValue(const std::string pName, const T &pValue) {
                
                std::string lName ( pName );
                boost::algorithm::to_lower(lName);                
                (*mValueMap)[lName] = pValue ;
                
                return true ;
            } ;

            void init() ;
            boost::shared_ptr<std::map<std::string, boost::any> > mValueMap ;
        };
        
    } // namespace input
} // namespace geotop

#endif // __INPUT_KEYWORDS__

