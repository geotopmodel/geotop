/** 
 * @file global_logger.h 
 * @Author Gianfranco Gallizia (skyglobe83@gmail.com) 
 * @copyright (C) 2014 eXact lab srl
 * 
 */ 
 
#ifndef GLOBAL_LOGGER_H_INCLUDED 
#define GLOBAL_LOGGER_H_INCLUDED 1 

#include "logger.h"
#include <fstream>
#include <string>

namespace geotop
{
    namespace logger
    {
        class GlobalLogger
        {
        public:
            /**
             * @brief Gives access to the Singleton
             * @return A pointer to the single instance
             *
             * See:
             * http://www.codeproject.com/Articles/1921/Singleton-Pattern-its-implementation-with-C
             */
            static GlobalLogger* getInstance();

            /**
             * @brief Destructor
             */
            ~GlobalLogger();

            /**
             * @brief Changes the minimum severity accepted
             * @param[in] minSeverity the minimum severity level accepted
             * 
             */
            void setSeverity(severity_levels minSeverity);

            /**
             * @brief Gives the full path of the log file
             * @return A constant pointer to the full path
             */
            std::string getLogFilePath();
            /**
             * @see Logger
             */
            void log(std::string const &message, severity_levels severity = Logger::DEFAULT_SEVERITY);
            void logf(const char* format, ...);
            void logsf(severity_levels severity, const char* format, ...);
            void writeAll(std::string const &message);
            void writefAll(const char* format, ...);
        private:
            GlobalLogger();
            static bool instanceFlag;
            static GlobalLogger* single;
            Logger mLogger;
            std::ofstream logfileStream;
            std::string append_path(std::string path, std::string filename);
        };
        
    }
}

#endif 

