/**
* @file logger.h
* @Author Gianfranco Gallizia (skyglobe83@gmail.com)
* @brief Log facility header
*
* Definitions header for geotop::logger
*/

#ifndef GEOTOP_LOGGER_H_INCLUDED
#define GEOTOP_LOGGER_H_INCLUDED 1

#include <string>
#include <iostream>
#include <cstdarg>
#include <cstdio>
#include <vector>

namespace geotop
{
    namespace logger
    {
        /**
         * @brief Severity levels
         *
         * Defines the severity levels for the logging facility
         */
        enum severity_levels
        {
            TRACE = 0,
            DEBUG = 1,
            NOTICE = 2,
            WARNING = 3,
            ERROR = 4,
            CRITICAL = 5
        };


        const std::string severity_labels[] = {
            "TRACE",
            "DEBUG",
            "NOTICE",
            "WARNING",
            "ERROR",
            "CRITICAL"
        };

        const int MAXMESSAGESIZE = 2048;


        /**
         * @brief The logging class
         */
        class Logger
        {
        public:
            /**
             * @brief Default constructor
             *
             * Uses std::clog as output stream
             */
            Logger();

            /**
             * @brief OStream constructor
             * @param[in] stream output stream where the log messages will be sent
             */
            Logger(std::ostream* stream);

            /**
             * @brief OStream constructor with minimum severity
             * @param[in] stream output stream where the log messages will be sent
             * @param[in] minSeverity the minimum level of severity accepted
             *
             */
            Logger(std::ostream* stream, severity_levels minSeverity);

            /**
             * @brief Destructor
             */
            virtual ~Logger();

            /**
             * @brief Default severity level
             */
            static const severity_levels DEFAULT_SEVERITY = NOTICE;

            /**
             * @brief Adds a new ostream to the logger
             * @param[in] stream The new stream to add
             * @param[in] minSeverity the minimum level of severity accepted
             *
             */
            void addOStream(std::ostream* stream, severity_levels minSeverity = DEFAULT_SEVERITY);

            /**
             * @brief Logs a message
             * @param[in] message The message to log
             * @param[in] severity The message's severity
             *
             * Takes message and sends it to an ostream
             */
            void log(std::string const &message, severity_levels severity = DEFAULT_SEVERITY);
           
            /**
             * @brief Logs a formatted message
             * @param[in] format a printf-style format string
             *
             * Takes a variable number of arguments and uses format to build a
             * log message then logs it to an ostream
             */
            void logf(const char* format, ...);

            /**
             * @brief Logs a formatted message with a user-selected severity
             * @param[in] severity The message's severity
             * @param[in] format a printf-style format string
             *
             * Like logf but with a user-defined severity level
             */
            void logsf(severity_levels severity, const char* format, ...);

        private:
            /**
             * @internal
             * @brief Defines a simple (ostream*, severity_levels) pair
             */
            struct LogStream
            {
                std::ostream* streamP;
                severity_levels severity;
            } ;
            std::vector<struct LogStream> mStreams;
        };
    }
}

#endif // GEOTOP_LOGGER_H_INCLUDED
