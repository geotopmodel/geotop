/**
 * @file logger.cc
 * @Author Gianfranco Gallizia (skyglobe83@gmail.com)
 * @brief Log facility class
 *
 * GEOtop logging facility
 */

#include "logger.h"

using namespace geotop::logger;

Logger::Logger(std::ostream* stream)
{
    mStreams.push_back(LogStream(stream,DEFAULT_SEVERITY));
}

Logger::Logger()
{
    mStreams.push_back(
                LogStream((std::ostream*)(&std::clog),DEFAULT_SEVERITY));
}

Logger::Logger(std::ostream* stream, severity_levels minSeverity)
{
    mStreams.push_back(LogStream(stream,minSeverity));
}

Logger::~Logger()
{
}

void Logger::addOStream(std::ostream* stream, severity_levels minSeverity)
{
    mStreams.push_back(LogStream(stream,minSeverity));
}

void Logger::log(std::string const &message, severity_levels severity)
{

    try
    {
        std::ostream* mStream;
        std::vector<LogStream>::iterator it = mStreams.begin();
        while(it != mStreams.end()){
            mStream = it->streamP;
            if (severity >= it->severity)
            *mStream << "[" << severity_labels[severity] << "]: "
                     << message << std::endl;
            mStream->flush();
            it++;
        }
    }
    catch(...)
    {
        //Catch and ignore any exception thrown
    }
}

void Logger::logf(const char* format, ...)
{
    char buffer[MAXMESSAGESIZE];
    va_list args;
    int chs;

    va_start(args, format);
    chs = vsprintf(buffer, format, args);

    if (chs != -1 && chs < MAXMESSAGESIZE)
    {
        std::string msg(buffer);
        log(msg, DEFAULT_SEVERITY);
    }
}

void Logger::logsf(severity_levels severity, const char* format, ...)
{
    char buffer[MAXMESSAGESIZE];
    va_list args;
    int chs;

    va_start(args, format);
    chs = vsprintf(buffer, format, args);

    if (chs != -1 && chs < MAXMESSAGESIZE)
    {
        std::string msg(buffer);
        log(msg, severity);
    }
}

