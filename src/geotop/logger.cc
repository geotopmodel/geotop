//
// Created by alberto on 4/18/18.
//

#include "logger.h"

// the default logger
Logger geolog;

Logger::Logger() :
        std_out{&std::cout},
        ofile{nullptr},
        _at_new_line{true}{
  push("geotop");
}

void Logger::push(const std::string &s) {
  std::string prefix;
  if (prefixes.size() > 0)
    prefix = prefixes.top();
  prefix += s + ":";
  prefixes.emplace(prefix);
}

std::string Logger::pop() {
  std::string s = prefixes.top();
  prefixes.pop();
  return s;
}

const std::string &Logger::prefix() const {
  static std::string s;
  if(prefixes.size() > 0)
    return prefixes.top();
  else
    return s;
}

Logger::Prefix::Prefix(const std::string &s, Logger &l) : log{&l} { this->log->push(s); }

Logger::Prefix::~Prefix() {
  log->pop();
}

Logger::Prefix::Prefix(const std::string &s): Prefix{s,geolog} {}
