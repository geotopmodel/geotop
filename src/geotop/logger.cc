//
// Created by alberto on 4/18/18.
//

#include "logger.h"
#include <limits>

// the default logger
Logger geolog;

Logger::Logger()
    : std_out{&std::cout}, ofile{nullptr}, _at_new_line{true},
      _console_level{std::numeric_limits<unsigned int>::max()},
      _file_level{std::numeric_limits<unsigned int>::max()} {
  push("geotop");
}

void Logger::pop() {
  if (!prefixes.empty()) {
    prefixes.pop();
  } else
    throw std::runtime_error{"You called pop() on and empty stack"};
}

const std::string &Logger::prefix() const {
  static std::string s;
  if (!prefixes.empty())
    return prefixes.top();
  else
    return s;
}
