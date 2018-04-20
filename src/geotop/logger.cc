//
// Created by alberto on 4/18/18.
//

#include "logger.h"
#include <limits>

// the default logger
Logger geolog;

Logger::Logger()
  : std_out{&std::cout},
    ofile{nullptr},
    _at_new_line{true},
    _console_level{std::numeric_limits<unsigned int>::max()},
    _file_level{std::numeric_limits<unsigned int>::max()} {
  push("geotop");
}

void Logger::push(const std::string& s) {
  std::string prefix;
  if (prefixes.size() > 0) prefix = prefixes.top();
  prefix += s + ":";
  prefixes.emplace(prefix);
}

std::string Logger::pop() {
  if (prefixes.size() > 0) {
    std::string s = prefixes.top();
    prefixes.pop();
    return s;
  }
  throw std::runtime_error{"You called pop() on and empty stack"};
}

const std::string& Logger::prefix() const {
  static std::string s;
  if (prefixes.size() > 0)
    return prefixes.top();
  else
    return s;
}

void Logger::attach_file_stream(std::ostream& o) {
  ofile = &o;
}

void Logger::detach_file_stream() {
  ofile = nullptr;
}

Logger::ScopedPrefix::ScopedPrefix(const std::string& s, Logger& l) : log{&l} {
  log->push(s);
}

Logger::ScopedPrefix::~ScopedPrefix() {
  log->pop();
}

Logger::ScopedPrefix::ScopedPrefix(const std::string& s)
  : ScopedPrefix{s, geolog} {}

Logger::ScopedConsoleLevel::ScopedConsoleLevel(const unsigned int level,
                                               Logger& l)
  : log{&l},
    console{[this](const unsigned int ll) { log->set_console_level(ll); },
            level,
    l.console_level()} {}

Logger::ScopedFileLevel::ScopedFileLevel(const unsigned int level, Logger& l)
  : log{&l},
    file{[this](const unsigned int ll) { log->set_file_level(ll); }, level,
    l.file_level()} {}

Logger::ScopedLevels::ScopedLevels(const unsigned int cl,
                                   const unsigned fl,
                                   Logger& l)
  : log{&l}, _fl{fl, l}, _cl{cl, l} {}

Logger::ScopedLevels::ScopedLevels(const unsigned int cl, Logger& l)
  : log{&l}, _fl{cl, l}, _cl{cl, l} {}
