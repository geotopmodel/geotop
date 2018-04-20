//
// Created by alberto on 4/18/18.
//

#ifndef GEOTOP_LOGGER_H
#define GEOTOP_LOGGER_H

#include <iostream>
#include <string>
#include <stack>
#include <sstream>
#include <iomanip>
#include <functional>

class Logger {
 public:
  /**
   * Utility class to add and remove a prefix to a Logger.
   * In the constructor it pushes the string to the Logger.
   * In the destructor it pops the just added prefix.
   *
   * The idea is that once you enter a scope (either a function, a part of a
   * function, etc.) you can create one ScopedPrefix object and you can forget
   * to call pop() since it automatically does.
   */
  class ScopedPrefix;

  /**
   * Helper class to set the console level of a logger inside a scope to a given
   * value. For example, this can be useful to either completely disable the
   * logging inside a scope or to enable a finer logging.
   */
  class ScopedConsoleLevel;

  /**
   * Helper class to set the file level of a logger inside a scope to a given
   * value. For example, this can be useful to either completely disable the
   * logging inside a scope or to enable a finer logging.
   */
  class ScopedFileLevel;

  /**
   * Adds "geotop" as first prefix. By default it prints to stdout and no file.
   * The depth levels for both console and file are not limited.
   */
  Logger();

  /**
   * Add @param s to the stack of prefixes. Note that each element of the stack
   * contains all the previous elements If first element is "first" and the
   * second element is "second", the prefix associated to the second is
   * "first:second:". While the prefix of the first element if "first:"
   */
  void push(const std::string& s);

  /**
   * Remove the last-added  and returns it
   */
  std::string pop();

  /**
   * @return the last prefix added
   */
  const std::string& prefix() const;

  /**
   * this is necessary to handle functions like
   * std::endl and std::flush
   *
   * Note that if you need to use std::setw, std::setprecision and the like,
   * you should adopt the following aproach
   *
   * std::ostringstream os;
   * os << std::setw(10) << var << std::setprecision(99) << othervar <<
   * std::endl; geolog << os.str();
   *
   * otherwise geolog will not see the effect of std::setw, std::setprecision
   * etc.
   */
  Logger& operator<<(std::ostream& (*p)(std::ostream&)) {
    std::ostringstream os;
    os << p;
    return *this << os.str();
  }

  /**
   * Attach @param o to the file stream @var ofile
   */
  void attach_file_stream(std::ostream& o);

  /**
   * Detach file stream. You may want to close the attached stream
   */
  void detach_file_stream();

  /**
   * return the depth level of the current prefix
   */
  unsigned int prefix_depth_level() const { return prefixes.size(); }

  unsigned int console_level() const { return _console_level; }

  void set_console_level(const unsigned int console_level) {
    _console_level = console_level;
  }

  unsigned int file_level() const { return _file_level; }

  void set_file_level(const unsigned int file_level) {
    _file_level = file_level;
  }

 private:
  /**
   * stack of prefixes
   */
  std::stack<std::string> prefixes;

  /**
   * pointer to a console channel. By default it points to std::cout
   */
  std::ostream* std_out;

  /**
   * pointer to a second channel. It makes sense to use it to write to a file
   */
  std::ostream* ofile;

  /**
   * If true, the last prefix is prepended to the string to print
   */
  bool _at_new_line;

  unsigned int _console_level;
  unsigned int _file_level;

  template <typename T>
  friend Logger& operator<<(Logger&, const T&);
};

template <typename T>
inline Logger& operator<<(Logger& l, const T& t) {
  std::ostringstream os;
  if (l._at_new_line) os << l.prefix();
  os << t;
  std::string s{os.str()};
  if (l.std_out && (l.prefix_depth_level() <= l.console_level()))
    *(l.std_out) << s;
  if (l.ofile && (l.prefix_depth_level() <= l.file_level())) *(l.ofile) << s;
  l._at_new_line = (s == "\n");
  return l;
}

class Logger::ScopedPrefix {
 public:
  /**
   * Push string @param s to logger @param l
   */
  ScopedPrefix(const std::string& s, Logger& l);

  /**
   * Push string @param s to the default logger geolog
   */
  ScopedPrefix(const std::string& s);

  /**
   * Call log->pop();
   */
  virtual ~ScopedPrefix();

 private:
  Logger* log;
};

class ScopedLevel {
public:
  ScopedLevel(std::function<unsigned int()> &logger_get_level,
                std::function<void(const unsigned int)> &logger_set_level, const unsigned int level,
                Logger &l);

private:
  std::function<void(const unsigned int)> set_level;
  std::function<unsigned int()> get_level;
  Logger *log;
  const unsigned int old_level;
};

class Logger::ScopedConsoleLevel{
 public:
  /**
   * set the console level of @param l to @param level
   */
  ScopedConsoleLevel(const unsigned int level, Logger& l);

  /**
   * set the console level of the default logger geolog to @param level
   */
  ScopedConsoleLevel(const unsigned int level);

  /**
   * restore the console level of the logger pointed to by @param log to it's
   * previous value
   */
  ~ScopedConsoleLevel();

 private:
  Logger* log;
  const unsigned int old_level;
};

class Logger::ScopedFileLevel {
public:
  /**
   * set the file level of @param l to @param level
   */
  ScopedFileLevel(const unsigned int level, Logger& l);

  /**
   * set the file level of the default logger geolog to @param level
   */
  ScopedFileLevel(const unsigned int level);

  /**
   * restore the file level of the logger pointed to by @param log to it's
   * previous value
   */
  ~ScopedFileLevel();

private:
  Logger* log;
  const unsigned int old_level;
};


/** the default logger */
extern Logger geolog;

#endif  // GEOTOP_LOGGER_H
