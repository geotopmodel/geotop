//
// Created by alberto on 4/18/18.
//

#ifndef GEOTOP_LOGGER_H
#define GEOTOP_LOGGER_H


#include <iostream>
#include <string>
#include <stack>
#include <sstream>

class Logger :public std::streambuf{
public:

  /**
   * Utility class to add and remove a prefix to a Logger.
   * In the constructor it pushes the string to the Logger.
   * In the destructor it pops the just added prefix.
   *
   * The idea is that once you enter a scope (either a function, a part of a function, etc.)
   * you can create one Prefix object and you can forget to call pop() since it automatically does.
   */
  class Prefix {
  public:

    /**
     * Push string @param s to logger @param l
     */
    Prefix(const std::string &s, Logger &l);

    /**
     * Push string @param s to the default logger geolog
     */
    Prefix(const std::string &s);

    /**
     * Call log->pop();
     */
    virtual ~Prefix();

  private:
    Logger *log;
  };

  /**
   * Adds "geotop" as first prefix. By default it prints to stdout and no file.
   */
  Logger();

  /**
   * Add @param s to the stack of prefixes
   */
  void push(const std::string &s);

  /**
   * Remove the last-added  and returns it
   */
  std::string pop();

  const std::string &prefix() const;

private:
  std::stack<std::string> prefixes;

  std::ostream *std_out;

  std::ostream *ofile;
  bool _at_new_line;
public:
  template <typename T> friend Logger& operator<<(Logger&, const T&);

  Logger& operator<<(std::ostream& (*p) (std::ostream &)){
    std::ostringstream os;
    os<<p;
    *this<<os.str();
    return *this;
  }
};

template <typename T>
inline
Logger& operator<<(Logger& l, const T& t){
  std::string pre;
  if(l._at_new_line) pre = l.prefix();
  std::ostringstream os;
  os << pre << t;
  if(l.std_out)
    *(l.std_out) << os.str();
  if(l.ofile)
    *(l.ofile) << os.str();
//  if (os.str().find('\n') != std::string::npos)
  if (os.str() == "\n")
    l._at_new_line = true;
  else
    l._at_new_line = false;
  return l;
}


extern Logger geolog;

#endif //GEOTOP_LOGGER_H
