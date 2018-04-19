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

class Logger{
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
   * Add @param s to the stack of prefixes. Note that each element of the stack contains all the previous elements
   * If first element is "first" and the second element is "second", the prefix associated to the second is
   * "first:second:". While the prefix of the first element if "first:"
   */
  void push(const std::string &s);

  /**
   * Remove the last-added  and returns it
   */
  std::string pop();

  /**
   * @return the last prefix added
   */
  const std::string &prefix() const;


  /**
   * this is necessary to handle functions like
   * std::endl and std::flush
   *
   * Note that if you need to use std::setw, std::setprecision and the like,
   * you should adopt the following aproach
   *
   * std::ostringstream os;
   * os << std::setw(10) << var << std::setprecision(99) << othervar << std::endl;
   * geolog << os.str();
   *
   * otherwise geolog will not see the effect of std::setw, std::setprecision etc.
   */
  Logger& operator<<(std::ostream& (*p) (std::ostream &)) {
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
private:

  /**
   * stack of prefixes
   */
  std::stack<std::string> prefixes;

  /**
   * pointer to a console channel. By default it points to std::cout
   */
  std::ostream *std_out;

  /**
   * pointer to a second channel. It makes sense to use it to write to a file
   */
  std::ostream *ofile;

  /**
   * If true, the last prefix is prepended to the string to print
   */
  bool _at_new_line;

public:
  template <typename T> friend Logger& operator<<(Logger&, const T&);
};

template <typename T>
inline
Logger& operator<<(Logger& l, const T& t){
  std::string pre;
  if(l._at_new_line) pre = l.prefix();
  std::ostringstream os;
  os << pre << t;
  std::string s{os.str()};
  if(l.std_out)
    *(l.std_out) << s;
  if(l.ofile)
    *(l.ofile) << s;
  l._at_new_line = (s == "\n");
  return l;
}


extern Logger geolog;

#endif //GEOTOP_LOGGER_H
