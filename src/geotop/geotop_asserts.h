//
// Created by alberto on 5/3/18.
//

#ifndef GEOTOP_GEOTOP_ASSERTS_H
#define GEOTOP_GEOTOP_ASSERTS_H

#include <iostream>
#include <sstream>
#include <string>

/**
 * Validity of pre- and post-conditions and any other requirements must be
 * properly checked. To this aim, Geotop offers a collection of GEO_ASSERT.
 * These assertions (and checks) are performed only when the code is compiled
 * in debug (or equivalent) mode. More precisely, they are skipped if the code
 * is compiled with the -DNDEBUG flag. Assertions are never enough. Put as many
 * assertions as you can without any worry for loss of performance (in release).
 *
 * Example of usage
 *
 * GEO_ASSERT(condition);
 *
 * GEO_ASSERT(condition) << "optional additional message" << std::endl;
 *
 * GEO_ASSERT_IN_RANGE(element,min,max); // check if element belongs to
 * [min,max]
 *
 * GEO_ASSERT_EQ(a,b); // check if a == b, do not use with floating point
 * numbers
 *
 * GEO_ASSERT_LT(a,b); // check if a < b
 *
 * GEO_ASSERT_LE(a,b); // check if a <= b
 *
 * GEO_ASSERT_GT(a,b); // check if a > b
 *
 * GEO_ASSERT_GE(a,b); // check if a >= b
 *
 * All the above, by default, throw std::runtime_error.
 *
 * If you whant to throw your own exception type you can use the macro
 * GEO_ASSERT as follows
 *
 * GEO_ASSERT(condition, exception_type);
 *
 * GEO_ASSERT(condition, exception_type) << "optional" << " message"
 *    << std::endl;
 *
 * The only constraint on the exception_type is that it must have a
 * constructor that takes a const std::string& or a const char * (as
 * the std::exceptions). For example
 *
 * struct my_exception: public std::runtime_error{
 *   using std::runtime_error::runtime_error; // using the same constructors
 *                                            // of the parent
 * };
 *
 * GEO_ASSERT(1>2,my_exception); // it will throw my_exception
 *
 * If you whant/need to create a specific type of assert with all the parameters
 * you want you can do as follows
 *
 * #define GEO_ASSERT_CUSTOM(a,b,c,d) \
 *   _GEO_ASSERT_( ( (a<b) && (c>d) ) || (b == c), std::runtime_error)  \
 *   << "all " << "what you whant to write to debug " << a << " " << b << " " \
 *   << c << " and " << d << std::endl;
 *
 * of course you are free to replace std::runtime_error to any
 * exception you like that can be constructed as explained before.
 *
 * If a condition must be always checked (i.e., also when the code is
 * compiled in release mode), use the GEO_ERROR interface
 *
 * GEO_ERROR(condition);  // throws an std::rutime_erorr
 *
 * GEO_ERROR(condition) << "optional" << " message"
 *    << std::endl;
 *
 * If you need to throw a particular exception, the syntax and the
 * requirements are the same for the assertions explained above.
 *
 * GEO_ERROR(condition, exception_type);
 * GEO_ERROR(condition, exception_type) << "optional" << " message"
 *    << std::endl;
 *
 * The user should use only the above interface. All the rest of this
 * file are technical details and for this reason they are put inside
 * an internal namespace
 */

namespace internal {

/**
 * Used to handle the optional message provided by the user.
 */
class MessageHandler {
public:
  MessageHandler() = default;
  MessageHandler(const MessageHandler &) = delete;

  template <typename T> inline MessageHandler &operator<<(const T &val) {
    _os << val;
    return *this;
  }

  template <typename T> inline MessageHandler &operator<<(T *const &p) {
    if (p == nullptr)
      _os << "nullptr";
    else
      _os << p;
    return *this;
  }

  inline MessageHandler &
  operator<<(std::ostream &(*basic_manipulator)(std::ostream &)) {
    _os << basic_manipulator;
    return *this;
  }

  inline MessageHandler &operator<<(const bool b) {
    return *this << (b ? "true" : "false");
  }

  const std::string get_string() const { return _os.str(); }

private:
  std::ostringstream _os;
};

/**
 * Helper class to manage the construction and throwing of the proper exception
 * type
 */
template <typename ET> struct AssertHelper {
  AssertHelper() = default;
  void operator=(const MessageHandler &m) { throw ET{m.get_string()}; }
};

/**
 * Used like /dev/null for the assertions when compiled in release mode
 */
class NullStream {
public:
  template <typename T> inline NullStream &operator<<(const T &) {
    return *this;
  }

  inline NullStream &operator<<(std::ostream &(*)(std::ostream &)) {
    return *this;
  }
};

} // end namespace internal

// now follows all the macros needed to provide the interface described at the
// beginning of the file

// trick to simulate the overloading of macro GEO_ASSERT depending on the number
// of arguments
#define SELECT_MACRO(_1, _2, NAME, ...) NAME

// when the condition is not satisfied, an exception is thrown

#define GEO_ERROR(...)                                                         \
  SELECT_MACRO(__VA_ARGS__, _GEO_ERROR2, _GEO_ERROR1, dummy)(__VA_ARGS__)

#define _GEO_ERROR2(cond, exception_type)                                      \
  if (!(cond))                                                                 \
    ::internal::AssertHelper<exception_type>() =			\
      internal::MessageHandler()                                               \
      << "\n\n"                                                                \
      << "------------------------------------------------------------------"  \
      << "---"                                                                 \
      << "\n"                                                                  \
      << "A runtime exception has been thrown\n\n"                             \
      << "       file: " << __FILE__ << '\n'                                   \
      << "   function: " << __PRETTY_FUNCTION__ << '\n'                        \
      << "       line: " << __LINE__ << '\n'

#define _GEO_ERROR1(cond) _GEO_ERROR2(cond, std::runtime_error)

#define GEO_ASSERT(...)                                                        \
  SELECT_MACRO(__VA_ARGS__, _GEO_ASSERT2, _GEO_ASSERT1, dummy)(__VA_ARGS__)

#ifndef NDEBUG
#define _GEO_ASSERT_(cond, exception_type) GEO_ERROR(cond, exception_type)
#else
#define _GEO_ASSERT_(cond, exception_type)                                     \
  internal::NullStream {}
#endif

#define _GEO_ASSERT(cond) _GEO_ASSERT_(cond, std::runtime_error)

#define _GEO_ASSERT2(cond, extype)                                             \
  _GEO_ASSERT_(cond, extype) << "  condition: " << #cond << " is not true\n\n"

#define _GEO_ASSERT1(cond) _GEO_ASSERT2(cond, std::runtime_error)

#define GEO_ASSERT_IN_RANGE(position, lower_bound, upper_bound)                \
  _GEO_ASSERT((position >= lower_bound) && (position <= upper_bound))          \
      << "Out of range: " << position << " is not in range [" << lower_bound   \
      << ", " << upper_bound << "]\n\n"

#define GEO_ASSERT_EQ(a, b)                                                    \
  _GEO_ASSERT((a == b)) << a << " was expected to be equal to " << b           \
                        << std::endl

#define GEO_ASSERT_LT(a, b)                                                    \
  _GEO_ASSERT((a < b)) << a << " was expected to be less than " << b           \
                       << std::endl

#define GEO_ASSERT_LE(a, b)                                                    \
  _GEO_ASSERT((a <= b)) << a << " was expected to be less or equal than " << b \
                        << std::endl

#define GEO_ASSERT_GT(a, b)                                                    \
  _GEO_ASSERT((a > b)) << a << " was expected to be greater than " << b        \
                       << std::endl

#define GEO_ASSERT_GE(a, b)                                                    \
  _GEO_ASSERT((a >= b)) << a << " was expected to be greater or equal than "   \
                        << b << std::endl

#define GEO_ERROR_IN_RANGE(position, lower_bound, upper_bound)                 \
  GEO_ERROR((position >= lower_bound) && (position <= upper_bound))            \
      << "Out of range: " << position << " is not in range [" << lower_bound   \
      << ", " << upper_bound << "]\n\n"

#define GEO_ERROR_EQ(a, b)                                                     \
  GEO_ERROR((a == b)) << a << " was expected to be equal to " << b << std::endl

#define GEO_ERROR_LT(a, b)                                                     \
  GEO_ERROR((a < b)) << a << " was expected to be less than " << b << std::endl

#define GEO_ERROR_LE(a, b)                                                     \
  GEO_ERROR((a <= b)) << a << " was expected to be less or equal than " << b   \
                      << std::endl

#define GEO_ERROR_GT(a, b)                                                     \
  GEO_ERROR((a > b)) << a << " was expected to be greater than " << b          \
                     << std::endl

#define GEO_ERROR_GE(a, b)                                                     \
  GEO_ERROR((a >= b)) << a << " was expected to be greater or equal than "     \
                      << b << std::endl
#endif // GEOTOP_GEOTOP_ASSERTS_H
