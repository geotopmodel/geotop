//
// Created by alberto on 5/3/18.
//

#ifndef GEOTOP_GEOTOP_ASSERTS_H
#define GEOTOP_GEOTOP_ASSERTS_H

#include <string>
#include <iostream>
#include <sstream>

/**
 * example of usage
 * GEO_ASSERT(condition);
 * GEO_ASSERT(condition) << "optional additional message" << std::endl;
 * GEO_ASSERT_IN_RANGE(element,min,max); // check if element belongs to [min,max]
 * GEO_ASSERT_EQ(a,b); // check if a == b, do not use with floating point numbers
 * GEO_ASSERT_LT(a,b); // check if a < b
 * GEO_ASSERT_LE(a,b); // check if a <= b
 * GEO_ASSERT_GT(a,b); // check if a > b
 * GEO_ASSERT_GE(a,b); // check if a >= b
 *
 * All the above, by default, throw std::runtime_error.
 *
 * If you whant to throw your own exception type you can use the macro GEO_ASSERT as follows
 *
 * GEO_ASSERT(condition, exception_type);
 * GEO_ASSERT(condition, exception_type) << "optional" << " message" << std::endl;
 *
 * The only constraint on the exception_type is that it must have a constructor that takes
 * a const std::string& or a const char * (as the std::exceptions). For example
 *
 * struct my_exception: public std::runtime_error{
 *   using std::runtime_error::runtime_error; // using the same constructors of the parent
 * };
 *
 * GEO_ASSERT(1>2,my_exception) << "example" << std::endl; // it will throw my_exception
 *
 * If you whant/need to create a specific type of assert with all the parameters you want
 * you can do as follows
 *
 * #define GEO_ASSERT_CUSTOM(a,b,c,d) \
 *   _GEO_ASSERT_( ( (a<b) && (c>d) ) || (b == c), std::runtime_error) << "all " \
 *   << "what you whant to write to debug " << a << " " << b << " " \
 *   << c << " and " << d << std::endl;
 *
 * of course you are free to replace std::runtime_error to any exception you like that can
 * be constructed as explained before.
 *
 * The user should use only the above interface. All the rest of this file are technical details
 * and for this reason they are put inside an internal namespace
 */

namespace internal {

  /**
   * Used to handle the optional message provided by the user.
   */
  class MessageHandler{
  public:
    MessageHandler() = default;
    MessageHandler(const MessageHandler&) = delete;

    template <typename T>
    inline MessageHandler& operator<< (const T& val){ _os << val; return *this;}

    template <typename  T>
    inline MessageHandler& operator<<(T* const& p){ if (p == nullptr) _os << "nullptr" ; else _os << p;
      return *this; }

    inline MessageHandler& operator<< (std::ostream& (*basic_manipulator)(std::ostream& )) { _os << basic_manipulator; return *this; }

    inline MessageHandler& operator<< (const bool b) {return *this << (b ? "true" : "false"); }

    const std::string get_string() const { return _os.str(); }
  private:
    std::ostringstream _os;
  };

  /**
   * Helper class to manage the construction and throwing of the proper exception type
   */
template <typename ET>
  struct AssertHelper{
    AssertHelper() = default;
    void operator=(const MessageHandler& m){
      throw ET(m.get_string());
    }
  };

  /**
   * Used like /dev/null for the assertions when compiled in release mode
   */
  class NullStream{
  public:
    template<typename T>
    inline NullStream& operator<<(const T&) {return *this;}

    inline NullStream& operator<<(std::ostream& (*) (std::ostream&)) {return *this;}
  };


} // end namespace internal

// now follows all the macros needed to provide the interface described at the beginning of the file

// trick to simulate the overloading of macro GEO_ASSERT depending on the number of arguments
#define SELECT_MACRO(_1,_2,NAME, ...) NAME

#define GEO_ASSERT(...) SELECT_MACRO(__VA_ARGS__, GEO_ASSERT2, GEO_ASSERT1, dummy) (__VA_ARGS__)

#ifndef NDEBUG
#define _GEO_ASSERT_(cond,exception_type)						\
if ( !(cond) ) \
    ::internal::AssertHelper<exception_type>() = internal::MessageHandler() << "\n\n"\
  << "------------------------------------------------------------------"\
  << "---"\
  << "\n"\
  << "A runtime exception has been thrown\n\n"\
  << "       file: " << __FILE__ << '\n'\
  << "   function: " << __PRETTY_FUNCTION__ << '\n'\
  << "       line: " << __LINE__ << "\n"
#else
#define _GEO_ASSERT_(cond, exception_type) internal::NullStream{}
#endif


#define GEO_ASSERT_(cond) _GEO_ASSERT_(cond,std::runtime_error)


#define GEO_ASSERT1(cond)						\
  GEO_ASSERT_(cond) << "  condition: " << #cond << " is not true\n\n"

#define GEO_ASSERT2(cond,extype)						\
  _GEO_ASSERT_(cond,extype) << "  condition: " << #cond << " is not true\n\n"

#define GEO_ASSERT_IN_RANGE(position,lower_bound,upper_bound)           \
  GEO_ASSERT_( (position >= lower_bound) && (position <= upper_bound) ) \
  << "Out of range: "							\
  << position << " is not in range [" << lower_bound <<", "		\
  << upper_bound << "]\n\n"

#define GEO_ASSERT_EQ(a,b) \
  GEO_ASSERT_( (a == b) ) << a << " is not equal to " << b << std::endl

#define GEO_ASSERT_LT(a,b) \
  GEO_ASSERT_( (a < b) ) << a << " was expected to be less than " << b << std::endl

#define GEO_ASSERT_LE(a,b) \
  GEO_ASSERT_( (a <= b) ) << a << " was expected to be less or equal than " << b << std::endl

#define GEO_ASSERT_GT(a,b) \
  GEO_ASSERT_( (a > b) ) << a << " was expected to be greater than " << b << std::endl

#define GEO_ASSERT_GE(a,b) \
  GEO_ASSERT_( (a >= b) ) << a << " was expected to be greater or equal than " << b << std::endl


#endif //GEOTOP_GEOTOP_ASSERTS_H
