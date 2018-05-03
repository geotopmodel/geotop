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
 * GEO_ASSERT_EQ(a,b); // check if a == b
 * GEO_ASSERT_LT(a,b); // check if a < b
 * GEO_ASSERT_LE(a,b); // check if a <= b
 * GEO_ASSERT_GT(a,b); // check if a > b
 * GEO_ASSERT_GE(a,b); // check if a >= b
 *
 * The user should use only the above interface. All the rest of this file are technical details
 * and for this reason they are put inside an internal namespace
 */

namespace internal {

  /**
   * Used like /dev/null for the assertions when compiled in release mode
   */
  class NullStream{
  public:
    template<typename T>
    inline NullStream& operator<<(const T&) {return *this;}

    inline NullStream& operator<<(std::ostream& (*) (std::ostream&)) {return *this;}
  };

  /**
   * Used to handle the optional message provided by the user
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


  struct AssertHelper{
    AssertHelper() = default;
    void operator=(const MessageHandler& m){
      _os << m.get_string();
      throw std::runtime_error{_os.str()};
    }
    std::ostringstream _os;
  };

#ifndef NDEBUG
#define GEO_ASSERT_(cond)						\
if ( !(cond) ) \
    ::internal::AssertHelper() = internal::MessageHandler() << "\n\n"\
  << "------------------------------------------------------------------"\
  << "---"\
  << "\n"\
  << "A runtime exception has been thrown\n\n"\
  << "       file: " << __FILE__ << '\n'\
  << "   function: " << __PRETTY_FUNCTION__ << '\n'\
  << "       line: " << __LINE__ << "\n"
#else
#define GEO_ASSERT_(cond) internal::NullStream{}
#endif

} // end namespace internal

  
#define GEO_ASSERT(cond)						\
  GEO_ASSERT_(cond) << "  condition: " << #cond << " is not true\n\n"
		

#define GEO_ASSERT_IN_RANGE(position,lower_bound,upper_bound)           \
  GEO_ASSERT_( (position >= lower_bound) && (position <= upper_bound) ) \
  << "Out of range: "							\
  << position << " is not in range [" << lower_bound <<", "		\
  << upper_bound << "]\n\n"

#define GEO_ASSERT_EQ(a,b) \
  GEO_ASSERT_( (a == b) ) << a << " is not equal to " << b << std::endl;

#define GEO_ASSERT_LT(a,b) \
  GEO_ASSERT_( (a < b) ) << a << " was expected to be less than " << b << std::endl

#define GEO_ASSERT_LE(a,b) \
  GEO_ASSERT_( (a <= b) ) << a << " was expected to be less or equal than " << b << std::endl;

#define GEO_ASSERT_GT(a,b) \
  GEO_ASSERT_( (a > b) ) << a << " was expected to be greater than " << b << std::endl

#define GEO_ASSERT_GE(a,b) \
  GEO_ASSERT_( (a >= b) ) << a << " was expected to be greater or equal than " << b << std::endl;


#endif //GEOTOP_GEOTOP_ASSERTS_H
