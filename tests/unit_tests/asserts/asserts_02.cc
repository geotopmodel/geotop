#include<gtest/gtest.h>
#include <geotop_asserts.h>

struct custom_exception : public std::runtime_error{
  using std::runtime_error::runtime_error;
};

TEST(CustomExceptions, 01){
#ifndef NDEBUG
  EXPECT_THROW(GEO_ASSERT(1>2,custom_exception), custom_exception);
#else
  EXPECT_NO_THROW(GEO_ASSERT(1>2,custom_exception));
#endif
}

TEST(CustomExceptions, 02){
  ::testing::internal::CaptureStderr();
try{
  GEO_ASSERT(1>2,custom_exception);
 } catch (custom_exception& e){
  std::cerr << e.what() << std::endl;
 }
#ifndef NDEBUG
 EXPECT_EQ("\n\n---------------------------------------------------------------------\nA runtime exception has been thrown\n\n       file: ../tests/unit_tests/asserts/asserts_02.cc\n   function: virtual void CustomExceptions_02_Test::TestBody()\n       line: 19\n  condition: 1>2 is not true\n\n\n", ::testing::internal::GetCapturedStderr());
#else
 EXPECT_EQ("", ::testing::internal::GetCapturedStderr());
#endif
}

// let us try to define a custom assert

#define GEO_ASSERT_CUSTOM(a,b,c,d) \
  _GEO_ASSERT_( ((a<b) && (c>d)) || (b == c), custom_exception) << "all " \
    << "what you whant to write to debug " << a << " " << b << " " \
    << c << " and " << d << std::endl;



TEST(CustomExceptions, 03){
  ::testing::internal::CaptureStderr();
try{
  GEO_ASSERT_CUSTOM(4,3,1,2);
 } catch (custom_exception& e){
  std::cerr << e.what() << std::endl;
 }
#ifndef NDEBUG
 EXPECT_EQ("\n\n---------------------------------------------------------------------\nA runtime exception has been thrown\n\n       file: ../tests/unit_tests/asserts/asserts_02.cc\n   function: virtual void CustomExceptions_03_Test::TestBody()\n       line: 42\nall what you whant to write to debug 4 3 1 and 2\n\n", ::testing::internal::GetCapturedStderr());
#else
 EXPECT_EQ("", ::testing::internal::GetCapturedStderr());
#endif
}


TEST(CustomExceptions, 04){
#ifndef NDEBUG
  EXPECT_THROW(GEO_ASSERT_CUSTOM(4,3,1,2), custom_exception);
#else
  EXPECT_NO_THROW(GEO_ASSERT(1>2,custom_exception));
#endif
}
