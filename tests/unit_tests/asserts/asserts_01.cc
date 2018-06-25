#include<gtest/gtest.h>
#include <geotop_asserts.h>

TEST(Assertions, 01){
  EXPECT_NO_THROW(GEO_ASSERT(1<2));
#ifndef NDEBUG
  EXPECT_THROW(GEO_ASSERT(1>2), std::runtime_error);
#else
  EXPECT_NO_THROW(GEO_ASSERT(1<2));
#endif
}

TEST(Assertions, 02){
  testing::internal::CaptureStderr();
  try
    {
      int a = 1;
      int b = 2;
      GEO_ASSERT(a>b) << "bla " << "blas " << std::endl << true << " lapack";
      (void)a;
      (void)b;
    }
  catch(std::exception& e){
    std::cerr << e.what() << std::endl;
  }
#ifndef NDEBUG
  EXPECT_EQ("\n\n---------------------------------------------------------------------\nA runtime exception has been thrown\n\n       file: ../tests/unit_tests/asserts/asserts_01.cc\n   function: virtual void Assertions_02_Test::TestBody()\n       line: 19\n  condition: a>b is not true\n\nbla blas \ntrue lapack\n", testing::internal::GetCapturedStderr());
#else
  EXPECT_EQ("", testing::internal::GetCapturedStderr());
#endif
}

TEST(Assertions, assert_in_range){
  int pos;
  int min = 2;
  int max = 10;

  pos = max + 8;
#ifndef NDEBUG
  EXPECT_THROW(GEO_ASSERT_IN_RANGE(pos, min, max), std::runtime_error);
#else
  EXPECT_NO_THROW(GEO_ASSERT_IN_RANGE(pos, min, max));
#endif
  pos = (max + min)/2 ;
  EXPECT_NO_THROW(GEO_ASSERT_IN_RANGE(pos, min, max));
}

TEST(Assertions, equality){
#ifndef NDEBUG
  EXPECT_THROW(GEO_ASSERT_EQ(2,3), std::runtime_error);
#else
  EXPECT_NO_THROW(GEO_ASSERT_EQ(2,3));
#endif
  EXPECT_NO_THROW(GEO_ASSERT_EQ(7,7));
}

TEST(Assertions, less){
#ifndef NDEBUG
  EXPECT_THROW(GEO_ASSERT_LT(3,2), std::runtime_error);
  EXPECT_THROW(GEO_ASSERT_LE(3,2), std::runtime_error);
#else
  EXPECT_NO_THROW(GEO_ASSERT_LT(3,2));
  EXPECT_NO_THROW(GEO_ASSERT_LE(3,2));
#endif
  EXPECT_NO_THROW(GEO_ASSERT_LT(2,3));
  EXPECT_NO_THROW(GEO_ASSERT_LE(2,2));
}

TEST(Assertions, greater){
#ifndef NDEBUG
  EXPECT_THROW(GEO_ASSERT_GT(2,5), std::runtime_error);
  EXPECT_THROW(GEO_ASSERT_GE(2,5), std::runtime_error);
#else
  EXPECT_NO_THROW(GEO_ASSERT_GT(2,5));
  EXPECT_NO_THROW(GEO_ASSERT_GE(2,5));
#endif
  EXPECT_NO_THROW(GEO_ASSERT_GT(5,2));
  EXPECT_NO_THROW(GEO_ASSERT_GE(2,2));
}
