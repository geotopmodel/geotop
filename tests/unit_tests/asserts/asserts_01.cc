#include<gtest/gtest.h>
#include <geotop_asserts.h>

TEST(MessageHandler, 01){
  internal::MessageHandler m{};
  EXPECT_EQ("",m.get_string());
}

TEST(MessageHandler, 02){
  internal::MessageHandler m{};
  m << "ok";
  EXPECT_EQ("ok",m.get_string());
}

TEST(MessageHandler, 03){
  internal::MessageHandler m{};
  int * p = nullptr;
  m << p;
  EXPECT_EQ("nullptr",m.get_string());
}

TEST(MessageHandler, 04){
  internal::MessageHandler m{};
  int * p = NULL;
  m << p;
  EXPECT_EQ("nullptr",m.get_string());
}

TEST(Assertions, 01){
  EXPECT_NO_THROW(GEO_ASSERT(1<2));
#ifndef NDEBUG
  EXPECT_ANY_THROW(GEO_ASSERT(1>2));
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
  EXPECT_EQ("\n\n---------------------------------------------------------------------\nA runtime exception has been thrown\n\n       file: ../tests/unit_tests/asserts/asserts_01.cc\n   function: virtual void Assertions_02_Test::TestBody()\n       line: 44\n  condition: a>b is not true\n\nbla blas \ntrue lapack\n", testing::internal::GetCapturedStderr());
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
  EXPECT_ANY_THROW(GEO_ASSERT_IN_RANGE(pos, min, max));
#else
  EXPECT_NO_THROW(GEO_ASSERT_IN_RANGE(pos, min, max));
#endif
  
  pos = (max + min)/2 ;
  EXPECT_NO_THROW(GEO_ASSERT_IN_RANGE(pos, min, max));
}

TEST(NullStream, test_1e8_assertions){
  size_t volatile i ;
  for(i=0; i < 1e8; ++ i )
    internal::NullStream{} << "ciao " << true << std::endl;
}
