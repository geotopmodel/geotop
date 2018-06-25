#include<gtest/gtest.h>
#include <geotop_asserts.h>

// NullStream does nothing. We check that doing nothing is fast
TEST(NullStream, test_1e6_assertions){
  size_t volatile i ;
  for(i=0; i < 1e6; ++ i )
    internal::NullStream{} << "this should never be printed " << true << std::endl;
}


