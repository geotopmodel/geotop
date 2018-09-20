#include <gtest/gtest.h>
#include <turtle.h>
#include <tensor.h>
#include <iterator> // std::prev

// tensor(d,r,c) => d = depth; r=row; c=column
TEST(Tensor, access_element){
  // DOUBLETENSOR pa{2,2,2};
  Tensor<int> t{2,2,2}; //  2x2x2 tensor

  double c{0.0};
  for (auto &x : t)
    x = ++c;
  
  EXPECT_EQ( 1.0, t(1,1,1) );
  EXPECT_EQ( 2.0, t(1,1,2) );
  EXPECT_EQ( 3.0, t(1,2,1) );
  EXPECT_EQ( 4.0, t(1,2,2) );
  
  EXPECT_DOUBLE_EQ( 5.0, t(2,1,1) );
  EXPECT_DOUBLE_EQ( 6.0, t(2,1,2) );
  EXPECT_DOUBLE_EQ( 7.0, t(2,2,1) );
  EXPECT_DOUBLE_EQ( 8.0, t(2,2,2) );
}
