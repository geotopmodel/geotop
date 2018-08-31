#include <gtest/gtest.h>
#include <matrix.h>

  TEST(Matrix, SquareBrackets){
  Matrix<int> m{2,2}; // 2x2 matrix
  int c{0};
   
  for (auto &x : m)
    x = ++c;

  EXPECT_EQ( 1, m.row(1)[1] );
  EXPECT_EQ( 2, m.row(1)[2] );
  EXPECT_EQ( 3, m.row(2)[1] );
  EXPECT_EQ( 4, m.row(2)[2] );
 
}

  TEST(Matrix, RoundBrackets){
  Matrix<int> m{2,2}; // 2x2 matrix
  int c{0};
   
  for (auto &x : m)
    x = ++c;

  EXPECT_EQ( 1, m.row(1)(1) );
  EXPECT_EQ( 2, m.row(1)(2) );
  EXPECT_EQ( 3, m.row(2)(1) );
  EXPECT_EQ( 4, m.row(2)(2) );
 
}
