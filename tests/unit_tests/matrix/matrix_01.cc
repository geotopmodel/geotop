#include <gtest/gtest.h>
#include <matrix.h>
#include <iterator> // std::prev

TEST(Matrix, constructor_4args){
  Matrix<int> m{3,1,5,2}; // 3x4 matrix
  EXPECT_EQ( std::size_t{3}, m.n_row );
  EXPECT_EQ( std::size_t{4}, m.n_col );
}

TEST(Matrix, constructor_2args){
  Matrix<int> m{3,5}; // 3x5 matrix
  EXPECT_EQ( std::size_t{3}, m.n_row );
  EXPECT_EQ( std::size_t{5}, m.n_col );
}

TEST(Matrix, initialization){
  Matrix<int> m{2,2}; // 2x2 matrix
  EXPECT_EQ( 0, m(1,1));
  EXPECT_EQ( 0, m(1,2));
  EXPECT_EQ( 0, m(2,1));
  EXPECT_EQ( 0, m(2,2));
  // EXPECT_EQ( 0, m(0,0)); // => expected to fail => OK
  // EXPECT_EQ( 0, m(3,3)); // => expected to fail => OK
}

TEST(Matrix, begin){
  Matrix<int> m{2,2}; // 2x2 matrix
  m(1,1) = 1;
  m(1,2) = 2;
  m(2,1) = 3;
  m(2,2) = 4;
  EXPECT_EQ( 1, *(m.begin()) );
}

TEST(Matrix, end){
  Matrix<int> m{2,2}; // 2x2 matrix
  m(1,1) = 1;
  m(1,2) = 2;
  m(2,1) = 3;
  m(2,2) = 4;
  EXPECT_EQ( 4, *std::prev(m.end()) );
}

TEST(Matrix, range_for){
  Matrix<double> m{2,2}; // 2x2 matrix
  double c{1.0};

  for (auto &x : m)
    x = ++c;

  EXPECT_DOUBLE_EQ( 2.0, m(1,1) );
  EXPECT_DOUBLE_EQ( 3.0, m(1,2) );
  EXPECT_DOUBLE_EQ( 4.0, m(2,1) );
  EXPECT_DOUBLE_EQ( 5.0, m(2,2) );

}

// TEST(Matrix, end){
//   Matrix<int> m(2,1,2,1); // 2x2 matrix
//   int c=0;
//   for (int i=0; i <3 ; ++i)
//     for (int j=0; j<3; ++j)
//       m[i][j]=c++;
//   testing::internal::CaptureStdout();
//   for (auto x : m)
//     std::cout << x << std::endl;
//   EXPECT_EQ("", testing::internal::GetCapturedStdout());
// }
