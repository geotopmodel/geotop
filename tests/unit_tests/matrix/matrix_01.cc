#include <gtest/gtest.h>
#include <matrix.h>
#include <iterator> // std::prev

TEST(Matrix, constructors){
  Matrix<int> m(3,1,5,2); // 3x4
  EXPECT_EQ( std::size_t{3}, m.n_row );
  EXPECT_EQ( std::size_t{4}, m.n_col );
}

TEST(Matrix, initialization){
  Matrix<int> m(3,2,3,2); // 2x2
  EXPECT_EQ( 0, m(0,0));
  EXPECT_EQ( 0, m(0,1));
  EXPECT_EQ( 0, m(1,0));
  EXPECT_EQ( 0, m(1,1));
}

TEST(Matrix, begin){
  Matrix<int> m(3,2,3,2); // 2x2
  m(0,0) = 1;
  m(0,1) = 2;
  m(1,0) = 3;
  m(1,1) = 4;
  EXPECT_EQ( 1, *(m.begin()) );
}

TEST(Matrix, end){
  Matrix<int> m(3,2,3,2); // 2x2
  m(0,0) = 1;
  m(0,1) = 2;
  m(1,0) = 3;
  m(1,1) = 4;
  EXPECT_EQ( 4, *std::prev(m.end()) );
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



// TEST(Matrix, end){
//   Matrix<double> m(2,1,2,1); // 2x2 matrix
//   m[0][0] = 1;
//   m[0][1] = 2;
//   m[1][0] = 3;
//   m[1][1] = 4;

//   EXPECT_EQ(*(m.end()), 0);
// }



// TEST(Matrix, range_for){
//   Matrix<double> m(3,3);
//   double c{1.0};
//   m[0][0] = -9999;
//   v[3][3] = -9999;

//   for (auto &x : m)
//     x = ++c;

//   EXPECT_DOUBLE_EQ(m[0][0], -9999);
//   EXPECT_DOUBLE_EQ(m[0][1], 2);
//   EXPECT_DOUBLE_EQ(m[0][2], 3);
//   EXPECT_DOUBLE_EQ(m[0][3], 4);

// }

// TEST(Matrix, range_for_zero){
//   Vector<double> v(3,0);
//   double c{0.0};
//   v[0] = -9999;
//   v[3] = -9999;
  
//   for (auto &x : v)
//     x = ++c;

//   EXPECT_DOUBLE_EQ(v[0], 1);
//   EXPECT_DOUBLE_EQ(v[1], 2);
//   EXPECT_DOUBLE_EQ(v[2], 3);
//   EXPECT_DOUBLE_EQ(v[3], 4);

// }

// TEST(Matrix, set_value){
//   Vector<double> v(3,0.);
//   double c{0.0};
  
//   for (auto &x : v)
//     x = ++c;

//   v = -9999.;

//   EXPECT_DOUBLE_EQ(v[0], -9999.);
//   EXPECT_DOUBLE_EQ(v[1], -9999.);
//   EXPECT_DOUBLE_EQ(v[2], -9999.);
//   EXPECT_DOUBLE_EQ(v[3], -9999.);
// }



// TEST(Matrix, copy_semantic){
//   Vector<double> v(3,0);
//   double c{0.0};
//   v[0] = -9999.;
//   v[3] = -9999.;
  
//   for (auto &x : v)
//     x = ++c;
  
//   Vector<double> v1{v};

//   v = -9999.0; // change v to check the vectors are not linked
//   EXPECT_DOUBLE_EQ(v1[0], 1);
//   EXPECT_DOUBLE_EQ(v1[1], 2);
//   EXPECT_DOUBLE_EQ(v1[2], 3);
//   EXPECT_DOUBLE_EQ(v1[3], 4);

//   Vector<double> v2{1};
//   v2[1] = -9999;

//   for (auto &x : v)
//     x = ++c;

//   v2=v1;
  
//   v = -9999.0; // change v1 to check the vectors are not linked

//   ASSERT_EQ(v1.size(), v2.size());
//   EXPECT_DOUBLE_EQ(v2[0], 1);
//   EXPECT_DOUBLE_EQ(v2[1], 2);
//   EXPECT_DOUBLE_EQ(v2[2], 3);
//   EXPECT_DOUBLE_EQ(v2[3], 4);
  
// }


// TEST(Matrix, out_of_range){
//   Matrix<double> v(3);
//   double c{1.0};
  
//   for (auto &x : v)
//     x = ++c;

//   EXPECT_NO_THROW(v.at(1));
//   EXPECT_NO_THROW(v.at(2));
//   EXPECT_NO_THROW(v.at(3));
  
//   EXPECT_ANY_THROW(v.at(0));
//   EXPECT_ANY_THROW(v.at(4));
//   EXPECT_ANY_THROW(v.at(-3));
//   EXPECT_ANY_THROW(v.at(5000));


//   EXPECT_NO_THROW(v(1));
//   EXPECT_NO_THROW(v(2));
//   EXPECT_NO_THROW(v(3));

// #ifndef NDEBUG
//   EXPECT_ANY_THROW(v(0));
//   EXPECT_ANY_THROW(v(4));
//   EXPECT_ANY_THROW(v(-3));
//   EXPECT_ANY_THROW(v(5000));
// #else
//   EXPECT_NO_THROW(v(0));
//   EXPECT_NO_THROW(v(4));
//   EXPECT_NO_THROW(v(-3));
//   EXPECT_NO_THROW(v(5000));
// #endif
  
// }

// TEST(Matrix, out_of_range_zero){
//   Vector<double> v(3,0);
//   double c{1.0};
  
//   for (auto &x : v)
//     x = ++c;

//   EXPECT_NO_THROW(v.at(0));
//   EXPECT_NO_THROW(v.at(1));
//   EXPECT_NO_THROW(v.at(2));
//   EXPECT_NO_THROW(v.at(3));
  
//   EXPECT_ANY_THROW(v.at(-3));
//   EXPECT_ANY_THROW(v.at(4));
//   EXPECT_ANY_THROW(v.at(5000));
  
// }

// TEST(Matrix, summation){
//   Matrix<double> v_0(3);
  
//   Matrix<double> v_1(3,0);

//   // vector lenght mismatch
// #ifndef NDEBUG
//   EXPECT_ANY_THROW(v_1 += v_0);
// #else
//   EXPECT_NO_THROW(v_1 += v_0);
// #endif
  
//   Matrix<double> v_2(3,0.);
//   double c{0.0};  
//   for (auto &x : v_1)
//     x = ++c;

//   for (auto &x : v_2)
//     x = ++c;

//   EXPECT_DOUBLE_EQ(v_2[0], 5.);
//   EXPECT_DOUBLE_EQ(v_2[1], 6.);
//   EXPECT_DOUBLE_EQ(v_2[2], 7.);
//   EXPECT_DOUBLE_EQ(v_2[3], 8.);

//   EXPECT_NO_THROW(v_2 += v_1);
  
//   EXPECT_DOUBLE_EQ(v_1[0], 1.);
//   EXPECT_DOUBLE_EQ(v_1[1], 2.);
//   EXPECT_DOUBLE_EQ(v_1[2], 3.);
//   EXPECT_DOUBLE_EQ(v_1[3], 4.);

//   EXPECT_DOUBLE_EQ(v_2[0], 6.);
//   EXPECT_DOUBLE_EQ(v_2[1], 8.);
//   EXPECT_DOUBLE_EQ(v_2[2], 10.);
//   EXPECT_DOUBLE_EQ(v_2[3], 12.);
// }
