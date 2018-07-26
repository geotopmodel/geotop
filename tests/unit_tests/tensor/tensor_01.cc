#include <gtest/gtest.h>
#include <tensor.h>
#include <iterator> // std::prev

TEST(Tensor, constructor_6args){
  Tensor<int> m{3,1,3,0,2,0}; // 3x4x3 tensor
  EXPECT_EQ( std::size_t{3}, m.n_row );
  EXPECT_EQ( std::size_t{4}, m.n_col );
  EXPECT_EQ( std::size_t{3}, m.n_dep );
}

TEST(Tensor, constructor_3args){
  Tensor<int> m{3,5,2}; // 3x5x2 tensor
  EXPECT_EQ( std::size_t{3}, m.n_row );
  EXPECT_EQ( std::size_t{5}, m.n_col );
  EXPECT_EQ( std::size_t{2}, m.n_dep );
}

TEST(Tensor, initialization){
  Tensor<int> m{2,2,2}; // 2x2x2 tensor
  EXPECT_EQ( 0, m(1,1,1) );
  EXPECT_EQ( 0, m(1,2,1) );
  EXPECT_EQ( 0, m(2,1,1) );
  EXPECT_EQ( 0, m(2,2,1) );
  
  EXPECT_EQ( 0, m(1,1,2) );
  EXPECT_EQ( 0, m(1,2,2) );
  EXPECT_EQ( 0, m(2,1,2) );
  EXPECT_EQ( 0, m(2,2,2) );
  // EXPECT_EQ( 0, m(0,0,0)); // => expected to fail => OK
  // EXPECT_EQ( 0, m(3,3,3)); // => expected to fail => OK
}

TEST(Tensor, begin){
  Tensor<int> m{2,2,2}; // 2x2x2 tensor
  m(1,1,1) = 1;
  m(1,2,1) = 2;
  m(2,1,1) = 3;
  m(2,2,1) = 4;
  
  m(1,1,2) = 5;
  m(1,2,2) = 6;
  m(2,1,2) = 7;
  m(2,2,2) = 8;
  EXPECT_EQ( 1, *(m.begin()) );
}

TEST(Tensor, end){
  Tensor<int> m{2,2,2}; // 2x2x2 tensor
  m(1,1,1) = 1;
  m(1,2,1) = 2;
  m(2,1,1) = 3;
  m(2,2,1) = 4;
  
  m(1,1,2) = 5;
  m(1,2,2) = 6;
  m(2,1,2) = 7;
  m(2,2,2) = 8;
  EXPECT_EQ( 8, *std::prev(m.end()) );
}

TEST(Tensor, range_for){
  Tensor<double> m{2,2,2}; // 2x2x2 tensor
  double c{1.0};
  m(1,1,1) = -9999;
  m(2,2,1) = -9999;

  for (auto &x : m)
    x = ++c;
  
  EXPECT_DOUBLE_EQ( 2.0, m(1,1,1) );
  EXPECT_DOUBLE_EQ( 3.0, m(1,2,1) );
  EXPECT_DOUBLE_EQ( 4.0, m(2,1,1) );
  EXPECT_DOUBLE_EQ( 5.0, m(2,2,1) );
  
  EXPECT_DOUBLE_EQ( 6.0, m(1,1,2) );
  EXPECT_DOUBLE_EQ( 7.0, m(1,2,2) );
  EXPECT_DOUBLE_EQ( 8.0, m(2,1,2) );
  EXPECT_DOUBLE_EQ( 9.0, m(2,2,2) );
}

TEST(Tensor, range_for_zero){
  Tensor<double> m{1,0,1,0,1,0}; // 2x2x2 tensor
  double c{0.0};
  m(0,0,0) = -9999;
  m(1,1,1) = -9999;

  for (auto &x : m)
    x = ++c;

  EXPECT_DOUBLE_EQ( 1.0, m(0,0,0) );
  EXPECT_DOUBLE_EQ( 2.0, m(0,1,0) );
  EXPECT_DOUBLE_EQ( 3.0, m(1,0,0) );
  EXPECT_DOUBLE_EQ( 4.0, m(1,1,0) );
  
  EXPECT_DOUBLE_EQ( 5.0, m(0,0,1) );
  EXPECT_DOUBLE_EQ( 6.0, m(0,1,1) );
  EXPECT_DOUBLE_EQ( 7.0, m(1,0,1) );
  EXPECT_DOUBLE_EQ( 8.0, m(1,1,1) );
  
}

TEST(Tensor, set_value){
  Tensor<double> m{1,0,1,0,1,0};
  double c{0.0};
  
  for (auto &x : m)
    x = ++c;

  m = -9999.;

  EXPECT_DOUBLE_EQ( m(0,0,0), -9999. );
  EXPECT_DOUBLE_EQ( m(0,1,0), -9999. );
  EXPECT_DOUBLE_EQ( m(1,0,0), -9999. );
  EXPECT_DOUBLE_EQ( m(1,1,0), -9999. );
  
  EXPECT_DOUBLE_EQ( m(0,0,1), -9999. );
  EXPECT_DOUBLE_EQ( m(0,1,1), -9999. );
  EXPECT_DOUBLE_EQ( m(1,0,1), -9999. );
  EXPECT_DOUBLE_EQ( m(1,1,1), -9999. );
}

TEST(Tensor, copy_constructor){
  Tensor<double> m{2,2,2};
  double c{0.0};
  
  for (auto &x : m)
    x = ++c;
  
  Tensor<double> m1{m};

  testing::internal::CaptureStdout();
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
   
  EXPECT_DOUBLE_EQ( m1(1,1,1), 1 );
  EXPECT_DOUBLE_EQ( m1(1,2,1), 2 );
  EXPECT_DOUBLE_EQ( m1(2,1,1), 3 );
  EXPECT_DOUBLE_EQ( m1(2,2,1), 4 );
  
  EXPECT_DOUBLE_EQ( m1(1,1,2), 5 );
  EXPECT_DOUBLE_EQ( m1(1,2,2), 6 );
  EXPECT_DOUBLE_EQ( m1(2,1,2), 7 );
  EXPECT_DOUBLE_EQ( m1(2,2,2), 8 );
 
}

TEST(Tensor, out_of_range){
  Tensor<double> m{3,3,3};

  // testing "at" (always throws errors if needed)
  EXPECT_NO_THROW( m.at(1,1,1) );
  EXPECT_NO_THROW( m.at(2,2,2) );
  EXPECT_NO_THROW( m.at(3,3,3) );
  
  EXPECT_ANY_THROW( m.at(0,1,1) );
  EXPECT_ANY_THROW( m.at(1,0,1) );
  EXPECT_ANY_THROW( m.at(1,1,0) );

  EXPECT_ANY_THROW( m.at(4,1,1) );
  EXPECT_ANY_THROW( m.at(1,4,1) );
  EXPECT_ANY_THROW( m.at(1,1,4) );
  

  // testing "()" (throws errors only if in DEBUG)
  EXPECT_NO_THROW( m.at(1,1,1) );
  EXPECT_NO_THROW( m.at(2,2,2) );
  EXPECT_NO_THROW( m.at(3,3,3) );

#ifndef NDEBUG 
  EXPECT_ANY_THROW( m.at(0,1,1) );
  EXPECT_ANY_THROW( m.at(1,0,1) );
  EXPECT_ANY_THROW( m.at(1,1,0) );

  EXPECT_ANY_THROW( m.at(4,1,1) );
  EXPECT_ANY_THROW( m.at(1,4,1) );
  EXPECT_ANY_THROW( m.at(1,1,4) );
  
#else
  EXPECT_ANY_THROW( m.at(0,1,1) );
  EXPECT_ANY_THROW( m.at(1,0,1) );
  EXPECT_ANY_THROW( m.at(1,1,0) );

  EXPECT_ANY_THROW( m.at(4,1,1) );
  EXPECT_ANY_THROW( m.at(1,4,1) );
  EXPECT_ANY_THROW( m.at(1,1,4) );
#endif
}

TEST(Tensor, out_of_range_zero){
  Tensor<double> m{3,0,3,0,1,0};
  
  EXPECT_NO_THROW( m.at(0,0,0) );
  EXPECT_NO_THROW( m.at(1,1,0) );
  EXPECT_NO_THROW( m.at(2,2,0) );
  EXPECT_NO_THROW( m.at(3,3,0) );

  EXPECT_NO_THROW( m.at(0,0,1) );
  EXPECT_NO_THROW( m.at(1,1,1) );
  EXPECT_NO_THROW( m.at(2,2,1) );
  EXPECT_NO_THROW( m.at(3,3,1) );
  
  EXPECT_ANY_THROW( m.at(-3,-3,-3) );
  EXPECT_ANY_THROW( m.at(4,4,4) );
  EXPECT_ANY_THROW( m.at(5000,5000,500) );
  
}


  // // print to check the resulting tensor
  // std::cout << "Original Tensor (m)" << std::endl;
  // for(int k=1; k<=2; k++){
  //   for(int i=1; i<=2; i++){
  //     for(int j=1; j<=2; j++){
  // 	std::cout << m(i,j,k) << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

  //  std::cout << "derived Tensor (m1)" << std::endl;
  // for(int k=1; k<=2; k++){
  //   for(int i=1; i<=2; i++){
  //     for(int j=1; j<=2; j++){
  // 	std::cout << m1(i,j,k) << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  // }
