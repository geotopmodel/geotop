#include <gtest/gtest.h>
#include <tensor.h>
#include <iterator> // std::prev

// tensor(d,r,c) => d = depth; r=row; c=column
TEST(Tensor, constructor_6args){
  Tensor<int> t{2,0,3,1,3,0}; // 3x4x3 tensor 
  EXPECT_EQ( std::size_t{3}, t.n_dep );
  EXPECT_EQ( std::size_t{3}, t.n_row );
  EXPECT_EQ( std::size_t{4}, t.n_col );
}

TEST(Tensor, constructor_3args){
  Tensor<int> t{2,3,5}; // 2x3x5 tensor
  EXPECT_EQ( std::size_t{2}, t.n_dep );
  EXPECT_EQ( std::size_t{3}, t.n_row );
  EXPECT_EQ( std::size_t{5}, t.n_col );
}

TEST(Tensor, initialization){
  Tensor<int> t{2,2,2}; // 2x2x2 tensor
  EXPECT_EQ( 0, t(1,1,1) );
  EXPECT_EQ( 0, t(1,1,2) );
  EXPECT_EQ( 0, t(1,2,1) );
  EXPECT_EQ( 0, t(1,2,2) );
  
  EXPECT_EQ( 0, t(2,1,1) );
  EXPECT_EQ( 0, t(2,1,2) );
  EXPECT_EQ( 0, t(2,2,1) );
  EXPECT_EQ( 0, t(2,2,2) );

#ifndef NDEBUG
  EXPECT_ANY_THROW( t(0,0,0)); 
  EXPECT_ANY_THROW( t(3,3,3));
#else
  EXPECT_NO_THROW( t(0,0,0) );
  EXPECT_NO_THROW( t(3,3,3) );
#endif
}

TEST(Tensor, begin){
  Tensor<int> t{2,2,2}; // 2x2x2 tensor
  t(1,1,1) = 1;
  t(1,1,2) = 2;
  t(1,2,1) = 3;
  t(1,2,2) = 4;
  
  t(2,1,1) = 5;
  t(2,1,2) = 6;
  t(2,2,1) = 7;
  t(2,2,2) = 8;
  EXPECT_EQ( 1, *(t.begin()) );
}

TEST(Tensor, end){
  Tensor<int> t{2,2,2}; // 2x2x2 tensor
  t(1,1,1) = 1;
  t(1,1,2) = 2;
  t(1,2,1) = 3;
  t(1,2,2) = 4;
  
  t(2,1,1) = 5;
  t(2,1,2) = 6;
  t(2,2,1) = 7;
  t(2,2,2) = 8;
  EXPECT_EQ( 8, *std::prev(t.end()) );
}

TEST(Tensor, range_for){
  Tensor<double> t{2,2,2}; // 2x2x2 tensor
  double c{1.0};
  t(1,1,1) = -9999;
  t(2,2,2) = -9999;

  for (auto &x : t)
    x = ++c;
  
  EXPECT_DOUBLE_EQ( 2.0, t(1,1,1) );
  EXPECT_DOUBLE_EQ( 3.0, t(1,1,2) );
  EXPECT_DOUBLE_EQ( 4.0, t(1,2,1) );
  EXPECT_DOUBLE_EQ( 5.0, t(1,2,2) );
  
  EXPECT_DOUBLE_EQ( 6.0, t(2,1,1) );
  EXPECT_DOUBLE_EQ( 7.0, t(2,1,2) );
  EXPECT_DOUBLE_EQ( 8.0, t(2,2,1) );
  EXPECT_DOUBLE_EQ( 9.0, t(2,2,2) );
}

TEST(Tensor, range_for_zero){
  Tensor<double> t{1,0,1,0,1,0}; // 2x2x2 tensor
  double c{0.0};
  t(0,0,0) = -9999;
  t(1,1,1) = -9999;

  for (auto &x : t)
    x = ++c;

  EXPECT_DOUBLE_EQ( 1.0, t(0,0,0) );
  EXPECT_DOUBLE_EQ( 2.0, t(0,0,1) );
  EXPECT_DOUBLE_EQ( 3.0, t(0,1,0) );
  EXPECT_DOUBLE_EQ( 4.0, t(0,1,1) );
  
  EXPECT_DOUBLE_EQ( 5.0, t(1,0,0) );
  EXPECT_DOUBLE_EQ( 6.0, t(1,0,1) );
  EXPECT_DOUBLE_EQ( 7.0, t(1,1,0) );
  EXPECT_DOUBLE_EQ( 8.0, t(1,1,1) );
  
}

TEST(Tensor, set_value){
  Tensor<double> t{1,0,1,0,1,0};
  double c{0.0};
  
  for (auto &x : t)
    x = ++c;

  t = -9999.;

  EXPECT_DOUBLE_EQ( t(0,0,0), -9999. );
  EXPECT_DOUBLE_EQ( t(0,0,1), -9999. );
  EXPECT_DOUBLE_EQ( t(0,1,0), -9999. );
  EXPECT_DOUBLE_EQ( t(0,1,1), -9999. );
  
  EXPECT_DOUBLE_EQ( t(1,0,0), -9999. );
  EXPECT_DOUBLE_EQ( t(1,0,1), -9999. );
  EXPECT_DOUBLE_EQ( t(1,1,0), -9999. );
  EXPECT_DOUBLE_EQ( t(1,1,1), -9999. );
}

TEST(Tensor, copy_constructor){
  Tensor<double> t{2,2,2};
  double c{0.0};
  
  for (auto &x : t)
    x = ++c;
  
  Tensor<double> t1{t};

  testing::internal::CaptureStdout();
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
   
  EXPECT_DOUBLE_EQ( t1(1,1,1), 1 );
  EXPECT_DOUBLE_EQ( t1(1,1,2), 2 );
  EXPECT_DOUBLE_EQ( t1(1,2,1), 3 );
  EXPECT_DOUBLE_EQ( t1(1,2,2), 4 );
  
  EXPECT_DOUBLE_EQ( t1(2,1,1), 5 );
  EXPECT_DOUBLE_EQ( t1(2,1,2), 6 );
  EXPECT_DOUBLE_EQ( t1(2,2,1), 7 );
  EXPECT_DOUBLE_EQ( t1(2,2,2), 8 );
 
}

TEST(Tensor, out_of_range){
  Tensor<double> t{3,3,3};

  // testing "at" (always throws errors if needed)
  EXPECT_NO_THROW( t.at(1,1,1) );
  EXPECT_NO_THROW( t.at(2,2,2) );
  EXPECT_NO_THROW( t.at(3,3,3) );
  
  EXPECT_ANY_THROW( t.at(0,1,1) );
  EXPECT_ANY_THROW( t.at(1,0,1) );
  EXPECT_ANY_THROW( t.at(1,1,0) );

  EXPECT_ANY_THROW( t.at(4,1,1) );
  EXPECT_ANY_THROW( t.at(1,4,1) );
  EXPECT_ANY_THROW( t.at(1,1,4) );
  

  // testing "()" (throws errors only if in DEBUG)
  EXPECT_NO_THROW( t.at(1,1,1) );
  EXPECT_NO_THROW( t.at(2,2,2) );
  EXPECT_NO_THROW( t.at(3,3,3) );

#ifndef NDEBUG 
  EXPECT_ANY_THROW( t.at(0,1,1) );
  EXPECT_ANY_THROW( t.at(1,0,1) );
  EXPECT_ANY_THROW( t.at(1,1,0) );

  EXPECT_ANY_THROW( t.at(4,1,1) );
  EXPECT_ANY_THROW( t.at(1,4,1) );
  EXPECT_ANY_THROW( t.at(1,1,4) );
  
#else
  EXPECT_ANY_THROW( t.at(0,1,1) );
  EXPECT_ANY_THROW( t.at(1,0,1) );
  EXPECT_ANY_THROW( t.at(1,1,0) );

  EXPECT_ANY_THROW( t.at(4,1,1) );
  EXPECT_ANY_THROW( t.at(1,4,1) );
  EXPECT_ANY_THROW( t.at(1,1,4) );
#endif
}

TEST(Tensor, out_of_range_zero){
  Tensor<double> t{1,0,3,0,3,0};
  
  EXPECT_NO_THROW( t.at(0,0,0) );
  EXPECT_NO_THROW( t.at(0,1,1) );
  EXPECT_NO_THROW( t.at(0,2,2) );
  EXPECT_NO_THROW( t.at(0,3,3) );

  EXPECT_NO_THROW( t.at(1,0,0) );
  EXPECT_NO_THROW( t.at(1,1,1) );
  EXPECT_NO_THROW( t.at(1,2,2) );
  EXPECT_NO_THROW( t.at(1,3,3) );
  
  EXPECT_ANY_THROW( t.at(-3,-3,-3) );
  EXPECT_ANY_THROW( t.at(4,4,4) );
  EXPECT_ANY_THROW( t.at(5000,5000,500) );
  
}


  // // print to check the resulting tensor
  // std::cout << "Original Tensor (t)" << std::endl;
  // for(int k=1; k<=2; k++){
  //   for(int i=1; i<=2; i++){
  //     for(int j=1; j<=2; j++){
  // 	std::cout << t(i,j,k) << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

  //  std::cout << "derived Tensor (t1)" << std::endl;
  // for(int k=1; k<=2; k++){
  //   for(int i=1; i<=2; i++){
  //     for(int j=1; j<=2; j++){
  // 	std::cout << t1(i,j,k) << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  // }
