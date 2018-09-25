#include <gtest/gtest.h>
#include <tensor.h>
#include <matrix.h>
#include <matrixview.h>

// tensor(d,r,c) => d = depth; r=row; c=column

TEST(Tensor, TensorElement){
  Tensor<int> t{2,2,2}; // 2x2x2 tensor
  int c{0};
   
  for (auto &x : t)
    x = ++c;

   testing::internal::CaptureStdout();
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  
  // (depth,row)[column] 
  EXPECT_EQ( 1, t[0] );
  EXPECT_EQ( 2, t[1] );
  EXPECT_EQ( 3, t[2] );
  EXPECT_EQ( 4, t[3] );

  EXPECT_EQ( 5, t[4] );
  EXPECT_EQ( 6, t[5] );
  EXPECT_EQ( 7, t[6] );
  EXPECT_EQ( 8, t[7] );
}

TEST(Tensor, RowView){
  Tensor<int> t{2,2,2}; // 2x2x2 tensor
  int c{0};
   
  for (auto &x : t)
    x = ++c;

   testing::internal::CaptureStdout();
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  
  // (depth,row)[column] 
  EXPECT_EQ( 1, t.row(1,1)[1] );
  EXPECT_EQ( 2, t.row(1,1)[2] );
  EXPECT_EQ( 3, t.row(1,2)[1] );
  EXPECT_EQ( 4, t.row(1,2)[2] );

  EXPECT_EQ( 5, t.row(2,1)[1] );
  EXPECT_EQ( 6, t.row(2,1)[2] );
  EXPECT_EQ( 7, t.row(2,2)[1] );
  EXPECT_EQ( 8, t.row(2,2)[2] );
}

TEST(Tensor, RowView_CheckIndex){
  Tensor<int> t{2,2,2}; // 2x2x2 tensor

   testing::internal::CaptureStdout();
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  
  EXPECT_NO_THROW( t.row(1,1) );
  EXPECT_NO_THROW( t.row(1,2) );
  EXPECT_NO_THROW( t.row(2,1) );
  EXPECT_NO_THROW( t.row(2,2) );

#ifndef NDEBUG
  EXPECT_ANY_THROW( t.row(0,0) );
  EXPECT_ANY_THROW( t.row(3,3) );

#else
  EXPECT_NO_THROW( t.row(0,0) );
  EXPECT_NO_THROW( t.row(3,3) );
#endif
}


TEST(Tensor, MatrixView){
  Tensor<int> t{2,2,2}; // 2x2x2 tensor
  int c{0};
   
  for (auto &x : t)
    x = ++c;

   testing::internal::CaptureStdout();
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  
  // (depth)[row][column] 
  EXPECT_EQ( 1, t.matrix(1)(1,1) );
  EXPECT_EQ( 2, t.matrix(1)(1,2) );
  EXPECT_EQ( 3, t.matrix(1)(2,1) );
  EXPECT_EQ( 4, t.matrix(1)(2,2) );

  EXPECT_EQ( 5, t.matrix(2)(1,1) );
  EXPECT_EQ( 6, t.matrix(2)(1,2) );
  EXPECT_EQ( 7, t.matrix(2)(2,1) );
  EXPECT_EQ( 8, t.matrix(2)(2,2) );
}


TEST(Tensor, TensorMatrix_CheckIndex){
  Tensor<int> t{2,2,2}; // 2x2x2 tensor

   testing::internal::CaptureStdout();
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  
  EXPECT_NO_THROW( t.matrix(1) );
  EXPECT_NO_THROW( t.matrix(2) );

#ifndef NDEBUG
  EXPECT_ANY_THROW( t.matrix(0) );
  EXPECT_ANY_THROW( t.matrix(3) );

#else
  EXPECT_NO_THROW( t.matrix(0) );
  EXPECT_NO_THROW( t.matrix(3) );
#endif
}
