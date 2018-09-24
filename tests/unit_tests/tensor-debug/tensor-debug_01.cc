#include <gtest/gtest.h>

#include <tensor.h> // for Tensor<double>
#include <iterator> // std::prev

#include <turtle.h> // for DOUBLETENSOR, includes:
// #include "t_alloc.h" // for new_doubletensor

#include <typeinfo>
#include <ostream>
// tensor(d,r,c) => d = depth; r=row; c=column
TEST(Tensor, Tensor_value){
  Tensor<double>* t = new Tensor<double>{2,2,2}; // 2x2x2 tensor

  double c{0.0};
  for (auto &x : *t) // need begin+end
    x = ++c;
  
  EXPECT_DOUBLE_EQ( 1.0, (*t)(1,1,1) );
  EXPECT_DOUBLE_EQ( 2.0, (*t)(1,1,2) );
  EXPECT_DOUBLE_EQ( 3.0, (*t)(1,2,1) );
  EXPECT_DOUBLE_EQ( 4.0, (*t)(1,2,2) );
  
  EXPECT_DOUBLE_EQ( 5.0, (*t)(2,1,1) );
  EXPECT_DOUBLE_EQ( 6.0, (*t)(2,1,2) );
  EXPECT_DOUBLE_EQ( 7.0, (*t)(2,2,1) );
  EXPECT_DOUBLE_EQ( 8.0, (*t)(2,2,2) );

  for(std::size_t k=1; k<=t->ndh; k++){
    for(std::size_t i=1; i<=t->nrh; i++){
      for(std::size_t j=1; j<=t->nch; j++){
  	std::cout << (*t)(k,i,j) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

TEST(Tensor, TENSOR_value){
  DOUBLETENSOR *pa; // 2x2x2 tensor
  
  pa = new_doubletensor(2,2,2); 

  int idx = 1;
  
  for(int k=1; k<=pa->ndh; k++){
    for(int i=1; i<=pa->nrh; i++){
      for(int j=1; j<=pa->nch; j++){
  	pa->co[k][i][j] = idx;
	idx ++;
      }
    }
  }
  
  for(int k=1; k<=pa->ndh; k++){
    for(int i=1; i<=pa->nrh; i++){
      for(int j=1; j<=pa->nch; j++){
  	std::cout << pa->co[k][i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  
  EXPECT_DOUBLE_EQ( 1.0, pa->co[1][1][1] );
  EXPECT_DOUBLE_EQ( 2.0, pa->co[1][1][2] );
  EXPECT_DOUBLE_EQ( 3.0, pa->co[1][2][1] );
  EXPECT_DOUBLE_EQ( 4.0, pa->co[1][2][2] );
  
  EXPECT_DOUBLE_EQ( 5.0, pa->co[2][1][1] );
  EXPECT_DOUBLE_EQ( 6.0, pa->co[2][1][2] );
  EXPECT_DOUBLE_EQ( 7.0, pa->co[2][2][1] );
  EXPECT_DOUBLE_EQ( 8.0, pa->co[2][2][2] );
}

TEST(Tensor, Compare_types){
  Tensor<double>* t = new Tensor<double>{2,2,2}; // 2x2x2 tensor

  DOUBLETENSOR *pa; // 2x2x2 tensor
  pa = new_doubletensor(2,2,2); 

  std::cout << std::endl;
  std::cout << "pa->co[1] is a " << typeid(pa->co[1]).name() << std::endl;
  std::cout << "pa->co[1] = " << pa->co[1] << std::endl;
  std::cout << "t->matrix(1) is a " << typeid(t->matrix(1)).name() << std::endl;
  std::cout << std::endl;
  std::cout << "pa->co[1][1] is a " << typeid(pa->co[1][1]).name() << std::endl;
  std::cout << "pa->co[1][1] = " << pa->co[1][1] << std::endl;
  std::cout << "t->row(1,1) is a " << typeid(t->row(1,1)).name() << std::endl;
  std::cout << std::endl;

  //std::cout << "t->matrix(1) = " << t->matrix(1) << std::endl;
  std::cout << std::endl;
  //std::cout << "t->row(1,1) = " << t->row(1,1) << std::endl;
  std::cout << std::endl;
}


TEST(Tensor, damny_references){
    Tensor<double>* t = new Tensor<double>{2,2,2}; // 2x2x2 tensor
    MatrixView<double> &&tr = t->matrix(1);
    tr(1,1) = 5;
    std::cout << tr(1,1) << std::endl;

    // DOUBLETENSOR *pa; // 2x2x2 tensor
    // pa = new_doubletensor(2,2,2);
    // double **par = pa->co[1];
    // par->co[1][1] = 5;
    // std::cout << par->co[1][1] << std::endl;
}
