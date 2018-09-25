#include <gtest/gtest.h>

#include <tensor.h> // for Tensor<double>
#include <iterator> // std::prev

#include <turtle.h> // for DOUBLETENSOR, includes:
// #include "t_alloc.h" // for new_doubletensor

#include <typeinfo>
#include <ostream>

// tensor(d,r,c) => d = depth; r=row; c=column

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


TEST(Tensor, Accessing_elements){
    Tensor<double>* t = new Tensor<double>{2,2,2}; // 2x2x2 tensor
    MatrixView<double> &&tr = t->matrix(1);
    tr(1,1) = 5;
    std::cout << std::endl;
    std::cout << "tr(1,1) = " << tr(1,1) << std::endl;
    std::cout << "tr(1,1) is a " << typeid(tr(1,1)).name() << std::endl;

    DOUBLETENSOR *pa; // 2x2x2 tensor
    pa = new_doubletensor(2,2,2);
    double **ppa = pa->co[1];
    ppa[1][1] = 5;
    std::cout << "ppa[1][1] = " << ppa[1][1] << std::endl;
    std::cout << "ppa[1][1] is a " << typeid(ppa[1][1]).name() << std::endl;
    std::cout << std::endl;
}
