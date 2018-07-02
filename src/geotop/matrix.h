//
// Created by elisa on 2018/07/02.
//

#ifndef GEOTOP_MATRIX_H
#define GEOTOP_MATRIX_H

#include "geotop_asserts.h"
#include <cassert>
#include <exception>
#include <memory>
#include <sstream>

template <typename T> class Matrix {
public:
    /**
  * @return the size of the matrix
  */
    std::size_t size() const noexcept { return _size; }

    /** pointer to the first element accessible element */
    T *begin() noexcept { return &co[nl]; }

    /** pointer to the one-past the last element */
    T *end() noexcept { return &co[nh + 1]; }

    /** const pointer to the first element accessible element */
    const T *begin() const noexcept { return &co[nl]; }

    /** const pointer to the one-past the last element */
    const T *end() const noexcept { return &co[nh + 1]; }

    /** subscripting operator */
    T* operator[](const std::size_t i) { return &co[ncol*i]; }

    /** destructor. default is fine */
    ~Matrix() = default;

    /** default constructor is deleted */
    Matrix() = delete;

/**
   * constructor
   * @param n_row number of rows
   * @param n_col number of columns
   */
    Matrix(std::size_t n_row, std::size_t n_col): nrow{n_row}, ncol{n_col} {}

    /** set all elements of the vector to @p v
  * this is useful to reinizialize all the elements of the vector to zero
  * my_matrix = 0
  */
    Matrix<T> &operator=(const T m) {
      for (auto &x : *this)
        x = m;
      return *this;
    }

    /** the actual data */
    std::unique_ptr<T[]> co;
    /** lower bound */
    std::size_t nl;
    /** upper bound */
    std::size_t nh;


private:
    std::size_t nrow, ncol;

};
#endif // GEOTOP_MATRIX_H
