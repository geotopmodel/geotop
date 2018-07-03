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
    /** pointer to the first accessible element */
    T *begin() noexcept { return &co[nrl][ncl]; }

    /** pointer to the one-past the last element */
    T *end() noexcept { return &co[nrh + 1][nch + 1]; }

    /** const pointer to the first element accessible element */
    const T *begin() const noexcept { return &co[nrl][ncl]; }

    /** const pointer to the one-past the last element */
    const T *end() const noexcept { return &co[nrh + 1][nch + 1]; }

    /** subscripting operator (non-checked) */
    T &operator[](const std::size_t i, const std::size_t j) noexcept { return co[i][j]; }

    /** subscripting operator (non-checked) */
    const T &operator[](const std::size_t i) const noexcept { return co[i][j]; }

// suggested by Alberto
 //   T *operator[](const std::size_t i) { return &co[ncol*i]; }

    /** destructor. default is fine */
    ~Matrix() = default;

    /** default constructor is deleted */
    Matrix() = delete;

/**
   * constructor
   * @param _nrl,_nrh lower and upper bound for rows
   * @param _ncl,_nch lower and upper bound for columns
   */
    Matrix(const std::size_t _nrl = 1, const std::size_t _nrh, const std::size_t _ncl = 1, const std::size_t _nch)
            : nrl{_nrl}, nrh{_nrh}, ncl{_ncl}, nch{_cch}, co{new T[nrh+1][nch+1]{}} {}

    /** set all elements of the vector to @p m
  * this is useful to reinizialize all the elements of the vector to zero
  * my_matrix = 0
  */
//    Matrix<T> &operator=(const T m) {
//      for (auto &x : *this)
//        x = m;
//      return *this;
//    } //    NON SO SE FUNZIONA....

    /** lower bound */
    std::size_t nrl, ncl;

    /** upper bound */
    std::size_t nrh, nch;

    /** the actual data */
    std::unique_ptr<T[]> co;

//private:
};
#endif // GEOTOP_MATRIX_H
