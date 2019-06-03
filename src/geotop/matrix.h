//
// Created by elisa on 2018/07/02.
//

#ifndef _GEOTOP_MATRIX_H
#define _GEOTOP_MATRIX_H

#include "geotop_asserts.h"
#include <memory> // to use std::unique_ptr
#include "rowview.h"

template <class T> class Matrix {

public:
    /** lower and upper bounds */
    std::size_t nrh, nrl; // rows
    std::size_t nch, ncl; // columns

    std::size_t n_row;
    std::size_t n_col;

    /** the actual data */
    std::unique_ptr<T[]> co;

    /** pointer to the first accessible element (needed by an iterator) */
    T *begin() noexcept { return &co[0]; }

    /** pointer to the one-past the last element (needed by an iterator)*/
    T *end() noexcept { return &co[n_row*n_col]; }

    /** const pointer to the first element accessible element */
    const T *begin() const noexcept { return &co[0]; }

    /** const pointer to the one-past the last element */
    const T *end() const noexcept { return &co[n_row*n_col]; }

    /** subscripting operator (non-checked) */
    T& operator[](const std::size_t i) noexcept {
        return co[i];
    }

    /** subscripting operator (non-checked) */
    const T& operator[](const std::size_t i) const noexcept {
        return co[i];
    }

    /** range-checked access operator */
    T &at(const std::size_t i, const std::size_t j) {
        GEO_ERROR_IN_RANGE(i, nrl, nrh);
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return (*this)[(i-nrl)*n_col+(j-ncl)];
    }

    /** range-checked access operator */
    const T &at(const std::size_t i, const std::size_t j) const {
        GEO_ERROR_IN_RANGE(i, nrl, nrh);
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return (*this)[(i-nrl)*n_col+(j-ncl)];
    }

    /**
    * access operator. When the code is compiled in debug mode, it performs
    * a range check. No check is done when the code is compiled in release mode.
    */
    T& operator()(const std::size_t i, const std::size_t j)
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(i,j);
#else
        return (*this)[(i-nrl)*n_col+(j-ncl)];
#endif
    }

    /**
        * access operator. When the code is compiled in debug mode, it performs
        * a range check. No check is done when the code is compiled in release mode.
        */
    const T& operator()(const std::size_t i, const std::size_t j) const
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(i,j);
#else
        return (*this)[(i-nrl)*n_col+(j-ncl)];
#endif
    }

    /** destructor. default is fine */
    ~Matrix() = default;

    /** default constructor is deleted */
    Matrix() = delete;

/**
   * constructor
   * @param _nrl,_nrh lower and upper bound for rows
   * @param _ncl,_nch lower and upper bound for columns
   */

    Matrix(const std::size_t _nrh, const std::size_t _nrl, const std::size_t _nch,  const std::size_t _ncl):
            nrh{_nrh}, nrl{_nrl}, nch{_nch}, ncl{_ncl},
            n_row{nrh-nrl+1}, n_col{nch-ncl+1}, co { new T[n_row*n_col]{} } {} // initialize all elements to 0

    Matrix(const std::size_t r, const std::size_t c): Matrix{r,1,c,1} {}

    /**
       * Copy constructor
       */
    Matrix(const Matrix<T> &m):
            nrh{m.nrh}, nrl{m.nrl}, nch{m.nch}, ncl{m.ncl},
            n_row{nrh-nrl+1}, n_col{nch-ncl+1}, co{new T[n_row*n_col]} {

        for (std::size_t i=0; i<n_row*n_col; ++i)
            (*this)[i] = m[i];
    }

    /** Move constructor */
    Matrix(Matrix<T> &&m) = default;

    /** Move assignment */
    Matrix<T>& operator=(Matrix<T> &&m) = default;

    /** Copy assignment */
    Matrix<T> &operator=(const Matrix<T> &m) {
        co.reset(); // release acquired memory
        *this = Matrix<T>{m}; // use move assignment and copy constructor
        return *this;
    }

    /** returns the i-th row of the matrix */
    RowView<T> row(const std::size_t i) {
        GEO_ASSERT_IN_RANGE(i, nrl, nrh);
        return RowView<T> { &co[(i-nrl)*n_col], nch, ncl }; // RowView<T> is constructed with the passed nch and ncl
    };

    /** set all elements of the matrix to @p m
       * this is useful to reinitialize all the elements of the matrix to zero
       * my_matrix = 0
       */
    Matrix<T> &operator=(const T m) {
        for (auto &x : *this)
            x = m;
        return *this;
    }

};

#endif // _GEOTOP_MATRIX_H
