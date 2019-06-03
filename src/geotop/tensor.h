//
// Created by elisa on 2018/07/26.
//

#ifndef _GEOTOP_TENSOR_H
#define _GEOTOP_TENSOR_H

#include "geotop_asserts.h"
#include <exception>
#include <memory> // to use std::unique_ptr
#include "matrix.h"
#include "matrixview.h"

template <class T> class RowView;
template <class T> class MatrixView;

template <class T> class Tensor {

public:
    /** lower and upper bounds */
    std::size_t ndh, ndl; // depth
    std::size_t nrh, nrl; // rows
    std::size_t nch, ncl; // columns

    std::size_t n_dep;
    std::size_t n_row;
    std::size_t n_col;

    /** the actual data */
    std::unique_ptr<T[]> co;

    /** pointer to the first accessible element (needed by an iterator) */
    T *begin() noexcept { return &co[0]; }

    /** pointer to the one-past the last element (needed by an iterator)*/
    T *end() noexcept { return &co[n_dep*n_row*n_col]; }

    /** const pointer to the first element accessible element */
    const T *begin() const noexcept { return &co[0]; }

    /** const pointer to the one-past the last element */
    const T *end() const noexcept { return &co[n_dep*n_row*n_col]; }

    /** destructor. default is fine */
    ~Tensor() = default;

    /** default constructor is deleted */
    Tensor() = delete;

    /** subscripting operator (non-checked) */
    T& operator[](const std::size_t i) noexcept {
        return co[i];
    }

    /** subscripting operator (non-checked) */
    const T& operator[](const std::size_t i) const noexcept {
        return co[i];
    }

    /** range-checked access operator */
    T &at(const std::size_t k, const std::size_t i, const std::size_t j) {
        GEO_ERROR_IN_RANGE(k, ndl, ndh);
        GEO_ERROR_IN_RANGE(i, nrl, nrh);
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return (*this)[(i-nrl)*n_col + (j-ncl) + (k-ndl)*(n_row*n_col)];
    }

    /** range-checked access operator */
    const T &at(const std::size_t k, const std::size_t i, const std::size_t j) const {
        GEO_ERROR_IN_RANGE(k, ndl, ndh);
        GEO_ERROR_IN_RANGE(i, nrl, nrh);
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return (*this)[(i-nrl)*n_col + (j-ncl) + (k-ndl)*(n_row*n_col)];
    }

    /**
    * access operator. When the code is compiled in debug mode, it performs
    * a range check. No check is done when the code is compiled in release mode.
    */
    T& operator()(const std::size_t k, const std::size_t i, const std::size_t j)
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(k,i,j);
#else
        return (*this)[(i-nrl)*n_col + (j-ncl) + (k-ndl)*(n_row*n_col)];

#endif
    }

    /**
        * access operator. When the code is compiled in debug mode, it performes
        * a range check. No check is done when the code is compiled in release mode.
        */
    const T& operator()(const std::size_t k, const std::size_t i, const std::size_t j) const
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(k,i,j);
#else
        return (*this)[(i-nrl)*n_col + (j-ncl) + (k-ndl)*(n_row*n_col)];
#endif
    }

    /**
   * constructor
   * @param _nrl,_nrh lower and upper bound for rows
   * @param _ncl,_nch lower and upper bound for columns
   * @param _ndl,_ndh lower and upper bound for depth
   */

    Tensor(const std::size_t _ndh,  const std::size_t _ndl,
           const std::size_t _nrh, const std::size_t _nrl, const std::size_t _nch,  const std::size_t _ncl):
            ndh{_ndh}, ndl{_ndl}, nrh{_nrh}, nrl{_nrl}, nch{_nch}, ncl{_ncl},
            n_dep{ndh-ndl+1}, n_row{nrh-nrl+1}, n_col{nch-ncl+1}, co { new T[(ndh-ndl+1)*(nrh-nrl+1)*(nch-ncl+1)]{} } {} // initialize all elements to 0

    Tensor(const std::size_t d, const std::size_t r, const std::size_t c): Tensor{d,1,r,1,c,1} {}

    /**
       * Copy constructor
       */
    Tensor(const Tensor<T> &t):
            ndh{t.ndh}, ndl{t.ndl}, nrh{t.nrh}, nrl{t.nrl}, nch{t.nch}, ncl{t.ncl},
            n_dep{t.n_dep}, n_row{t.n_row}, n_col{t.n_col},
            co{new T[n_dep*n_row*n_col]} {

        for (std::size_t i=0; i<n_dep*n_row*n_col; ++i)
            (*this)[i] = t[i];
    }

    /** Move constructor */
    Tensor(Tensor<T> &&t) = default;

    /** Move assignment */
    Tensor<T>& operator=(Tensor<T> &&t) = default;

    /** Copy assignment */
    Tensor<T> &operator=(const Tensor<T> &t) {
        co.reset(); // release acquired memory
        *this = Tensor<T>{t}; // use move assignment and copy constructor
        return *this;
    }

    /** set all elements of the tensor to @p t
       * this is useful to reinitialize all the elements of the tensor to zero
       * my_tensor = 0
       */
    Tensor<T> &operator=(const T t) {
        for (auto &x : *this)
            x = t;
        return *this;
    }

    /** returns the i-th row of the k-th layer of the tensor */
    RowView<T> row(const std::size_t k, const std::size_t i) {
        GEO_ASSERT_IN_RANGE(k, ndl, ndh);
        GEO_ASSERT_IN_RANGE(i, nrl, nrh);
        return RowView<T> { &co[(i-nrl)*n_col + (k-ndl)*(n_row*n_col)], nch, ncl}; // initialization
    };

    MatrixView<T> matrix(const std::size_t k) {
        GEO_ASSERT_IN_RANGE(k, ndl, ndh);
        return MatrixView<T> { &co[(k-ndl)*(n_row*n_col)], nrh, nrl, nch, ncl}; // initialization
    };

};

#endif // _GEOTOP_TENSOR_H
