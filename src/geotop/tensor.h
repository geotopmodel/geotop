//
// Created by elisa on 2018/07/26.
//

#ifndef GEOTOP_TENSOR_H
#define GEOTOP_TENSOR_H

#include "geotop_asserts.h"
#include <exception>
#include <memory> // to use std::unique_ptr
#include "matrix.h"

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
    std::unique_ptr<T[]> co; // it has to be written after the lower and upper limits (otherwise: sigfault)

    /** pointer to the first accessible element (needed by an iterator) */
    T *begin() noexcept { return &co[0]; }

    /** pointer to the one-past the last element (needed by an iterator)*/
    T *end() noexcept { return &co[(ndh-ndl+1)*(nrh-nrl+1)*(nch-ncl+1)]; }

    /** const pointer to the first element accessible element */
    const T *begin() const noexcept { return &co[0]; }

    /** const pointer to the one-past the last element */
    const T *end() const noexcept { return &co[(ndh-ndl+1)*(nrh-nrl+1)*(nch-ncl+1)]; }


    /** destructor. default is fine */
    ~Tensor() = default;

    /** default constructor is deleted */
    Tensor() = delete;

/**
   * constructor
   * @param _nrl,_nrh lower and upper bound for rows
   * @param _ncl,_nch lower and upper bound for columns
   * @param _ndl,_ndh lower and upper bound for depth
   */

    Tensor(const std::size_t _ndh,  const std::size_t _ndl,
           const std::size_t _nrh, const std::size_t _nrl, const std::size_t _nch,  const std::size_t _ncl):
            ndh{_ndh}, ndl{_ndl}, nrh{_nrh}, nrl{_nrl}, nch{_nch}, ncl{_ncl},
            n_dep{ndh-ndl+1}, n_row{nrh-nrl+1}, n_col{nch-ncl+1}, co { new T[(ndh-ndl+1)*(nrh-nrl+1)*(nch-ncl+1)]{} } {}
    // .................................................................................|-->initialize all elements to 0

    Tensor(const std::size_t d, const std::size_t r, const std::size_t c): Tensor{d,1,r,1,c,1} {}

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

    T& operator[](const std::size_t i) noexcept {
        return co[i];
    }

    const T& operator[](const std::size_t i) const noexcept {
        return co[i];
    }

    /**
    * access operator. When the code is compiled in debug mode, it performes
    * a range check. No check is done when the code is compiled in release mode.
    */
    T& operator()(const std::size_t k, const std::size_t i, const std::size_t j) {
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

    /** set all elements of the vector to @p v
        * this is useful to reinizialize all the elements of the vector to zero
        * my_Tensor = 0
        */
    Tensor<T> &operator=(const T v) {
        for (auto &x : *this)
            x = v;
        return *this;
    }

    /**
       * Copy constructor
       */
    Tensor(const Tensor<T> &t)
            : ndl{t.ndl}, ndh{t.ndh}, nrl{t.nrl}, nrh{t.nrh}, ncl{t.ncl}, nch{t.nch},
              n_dep{ndh-ndl+1}, n_row{nrh-nrl+1}, n_col{nch-ncl+1},
              co{new T[(ndh-ndl+1)*(nrh-nrl+1)*(nch-ncl+1)]} {
        for (std::size_t i=0; i<(ndh-ndl+1)*(nrh-nrl+1)*(nch-ncl+1); ++i)
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

    /** Given a layer (k) and a row (i), it gives all the elements corrisponding to different columns */
    RowView<T> row(const std::size_t k, const std::size_t i) {
        GEO_ASSERT_IN_RANGE(k, ndl, ndh);
        GEO_ASSERT_IN_RANGE(i, nrl, nrh);
        return RowView<T> { &co[(i-nrl)*n_col + (k-ndl)*(n_row*n_col)], nch, ncl}; // intialization
    };

    MatrixView<T> matrix(const std::size_t k) {
        GEO_ASSERT_IN_RANGE(k, ndl, ndh);
        return MatrixView<T> { &co[(k-ndl)*(n_row*n_col)], nrh, nrl, nch, ncl}; // intialization
    };

};
// ----------------------------------------------------------------------------------------------------------------
#endif // GEOTOP_Tensor_H
