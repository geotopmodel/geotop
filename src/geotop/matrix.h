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

template <class T> class Matrix {
public:
    /** pointer to the first accessible element (needed by an iterator) */
    T *begin() noexcept { return &co[0]; }

    /** pointer to the one-past the last element (needed by an iterator)*/
    T *end() noexcept { return &co[(nrh-nrl+1)*(nch-ncl+1)]; }

    /** const pointer to the first element accessible element */
    const T *begin() const noexcept { return &co[0]; }

    /** const pointer to the one-past the last element */
    const T *end() const noexcept { return &co[(nrh-nrl+1)*(nch-ncl+1)]; }

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
            n_row{nrh-nrl+1}, n_col{nch-ncl+1}, co { new T[(nrh-nrl+1)*(nch-ncl+1)]{} } {}

    Matrix(const std::size_t r, const std::size_t c): Matrix{r,1,c,1} {}

    /** range-checked access operator */
    T &at(const std::size_t i, const std::size_t j) {
        GEO_ERROR_IN_RANGE(i, nrl, nrh);
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return co[(i-nrl)*n_col+(j-ncl)];
    }

    /** range-checked access operator */
    const T &at(const std::size_t i, const std::size_t j) const {
        GEO_ERROR_IN_RANGE(i, nrl, nrh);
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return co[(i-nrl)*n_col+(j-ncl)];
    }

    /**
    * access operator. When the code is compiled in debug mode, it performes
    * a range check. No check is done when the code is compiled in release mode.
    */
    T& operator()(const std::size_t i, const std::size_t j) {
#ifndef NDEBUG
        return at(i,j);
#else
        return co[(i-nrl)*n_col+(j-ncl)];
#endif
    }
    /**
        * access operator. When the code is compiled in debug mode, it performes
        * a range check. No check is done when the code is compiled in release mode.
        */
    const T& operator()(const std::size_t i, const std::size_t j) const {
#ifndef NDEBUG
        return at(i,j);
#else
        return co[(i-nrl)*n_col+(j-ncl)];
#endif
    }
//    T& operator()(const std::size_t i, const std::size_t j) noexcept {return co[(i-nrl)*n_col+(j-ncl)]; }
// SE SCRIVEVO   EXPECT_NO_THROW(m.at(2,2)); MI DICEVA CHE IL TEST ERA GIUSTO!!!

    /** set all elements of the vector to @p v
        * this is useful to reinizialize all the elements of the vector to zero
        * my_matrix = 0
        */
    Matrix<T> &operator=(const T v) {
        for (auto &x : *this)
            x = v;
        return *this;
    }

    std::size_t nrh, nrl; // rows
    std::size_t nch, ncl; // columns

    std::size_t n_row;
    std::size_t n_col;

    std::unique_ptr<T[]> co;

//private:

};
#endif // GEOTOP_MATRIX_H
