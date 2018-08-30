//
// Created by elisa on 2018/07/02.
//

#ifndef GEOTOP_MATRIX_H
#define GEOTOP_MATRIX_H

#include "geotop_asserts.h"
#include <memory>
#include "rowview.h"

template <class T> class Matrix {
    //private:
    // ...

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
    // ............................................................................|-->initialize all elements to 0

    Matrix(const std::size_t r, const std::size_t c): Matrix{r,1,c,1} {}

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
        * access operator. When the code is compiled in debug mode, it performes
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

    /** set all elements of the vector to @p v
        * this is useful to reinizialize all the elements of the vector to zero
        * my_matrix = 0
        */
    Matrix<T> &operator=(const T v) {
        for (auto &x : *this)
            x = v;
        return *this;
    }

    /**
       * Copy constructor
       */
    Matrix(const Matrix<T> &m)
            : nrl{m.nrl}, nrh{m.nrh}, ncl{m.ncl}, nch{m.nch}, n_col{nch-ncl+1}, co{new T[(nrh-nrl+1)*(nch-ncl+1)]} {
        for (std::size_t i=0; i<(nrh-nrl+1)*(nch-ncl+1); ++i)
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


    RowView<T> row(const std::size_t i) {
        GEO_ASSERT_IN_RANGE(i, nrl, nrh);
        return RowView<T> { &co[(i-nrl)*n_col], nch, ncl }; // RowView<T> is constructed with the passed nch and ncl
    };

};
// ----------------------------------------------------------------------------------------------------------------
#endif // GEOTOP_MATRIX_H
