//
// Created by elisa on 2018/07/27.
//

#ifndef GEOTOP_MATRIXVIEW_H
#define GEOTOP_MATRIXVIEW_H

#include "geotop_asserts.h"
#include "rowview.h"

template <class T> class MatrixView {
    /** class used to access a matrix of a Tensor<T> */
public:
    /** the actual data */
    T *elem;

    /** lower and upper bounds */
    std::size_t nrh, nrl; // rows
    std::size_t nch, ncl; // columns

    /** subscripting operator (non-checked) */
    RowView<T> operator[](const std::size_t i) noexcept {
        return RowView<T>{elem, nch, ncl};
    }

    /** subscripting operator (non-checked) */
    const RowView<T> operator[](const std::size_t i) const noexcept {
        return RowView<T>{elem, nch, ncl};
    }

    /** range-checked access operator */
    T &at(const std::size_t i, const std::size_t j) {
        GEO_ERROR_IN_RANGE(i, nrl, nrh);
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return elem[(i-nrl)*(nch-ncl+1)+(j-ncl)];
    }

    /** range-checked access operator */
    const T &at(const std::size_t i, const std::size_t j) const {
        GEO_ERROR_IN_RANGE(i, nrl, nrh);
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return elem[(i-nrl)*(nch-ncl+1)+(j-ncl)];
    }

    /**
* access operator. When the code is compiled in debug mode, it performs
* a range check. No check is done when the code is compiled in release mode.
*/
    T &operator()(const std::size_t i, const std::size_t j)
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(i,j);
#else
        return elem[(i-nrl)*(nch-ncl+1)+(j-ncl)];
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
        return elem[(i-nrl)*(nch-ncl+1)+(j-ncl)];
#endif
    }

};

#endif // GEOTOP_MATRIXVIEW_H
