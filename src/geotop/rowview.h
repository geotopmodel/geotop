//
// Created by elisa on 2018/07/27.
//

#ifndef GEOTOP_ROWVIEW_H
#define GEOTOP_ROWVIEW_H

#include "geotop_asserts.h"

template <class T> class RowView {
    /** This class is used to access a row of a Matrix<T> */
public:

    /** the actual data */
    T *elem;

    /** lower and upper bounds */
    std::size_t nch;
    std::size_t ncl;

    /** range-checked access operator */
    T &at(const std::size_t j) {
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return (*this)[j];
    }

    /** range-checked access operator */
    const T &at(const std::size_t j) const {
        GEO_ERROR_IN_RANGE(j, ncl, nch);
        return (*this)[j];
    }

    T &operator[](const std::size_t j) noexcept {
        return elem[j - ncl];
    }

    const T &operator[](const std::size_t j) const noexcept {
        return elem[j - ncl];
    }

    /**
* access operator. When the code is compiled in debug mode, it performs
* a range check. No check is done when the code is compiled in release mode.
*/
    T& operator()(const std::size_t j)
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(j);
#else
        return (*this)[j];

#endif
    }

    /**
        * access operator. When the code is compiled in debug mode, it performs
        * a range check. No check is done when the code is compiled in release mode.
        */
    const T& operator()(const std::size_t j) const
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(j);
#else
        return (*this)[j];
#endif
    }

};
#endif // GEOTOP_ROWVIEW_H
