//
// Created by elisa on 2018/07/27.
//

#ifndef GEOTOP_MATRIXVIEW_H
#define GEOTOP_MATRIXVIEW_H

#include "geotop_asserts.h"
#include <cassert>
#include <exception>
#include <memory>
#include <sstream>

#include "rowview.h"

template <class T> class MatrixView {
public:
    T *elem;
    std::size_t nrh;
    std::size_t nrl;
    std::size_t nch;
    std::size_t ncl;

    RowView<T> operator[](const std::size_t i) noexcept {
        return RowView<T>{elem, nch, ncl};
    }
    const RowView<T> operator[](const std::size_t i) const noexcept {
        return RowView<T>{elem, nch, ncl};
    }

    T &operator()(const std::size_t i, const std::size_t j) noexcept {
        return elem[(i-nrl)*(nch-ncl+1) + (j-ncl)];
    }

    const T &operator()(const std::size_t i, const std::size_t j) const noexcept {
        return elem[(i-nrl)*(nch-ncl+1) + (j-ncl)];
    }
};

#endif // GEOTOP_MATRIXVIEW_H
