//
// Created by elisa on 2018/07/27.
//

#ifndef GEOTOP_ROWVIEW_H
#define GEOTOP_ROWVIEW_H

#include "geotop_asserts.h"
#include <cassert>
#include <exception>
#include <memory>
#include <sstream>

template <class T> class RowView {
public:
    T *elem;
    std::size_t nch;
    std::size_t ncl;

    T &operator[](const std::size_t j) noexcept {
        return elem[j - ncl];
    }

    const T &operator[](const std::size_t j) const noexcept {
        return elem[j - ncl];
    }

};
#endif // GEOTOP_ROWVIEW_H
