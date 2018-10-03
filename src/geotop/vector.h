//
// Created by alberto on 4/19/18.
//

#ifndef GEOTOP_VECTOR_H
#define GEOTOP_VECTOR_H

#include "geotop_asserts.h"
#include <memory> // to use std::unique_ptr

template <typename T> class Vector {

private:
    /** size of the array */
    std::size_t _size;

public:
    /** lower bound */
    std::size_t nl;
    /** upper bound */
    std::size_t nh;

    /** the actual data */
    std::unique_ptr<T[]> co;

    /**
     * @return the size of the vector
     */
    std::size_t size() const noexcept { return _size; }

    /** pointer to the first accessible element */
    T *begin() noexcept { return &co[nl]; }

    /** pointer to the one-past the last element */
    T *end() noexcept { return &co[nh + 1]; }

    /** const pointer to the first element accessible element */
    const T *begin() const noexcept { return &co[nl]; }

    /** const pointer to the one-past the last element */
    const T *end() const noexcept { return &co[nh + 1]; }

    /** subscripting operator (non-checked) */
    T &operator[](const std::size_t i) noexcept { return co[i]; }

    /** subscripting operator (non-checked) */
    const T &operator[](const std::size_t i) const noexcept { return co[i]; }

    /** range-checked access operator */
    T &at(const std::size_t i) {
        GEO_ERROR_IN_RANGE(i, nl, nh);
        return (*this)[i];
    }

    /** range-checked access operator */
    const T &at(const std::size_t i) const {
        GEO_ERROR_IN_RANGE(i, nl, nh);
        return (*this)[i];
    }

    /**
     * access operator. When the code is compiled in debug mode, it performs
     * a range check. No check is done when the code is compiled in release mode.
     */
    T &operator()(const std::size_t i)
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(i);
#else
        return (*this)[i];
#endif
    }

    /**
     * access operator. When the code is compiled in debug mode, it performes
     * a range check. No check is done when the code is compiled in release mode.
     */
    const T &operator()(const std::size_t i) const
#ifdef NDEBUG
    noexcept
#endif
    {
#ifndef NDEBUG
        return at(i);
#else
        return (*this)[i];
#endif
    }

    /** destructor. default is fine */
    ~Vector() = default;

    /** default constructor is deleted */
    Vector() = delete;

    /**
     * constructor
     * @param ub upper bound of the range of accessible elements
     * @param lb lower bound of the range of accessible elememts
     *
     * you can access elements in the range [l,n] boundaries included
     */
    explicit Vector(const std::size_t ub, const std::size_t lb = 1)
            : _size{ub - lb + 1}, nl{lb}, nh{ub}, co{new T[nh + 1]{}} {} // initialize all elements to 0

    /**
     * Copy constructor
     */
    Vector(const Vector<T> &v)
            : _size{v._size}, nl{v.nl}, nh{v.nh}, co{new T[nh + 1]} {
        for (auto i = nl; i <= nh; ++i)
            (*this)[i] = v[i];
    }

    /** Move constructor */
    Vector(Vector<T> &&v) = default;

    /** Move assignment */
    Vector<T>& operator=(Vector<T> &&v) = default;


    /** Copy assignment */
    Vector<T> &operator=(const Vector<T> &v) {
        co.reset(); // release acquired memory
        *this = Vector<T>{v}; // use move assignment and copy constructor
        return *this;
    }

    /** set all elements of the vector to @p v
     * this is useful to reinizialize all the elements of the vector to zero
     * my_vector = 0
     */
    Vector<T> &operator=(const T v) {
        for (auto &x : *this)
            x = v;
        return *this;
    }

    /**
     * Component-wise summation
     */
    Vector<T> &operator+=(const Vector<T> &v) {

        GEO_ASSERT_EQ(nl, v.nl) << "vector length mismatch\n";
        GEO_ASSERT_EQ(nh, v.nh) << "vector length mismatch\n";

        for (auto i = nl; i <= nh; ++i)
            (*this)[i] += v[i];
        return *this;
    }


};

#endif // GEOTOP_VECTOR_H
