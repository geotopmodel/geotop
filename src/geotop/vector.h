//
// Created by alberto on 4/19/18.
//

#ifndef GEOTOP_VECTOR_H
#define GEOTOP_VECTOR_H

#include "geotop_asserts.h"
#include <cassert>
#include <exception>
#include <memory>
#include <sstream>

template <typename T> struct Vector {

  /**
   * @return the size of the vector
   */
  std::size_t size() const noexcept { return _size; }

  /** pointer to the first element accessible element */
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
    check(i);
    return co[i];
  }

  /** range-checked access operator */
  const T &at(const std::size_t i) const {
    check(i);
    return co[i];
  }

  /** destructor. default is fine */
  ~Vector() = default;

  /** default constructor is deleted */
  Vector() = delete;

  /**
   * constructor
   * @param n upper bound of the range of accessible elements
   * @param l lower bound of the range of accessible elememts
   *
   * you can access elements in the range [l,n] boundaries included
   */
  explicit Vector(const std::size_t n, const std::size_t l = 1)
      : nl{l}, nh{n}, co{new T[nh + 1]{}}, _size{nh - nl + 1} {}

  /**
   * Copy constructor
   * @param v
   */
  Vector(const Vector<T> &v)
      : nl{v.nl}, nh{v.nh}, co{new T[nh + 1]}, _size{v._size} {
    for (auto i = nl; i <= nh; ++i)
      co[i] = v.co[i];
  }

  /** Copy assignment */
  Vector<T> &operator=(const Vector<T> &v) {
    nl = v.nl;
    nh = v.nh;
    co.reset(new T[nh + 1]);
    for (std::size_t i = nl; i <= nh; ++i)
      co[i] = v.co[i];
    return *this;
  }

  /** set all elements of the vector to @p v
   * this is useful to reinizialize all the elements of the vector to zero
   * my_vector = 0
   */
  Vector<T> &operator=(const T v) {
    for (auto i = nl; i <= nh; ++i)
      co[i] = v;
    return *this;
  }

  /**
   * Component-wise summation
   */
  Vector<T> &operator+=(const Vector<T> &v) {

    GEO_ASSERT_EQ(nl, v.nl) << "vector length mismatch\n";
    GEO_ASSERT_EQ(nh, v.nh) << "vector length mismatch\n";

    for (auto i = nl; i <= nh; ++i)
      co[i] += v.co[i];
    return *this;
  }

  /** lower bound */
  std::size_t nl;
  /** upper bound */
  std::size_t nh;

  /** the actual data */
  std::unique_ptr<T[]> co;

private:
  /** size of the array */
  std::size_t _size;

  /**
   * helper function used to check if @param i is within the valid range
   */
  void check(const std::size_t i) const {
    if (i < nl || i > nh) {
      std::ostringstream os;
      os << i << " does not belong to range [" << nl << ", " << nh << "]."
         << std::endl;
      throw std::out_of_range{os.str()};
    }
  }
};

#endif // GEOTOP_VECTOR_H
