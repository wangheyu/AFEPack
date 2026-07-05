#ifndef _MATRIXBASE_H_
#define _MATRIXBASE_H_

#include <vector>
#include <cblas.h>
#include <lapacke.h>
#include <stdexcept>
#include <algorithm>
#include <concepts>

// Concepts for constraints
template <typename T>
concept FloatingPoint = std::floating_point<T>;

template <typename T>
concept Integral = std::integral<T>;

// Base class using CRTP
template <typename Derived, FloatingPoint T>
class MatrixBase {
public:
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<Derived&>(*this); }

  size_t rows() const { return derived().rows(); }
  size_t cols() const { return derived().cols(); }
  T operator()(size_t i, size_t j) const { return derived()(i, j); }
  T& operator()(size_t i, size_t j) { return derived()(i, j); }

  Derived operator+(const Derived& other) const { return derived().add(other); }
  Derived operator*(const Derived& other) const { return derived().multiply(other); }
};

#endif//_MATRIXBASE_H_
