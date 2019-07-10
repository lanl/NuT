// Vector.h
// T. M. Kelley
// Dec 21, 2018
// (c) Copyright 2018 LANSLLC, all rights reserved

#pragma once

#include "soft_equiv.hh"
#include <array>
#include <cmath>
#include <iostream>
#include <type_traits>

namespace nut {

// fwd declarations...
struct Vector;
struct Vector2;
struct Vector1;
// ... enable these declarations
inline Vector
operator-(Vector const & l, Vector const & r);
inline Vector2
operator-(Vector2 const & l, Vector2 const & r);
inline Vector1
operator-(Vector1 const & l, Vector1 const & r);

/**\brief A 3D Cartesian vector */
struct Vector {
  static constexpr size_t dim{3};

  using value_type = double;

  using iterator = typename std::array<value_type, dim>::iterator;

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  // interface
  constexpr static size_t size() {
    return dim;
  }

  double norm() const { return sqrt(this->dot(*this)); }

  bool operator==(Vector const & r) const
  {
    return data_[0] == r.data_[0] && data_[1] == r.data_[1] &&
           data_[2] == r.data_[2];
  }

  bool operator!=(Vector const & r) const { return !(*this == r); }

  double dot(Vector const & other) const
  {
    double const d = data_[0] * other.data_[0] + data_[1] * other.data_[1] +
                     data_[2] * other.data_[2];
    return d;
  }

  /**\brief New vector is old vector scaled by d */
  Vector scale(double const d) const
  {
    return {d * data_[0], d * data_[1], d * data_[2]};
  }

  /**\brief Scale this vector (mutate) */
  void scale_this(double const d)
  {
    data_[0] *= d;
    data_[1] *= d;
    data_[2] *= d;
    return;
  }

  Vector operator*(double const d) const
  {
    Vector t = this->scale(d);
    return t;
  }

  double const & operator[](size_t i) const
  {
    // Require( i < 3);
    return data_[i];
  }

  double & operator[](size_t i)
  {
    // Require( i < 3);
    return data_[i];
  }

  /**\brief Extract component parallel to u */
  inline Vector parallel(Vector const & u) const
  {
    double norm_u_sq{u.dot(u)};
    if(nut::soft_equiv(0.0, norm_u_sq)) { return Vector{0.0, 0.0, 0.0}; }
    double magnitude{this->dot(u)};
    Vector parallel{u.scale(magnitude / norm_u_sq)};
    return parallel;
  }  // parallel

  /**\brief Extract component perpendicular to u */
  inline Vector perpendicular(Vector const & u) const
  {
    Vector const parl{this->parallel(u)};
    Vector const perp{*this - parl};
    return perp;
  }

  /**\brief Reflect vector from plane described by normal */
  inline Vector reflect(Vector const & n) const
  {
    Vector parl{this->parallel(n)};
    Vector perp{this->perpendicular(n)};
    // perpendicular component unchanged, parallel component reversed
    return perp - parl;
  }

  // data
  std::array<double, 3> data_;

  friend std::ostream & operator<<(std::ostream & s, Vector const & v)
  {
    s << "{";
    for(uint32_t d = 0; d < 2; ++d) { s << v.data_[d] << ","; }
    s << v.data_[dim - 1] << "}";
    return s;
  }
};  // struct Vector

inline Vector operator*(double const d, Vector const & v)
{
  Vector t = v.scale(d);
  return t;
}

inline Vector
operator-(Vector const & l, Vector const & r)
{
  return {l.data_[0] - r.data_[0], l.data_[1] - r.data_[1],
          l.data_[2] - r.data_[2]};
}

inline Vector
operator+(Vector const & l, Vector const & r)
{
  return {l.data_[0] + r.data_[0], l.data_[1] + r.data_[1],
          l.data_[2] + r.data_[2]};
}

// ------------------------------- end of Vector -------------------------------

/**\brief A 4D vector */
template <typename FP_T>
struct Vector4 {
  static_assert(std::is_floating_point_v<FP_T>,
                "Template parameter FP_T must be a floating point type");

  static constexpr size_t dim{4};

  using value_type = double;

  using iterator = typename std::array<value_type, dim>::iterator;

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  // interface
  constexpr static size_t size() {
    return dim;
  }

  FP_T norm() const { return sqrt(this->dot(*this)); }

  bool operator==(Vector4 const & r) const
  {
    return data_[0] == r.data_[0] && data_[1] == r.data_[1] &&
           data_[2] == r.data_[2] && data_[3] == r.data_[3];
  }

  bool operator!=(Vector4 const & r) const { return !(*this == r); }

  FP_T dot(Vector4 const & other) const
  {
    FP_T const d = data_[0] * other.data_[0] + data_[1] * other.data_[1] +
                   data_[2] * other.data_[2] + data_[3] * other.data_[3];
    return d;
  }

  /**\brief New vector is old vector scaled by d */
  Vector4 scale(FP_T const d) const
  {
    return {d * data_[0], d * data_[1], d * data_[2], d * data_[3]};
  }

  /**\brief Scale this vector (mutate) */
  void scale_this(FP_T const d)
  {
    data_[0] *= d;
    data_[1] *= d;
    data_[2] *= d;
    data_[3] *= d;
    return;
  }

  Vector4 operator*(FP_T const d) const
  {
    Vector4 t = this->scale(d);
    return t;
  }

  FP_T const & operator[](size_t i) const
  {
    // Require( i < 4);
    return data_[i];
  }

  FP_T & operator[](size_t i)
  {
    // Require( i < 4);
    return data_[i];
  }

  /**\brief Extract component parallel to u */
  inline Vector4 parallel(Vector4 const & u) const
  {
    FP_T norm_u_sq{u.dot(u)};
    if(nut::soft_equiv(0.0, norm_u_sq)) { return Vector4{0.0, 0.0, 0.0, 0.0}; }
    FP_T magnitude{this->dot(u)};
    Vector4 parallel_comp{u.scale(magnitude / norm_u_sq)};
    return parallel_comp;
  }  // parallel

  /**\brief Extract component perpendicular to u */
  inline Vector4 perpendicular(Vector4 const & u) const
  {
    Vector4 const parl{this->parallel(u)};
    Vector4 const perp{*this - parl};
    return perp;
  }

  /**\brief Reflect vector from plane described by normal */
  inline Vector4 reflect(Vector4 const & n) const
  {
    Vector4 parl{this->parallel(n)};
    Vector4 perp{this->perpendicular(n)};
    // perpendicular component unchanged, parallel component reversed
    return perp - parl;
  }

  // data
  std::array<FP_T, 4> data_;

  friend std::ostream & operator<<(std::ostream & s, Vector4 const & v)
  {
    s << "{";
    for(uint32_t d = 0; d < 3; ++d) { s << v.data_[d] << ","; }
    s << v.data_[dim - 1] << "}";
    return s;
  }
};  // struct Vector4

template <typename FP_T>
inline Vector4<FP_T> operator*(FP_T const d, Vector4<FP_T> const & v)
{
  Vector4<FP_T> t = v.scale(d);
  return t;
}

template <typename FP_T>
inline Vector4<FP_T>
operator-(Vector4<FP_T> const & l, Vector4<FP_T> const & r)
{
  return {l.data_[0] - r.data_[0], l.data_[1] - r.data_[1],
          l.data_[2] - r.data_[2], l.data_[3] - r.data_[3]};
}

template <typename FP_T>
inline Vector4<FP_T>
operator+(Vector4<FP_T> const & l, Vector4<FP_T> const & r)
{
  return {l.data_[0] + r.data_[0], l.data_[1] + r.data_[1],
          l.data_[2] + r.data_[2], l.data_[3] + r.data_[3]};
}
// -------------------------------- end Vector4 --------------------------------

/**\brief A 2D vector */
struct Vector2 {
  static constexpr size_t dim{2};

  using value_type = double;

  using iterator = typename std::array<value_type, dim>::iterator;

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  // interface
  constexpr static size_t size() {
    return dim;
  }

  double norm() const { return sqrt(this->dot(*this)); }

  bool operator==(Vector2 const & r) const
  {
    return data_[0] == r.data_[0] && data_[1] == r.data_[1];
  }

  bool operator!=(Vector2 const & r) const { return !(*this == r); }

  double dot(Vector2 const & other) const
  {
    double const d = data_[0] * other.data_[0] + data_[1] * other.data_[1];
    return d;
  }

  /**\brief New vector is old vector scaled by d */
  Vector2 scale(double const d) const { return {d * data_[0], d * data_[1]}; }

  /**\brief Scale this vector (mutate) */
  void scale_this(double const d)
  {
    data_[0] *= d;
    data_[1] *= d;
    return;
  }

  Vector2 operator*(double const d) const
  {
    Vector2 t = this->scale(d);
    return t;
  }

  double const & operator[](size_t i) const
  {
    // Require( i < dim);
    return data_[i];
  }

  double & operator[](size_t i)
  {
    // Require( i < dim);
    return data_[i];
  }

  /**\brief Extract component parallel to u */
  inline Vector2 parallel(Vector2 const & u) const
  {
    double norm_u_sq{u.dot(u)};
    if(nut::soft_equiv(0.0, norm_u_sq)) { return Vector2{0.0, 0.0}; }
    double magnitude{this->dot(u)};
    Vector2 parallel{u.scale(magnitude / norm_u_sq)};
    return parallel;
  }  // parallel

  /**\brief Extract component perpendicular to u */
  inline Vector2 perpendicular(Vector2 const & u) const
  {
    Vector2 const parl{this->parallel(u)};
    Vector2 const perp{*this - parl};
    return perp;
  }

  /**\brief Reflect vector from plane described by normal */
  inline Vector2 reflect(Vector2 const & n) const
  {
    Vector2 parl{this->parallel(n)};
    Vector2 perp{this->perpendicular(n)};
    // perpendicular component unchanged, parallel component reversed
    return perp - parl;
  }

  // data
  std::array<double, dim> data_;
  friend std::ostream & operator<<(std::ostream & s, Vector2 const & v)
  {
    s << "{" << v.data_[0] << "," << v.data_[1] << "}";
    return s;
  }
};  // struct Vector2

inline Vector2 operator*(double const d, Vector2 const & v)
{
  Vector2 t = v.scale(d);
  return t;
}

inline Vector2
operator+(Vector2 const & l, Vector2 const & r)
{
  return {l.data_[0] + r.data_[0], l.data_[1] + r.data_[1]};
}

inline Vector2
operator-(Vector2 const & l, Vector2 const & r)
{
  return {l.data_[0] - r.data_[0], l.data_[1] - r.data_[1]};
}

/**\brief A 1D Cartesian vector */
struct Vector1 {
  static constexpr size_t dim{1};

  using value_type = double;

  using iterator = typename std::array<value_type, dim>::iterator;

  using const_iterator = typename std::array<value_type, dim>::const_iterator;

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  const_iterator begin() const { return data_.begin(); }

  const_iterator end() const { return data_.end(); }

  // interface
  constexpr static size_t size() {
    return dim;
  }

  double norm() const { return sqrt(this->dot(*this)); }

  bool operator==(Vector1 const & r) const { return data_[0] == r.data_[0]; }

  bool operator!=(Vector1 const & r) const { return !(*this == r); }

  double dot(Vector1 const & other) const
  {
    double const d = data_[0] * other.data_[0];
    return d;
  }

  /**\brief New vector is old vector scaled by d */
  Vector1 scale(double const d) const { return {d * data_[0]}; }

  /**\brief Scale this vector (mutate) */
  void scale_this(double const d)
  {
    data_[0] *= d;
    return;
  }

  Vector1 operator*(double const d) const
  {
    Vector1 t = this->scale(d);
    return t;
  }

  double const & operator[](size_t i) const
  {
    // Require( i < dim);
    return data_[i];
  }

  double & operator[](size_t i)
  {
    // Require( i < dim);
    return data_[i];
  }

  /**\brief We're all parallel all the time */
  inline Vector1 parallel(Vector1 const & u) const
  {
    double norm_u_sq{u.dot(u)};
    if(nut::soft_equiv(0.0, norm_u_sq)) { return Vector1{0.0}; }
    return *this;
  }  // parallel

  /**\brief Extract component perpendicular to u */
  inline Vector1 perpendicular(Vector1 const & /*u*/) const
  {
    return Vector1{0.0};
  }

  /**\brief Reflect vector from plane described by normal */
  inline Vector1 reflect(Vector1 const & /*n*/) const
  {
    return Vector1{-1.0 * data_[0]};
  }

  // data
  std::array<double, dim> data_;

  friend std::ostream & operator<<(std::ostream & s, Vector1 const & v)
  {
    s << "{" << v.data_[0] << "}";
    return s;
  }
};  // struct Vector1

inline Vector1 operator*(double const d, Vector1 const & v)
{
  Vector1 t = v.scale(d);
  return t;
}

inline Vector1
operator+(Vector1 const & l, Vector1 const & r)
{
  return {l.data_[0] + r.data_[0]};
}

inline Vector1
operator-(Vector1 const & l, Vector1 const & r)
{
  return {l.data_[0] - r.data_[0]};
}

}  // namespace nut

// End of file
