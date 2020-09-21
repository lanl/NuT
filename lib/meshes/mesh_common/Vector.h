// Vector.h
// T. M. Kelley
// Dec 21, 2018
// (c) Copyright 2018 LANSLLC, all rights reserved

#pragma once

#include "base/soft_equiv.h"
#include <array>
#include <cmath>
#include <type_traits>

namespace nut_mesh {

// fwd declarations...
struct Vector;
// struct Vector2;
struct Vector1;

inline Vector
operator-(Vector const & l, Vector const & r);
// inline Vector2
// operator-(Vector2 const & l, Vector2 const & r);
inline Vector1
operator-(Vector1 const & l, Vector1 const & r);
// template<typename FP_T>
// inline Vector4<FP_T> operator-(Vector4<FP_T> const & l,
//   Vector4<FP_T> const & r);

/**\brief A 3D Cartesian vector */
struct Vector {
  using value_type = double;

  static constexpr size_t dim{3};

  // interface
  constexpr static size_t size() { return dim; }

  double norm() const { return sqrt(this->dot(*this)); }

  Vector & operator=(Vector const &) = default;

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
    if(nut_mesh::soft_equiv(0.0, norm_u_sq)) { return Vector{0.0, 0.0, 0.0}; }
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

  friend std::ostream & operator<<(std::ostream & s, Vector const & v)
  {
    s << "{";
    for(uint32_t d = 0; d < (dim - 1); ++d) { s << v.data_[d] << ","; }
    s << v.data_[dim - 1] << "}";
    return s;
  }

  friend std::istream & operator>>(std::istream & s, Vector & v)
  {
    // read two possible formats
    // <w>vx<w>vy<w>vz
    //       or
    // {vx,vy,vz}
    char c;
    s >> c;
    if(c == '{') {
      for(uint32_t d = 0; d < dim; ++d) { s >> v.data_[d] >> c; }
    }
    else {
      s.putback(c);
      for(uint32_t d = 0; d < dim; ++d) { s >> v.data_[d]; }
    }
    return s;
  }

  auto begin() const { return data_.begin(); }

  auto end() const { return data_.end(); }

  // data
  std::array<double, 3> data_;
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

  using value_type = FP_T;

  static constexpr size_t dim{4};

  // interface
  constexpr static size_t size() { return dim; }

  FP_T norm() const { return sqrt(this->dot(*this)); }

  bool operator==(Vector4 const & r) const
  {
    return data_[0] == r.data_[0] && data_[1] == r.data_[1] &&
           data_[2] == r.data_[2] && data_[3] == r.data_[3];
  }

  bool operator!=(Vector4 const & r) const { return !(*this == r); }

  Vector4 & operator=(Vector4 const &) = default;

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
    if(nut_mesh::soft_equiv(0.0, norm_u_sq)) {
      return Vector4{0.0, 0.0, 0.0, 0.0};
    }
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

  friend std::ostream & operator<<(std::ostream & s, Vector4 const & v)
  {
    s << "{";
    for(uint32_t d = 0; d < (dim - 1); ++d) { s << v.data_[d] << ","; }
    s << v.data_[dim - 1] << "}";
    return s;
  }

  friend std::istream & operator>>(std::istream & s, Vector4 & v)
  {
    // read two possible formats
    // <w>vx<w>vy<w>vz
    //       or
    // {vx,vy,vz}
    char c;
    s >> c;
    if(c == '{') {
      for(uint32_t d = 0; d < dim; ++d) { s >> v.data_[d] >> c; }
    }
    else {
      s.putback(c);
      for(uint32_t d = 0; d < dim; ++d) { s >> v.data_[d]; }
    }
    return s;
  }

  auto begin() const { return data_.begin(); }

  auto end() const { return data_.end(); }

  // data
  std::array<FP_T, 4> data_;
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
template <typename FP_T>
struct Vector2 {
  static_assert(std::is_floating_point_v<FP_T>,
                "Template parameter FP_T must be a floating point type");

  using value_type = FP_T;

  static constexpr size_t dim{2};

  // ctors
  Vector2(FP_T d1, FP_T d2) {
    data_[0] = d1;
    data_[1] = d2;
    return;
  }

  // interface
  constexpr static size_t size() { return dim; }

  FP_T norm() const { return sqrt(this->dot(*this)); }

  bool operator==(Vector2 const & r) const
  {
    return data_[0] == r.data_[0] && data_[1] == r.data_[1];
  }

  bool operator!=(Vector2 const & r) const { return !(*this == r); }

  Vector2 & operator=(Vector2 const &) = default;

  FP_T dot(Vector2 const & other) const
  {
    FP_T const d = data_[0] * other.data_[0] + data_[1] * other.data_[1];
    return d;
  }

  /**\brief New vector is old vector scaled by d */
  Vector2 scale(FP_T const d) const { return {d * data_[0], d * data_[1]}; }

  /**\brief Scale this vector (mutate) */
  void scale_this(FP_T const d)
  {
    data_[0] *= d;
    data_[1] *= d;
    return;
  }

  Vector2 operator*(FP_T const d) const
  {
    Vector2 t = this->scale(d);
    return t;
  }

  FP_T const & operator[](size_t i) const
  {
    // Require( i < dim);
    return data_[i];
  }

  FP_T & operator[](size_t i)
  {
    // Require( i < dim);
    return data_[i];
  }

  /**\brief Extract component parallel to u */
  inline Vector2 parallel(Vector2 const & u) const
  {
    FP_T norm_u_sq{u.dot(u)};
    if(nut_mesh::soft_equiv(0.0, norm_u_sq)) { return Vector2{0.0, 0.0}; }
    FP_T magnitude{this->dot(u)};
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

  friend std::ostream & operator<<(std::ostream & s, Vector2 const & v)
  {
    s << "{";
    for(uint32_t d = 0; d < (dim - 1); ++d) { s << v.data_[d] << ","; }
    s << v.data_[dim - 1] << "}";
    return s;
  }

  friend std::istream & operator>>(std::istream & s, Vector2 & v)
  {
    // read two possible formats
    // <w>vx<w>vy<w>vz
    //       or
    // {vx,vy,vz}
    char c;
    s >> c;
    if(c == '{') {
      for(uint32_t d = 0; d < dim; ++d) { s >> v.data_[d] >> c; }
    }
    else {
      s.putback(c);
      for(uint32_t d = 0; d < dim; ++d) { s >> v.data_[d]; }
    }
    return s;
  }

  auto begin() const { return data_.begin(); }

  auto end() const { return data_.end(); }

  // data
  std::array<FP_T, dim> data_;
};  // struct Vector2

template <typename FP_T>
inline Vector2<FP_T> operator*(double const d, Vector2<FP_T> const & v)
{
  Vector2 t = v.scale(d);
  return t;
}

template <typename FP_T>
inline Vector2<FP_T>
operator+(Vector2<FP_T> const & l, Vector2<FP_T> const & r)
{
  return {l.data_[0] + r.data_[0], l.data_[1] + r.data_[1]};
}

template <typename FP_T>
inline Vector2<FP_T>
operator-(Vector2<FP_T> const & l, Vector2<FP_T> const & r)
{
  return {l.data_[0] - r.data_[0], l.data_[1] - r.data_[1]};
}

// -----------------------------------------------------------------------------
//
// -------------------------------- end Vector2 --------------------------------

/**\brief A 1D Cartesian vector */
struct Vector1 {
  using value_type = double;

  static constexpr size_t dim{1};

  // interface
  constexpr static size_t size() { return dim; }

  double norm() const { return sqrt(this->dot(*this)); }

  bool operator==(Vector1 const & r) const { return data_[0] == r.data_[0]; }

  bool operator!=(Vector1 const & r) const { return !(*this == r); }

  Vector1 & operator=(Vector1 const &) = default;

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
    if(nut_mesh::soft_equiv(0.0, norm_u_sq)) { return Vector1{0.0}; }
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

  friend std::ostream & operator<<(std::ostream & s, Vector1 const & v)
  {
    s << "{";
    for(uint32_t d = 0; d < (dim - 1); ++d) { s << v.data_[d] << ","; }
    s << v.data_[dim - 1] << "}";
    return s;
  }

  friend std::istream & operator>>(std::istream & s, Vector1 & v)
  {
    // read two possible formats
    // <w>vx<w>vy<w>vz
    //       or
    // {vx,vy,vz}
    char c;
    s >> c;
    if(c == '{') {
      for(uint32_t d = 0; d < dim; ++d) { s >> v.data_[d] >> c; }
    }
    else {
      s.putback(c);
      for(uint32_t d = 0; d < dim; ++d) { s >> v.data_[d]; }
    }
    return s;
  }

  auto begin() const { return data_.begin(); }

  auto end() const { return data_.end(); }

  // data
  std::array<double, dim> data_;
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

}  // namespace nut_mesh

// End of file
