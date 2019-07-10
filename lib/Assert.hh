// Assert.hh
// T. M. Kelley
// Jan. 10, 2010
// Header for Assert
// (c) Copyright 2011 LANSLLC all rights reserved.

#pragma once

#include <sstream>
#include <stdexcept>

namespace nut {
/** call resize in a way that, if the alloc fails, makes some noise on
 * the way out. */
template <typename vec_t>
static inline void
resize_vec(vec_t & v,
           size_t const s,
           char const * const fname,
           char const * const elemName)
{
  try {
    v.resize(s);
  } catch(std::exception & exc) {
    std::stringstream errstr;
    errstr << fname << ": failed to allocate memory for " << s << " "
           << elemName << "s.\n";
    throw exc;
  }
}

static inline void
Insist(bool const cond, char const * const msg)
{
  if(!cond) { throw std::runtime_error(msg); }
}

#ifdef REQUIRE_ON

/** Require: a condition that, if not true, is irrecoverable.
 * Uses char * in the interface to maintain performance when
 * assertions are turned off. (If we used std::string, the caller
 * would still construct that string, even if the body of the
 * function is empty.) */
static inline void
Require(bool const cond, char const * const msg)
{
  if(!cond) { throw std::runtime_error(msg); }
}

/** Check that something is true. */
static inline void
Check(bool const cond, char const * const msg)
{
  if(!cond) { throw std::runtime_error(msg); }
}

/** Two named things are equal. */
template <typename eq_t>
inline void
Equal(eq_t const & a,
      eq_t const & b,
      char const * const a_name_cstr,
      char const * const b_name_cstr)
{
  std::string const a_name(a_name_cstr);
  std::string const b_name(b_name_cstr);
  std::string const errstr(a_name + " != " + b_name);
  Require(a == b, errstr.c_str());
}

/** min < x < max */
template <typename comp_t>
inline void
InOpenRange(comp_t const & x,
            comp_t const & min,
            comp_t const & max,
            char const * const name)
{
  std::stringstream errstr;
  bool cond = x > min && x < max;
  if(!cond) {
    std::stringstream errstr;
    errstr << name << " (" << x << ") was not in range (" << min << "," << max
           << ")";
    Require(false, errstr.str().c_str());
  }
}  // InOpenRange

/** x > min */
template <typename gt_t>
inline void
GreaterThan(gt_t const & x, gt_t const & min, char const * const name)
{
  bool cond = x > min;
  if(!cond) {
    std::stringstream errstr;
    errstr << name << " (" << x << ") <= " << min;
    Require(false, errstr.str().c_str());
  }
}  // GreaterThan

/** x >= min */
template <typename gt_t>
inline void
GreaterEqual(gt_t const & x, gt_t const & min, char const * const name)
{
  bool cond = x >= min;
  if(!cond) {
    std::stringstream errstr;
    errstr << name << " (" << x << ") <= " << min;
    Require(false, errstr.str().c_str());
  }
}  // GreaterEqual

/** x < max */
template <typename gt_t>
inline void
LessThan(gt_t const & x, gt_t const & max, char const * const name)
{
  bool cond = x < max;
  if(!cond) {
    std::stringstream errstr;
    errstr << name << " (" << x << ") >= " << max;
    Require(false, errstr.str().c_str());
  }
}  // LessThan

/** x < max */
template <typename gt_t>
inline void
LessThan(gt_t const & x,
         gt_t const & max,
         char const * const name,
         char const * const max_name)
{
  bool cond = x < max;
  if(!cond) {
    std::stringstream errstr;
    errstr << name << " (" << x << ") >= " << max_name << " (" << max << ").";
    Require(false, errstr.str().c_str());
  }
}  // LessThan

#else

#define Require Require_off
#define Check Check_off
#define Equal Equal_off
#define InOpenRange InOpenRange_off
#define GreaterThan GreaterThan_off
#define GreaterEqual GreaterEqual_off
#define LessThan LessThan_off

static inline void
Require(bool const, char const * const)
{
}

static inline void
Check(bool const, char const * const)
{
}

template <typename eq_t>
void
Equal(eq_t const &, eq_t const &, char const * const, char const * const)
{
}

template <typename comp_t>
void
InOpenRange(comp_t const &, comp_t const &, comp_t const &, char const * const)
{
}

template <typename gt_t>
void
GreaterThan(gt_t const &, gt_t const &, char const * const)
{
}

template <typename gt_t>
inline void
LessThan(gt_t const &, gt_t const &, char const * const)
{
}

template <typename gt_t>
inline void
LessThan(gt_t const &, gt_t const &, char const * const, char const * const)
{
}

template <typename gt_t>
void
GreaterEqual(gt_t const &, gt_t const &, char const * const)
{
}

#endif /* #ifdef REQUIRE_ON */

}  // namespace nut
