// RNG.hh
// T. M. Kelley
// Dec 21, 2010
// Header for RNG: classes that provide loverly random numbers
// (c) Copyright 2010 LANSLLC all rights reserved.

#pragma once

#include "Random123/philox.h"
#include <iostream>
#include <iterator>
#include <vector>
#include <stdint.h>

namespace nut {

/*! Uses the Random123 Philox counter based RNG. Initialize with a
 * key and a particle index. Key might consist of one word seed and
 * a thread index, for instance. This class will use the 10 round
 * version, for comparability with the Haskell McPhD code. */
class Philox4x32_RNG {
public:
  typedef r123::Philox4x32 gen_t;
  typedef gen_t::ctr_type ctr_t;
  typedef gen_t::key_type key_t;
  typedef uint32_t seed_t;
  typedef seed_t const csd_t;

  Philox4x32_RNG(ctr_t const c, key_t const k) : m_c(c), m_k(k), m_idx(0)
  {
    this->advance();
  }

  /** For compatibility with STL. */
  Philox4x32_RNG() {}

  double random()
  {
    double ret = m_bank[m_idx];
    if(m_idx == 0)
      m_idx++;
    else {
      m_idx = 0;
      this->advance();
    }
    return ret;
  }

  static key_t make_key(uint32_t lo, uint32_t hi)
  {
    key_t key;
    key.v[0] = lo;
    key.v[1] = hi;
    return key;
  }

  static ctr_t make_ctr(csd_t w4, csd_t w3, csd_t w2, csd_t w1)
  {
    ctr_t c;
    c.v[0] = w4;
    c.v[1] = w3;
    c.v[2] = w2;
    c.v[3] = w1;
    return c;
  }

  friend std::ostream & operator<<(std::ostream &, Philox4x32_RNG const &);

private:
  void advance()
  {
    ctr_t out = gen(m_c, m_k);
    this->incr_ctr();
    this->bank_ctr(out);
  }

  void incr_ctr() { m_c[3]++; }

  void bank_ctr(ctr_t const ctr)
  {
    m_bank[0] = this->words_to_double(ctr[0], ctr[1]);
    m_bank[1] = this->words_to_double(ctr[2], ctr[3]);
  }

  double words_to_double(uint32_t const w1, uint32_t const w2)
  {
    int32_t const u1(w1);
    int32_t const u2(w2);
    double const m_inv_52 = 2.220446049250313080847263336181640625e-16;
    double const m_inv_53 = 1.1102230246251565404236316680908203125e-16;
    double const m_inv_32 = 2.3283064365386962890625e-10;
    return u1 * m_inv_32 + (0.5 + m_inv_53) + (u2 & 0xFFFFF) * m_inv_52;
  }

  // state
  gen_t gen;
  ctr_t m_c;
  key_t m_k;
  double m_bank[2];
  uint32_t m_idx;
};  // class Philox4x32_RNG

/* ---------------- Other RNG's -------------------*/

/** This fake RNG is used for testing in certain limited circumstances.
    It encapsulates a buffer of "random" seeds and just feeds them out
    in sequence. This is useful for testing certain situations with a
    stacked deck, for example, allowing us to make sure that events
    are correctly selected.*/
template <typename fp_t>
class Buffer_RNG {
public:
  typedef std::pair<Buffer_RNG<fp_t>, Buffer_RNG<fp_t> > new_gens;

public:
  /* get a random number */
  fp_t random();

  Buffer_RNG(fp_t const * const rns, size_t const n_rns, bool silent = false)
      : m_rns(rns, rns + n_rns),
        m_idx(0),
        m_n_issued(0),
        m_old_idx(0),
        m_n_warns(0),
        m_max_warns(2),
        m_warn_next(false),
        m_silent(silent)
  {
  }

  bool operator==(Buffer_RNG const & r) const
  {
    bool same = m_rns == r.m_rns && m_n_issued == r.m_n_issued &&
                m_old_idx == r.m_old_idx && m_n_warns == r.m_n_warns &&
                m_max_warns == r.m_max_warns && m_warn_next == r.m_warn_next;
    return same;
  }

  new_gens split() const { return new_gens(*this, *this); }

private:
  std::vector<fp_t> m_rns;  // shd be const, but can't use push_back if so

  size_t m_idx;

  size_t m_n_issued;

  // counters to track rolling over and warning about it
  size_t m_old_idx;

  uint32_t m_n_warns;

  uint32_t m_max_warns;

  bool m_warn_next;

  bool m_silent;

private:
};  // Buffer_RNG

template <typename fp_t>
fp_t
Buffer_RNG<fp_t>::random()
{
  // warn if stream rolled over on the preceding call (see below)
  if(m_warn_next && !m_silent) {
    std::cerr << "Buffer_RNG rolled over from " << m_old_idx << " to " << m_idx
              << std::endl;
    m_warn_next = false;
    m_n_warns += 1;
  }

  fp_t val = m_rns[m_idx];
  size_t const new_idx = (m_idx + 1) % m_rns.size();

  // if we roll over, then the next time we're asked for a random seed
  // (i.e. when we are recycling the stream) issue a warning
  if(new_idx < m_idx) {
    m_warn_next = true;
    m_old_idx = m_idx;
    if(m_n_warns == m_max_warns) {
      // do something effective
    }
  }

  // stream state
  m_idx = new_idx;
  m_n_issued++;

  return val;
}  // random

/** Uses L'Ecuyer's Multiplicative Linear Congruential Generator
 * {Comm. ACM, v. 31, p. 742 (1988)}, as implemented in the Haskell
 * System.Random library, to permit direct comparison. Outdated after
 * Nov-Dec 2011. */
class MLCG {
public:
  typedef std::pair<MLCG, MLCG> new_gens;

public:
  explicit MLCG(int32_t s)
  {
    if(s < 0) { s = -s; }
    int32_t q = s / 2147483562;
    int32_t s1 = s % 2147483562;
    int32_t s2 = q % 2147483398;
    this->m_s1 = s1 + 1;
    this->m_s2 = s2 + 1;
    return;
  }  // ctor

  MLCG(int32_t s1, int32_t s2) : m_s1(s1), m_s2(s2) {}

  /*!\brief form two new RNGs from this one. NB: a comment in the Haskell*
   * source says, "no statistical foundation for this!". */
  new_gens split()
  {
    int32_t new_s1 = m_s1 == 2147483562 ? 1 : m_s1 + 1;
    int32_t new_s2 = m_s2 == 1 ? 2147483398 : m_s2 - 1;
    // advance state
    this->next();
    int32_t t1 = m_s1;
    int32_t t2 = m_s2;
    return new_gens(MLCG(new_s1, t2), MLCG(t1, new_s2));
  }

  /*!\brief uniform random deviate in [0,1). */
  double random()
  {
    double x = (double)this->random_int();
    int64_t const h = int64_t(2147483647);
    int64_t const l = int64_t(-2147483648LL);
    int64_t const k = h - l + 1;
    // int64_t const k = int64_t(4294967295LL + 1LL);
    double scaled_x = 0.5 + 1.0 / (double)k * x;
    return scaled_x;
  }

  void dump_state(std::ostream & o) const
  {
    o << m_s1 << "," << m_s2;
    return;
  }

private:
  /*!\brief pull the next integer from the stream, and advance the
   * state of this object. */
  int32_t next()
  {
    int32_t const c1 = 53668;
    int32_t const d1 = 12211;
    int32_t const e1 = 40014;

    int32_t const c2 = 52774;
    int32_t const d2 = 3791;
    int32_t const e2 = 40692;

    int32_t const k1 = m_s1 / c1;
    int32_t const k2 = m_s2 / c2;
    int32_t const s1a = e1 * (m_s1 - k1 * c1) - k1 * d1;
    int32_t const s2a = e2 * (m_s2 - k2 * c2) - k2 * d2;
    int32_t const s1b = s1a < 0 ? (s1a + 2147483563) : s1a;
    int32_t const s2b = s2a < 0 ? (s2a + 2147483399) : s2a;

    this->m_s1 = s1b;
    this->m_s2 = s2b;

    int32_t const z = s1b - s2b;
    int32_t const z1 = z < 1 ? (z + 2147483562) : z;

    return z1;
  }  // next

  /* 64b to match Haskell's internal use of arbitrary-size Integer */
  int64_t random_int()
  {
    int64_t const b = 2147483561;
    int64_t const h = 2147483647LL;
    int64_t const l = -2147483648LL;
    int64_t const k = h - l + 1;
    int64_t const n = iLogBase(b, k);
    int64_t const v = this->f(n, 1);
    return l + (v % k);
  }

  int64_t iLogBase(int64_t b, int64_t i)
  {
    int64_t sum(1);
    while(i > b) {
      i = i / b;
      sum++;
    }
    return sum;
  }

  int64_t f(int64_t n, int64_t acc)
  {
    while(n > 0) {
      int64_t x = this->next();
      acc = x + acc * 2147483561;
      n--;
    }
    return acc;
  }

  int32_t m_s1, m_s2;
};  // class MLCG

/* Uses Numerical Recipe's LCG RNG */
class LCG_RNG {
  /* Long period (> 2 × 1018) random number generator of L’Ecuyer with
     Bays-Durham shuffle and added safeguards. Returns a uniform random
     deviate between 0.0 and 1.0 (exclusive of the endpoint values).
     Call with idum a negative integer to initialize; thereafter, do not
     alter idum between successive deviates in a sequence. RNMX should
     approximate the largest floating value that is less than 1.*/
#define NTAB 32

  int32_t const IM1;
  int32_t const IM2;
  int32_t const IMM1;
  int32_t const IA1;
  int32_t const IA2;
  int32_t const IQ1;
  int32_t const IQ2;
  int32_t const IR1;
  int32_t const IR2;
  float const AM;
  float const NDIV;
  float const EPS;
  float const RNMX;
  double const AM_d;

  int32_t m_idum;
  int32_t m_idum2;
  int32_t m_iy;
  std::vector<int32_t> m_iv;

public:
  void dump_state(std::ostream & outstr) const
  {
    std::copy(&m_iv[0], &m_iv[NTAB],
              std::ostream_iterator<int32_t>(outstr, ","));
  }

  int32_t const * get_state() { return &m_iv[0]; }

  size_t get_state_size() const { return NTAB; }

  // call with negative seed
  explicit LCG_RNG(int32_t const idum_init)
      : IM1(2147483563),
        IM2(2147483399),
        IMM1(IM1 - 1),
        IA1(40014),
        IA2(40692),
        IQ1(53668),
        IQ2(52774),
        IR1(12211),
        IR2(3791),
        // NTAB(32),
        AM(1.0 / IM1),
        NDIV(1 + IMM1 / NTAB),
        EPS(1.2e-7),
        RNMX(1.0 - EPS),
        AM_d(1.0 / IM1),
        m_idum(idum_init),
        m_idum2(123456789),
        m_iy(0),
        m_iv(NTAB, 0)

  {
    if(-(m_idum) < 1) {
      m_idum = 1;  // Be sure to prevent idum = 0.
    }
    else {
      m_idum = -(m_idum);
    }
    m_idum2 = (m_idum);
    for(size_t i = 0; i <= NTAB + 7;
        i++) {  // Load the shuffle table (after 8 warm-ups).
      size_t const j = NTAB + 7 - i;
      int32_t k = m_idum / IQ1;
      m_idum = IA1 * (m_idum - k * IQ1) - k * IR1;
      if(m_idum < 0) { m_idum += IM1; }
      if(j < NTAB) { m_iv[j] = m_idum; }
    }
    m_iy = m_iv[0];
    return;
  }  // ctor

  // copy ctor
  LCG_RNG(LCG_RNG const & r)
      : IM1(r.IM1),
        IM2(r.IM2),
        IMM1(r.IMM1),
        IA1(r.IA1),
        IA2(r.IA2),
        IQ1(r.IQ1),
        IQ2(r.IQ2),
        IR1(r.IR1),
        IR2(r.IR2),
        AM(r.AM),
        NDIV(r.NDIV),
        EPS(r.EPS),
        RNMX(r.RNMX),
        AM_d(r.AM_d),
        m_idum(r.m_idum),
        m_idum2(r.m_idum2),
        m_iy(r.m_iy),
        m_iv(r.m_iv)
  {
  }

  // hide assignment operator
  LCG_RNG & operator=(LCG_RNG const & l)
  {
    this->m_idum = l.m_idum;
    this->m_idum2 = l.m_idum2;
    this->m_iy = l.m_iy;
    this->m_iv = l.m_iv;
    return *this;
  }

  float ranf()
  {
    float const temp = AM * ranf_core();
    if(temp > RNMX)
      return RNMX;  // Because users don’t expect endpoint values.
    else
      return temp;
  }

  // use two floats to generate a double. Meh.
  double random()
  {
    // double const FLT_MULT = 2.3283064365386962891e-10;
    double f1 = AM_d * (double)this->ranf_core();
    double f2 = AM_d * (double)this->ranf_core();
    return f1 + AM_d * f2;
  }

private:
  // this is the core of the generator--it returns an unscaled number.
  float ranf_core()
  {
    int32_t j;
    int32_t k;

    k = m_idum / IQ1;  //  Start here when not initializing.
    m_idum = IA1 * (m_idum - k * IQ1) - k * IR1;
    if(m_idum < 0) {  // Compute idum=(IA1*idum) % IM1 without
      m_idum += IM1;  // overflows by Schrage’s method.
    }
    k = m_idum2 / IQ2;
    m_idum2 = IA2 * (m_idum2 - k * IQ2) -
              k * IR2;  // Compute idum2=(IA2*idum) % IM2 likewise.
    if(m_idum2 < 0) { m_idum2 += IM2; }
    j = m_iy / NDIV;  // Will be in the range 0..NTAB-1.
    if(j >= NTAB) {
      std::cerr << "index j out of range! j = " << j << ", NTAB = " << NTAB
                << std::endl;
    }
    m_iy = m_iv.at(j) - m_idum2;  // Here idum is shuffled, idum and idum2 are
    m_iv.at(j) = m_idum;          // combined to generate output.
    if(m_iy < 1) { m_iy += IMM1; }
    return m_iy;
  }

#undef NTAB  // Try to limit the scope of this monstrosity

};  // LCG

}  // namespace nut

// End of file
