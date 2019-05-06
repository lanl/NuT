// Luminosity.hh
// T. M. Kelley
// Jun 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#include <vector>

namespace nut {
template <typename fp_t>
struct Luminosity {
  typedef std::vector<fp_t> vf;

  vf nue;
  vf nueb;
  vf nux;

  explicit Luminosity(size_t n)
      : nue(n, fp_t(0)), nueb(n, fp_t(0)), nux(n, fp_t(0))
  {
  }

  size_t size() const { return nue.size(); }

};  // Luminosity

template <typename fp_t>
struct lum_t {
  fp_t nue;
  fp_t nueb;
  fp_t nux;
};  // lum

}  // namespace nut

// version
// $Id$

// End of file
