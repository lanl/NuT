// utilities_io.hh
// T. M. Kelley
// Oct 07, 2015
// (c) Copyright 2015 LANSLLC, all rights reserved

#ifndef UTILITIES_IO_HH
#define UTILITIES_IO_HH

#include "Assert.hh"
#include "types.hh"
#include <ostream>
#include <string>

/**!\file commonly used assertions, checks, exceptions, etc. */

namespace nut {
/** \brief Print arrays in Python style. */
template <typename T, typename ItT>
void
print_per_cell_field(std::string const & name,
                     ItT begin,
                     ItT end,
                     std::ostream & o)
{
  o << name << " = [";
  std::copy(begin, end, std::ostream_iterator<T>(o, ","));
  o << "]\n";
}

template <typename esc, typename vescIt>
void
print_esc_spectrum(std::string const & name,
                   vescIt const begin,
                   vescIt const end,
                   std::ostream & o)
{
  o << name << " = [";
  for(vescIt curr = begin; curr != end; ++curr) {
    o << "(" << (*curr).first << "," << (*curr).second << "),";
  }
  o << "]\n";
}

/**! Assert that the cell c > 0 and c <= n_cells */
inline void
cellOK(cell_t const c, cell_t const n_cells)
{
  InOpenRange(c, cell_t(0), n_cells + 1, "cell");
}

/**! Assert that cell c is ok (above) and return index into arrays */
inline cell_t
make_idx(cell_t const c, cell_t const n_cells)
{
  cellOK(c, n_cells);
  return c - 1;
}

// throw an exception when an unknown species is encountered
inline void
unhandled_species_exc(Species const s, std::string const & method)
{
  std::stringstream errstr;
  errstr << method << ": Unhandled neutrino species " << s
         << ", better known as " << species_name(s);
  throw arg_error(errstr.str());
}

template <typename fp_t>
inline void
nrgOK(fp_t const e, char const * const s)
{
  GreaterThan(e, fp_t(0), s);
}

template <typename fp_t>
inline void
nrgOK(fp_t const e)
{
  char const * const s = "neutrino energy";
  GreaterThan(e, fp_t(0), s);
}

}  // namespace nut

#endif  // include guard

// End of file
