// Cell_Data.hh
// T. M. Kelley
// Mar 04, 2020
// (c) Copyright 2020 Triad National Security, all rights reserved

#pragma once

#include <iterator>
#include <vector>

namespace nut {

template <typename flp_t, typename vector_t>
struct Cell_Data {
  // types
  using fp_t = flp_t;

  // Densities
  fp_t rho_p;        //!< density of protons
  fp_t rho_e_minus;  //!< density of electrons
  fp_t rho_e_plus;   //!< density of positrons
  fp_t rho_A;        //!< density of nucleona
  fp_t y_e;          //!< electron number density fraction n_e/(n_p + n_n)
  fp_t abar;         //!< Average nucleon mass

  // Temperatures
  fp_t T_p;        //!< temperature of protons
  fp_t T_e_minus;  //!< temperature of electrons
  fp_t T_e_plus;   //!< temperature of positrons

  // Luminosities
  fp_t l_nue;   //!< luminosity of nu_e
  fp_t l_nueb;  //!< luminosity of (nu_e)bar
  fp_t l_nux;   //!< luminosity of nu_x

  // velocity
  vector_t velocity;

  fp_t rho_nuc() const { return rho_p + rho_A; }

};  // Cell_Data

/**\brief Get an iterator over nu_e luminosity data given an iterator to a
 * vector of cell_data_t
 *
 * \param it: iterator on a vector<cell_data_t>
 * \return a Cell_Data_Iterator that selects the nu_e luminosity field
 */
template <typename cell_data_t>
auto
make_nu_e_iterator(typename std::vector<cell_data_t>::const_iterator && it);

/**\brief Get an iterator over nu_ebar luminosity data given an iterator to a
 * vector of cell_data_t
 *
 * \param it: iterator on a vector<cell_data_t>
 * \return a Cell_Data_Iterator that selects the nu_ebar luminosity field
 */
template <typename cell_data_t>
auto
make_nu_ebar_iterator(typename std::vector<cell_data_t>::const_iterator && it);

/**\brief Get an iterator over nu_x luminosity data given an iterator to a
 * vector of cell_data_t
 *
 * \param it: iterator on a vector<cell_data_t>
 * \return a Cell_Data_Iterator that selects the nu_x luminosity field
 */
template <typename cell_data_t>
auto
make_nu_x_iterator(typename std::vector<cell_data_t>::const_iterator && it);

/**\brief Class to iterate over a single field of a Cell_Data struct.
 *
 * It's a bit ugly, but this should be simplfied by using / extending
 * the make_<field>_iterator functions above. */
template <typename Cell_Data_T,
          typename Cell_Data_Projector,
          typename Return_T = decltype(
              std::declval<Cell_Data_Projector>()(std::declval<Cell_Data_T>()))>
struct Cell_Data_Iterator
    : public std::iterator<std::input_iterator_tag, Return_T> {
  // types
  using value_type = Return_T;
  using vector_t = std::vector<Cell_Data_T>;
  using vector_it_t = typename vector_t::const_iterator;

  Cell_Data_Iterator(vector_it_t & it, Cell_Data_Projector p)
      : m_current(it), m_proj(p)
  {
  }

  Cell_Data_Iterator(vector_it_t && it, Cell_Data_Projector p)
      : m_current(it), m_proj(p)
  {
  }

  bool operator==(Cell_Data_Iterator const & other) const
  {
    return m_current == other.m_current;
  }

  bool operator!=(Cell_Data_Iterator const & other) const
  {
    return !(*this == other);
  }

  Cell_Data_Iterator & operator++()
  {
    m_current++;
    return *this;
  }

  Cell_Data_Iterator operator++(int)
  {
    Cell_Data_Iterator retval = *this;
    ++(*this);
    return retval;
  }

  value_type operator*() const { return m_proj(*m_current); }

  vector_it_t m_current;

  Cell_Data_Projector m_proj;
};

template <typename cell_data_t>
auto
make_nu_e_iterator(typename std::vector<cell_data_t>::const_iterator && it)
{
  auto fetch_nue = [](cell_data_t const & d) { return d.l_nue; };
  return Cell_Data_Iterator<cell_data_t, decltype(fetch_nue)>(it, fetch_nue);
}

template <typename cell_data_t>
auto
make_nu_ebar_iterator(typename std::vector<cell_data_t>::const_iterator && it)
{
  auto fetch_nueb = [](cell_data_t const & d) { return d.l_nueb; };
  return Cell_Data_Iterator<cell_data_t, decltype(fetch_nueb)>(it, fetch_nueb);
}

template <typename cell_data_t>
auto
make_nu_x_iterator(typename std::vector<cell_data_t>::const_iterator && it)
{
  auto fetch_nux = [](cell_data_t const & d) { return d.l_nux; };
  return Cell_Data_Iterator<cell_data_t, decltype(fetch_nux)>(it, fetch_nux);
}

}  // namespace nut

// End of file
