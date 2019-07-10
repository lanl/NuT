// Tally.hh
// T. M. Kelley
// Jan. 10, 2010
// Header for Tally
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef TALLY_H
#define TALLY_H

#include "Assert.hh"
#include "detail/Vector.h"
#include "types.hh"
#include "utilities.hh"
#include "utilities_io.hh"
#include <algorithm>  // transform
#include <numeric>
#include <sstream>
#include <vector>

namespace nut {

/*! \brief arrays to keep track of tally in each mesh cell
 * \tparam <fp_t> {floating point type}
 */
template <
    typename fpt,
    size_t Dim,
    /* Selects Vector if Dim == 3, otherwise Vector1 */
    typename vector_ty = std::conditional_t<
        Dim == static_cast<size_t>(3u),
        nut::Vector,
        std::conditional_t<Dim == static_cast<size_t>(1u), nut::Vector1, void>>>
struct Tally {
  static_assert(Dim == 1u || Dim == 3u,
                "Tally implemented only for one or three dimensional meshes");

  static constexpr size_t dim = Dim;

  using cntr_t = uint32_t;  // counter type
  using fp_t = fpt;
  using FP_T = fp_t;
  using vector_t = vector_ty;
  using vec_cr = vector_t const &;
  using vc = std::vector<cntr_t>;
  using vf = std::vector<fp_t>;
  using vv = std::vector<vector_t>;
  using fp_cr = fp_t const &;
  using esc = std::pair<fp_t, fp_t>;
  using vesc = std::vector<esc>;

  vf energy;
  vv momentum;

  // particle counters: number & energy weight of each
  // species created in each cell (can be negative)
  vc n_n;
  vc n_p;
  vc n_e_minus;
  vc n_e_plus;

  vf ew_n;
  vf ew_p;
  vf ew_e_minus;
  vf ew_e_plus;

  // cell-wise event counters
  vc n_escape;
  vc n_reflect;
  vc n_cell_bdy;
  vc n_cutoff;
  vc n_nucl_el_scat;  // elastic scatter from nucleon

  vc n_nu_e_el_scat;       // nu_e -- electron scatter
  vc n_nu_e_bar_pos_scat;  // nu_e_bar -- positron scatter
  vc n_nu_x_el_scat;
  vc n_nu_x_bar_pos_scat;

  vf ew_escaped;

  // for nucleon absorption, it's essential to track
  // momentum and energy dep by species
  vc n_nu_e_nucl_abs;
  vc n_nu_e_bar_nucl_abs;
  vc n_nu_x_nucl_abs;
  vf ew_nu_e_nucl_abs;
  vf ew_nu_e_bar_nucl_abs;
  vf ew_nu_x_nucl_abs;

  // census counters
  vc n_census_nu_e;
  vc n_census_nu_e_bar;
  vc n_census_nu_x;
  vc n_census_nu_x_bar;

  vf ew_census_nu_e;
  vf ew_census_nu_e_bar;
  vf ew_census_nu_x;
  vf ew_census_nu_x_bar;

  vesc escape_spectrum;

  fp_t path_length;

  // interface

  /*!\brief structure of arrays to store tally data.
   */
  explicit Tally(size_t const ncells = 0)
      : energy(ncells, fp_t(0)),
        momentum(ncells, {fp_t(0)}),
        n_n(ncells, cntr_t(0)),
        n_p(ncells, cntr_t(0)),
        n_e_minus(ncells, cntr_t(0)),
        n_e_plus(ncells, cntr_t(0)),
        ew_n(ncells, fp_t(0)),
        ew_p(ncells, fp_t(0)),
        ew_e_minus(ncells, fp_t(0)),
        ew_e_plus(ncells, fp_t(0)),
        n_escape(ncells, cntr_t(0)),
        n_reflect(ncells, cntr_t(0)),
        n_cell_bdy(ncells, cntr_t(0)),
        n_cutoff(ncells, cntr_t(0)),
        n_nucl_el_scat(ncells, cntr_t(0)),
        n_nu_e_el_scat(ncells, cntr_t(0)),
        n_nu_e_bar_pos_scat(ncells, cntr_t(0)),
        n_nu_x_el_scat(ncells, cntr_t(0)),
        n_nu_x_bar_pos_scat(ncells, cntr_t(0)),
        ew_escaped(ncells, fp_t(0)),
        n_nu_e_nucl_abs(ncells, cntr_t(0)),
        n_nu_e_bar_nucl_abs(ncells, cntr_t(0)),
        n_nu_x_nucl_abs(ncells, cntr_t(0)),
        ew_nu_e_nucl_abs(ncells, fp_t(0)),
        ew_nu_e_bar_nucl_abs(ncells, fp_t(0)),
        ew_nu_x_nucl_abs(ncells, fp_t(0)),
        n_census_nu_e(ncells, cntr_t(0)),
        n_census_nu_e_bar(ncells, cntr_t(0)),
        n_census_nu_x(ncells, cntr_t(0)),
        n_census_nu_x_bar(ncells, cntr_t(0)),
        ew_census_nu_e(ncells, fp_t(0)),
        ew_census_nu_e_bar(ncells, fp_t(0)),
        ew_census_nu_x(ncells, fp_t(0)),
        ew_census_nu_x_bar(ncells, fp_t(0)),
        path_length(fp_t(0)),
        m_n_cells(ncells),
        m_initialized(0 == ncells ? false : true)
  {
  }

  /** \brief Resize tally: only callable once. */
  void resize(size_t const ncells)
  {
    Require(!m_initialized, "Tally already initialized");
    m_initialized = true;

    energy.resize(ncells);
    momentum.resize(ncells);
    n_n.resize(ncells);
    n_p.resize(ncells);
    n_e_minus.resize(ncells);

    n_e_plus.resize(ncells);
    ew_n.resize(ncells);
    ew_p.resize(ncells);
    ew_e_minus.resize(ncells);
    ew_e_plus.resize(ncells);

    n_escape.resize(ncells);
    n_reflect.resize(ncells);
    n_cell_bdy.resize(ncells);
    n_cutoff.resize(ncells);
    n_nucl_el_scat.resize(ncells);

    n_nu_e_el_scat.resize(ncells);
    n_nu_e_bar_pos_scat.resize(ncells);
    n_nu_x_el_scat.resize(ncells);
    n_nu_x_bar_pos_scat.resize(ncells);
    ew_escaped.resize(ncells);

    n_nu_e_nucl_abs.resize(ncells);
    n_nu_e_bar_nucl_abs.resize(ncells);
    n_nu_x_nucl_abs.resize(ncells);
    ew_nu_e_nucl_abs.resize(ncells);
    ew_nu_e_bar_nucl_abs.resize(ncells);

    ew_nu_x_nucl_abs.resize(ncells);
    n_census_nu_e.resize(ncells);
    n_census_nu_e_bar.resize(ncells);
    n_census_nu_x.resize(ncells);
    n_census_nu_x_bar.resize(ncells);

    ew_census_nu_e.resize(ncells);
    ew_census_nu_e_bar.resize(ncells);
    ew_census_nu_x.resize(ncells);
    ew_census_nu_x_bar.resize(ncells);

    m_n_cells = ncells;
    return;
  }

  void merge(Tally<fp_t, dim> const & other)
  {
    Require(other.n_cells() == this->n_cells(),
            "Cannot merge tallies with different sizes");

    merge_vectors(other.energy, energy);
    merge_vectors(other.momentum, momentum);
    merge_vectors(other.n_n, n_n);
    merge_vectors(other.n_p, n_p);
    merge_vectors(other.n_e_minus, n_e_minus);
    merge_vectors(other.n_e_plus, n_e_plus);
    merge_vectors(other.ew_n, ew_n);
    merge_vectors(other.ew_p, ew_p);
    merge_vectors(other.ew_e_minus, ew_e_minus);
    merge_vectors(other.ew_e_plus, ew_e_plus);
    merge_vectors(other.n_escape, n_escape);
    merge_vectors(other.n_reflect, n_reflect);
    merge_vectors(other.n_cell_bdy, n_cell_bdy);
    merge_vectors(other.n_cutoff, n_cutoff);
    merge_vectors(other.n_nucl_el_scat, n_nucl_el_scat);
    merge_vectors(other.n_nu_e_el_scat, n_nu_e_el_scat);
    merge_vectors(other.n_nu_e_bar_pos_scat, n_nu_e_bar_pos_scat);
    merge_vectors(other.n_nu_x_el_scat, n_nu_x_el_scat);
    merge_vectors(other.n_nu_x_bar_pos_scat, n_nu_x_bar_pos_scat);
    merge_vectors(other.ew_escaped, ew_escaped);
    merge_vectors(other.n_nu_e_nucl_abs, n_nu_e_nucl_abs);
    merge_vectors(other.n_nu_e_bar_nucl_abs, n_nu_e_bar_nucl_abs);
    merge_vectors(other.n_nu_x_nucl_abs, n_nu_x_nucl_abs);
    merge_vectors(other.ew_nu_e_nucl_abs, ew_nu_e_nucl_abs);
    merge_vectors(other.ew_nu_e_bar_nucl_abs, ew_nu_e_bar_nucl_abs);
    merge_vectors(other.ew_nu_x_nucl_abs, ew_nu_x_nucl_abs);
    merge_vectors(other.n_census_nu_e, n_census_nu_e);
    merge_vectors(other.n_census_nu_e_bar, n_census_nu_e_bar);
    merge_vectors(other.n_census_nu_x, n_census_nu_x);
    merge_vectors(other.n_census_nu_x_bar, n_census_nu_x_bar);
    merge_vectors(other.ew_census_nu_e, ew_census_nu_e);
    merge_vectors(other.ew_census_nu_e_bar, ew_census_nu_e_bar);
    merge_vectors(other.ew_census_nu_x, ew_census_nu_x);
    merge_vectors(other.ew_census_nu_x_bar, ew_census_nu_x_bar);

    path_length += other.path_length;
    return;
  }  // merge

  uint32_t n_cells() const { return m_n_cells; }

  void accum_pl(fp_t const pl)
  {
    Require(pl >= fp_t(0), "path length < 0!");
    path_length += pl;
  }

  uint32_t total_mc_steps() const
  {
    return std::accumulate(n_escape.begin(), n_escape.end(), 0) +
           std::accumulate(n_reflect.begin(), n_reflect.end(), 0) +
           std::accumulate(n_cell_bdy.begin(), n_cell_bdy.end(), 0) +
           std::accumulate(n_cutoff.begin(), n_cutoff.end(), 0) +
           std::accumulate(n_nucl_el_scat.begin(), n_nucl_el_scat.end(), 0) +
           std::accumulate(n_nu_e_el_scat.begin(), n_nu_e_el_scat.end(), 0) +
           std::accumulate(n_nu_e_bar_pos_scat.begin(),
                           n_nu_e_bar_pos_scat.end(), 0) +
           std::accumulate(n_nu_x_el_scat.begin(), n_nu_x_el_scat.end(), 0) +
           std::accumulate(n_nu_x_bar_pos_scat.begin(),
                           n_nu_x_bar_pos_scat.end(), 0) +
           std::accumulate(n_nu_e_nucl_abs.begin(), n_nu_e_nucl_abs.end(), 0) +
           std::accumulate(n_nu_e_bar_nucl_abs.begin(),
                           n_nu_e_bar_nucl_abs.end(), 0) +
           std::accumulate(n_nu_x_nucl_abs.begin(), n_nu_x_nucl_abs.end(), 0) +
           std::accumulate(n_census_nu_e.begin(), n_census_nu_e.end(), 0) +
           std::accumulate(n_census_nu_e_bar.begin(), n_census_nu_e_bar.end(),
                           0) +
           std::accumulate(n_census_nu_x.begin(), n_census_nu_x.end(), 0) +
           std::accumulate(n_census_nu_x_bar.begin(), n_census_nu_x_bar.end(),
                           0);
  }  // total_mc_steps

  uint32_t total_n_census() const
  {
    return std::accumulate(n_census_nu_e.begin(), n_census_nu_e.end(), 0) +
           std::accumulate(n_census_nu_e_bar.begin(), n_census_nu_e_bar.end(),
                           0) +
           std::accumulate(n_census_nu_x.begin(), n_census_nu_x.end(), 0) +
           std::accumulate(n_census_nu_x_bar.begin(), n_census_nu_x_bar.end(),
                           0);
  }  // total_mc_steps

  // template <typename vector_t>
  void deposit_inelastic_scat(cell_t const cell,
                              fp_cr e_i,
                              fp_cr e_f,
                              vector_t const & omega_i,
                              vector_t const & omega_f,
                              fp_cr wt,
                              nut::Species const species)
  {
    // We tally omega * ew, need to divide by c at end of time step
    // Do we need to track momentum and energy seperately by species?
    static_assert(vector_t::dim == Dim,
                  "template vector dimenstion must equal tally dimension");
    nrgOK(e_i, "initial energy");
    nrgOK(e_f, "final energy");
    cell_t index = make_idx(cell, m_n_cells);

    vector_t mom_dep = (e_i * omega_i) - (e_f * omega_f);

    for(size_t i = 0; i < vector_t::dim; ++i) {
      momentum[index][i] += wt * mom_dep[i];
    }
    energy[index] += wt * (e_i - e_f);
    return;
  }  // deposit_inelastic_el_scat

  void deposit_energy(cell_t const cell, fp_cr wt, fp_cr e)
  {
    // particle could take energy from medium, thus no nrgOK
    cell_t const index = make_idx(cell, m_n_cells);
    energy[index] += e * wt;
    return;
  }  // deposit_energy

  // template <typename vector_t>
  void deposit_momentum_elastic(cell_t const cell,
                                vector_t const & omega,
                                fp_cr e,
                                fp_cr wt)
  {
    // We tally omega * ew, need to divide by c at end of time step
    Equal(vector_t::dim, Dim, "template vector dimension", "tally dimension");
    nrgOK(e, "energy");
    cell_t const index = make_idx(cell, m_n_cells);
    vector_t momm = wt * omega * e;
    for(size_t i = 0; i < dim; ++i) { momentum[index][i] += momm[i]; }
    return;
  }  // deposit_momentum_elastic

  void count_lepton_scatter(cell_t const cell,
                            Species const species,
                            size_t const n = 1)
  {
    // GreaterThan(n,0,"# nucleon absorptions");
    cell_t const index = make_idx(cell, m_n_cells);
    switch(species) {
      case nu_e: this->n_nu_e_el_scat[index] += 1; break;
      case nu_e_bar: this->n_nu_e_bar_pos_scat[index] += 1; break;
      case nu_x:
      case nu_mu:
      case nu_tau:  // yes--fall through
        this->n_nu_x_el_scat[index] += 1;
        break;
      case nu_mu_bar:
      case nu_tau_bar:  // yes--fall through
        this->n_nu_x_bar_pos_scat[index] += 1;
        break;
      default: unhandled_species_exc(species, "count_electron_scatter");
    }
    return;
  }

  void count_nucleon_abs(cell_t const cell,
                         Species const species,
                         fp_cr ew,
                         cntr_t n = 1)
  {
    GreaterThan(n, 0u, "# nucleon absorptions");
    cell_t const index = make_idx(cell, m_n_cells);
    // tally ew absorbed by species
    if(species == nu_e) {
      n_nu_e_nucl_abs[index] += n;
      ew_nu_e_nucl_abs[index] += ew;
    }
    else if(species == nu_e_bar) {
      n_nu_e_bar_nucl_abs[index] += n;
      ew_nu_e_bar_nucl_abs[index] += ew;
    }
    else if(species == nu_x) {
      n_nu_x_nucl_abs[index] += n;
      ew_nu_x_nucl_abs[index] += ew;
    }
    else {
      std::stringstream errstr;
      errstr << "Tally::count_nucleon_abs: unknown/unhandled species "
             << species << ".";
      throw(std::runtime_error(errstr.str()));
    }
    return;
  }

  void count_nucleon_elastic_scatter(cell_t const cell, cntr_t const n = 1)
  {
    GreaterThan(n, 0u, "# nucleon elastic scatter");
    cell_t const index = make_idx(cell, m_n_cells);
    n_nucl_el_scat[index] += n;
    return;
  }

  void count_escape(cell_t const cell, fp_cr ew, fp_cr e, cntr_t const n = 1)
  {
    GreaterThan(n, 0u, "# escapes");
    cell_t const index = make_idx(cell, m_n_cells);
    n_escape[index] += n;
    ew_escaped[index] += ew;
    escape_spectrum.push_back(esc(e, ew));
    return;
  }

  void count_reflect(cell_t const cell, cntr_t const n = 1)
  {
    GreaterThan(n, 0u, "# reflects");
    cell_t const index = make_idx(cell, m_n_cells);
    n_reflect[index] += n;
    return;
  }

  void count_cell_bdy(cell_t const cell, cntr_t const n = 1)
  {
    GreaterThan(n, 0u, "# cell bdys");
    cell_t const index = make_idx(cell, m_n_cells);
    n_cell_bdy[index] += n;
    return;
  }

  void count_cutoff(cell_t const cell, cntr_t const n = 1)
  {
    GreaterThan(n, 0u, "# cutoffs");
    cell_t const index = make_idx(cell, m_n_cells);
    n_cutoff[index] += n;
    return;
  }

  void count_census(cell_t const cell,
                    fp_cr energy_weight,
                    Species const & species,
                    cntr_t const n = 1)
  {
    GreaterThan(energy_weight, fp_t(0.0), "energy weight");
    cell_t const index = make_idx(cell, m_n_cells);
    switch(species) {
      case nu_e:
        n_census_nu_e[index] += n;
        ew_census_nu_e[index] += energy_weight;
        break;
      case nu_e_bar:
        n_census_nu_e_bar[index] += n;
        ew_census_nu_e_bar[index] += energy_weight;
        break;
      case nu_mu:
      case nu_tau:
        n_census_nu_x[index] += n;
        ew_census_nu_x[index] += energy_weight;
        break;
      case nu_mu_bar:
      case nu_tau_bar:
        n_census_nu_x_bar[index] += n;
        ew_census_nu_x_bar[index] += energy_weight;
        break;
      default: unhandled_species_exc(species, "Tally::count_census");
    }
    return;
  }

private:
  cell_t m_n_cells;

  bool m_initialized;

};  // Tally

}  // namespace nut

#endif

// handy list of all vectors for code generation:
// ["energy","momentum","n_n","n_p","n_e_minus","n_e_plus","ew_n","ew_p","ew_e_minus","ew_e_plus","n_escape","n_reflect","n_cell_bdy","n_cutoff","n_nucl_el_scat","n_nu_e_el_scat","n_nu_e_bar_pos_scat","n_nu_x_el_scat","n_nu_x_bar_pos_scat","ew_escaped","n_nu_e_nucl_abs","n_nu_e_bar_nucl_abs","n_nu_x_nucl_abs","ew_nu_e_nucl_abs","ew_nu_e_bar_nucl_abs","ew_nu_x_nucl_abs","n_census_nu_e","n_census_nu_e_bar","n_census_nu_x","n_census_nu_x_bar","ew_census_nu_e","ew_census_nu_e_bar","ew_census_nu_x","ew_census_nu_x_bar"]

// version
// $Id$

// End of file
