// Opacity.hh
// T. M. Kelley
// Jan 12, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifdef USING_HFBC_SIGMAS
#include "sigmas_HBFC.hh"
#else
#error "Undefined cross sections!!!"
#endif

#include "Cell_Data.hh"
#include "Planck.hh"
#include "types.hh"
#include "utilities_io.hh"

namespace nut {

template <typename fp_t, typename vector_t>
class Opacity {
public:
  // types
  using Cell_Data_T = Cell_Data<fp_t, vector_t>;
  using fp_cr = fp_t const &;
  using vec_cell_data_t = std::vector<Cell_Data_T>;

  /** Get pair of begin and end iterators of nu_e luminance data.
   *
   * These iterators march over the array of structs picking out just the
   * nu_e luminance field. It a simple thing, and difficult to say that type.
   */
  auto get_nue_luminance_data() const {
    return make_pair(make_nu_e_iterator<Cell_Data_T>(m_cell_data.begin()),
                     make_nu_e_iterator<Cell_Data_T>(m_cell_data.end()));
  }

  /** Get pair of begin and end iterators of nu_ebar luminance data.
   *
   * These iterators march over the array of structs picking out just the
   * nu_ebar luminance field. It a simple thing, and difficult to say that type.
   */
  auto get_nuebar_luminance_data() const {
    return make_pair(make_nu_ebar_iterator<Cell_Data_T>(m_cell_data.begin()),
                     make_nu_ebar_iterator<Cell_Data_T>(m_cell_data.end()));
  }

  /** Get pair of begin and end iterators of nu_x luminance data.
   *
   * These iterators march over the array of structs picking out just the
   * nu_x luminance field. It a simple thing, and difficult to say that type.
   */
  auto get_nux_luminance_data() const {
    return make_pair(make_nu_x_iterator<Cell_Data_T>(m_cell_data.begin()),
                     make_nu_x_iterator<Cell_Data_T>(m_cell_data.end()));
  }

  // ctor
  explicit Opacity(vec_cell_data_t && data, bool have_composition_ = false)
      : m_cell_data(std::forward<vec_cell_data_t>(data)),
        m_n_cells(m_cell_data.size()),
        m_have_composition(have_composition_)
  {
  }

  size_t const num_cells() const { return m_n_cells; }

  fp_cr T_p(cell_t const cidx) const
  {
    cellOK(cidx, m_n_cells + 1);
    return m_cell_data.at(cidx - 1).T_p;
  }

  Cell_Data_T const &cell_data(cell_t const cidx) const
  {
    cellOK(cidx, m_n_cells + 1);
    return m_cell_data.at(cidx - 1);
  }

  vector_t const & velocity(cell_t const cidx) const
  {
    cellOK(cidx, m_n_cells + 1);
    return m_cell_data.at(cidx - 1).velocity;
  }

  /* total collision opacity */
  fp_t sigma_collide(cell_t const & cell,
                     fp_t const e_nu,
                     Species const nu_species) const
  {
    cellOK(cell, m_n_cells);
    nrgOK(e_nu);
    fp_t s_coll =
        sigma_N_total(cell, e_nu) + sigma_lepton(cell, e_nu, nu_species);
    return s_coll;
  }  // sigma_collide

  /* total nucleon opacity */
  fp_t sigma_N_total(cell_t const cell, fp_t const e_nu) const
  {
    cellOK(cell, m_n_cells);
    nrgOK(e_nu);
    fp_t const sig_N_abs = sigma_N_abs(cell, e_nu);
    fp_t const sig_N_el = sigma_N_elastic(cell, e_nu);
    return sig_N_abs + sig_N_el;
  }

  /** total lepton opacity */
  fp_t sigma_lepton(cell_t const cell,
                    fp_t const e_nu,
                    Species const nu_species) const
  {
    // compute total interaction cross-section for lepton interaction
    nrgOK(e_nu);
    fp_t const sig_e_m = sigma_nu_e_minus(cell, e_nu, nu_species);
    fp_t const sig_e_p = sigma_nu_e_plus(cell, e_nu, nu_species);
    fp_t const sigma = sig_e_m + sig_e_p;
    return sigma;
  }  // sigma_lepton

  /**nucleon absorption */
  fp_t sigma_N_abs(cell_t const cell, fp_t const e_nu) const
  {
    nrgOK(e_nu);
    fp_t const sig_N_abs =
        m_cell_data.at(cell - 1).rho_nuc() * sigmas::nu_N_abs(e_nu);
    return sig_N_abs;
  }  // sigma_N_abs

  /** nucleon elastic scatter */
  fp_t sigma_N_elastic(cell_t const cell, fp_t const e_nu) const
  {
    nrgOK(e_nu);
    size_t const index = make_idx(cell, m_n_cells);
    // we treat protons and neutrons distinctly from the nuclei
    // for elastic scattering *if* we have composition information.
    fp_t const sig_N_el =
        sigmas::nu_N_elastic(e_nu) * (m_cell_data[index].rho_p);
    fp_t sigma(sig_N_el);
    if(m_have_composition == true) {
      fp_t const abar = m_cell_data[index].abar;
      sigma += m_cell_data[index].rho_A * sigmas::nu_A_elastic(e_nu, abar);
    }
    return sigma;
  }  // sigma_N_elastic
  // nu-lepton cross sections, using mean temp as lepton energy

  fp_t sigma_nu_e_minus(cell_t const cell,
                        fp_cr e_nu,
                        Species const nu_species) const
  {
    // compute opacity for e- interaction
    cellOK(cell, m_n_cells);
    nrgOK(e_nu);
    fp_t sigma(0);

    cell_t const index = make_idx(cell, m_n_cells);
    fp_t const e_e_minus = m_cell_data[index].T_e_minus;
    fp_t const rho = m_cell_data[index].rho_e_minus;

    switch(nu_species) {
      case nu_e: sigma = sigma_nu_e_e_minus(e_nu, e_e_minus, rho); break;
      case nu_e_bar:
        sigma = sigma_nu_e_bar_e_minus(e_nu, e_e_minus, rho);
        break;
      case nu_x:
      case nu_mu:
      case nu_tau: sigma = sigma_nu_x_e_minus(e_nu, e_e_minus, rho); break;
      case nu_mu_bar:
      case nu_tau_bar:
        sigma = sigma_nu_x_bar_e_minus(e_nu, e_e_minus, rho);
        break;
      default:
        // should not get here
        unhandled_species_exc(nu_species, "Opacity::sigma_nu_e_minus");
    }  // switch
    return sigma;
  }

  fp_t sigma_nu_e_plus(cell_t const cell,
                       fp_cr e_nu,
                       Species const nu_species) const
  {
    cellOK(cell, m_n_cells);
    nrgOK(e_nu);
    fp_t sigma(0);
    cell_t const index = make_idx(cell, m_n_cells);

    fp_t const e_e_plus = m_cell_data[index].T_e_plus;
    fp_t const rho = m_cell_data[index].rho_e_plus;

    switch(nu_species) {
      case nu_e: sigma = sigma_nu_e_e_plus(e_nu, e_e_plus, rho); break;
      case nu_e_bar: sigma = sigma_nu_e_bar_e_plus(e_nu, e_e_plus, rho); break;
      case nu_x:
      case nu_mu:
      case nu_tau: sigma = sigma_nu_x_e_plus(e_nu, e_e_plus, rho); break;
      case nu_mu_bar:
      case nu_tau_bar:
        sigma = sigma_nu_x_bar_e_plus(e_nu, e_e_plus, rho);
        break;
      default:
        // should not get here
        unhandled_species_exc(nu_species, "Opacity::sigma_nu_e_plus");
    }
    return sigma;
  }  // sigma_nu_e_plus

  fp_t sigma_nu_e_e_minus(fp_cr e_nu, fp_cr e_e_minus, fp_cr rho) const
  {
    fp_t sigma = sigmas::nu_e_e_minus(e_nu, e_e_minus) * rho;
    return sigma;
  }

  fp_t sigma_nu_e_bar_e_minus(fp_cr e_nu, fp_cr e_e_minus, fp_cr rho) const
  {
    fp_t sigma = sigmas::nu_e_bar_e_minus(e_nu, e_e_minus) * rho;
    return sigma;
  }

  fp_t sigma_nu_x_e_minus(fp_cr e_nu, fp_cr e_e_minus, fp_cr rho) const
  {
    fp_t sigma = sigmas::nu_x_e_minus(e_nu, e_e_minus) * rho;
    return sigma;
  }

  fp_t sigma_nu_x_bar_e_minus(fp_cr e_nu, fp_cr e_e_minus, fp_cr rho) const
  {
    fp_t sigma = sigmas::nu_x_bar_e_minus(e_nu, e_e_minus) * rho;
    return sigma;
  }

  fp_t sigma_nu_e_e_plus(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho) const
  {
    fp_t sigma = sigmas::nu_e_e_plus(e_nu, e_e_plus) * rho;
    return sigma;
  }

  fp_t sigma_nu_e_bar_e_plus(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho) const
  {
    fp_t sigma = sigmas::nu_e_bar_e_plus(e_nu, e_e_plus) * rho;
    return sigma;
  }

  fp_t sigma_nu_x_e_plus(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho) const
  {
    fp_t sigma = sigmas::nu_x_e_plus(e_nu, e_e_plus) * rho;
    return sigma;
  }

  fp_t sigma_nu_x_bar_e_plus(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho) const
  {
    fp_t sigma = sigmas::nu_x_bar_e_plus(e_nu, e_e_plus) * rho;
    return sigma;
  }

  // state:
  vec_cell_data_t m_cell_data;

  cell_t const m_n_cells;

  bool m_have_composition;
};  // class Opacity

}  // namespace nut

// End of file
