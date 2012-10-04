// Opacity.hh
// T. M. Kelley
// Jan 12, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifdef USING_HFBC_SIGMAS
#include "sigmas_HBFC.hh"
#else
#error "Undefined cross sections!!!"
#endif

#include "Density.hh"
#include "Temperature.hh"
#include "Planck.hh"
#include "types.hh"
#include "utilities.hh"


namespace nut{

    template <typename fp_t,
              // typename rng_t,
              typename den_t = Density<fp_t>,
              typename temp_t = Temperature<fp_t> >
    class Opacity
    {
    public:
        // types
        typedef den_t density_t;
        typedef temp_t temperature_t;
        typedef fp_t const & fp_cr;

        // ctor
        Opacity(density_t const & rho_, temperature_t const & T_)
            : m_n_cells(rho_.n_cells),m_rho(rho_),m_T(T_) {}

        fp_cr temp(cell_t const cidx) const
        {
            cellOK(cidx,m_n_cells+1);
            return m_T.T_p.at(cidx-1);
        }

        /* total collision opacity */
        fp_t sigma_collide( cell_t const & cell,
                            fp_t const e_nu,
                            Species const nu_species) const {
            cellOK(cell,m_n_cells);
            nrgOK(e_nu);
            fp_t s_coll =
                sigma_N_total(cell, e_nu) +
                sigma_lepton(cell,e_nu, nu_species);
            return s_coll;
        } // sigma_collide

        /* total nucleon opacity */
        fp_t sigma_N_total(cell_t const cell, fp_t const e_nu) const {
            cellOK(cell,m_n_cells);
            nrgOK(e_nu);
            fp_t const sig_N_abs = sigma_N_abs(cell,e_nu);
            fp_t const sig_N_el  = sigma_N_elastic(cell,e_nu);
            return  sig_N_abs + sig_N_el;
        }


        /** total lepton opacity */
        fp_t sigma_lepton(cell_t const cell, fp_t const e_nu,
                          Species const nu_species) const {
            // compute total interaction cross-section for lepton interaction
            nrgOK(e_nu);
            fp_t const sig_e_m = sigma_nu_e_minus(cell, e_nu, nu_species);
            fp_t const sig_e_p = sigma_nu_e_plus( cell, e_nu, nu_species);
            fp_t const sigma =  sig_e_m + sig_e_p;
            return sigma;
        } // sigma_lepton


        /**nucleon absorption */
        fp_t sigma_N_abs(cell_t const cell, fp_t const e_nu) const {
            nrgOK(e_nu);
            fp_t const sig_N_abs = m_rho.rho_nuc(cell) * sigmas::nu_N_abs(e_nu);
            return sig_N_abs;
        } // sigma_N_abs

        /** nucleon elastic scatter */
        fp_t sigma_N_elastic(cell_t const cell, fp_t const e_nu) const {
            nrgOK(e_nu);
            size_t const index = make_idx(cell,m_n_cells);
            // we treat protons and neutrons distinctly from the nuclei
            // for elastic scattering *if* we have composition information.
            fp_t const sig_N_el =
                sigmas::nu_N_elastic(e_nu) * (m_rho.rho_p[index]);
            fp_t sigma(sig_N_el);
            if( m_rho.have_composition == true)
            {
                fp_t const abar = m_rho.abar[index];
                sigma += m_rho.rho_A[index] * sigmas::nu_A_elastic(e_nu,abar);
            }
            return sigma;
        } // sigma_N_elastic


// nu-lepton cross sections, sampling a distribution via a power-law approx.
// to determine the lepton energy

        // template <typename rng_t>
        // fp_t sample_sigma_nu_e_minus(cell_t const cell,fp_t const e_nu,
        //                              Species const nu_species,
        //                              fp_t const alpha, rng_t & rng){
        //     //compute opacity for e- interaction
        //     cellOK(cell,m_n_cells);
        //     nrgOK(e_nu);
        //     fp_t sigma(0);

        //     cell_t const index = make_idx(cell,m_n_cells);
        //     fp_t const ebar = m_T.T_e_minus[index];
        //     nrgOK(ebar,"ebar");
        //     fp_t const rho = m_rho.rho_e_minus[index];
        //     fp_t e_e_minus = gen_power_law_energy(alpha,ebar,rng);
        //     nrgOK(e_e_minus,"energy e_minus");

        //     switch(nu_species)
        //     {
        //     case nu_e:
        //         sigma = sigma_nu_e_e_minus(e_nu, e_e_minus, rho);
        //         break;
        //     case nu_e_bar:
        //         sigma = sigma_nu_e_bar_e_minus(e_nu, e_e_minus, rho);
        //         break;
        //     case nu_x:
        //     case nu_mu:
        //     case nu_tau:
        //         sigma = sigma_nu_x_e_minus(e_nu, e_e_minus, rho);
        //         break;
        //     case nu_mu_bar:
        //     case nu_tau_bar:
        //         sigma = sigma_nu_x_bar_e_minus(e_nu, e_e_minus, rho);
        //         break;
        //     default:
        //         // should not get here
        //         unhandled_species_exc(nu_species,
        //                               "Opacity::sample_sigma_nu_e_minus");
        //     } // switch
        //     return sigma;
        // } // sample_sigma_nu_e_minus


        // template <typename rng_t>
        // fp_t sample_sigma_nu_e_plus(cell_t const cell,fp_cr e_nu,
        //                             Species const nu_species,
        //                             fp_cr alpha,rng_t & rng) {
        //     //compute opacity for e- interaction
        //     cellOK(cell,m_n_cells);
        //     nrgOK(e_nu);
        //     fp_t sigma(0);

        //     cell_t const index = make_idx(cell,m_n_cells);
        //     fp_t ebar = m_T.T_e_plus[index];
        //     nrgOK(ebar,"ebar");
        //     fp_t e_e_plus = gen_power_law_energy(alpha,ebar,rng);
        //     nrgOK(e_e_plus, "energy e_plus");
        //     fp_t const rho = m_rho.rho_e_plus[index];

        //     switch( nu_species)
        //     {
        //     case nu_e:
        //         sigma = sigma_nu_e_e_plus(cell,e_nu,alpha,rng, e_e_plus, rho);
        //         break;
        //     case nu_e_bar:
        //         sigma = sigma_nu_e_bar_e_plus(cell,e_nu,alpha,rng, e_e_plus, rho);
        //         break;
        //     case nu_mu:
        //     case nu_tau:
        //         sigma = sigma_nu_x_e_plus(cell,e_nu,alpha,rng, e_e_plus, rho);
        //         break;
        //     case nu_mu_bar:
        //     case nu_tau_bar:
        //         sigma = sigma_nu_x_bar_e_plus(cell,e_nu,alpha,rng, e_e_plus, rho);
        //         break;
        //     default:
        //         // should not get here
        //         unhandled_species_exc(nu_species,
        //                               "Opacity::sample_sigma_nu_e_plus");
        //     }
        //     return sigma;
        // } // sample_sigma_nu_e_plus


// nu-lepton cross sections, using mean temp as lepton energy

        fp_t sigma_nu_e_minus(cell_t const cell,fp_cr e_nu,
                              Species const nu_species) const {
            //compute opacity for e- interaction
            cellOK(cell,m_n_cells);
            nrgOK(e_nu);
            fp_t sigma(0);

            cell_t const index   = make_idx(cell,m_n_cells);
            fp_t const e_e_minus = m_T.T_e_minus[index];
            fp_t const rho       = m_rho.rho_e_minus[index];

            switch(nu_species)
            {
            case nu_e:
                sigma = sigma_nu_e_e_minus(e_nu,e_e_minus,rho);
                break;
            case nu_e_bar:
                sigma = sigma_nu_e_bar_e_minus(e_nu,e_e_minus,rho);
                break;
            case nu_x:
            case nu_mu:
            case nu_tau:
                sigma = sigma_nu_x_e_minus(e_nu,e_e_minus,rho);
                break;
            case nu_mu_bar:
            case nu_tau_bar:
                sigma = sigma_nu_x_bar_e_minus(e_nu,e_e_minus,rho);
                break;
            default:
                // should not get here
                unhandled_species_exc(nu_species,"Opacity::sigma_nu_e_minus");
            } // switch
            return sigma;
        }


        fp_t sigma_nu_e_plus(cell_t const cell,fp_cr e_nu,
                             Species const nu_species) const {
            cellOK(cell,m_n_cells);
            nrgOK(e_nu);
            fp_t sigma(0);
            cell_t const index = make_idx(cell,m_n_cells);

            fp_t const e_e_plus = m_T.T_e_plus[index];
            fp_t const rho = m_rho.rho_e_plus[index];

            switch( nu_species)
            {
            case nu_e:
                sigma = sigma_nu_e_e_plus(e_nu, e_e_plus, rho);
                break;
            case nu_e_bar:
                sigma = sigma_nu_e_bar_e_plus(e_nu, e_e_plus, rho);
                break;
            case nu_x:
            case nu_mu:
            case nu_tau:
                sigma = sigma_nu_x_e_plus(e_nu, e_e_plus, rho);
                break;
            case nu_mu_bar:
            case nu_tau_bar:
                sigma = sigma_nu_x_bar_e_plus(e_nu, e_e_plus, rho);
                break;
            default:
                // should not get here
                unhandled_species_exc(nu_species,"Opacity::sigma_nu_e_plus");
            }
            return sigma;
        } // sigma_nu_e_plus


        fp_t sigma_nu_e_e_minus(fp_cr e_nu, fp_cr e_e_minus, fp_cr rho) const {
            fp_t sigma = sigmas::nu_e_e_minus(e_nu,e_e_minus) * rho;
            return sigma;
        }


        fp_t sigma_nu_e_bar_e_minus(fp_cr e_nu, fp_cr e_e_minus, fp_cr rho) const {
            fp_t sigma = sigmas::nu_e_bar_e_minus(e_nu,e_e_minus) * rho;
            return sigma;
        }


        fp_t sigma_nu_x_e_minus(fp_cr e_nu, fp_cr e_e_minus, fp_cr rho) const {
            fp_t sigma = sigmas::nu_x_e_minus(e_nu,e_e_minus) * rho;
            return sigma;
        }


        fp_t sigma_nu_x_bar_e_minus(fp_cr e_nu, fp_cr e_e_minus, fp_cr rho) const {
            fp_t sigma = sigmas::nu_x_bar_e_minus(e_nu,e_e_minus) * rho;
            return sigma;
        }


        fp_t sigma_nu_e_e_plus(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho) const {
            fp_t sigma = sigmas::nu_e_e_plus(e_nu,e_e_plus) * rho;
            return sigma;
        }


        fp_t sigma_nu_e_bar_e_plus(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho) const {
            fp_t sigma = sigmas::nu_e_bar_e_plus(e_nu,e_e_plus) * rho;
            return sigma;
        }


        fp_t sigma_nu_x_e_plus(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho) const {
            fp_t sigma = sigmas::nu_x_e_plus(e_nu,e_e_plus) * rho;
            return sigma;
        }


        fp_t sigma_nu_x_bar_e_plus(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho) const {
            fp_t sigma = sigmas::nu_x_bar_e_plus(e_nu,e_e_plus) * rho;
            return sigma;
        }

        // to do: construct a template member function that takes
        // a functor as parameter. Replace all the above with one call
        // variously instantiated with calls to the sigmas.
        // template <typename x_sec>
        // static
        // fp_t sigma(fp_cr e_nu, fp_cr e_e_plus, fp_cr rho, x_sec xs){
        //     fp_t sigma = xs(e_nu,e_e_plus) * rho;
        //     return sigma;
        // } // sigma


        // state:
        cell_t const m_n_cells;

        density_t const & m_rho;

        temperature_t const & m_T;

    }; // class Opacity


} // nut::

// End of file
