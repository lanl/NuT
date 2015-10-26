// Density.hh
// T. M. Kelley
// Jan. 10, 2010
// Header for Density
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef DENSITY_H
#define DENSITY_H

#include "Assert.hh"
#include "utilities_io.hh"
#include <vector>
#include <iostream>


namespace nut
{
    /*! \brief arrays to keep track of density in each mesh cell
     * \tparam <fp_t> {floating point type}
     */
    template <typename fp_t>
    struct Density
    {
        typedef std::vector<fp_t> vf;

        vf rho_p;
        vf rho_e_minus;
        vf rho_e_plus;
        vf rho_A;
        vf y_e;
        vf abar;

        cell_t const n_cells;

        bool have_composition;

        /*!\brief structure of arrays to store density and composition data.
         *
         * Densities should be number densities--sadly, do not enforce
         * that at present.

         * Densities are stored for the various quantities, as well as
         * (optionally) the electron number fraction y_e = n_e/(n_p+n_n),
         * and the mean atomic number. The density vectors must all be the same
         * size (one element per mesh cell).  If y_e and abar are provided,
         * they must also be the same size (n_cells); if not, they must have
         * size 0.
         *
         * Also note that y_e and abar are sketchy at present: they're defined
         * here, but likely not used properly elsewhere.
         */
        Density(vf const & rho_p_,
                vf const & rho_e_minus_,
                vf const & rho_e_plus_,
                vf const & rho_A_,
                vf const & y_e_,
                vf const & abar_
            )
            : rho_p(rho_p_),
              rho_e_minus(rho_e_minus_),
              rho_e_plus(rho_e_plus_),
              rho_A(rho_A_),
              n_cells(rho_p.size())
            {
                if( (y_e.empty() and !abar.empty() ) ||
                    (!y_e.empty() and abar.empty() )  )
                {
                    std::cerr << "Density ctor: abar and y_e must"
                              << " either both be full or both empty"
                              << std::endl;
                }
                y_e.assign(y_e_.begin(),y_e_.end());
                abar.assign(abar_.begin(),abar_.end());
                have_composition = !abar.empty();
                this->check_sizes();
                return;
            }

        // deep copy ctor
        explicit Density( Density const & d)
            : rho_p(d.rho_p),
              rho_e_minus(d.rho_e_minus),
              rho_e_plus(d.rho_e_plus),
              rho_A(d.rho_A),
              y_e(d.y_e),
              abar(d.abar),
              n_cells(d.n_cells),
              have_composition(d.have_composition)
            {}

        // size_t ctor
        Density(size_t const & sz,
                bool have_comp)
            : rho_p(sz),
              rho_e_minus(sz),
              rho_e_plus(sz),
              rho_A(sz),
              y_e(sz),
              abar(sz),
              n_cells(sz),
              have_composition(have_comp)
            {}

        /**Get neutron, proton, and A density summed. i in (0,n_cells]. */
        fp_t rho_nuc(cell_t const c) const {
            cell_t const i = make_idx(c,n_cells);
            return rho_p[i] + rho_A[i]; // rho_n[i];
        }


        /*!\brief check sizes of arrays for consistency at construction */
        void check_sizes(){
            size_t n_cells = rho_p.size();
            dbc::Equal(rho_e_minus.size(), n_cells,"rho_e_minus","n_cells");
            dbc::Equal(rho_e_plus.size(),  n_cells,"rho_e_plus","n_cells");
            dbc::Equal(rho_A.size(),       n_cells,"rho_A","n_cells");
            if(this->have_composition)
            {
                dbc::Equal(y_e.size(), n_cells,"y_e","n_cells");
                dbc::Equal(abar.size(),n_cells,"abar","n_cells");
            }
            return;
        } // check_sizes

        size_t size() const {return rho_p.size();}

    }; // Density

    template <typename fp_t>
    struct dens_t
    {
        fp_t nucleon; // unbound nucleons
        fp_t e_minus;
        fp_t e_plus;
    }; //density_t

} // nut::

#endif



// version
// $Id$

// End of file
