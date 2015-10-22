// MatState.hh
// T. M. Kelley
// Jul 13, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

/* Intent is to aggregate object creation between fileio and app */

#ifndef MATSTATE_HH
#define MATSTATE_HH

#include "Density.hh"
#include "Luminosity.hh"
#include "Temperature.hh"
#include "Velocity.hh"
#include "fileio.hh"
#include "constants.hh"
#include "Assert.hh"
#include <algorithm>

namespace nut
{
    namespace P // transformations specific to the p.1-p.7 file format
    {
        template <typename fp_t,size_t dim = 1>
        struct state_entry_t
        {
            dens_t<fp_t>  d;
            temp_t<fp_t>  t;
            lum_t<fp_t>   l;
            vel_t<dim> v;
        };

        template <typename fp_t, size_t dim = 1>
        state_entry_t<fp_t, dim> row_to_state_entry(MatStateRowP<fp_t,dim> const & row);

    } // P::

    template <typename fp_t,size_t dim = 1>
    struct MatState
    {
        typedef Density<fp_t>      Density_T;
        typedef Luminosity<fp_t>   Luminosity_T;
        typedef Temperature<fp_t>  Temperature_T;
        typedef Velocity<fp_t,dim> Velocity_T;

        Density_T      density;
        Luminosity_T   luminosity;
        Temperature_T  temperature;
        Velocity_T     velocity;

        explicit MatState(std::vector< MatStateRowP<fp_t,dim> > const & rows);

    private:
        void add_entry( P::state_entry_t<fp_t,dim> const & entry, size_t const pos);
    }; // MatState

    // ctor

    template <typename fp_t,size_t dim>
    MatState<fp_t,dim>::MatState(std::vector< MatStateRowP<fp_t,dim> > const & rows)
        : density(rows.size(),false)
        , luminosity(rows.size())
        , temperature(rows.size())
        , velocity(rows.size())
    {
        // for each line,
        for(size_t i = 0; i < rows.size(); ++i)
        {
            // 1. extract & convert interesting quantities,
            P::state_entry_t<fp_t,dim> entry = P::row_to_state_entry(rows[i]);
            // 2. push back onto relevant arrays
            add_entry(entry, i);
        }
        return;
    } // MatState :: ctor

    template <typename fp_t,size_t dim>
    void
    MatState<fp_t,dim>::add_entry(P::state_entry_t<fp_t,dim> const & entry,
                              size_t const i
        )
    {
        density.rho_p[i]       = entry.d.nucleon;
        density.rho_e_minus[i] = entry.d.e_minus;
        density.rho_e_plus[i]  = entry.d.e_plus;

        luminosity.nue[i]  = entry.l.nue;
        luminosity.nueb[i] = entry.l.nueb;
        luminosity.nux[i]  = entry.l.nux;

        temperature.T_p[i]       = entry.t.t_p;
        temperature.T_e_minus[i] = entry.t.t_e_m;
        temperature.T_e_plus[i]  = entry.t.t_e_p;

        velocity.vs[i] = entry.v.v;
    }

    namespace
    {
        template <typename fp_t,size_t dim>
        fp_t extract_radius(MatStateRowP<fp_t,dim> const & row){return row.radius;}
    }

    /*!\brief: create a mesh within the given limits; identify the limiting
    * indices in the rows vector. The limiting indices use STL begin,end
    * convention. */
    template <typename mesh_t, typename fp_t, size_t dim = 1>
    mesh_t rows_to_mesh(std::vector<MatStateRowP<fp_t,dim> > const & rows,
                        fp_t const llimit,
                        fp_t const ulimit,
                        size_t & llimitIdx,
                        size_t & ulimitIdx
        )
    {
        // bool const  lims_ok = (ulimit > llimit);
        // Require(lims_ok , "rows_to_mesh: lower limit >= upper limit");
        size_t nrows = rows.size();
        std::vector<typename mesh_t::geom_t> bndsTmp(nrows + 1);
        size_t lIdx = 0;
        size_t uIdx = nrows;
        // get radii from rows
        std::vector<fp_t> rads(nrows);
        std::transform(rows.begin(),rows.end(),rads.begin(),extract_radius<fp_t,dim>);
        bndsTmp[0] = rads[0] - (rads[1] - rads[0])/2;
        for(size_t i = 0; i < nrows-1; ++i)
        {
            bndsTmp[i+1] = rads[i] + (rads[i+1] - rads[i])/2;
        }
        bndsTmp[nrows] = rads[nrows-1] + (rads[nrows-1] - rads[nrows-2])/2;

        // would use STL algs here, but need a actual indices
        // find limiting indices: want the greatest index such that
        // bndsTmp[i] < llimit
        while(bndsTmp[lIdx+1] < llimit && lIdx < nrows){lIdx++;}
        if(lIdx == nrows)
        {
            std::stringstream errstr;
            errstr << "rows_to_mesh: lower bound (" << llimit
                   << ") not within range of input (max input radius: "
                   << rads[nrows-1] << ", which comes out to upper bdy of "
                   << bndsTmp[nrows] << std::endl;
            throw(std::runtime_error(errstr.str()));
        }
        // now look for the least upper index such that
        // bndsTmp[uIdx] > ulimit
        uIdx = lIdx;
        while(bndsTmp[uIdx] < ulimit && uIdx < nrows){uIdx++;}

        // outputs
        ulimitIdx = uIdx;
        llimitIdx = lIdx;
        size_t const ncells = uIdx - lIdx;
        std::vector<typename mesh_t::geom_t> bounds(ncells + 1);
        std::vector<typename mesh_t::bdy_desc_t> descs(ncells + 1);

        std::copy(&bndsTmp[lIdx],&bndsTmp[uIdx+1],bounds.begin());
        descs[0] = bdy_types::descriptor::R;
        descs[ncells] = bdy_types::descriptor::V;
        for(size_t i = 1; i < ncells; ++i){descs[i] = bdy_types::descriptor::T;}
        return mesh_t(bounds,descs);
    } // rows_to_mesh



    namespace P
    {
        /*! functions for extracting and converting particular quantities
         * from a MatStateRow */
        template <typename fp_t,size_t dim = 1>
        dens_t<fp_t> row_to_rho(MatStateRowP<fp_t,dim> const & row)
        {
            dens_t<fp_t> r;
            r.nucleon = row.density/pmg;
            r.e_minus = row.ye * r.nucleon;
            r.e_plus  = fp_t(0); // no positrons in this model (???)
            return r;
        }

        template <typename fp_t,size_t dim = 1>
        temp_t<fp_t> row_to_temp(MatStateRowP<fp_t,dim> const & row)
        {
            temp_t<fp_t> t;
            // convert to MeV: T's in file are in 10^9 K, k_B in MeV/K
            t.t_e_m = t.t_e_p = t.t_p = k_B * fp_t(1e9) * row.temperature;
            return t;
        }

        template <typename fp_t,size_t dim = 1>
        lum_t<fp_t> row_to_lum(MatStateRowP<fp_t,dim> const & row)
        {
            lum_t<fp_t> l;
            l.nue  = 1e51*(row.lnue_capture  + row.lnue_pair);
            l.nueb = 1e51*(row.lnueb_capture + row.lnueb_pair);
            l.nux  = 1e51*(row.lnux_pair);
            return l;
        }

        template <typename fp_t,size_t dim = 1>
        vel_t<dim> row_to_vel(MatStateRowP<fp_t,dim> const & row)
        {
            vel_t<dim> v;
            v.v = row.velocity;
            return v;
        }

        template <typename fp_t,size_t dim>
        state_entry_t<fp_t,dim> row_to_state_entry(MatStateRowP<fp_t,dim> const & row)
        {
            state_entry_t<fp_t,dim> entry;
            entry.d = row_to_rho(row);
            entry.t = row_to_temp(row);
            entry.l = row_to_lum(row);
            entry.v = row_to_vel(row);
            return entry;
        } // row_to_state_entry

    } // P::

} // nut::

#endif // include guard

// version
// $Id$

// End of file
