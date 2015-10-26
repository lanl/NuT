// sourcery.hh
// T. M. Kelley
// Mar 23, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef SOURCERY_HH
#define SOURCERY_HH

#include "Assert.hh"
#include "lorentz.hh"
#include "Planck.hh"
#include "utilities.hh"
#include <algorithm>
#include <numeric>
#include <ostream>

namespace nut
{
    /*!\brief for each cell, compute number of particles, energy weight of
     * each, and the total energy per cell. */
    template <typename FP_T, typename SZ_T>
    struct src_stats_t
    {
        typedef FP_T fp_t;
        typedef SZ_T sz_t;
        typedef std::vector<sz_t> vsz;
        typedef std::vector<fp_t> vfp;
        vsz ns,cidxs;
        vfp ews, es;
        explicit src_stats_t(size_t const n)
            : ns(n, sz_t(0) ),
              cidxs(n, sz_t(0) ),
              ews(n,fp_t(0) ),
              es(n ,fp_t(0) )
        {}

        size_t
        size() const
        {
            dbc::Require(ns.size() == ews.size()
                    && ns.size() == es.size()
                    && ns.size() == cidxs.size(),
                    "src_stats_t.size: sizes disagree");
            return ns.size();
        } // size

        void
        push_back(sz_t const n, sz_t const c, fp_t const e, fp_t const w)
        {
            ns.push_back(n);
            cidxs.push_back(c);
            es.push_back(e);
            ews.push_back(w);
            return;
        }

        void
        resize(size_t const s)
        {
            ns.resize(s);
            cidxs.resize(s);
            es.resize(s);
            ews.resize(s);
            return;
        }

        void
        reserve(size_t const s)
        {
            ns.reserve(s);
            cidxs.reserve(s);
            es.reserve(s);
            ews.reserve(s);
            return;
        }

        size_t
        n_particles() const
        {
            return std::accumulate(ns.begin(),ns.end(),sz_t(0));
        } // n_particles


    }; // src_stats_t


    template <typename fp_t, typename sz_t>
    void calc_src_stats_lum(
        std::vector<fp_t> const & lums,
        std::vector<sz_t> const & cidxs,
        src_stats_t<fp_t,sz_t> & stats,  // output
        fp_t const dt,
        size_t const ntot,
        size_t const ncen
        )
    {
        dbc::Require(ntot >= ncen,
                "calc_src_stats_lum: ntot must be >= ncen");
        dbc::Equal(stats.size(), lums.size(),"stats.size()", "lums.size()");
        dbc::Equal(stats.size(),cidxs.size(),"stats.size()","idxs.size()");
        // populate energies
        std::transform(lums.begin(),lums.end(),
                       stats.es.begin(),
                       mult<fp_t>(dt));
        std::copy(cidxs.begin(),cidxs.end(),stats.cidxs.begin());
        // compute # per cell
        fp_t const ev_tot = sum(stats.es);
        dbc::GreaterThan(ev_tot,fp_t(0),"total energy must be > 0");

        std::vector<fp_t> fracs(lums.size(),fp_t(0));
        std::transform(stats.es.begin(),stats.es.end(),
                       fracs.begin(),
                       div_by<fp_t>(ev_tot));

        size_t const nv_tot(ntot - ncen);
        dbc::GreaterThan(nv_tot,size_t(0),
                    "total # particles must be >= # census particles");
        std::transform(fracs.begin(),fracs.end(),
                       stats.ns.begin(),
                       int_mult<fp_t,size_t>( fp_t(nv_tot) ));
        // compute energy weight per cell
        std::transform(stats.es.begin(),stats.es.end(),
                       stats.ns.begin(),
                       stats.ews.begin(),
                       div_if_ne0<fp_t>);

        return;
    } // calc_src_stats_lum

    /*!\brief generate a single particle
     * */
    template <typename mesh_t, typename geom_t, typename rng_t, typename part_t>
    part_t
    gen_init_particle(mesh_t const & mesh,
                      cell_t const cell,
                      geom_t const dt,
                      geom_t const alpha,
                      Species const s,
                      geom_t const ew, // e wt
                      geom_t const T,  // temp/cell
                      vec_t<part_t::dim> const v,
                      rng_t & rng
        )
    {
        uint32_t const dim = part_t::dim;

        // Mean energy of Fermionic Planckian is 7 pi^4/(180 zeta(3) * (k_B T).
        // Prefactor ~3.15137
        geom_t const ebar = 3.15137 * T;
        geom_t const urd  = rng.random();
        geom_t const r    = mesh.sample_position(urd,cell);

        geom_t const osd  = rng.random();

        geom_t const oc   = 2.0 * osd - 1.0;
        geom_t const ec   = gen_power_law_energy(alpha, ebar, rng);
        EandOmega<dim> enol    = mesh_t::LT_to_lab(v, ec, oc);
        geom_t const e    = enol.first;
        vec_t<dim> const o    = enol.second;

        return part_t(r,o,e,dt,ew,cell,rng,s);
    } // gen_init_particle

} // nut::


#endif //SOURCERY_HH

// version
// $Id$

// End of file
