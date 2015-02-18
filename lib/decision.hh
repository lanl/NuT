// decision.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef DECISION_HH
#define DECISION_HH

#include "types.hh"
#include "constants.hh"
#include <utility> // pair
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "lorentz.hh"
// #include <iomanip>

namespace nut
{

    namespace
    {
        /** throw an exception when code fails to decide an event */
        void
        unresolved_event(size_t lineno, std::string const & eventstr);

        /** Compare std::pairs on the second element. */
        template <typename pair> bool
        pair_min_2nd(pair const & p1, pair const & p2){
            return p1.second < p2.second;}

    }

    template <typename particle_t, typename mesh_t, typename opacity_t,
              typename velocity_t>
    event_n_dist
    decide_event( particle_t & p, // non-const b/c of RNG
                  mesh_t const & mesh,
                  opacity_t const & opacity,
                  velocity_t const & velocity)
    {
        typedef typename particle_t::fp_t fp_t;
        // typedef std::pair<events::Event,geom_t> event_n_dist;
        typedef std::vector<event_n_dist> vend_t;
        typedef vend_t::iterator vend_it;
        typedef typename mesh_t::d_to_b_t d_to_b_t;

        vend_t e_n_ds(0);
        e_n_ds.reserve(3);

        cell_t const cell  = p.cell;
        fp_t const tleft   = p.t;
        fp_t const x       = p.x;
        Species const species = p.species;

        // compute distance to events, push onto vector
        // compute cross-section in comoving frame

        geom_t v = velocity.v(cell);
        geom_t const eli  = p.e;
        geom_t const oli  = p.omega;
        // LT to comoving frame (compute interaction comoving).
        EandOmega eno_cmi = LT_to_comoving_sphere1D(v,eli,oli);
        geom_t const eci  = eno_cmi.first;

        fp_t const sig_coll   = opacity.sigma_collide(cell,eci,species);

        fp_t const random_dev = p.rng.random();
        fp_t const ignored    = p.rng.random(); // to keep pace with McPhD

        geom_t const d_coll     = (sig_coll != fp_t(0)) ?
            -std::log(random_dev)/sig_coll : huge;
        e_n_ds.push_back(event_n_dist(events::collision,d_coll));

        d_to_b_t dnf = mesh.distance_to_bdy(x,oli,cell);
        geom_t const d_bdy = dnf.d;
        e_n_ds.push_back(event_n_dist(events::boundary,d_bdy));

        geom_t const d_step_end = c * tleft;
        e_n_ds.push_back(event_n_dist(events::step_end,d_step_end));

        // pick (event,dist) pair with shortest distance
        event_n_dist closest = *(std::min_element(e_n_ds.begin(),e_n_ds.end(),
                                                  pair_min_2nd<event_n_dist>));

        // if needed, further resolve events
        if(closest.first == events::collision)
        {
            // need to compute the comoving energy at the scattering site,
            // in order to compute the different cross sections there.

            typename mesh_t::coord_t const scat_site(
                mesh.new_coordinate(x,oli,d_coll));
            geom_t const o_sct = scat_site.omega;
            EandOmega const eno_scat = LT_to_comoving_sphere1D(v,eli,o_sct);
            fp_t const ecscat = eno_scat.first;
            closest.first =
                decide_scatter_event(p.rng,ecscat,cell,opacity,species);
        }
        else
        {
            if(closest.first == events::boundary)
            {
                closest.first = decide_boundary_event(mesh,cell,dnf.face);
            }
            // compensate RNG if decide_scatter_event not called
            p.rng.random();
        }
        return closest;
    } // decide_event


    template <typename MeshT>
    events::Event
    decide_boundary_event( MeshT const & mesh, cell_t const cell,
                           cell_t const face)
    {
        using namespace bdy_types;
        using namespace events;
        Event event(null);
        bdy_types::descriptor b_type( mesh.bdy_type(cell,face));
        switch(b_type)
        {
        case V:
            event = escape;
            break;
        case R:
            event = reflect;
            break;
        case T:
            // 1D specific
            event = face == 0 ? cell_low_x_boundary : cell_high_x_boundary;
            break;
        default:
            unresolved_event(__LINE__,"decide_boundary_event");
        };
        return event;
    } // decide_boundary_event


    template <typename rng_t, typename fp_t, typename opacity_t>
    events::Event
    decide_scatter_event(rng_t & rng, fp_t const nu_nrg, cell_t const cell,
                         opacity_t const & op, Species const s)
    {
        typedef const fp_t fp_c;
        // cell_t const idx = cell - 1;
        fp_c sig_collide = op.sigma_collide(cell,nu_nrg,s);

        // selectors
        fp_c p_type_sel = rng.random();  // particle interactor
        events::Event event(events::null);

        fp_c N_total = op.sigma_N_total(cell,nu_nrg);    // nucleon
        fp_c prob_nucleon = N_total/sig_collide;

        // other top-level cross-sections here

        fp_c prob_abs      = op.sigma_N_abs(cell,nu_nrg)/sig_collide;
        fp_c prob_elastic  = op.sigma_N_elastic(cell,nu_nrg)/sig_collide;
        fp_c prob_electron = op.sigma_nu_e_minus(cell,nu_nrg,s)/sig_collide;
        fp_c prob_positron = op.sigma_nu_e_plus( cell,nu_nrg,s)/sig_collide;

        if(p_type_sel < prob_abs)
        {
            event = events::nucleon_abs;
        }
        else if(p_type_sel < prob_abs + prob_elastic)
        {
            event = events::nucleon_elastic_scatter;
        }
        else if(p_type_sel < prob_nucleon + prob_electron)
        {
            event = events::electron_scatter;
        }
        else if(p_type_sel < prob_nucleon + prob_electron + prob_positron)
        {
            event = events::positron_scatter;
        }
        else
        {
            unresolved_event(__LINE__,"scatter");
        }
        return event;
    } // decide_scatter_event


    namespace
    {
        void
        unresolved_event(size_t lineno, std::string const & eventstr)
        {
            std::stringstream msg;
            msg << "decision.hh:" << lineno << ": unresolved event: "
                << eventstr;
            throw std::domain_error(msg.str());
        } // unresolved_event
    } // anonymous::


} // nut::


#endif // include guard

// version
// $Id$

// End of file
