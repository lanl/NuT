// application.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef APPLY_EVENT_HH
#define APPLY_EVENT_HH

/**\file apply events to particles */

#include "Event.hh"
#include "Fates.hh"
#include "Planck.hh"
#include "Tally.hh"
#include "Velocity.hh"
#include "constants.hh"
#include "lorentz.hh"
#include "types.hh"
#include <stdexcept>
#include <sstream>
// #include <iomanip>


namespace nut
{
    template <typename p_t, typename tally_t,
              typename vel_t, typename Mesh_T>
    void apply_nucleon_elastic_scatter(p_t & p, tally_t & t, vel_t const & vel);

    template <typename p_t, typename tally_t,
              typename velocity_t, typename Mesh_T>
    void apply_lepton_scatter(p_t & p, tally_t & t, typename tally_t::FP_T const e_lep,
                              velocity_t const & vel);

    template <typename Mesh_T, typename Particle_T, typename fp_t>
    void stream_particle(Particle_T & p, fp_t const d, Mesh_T const & mesh)
    {
        typename Mesh_T::coord_t newcoord =
            mesh.new_coordinate(p.x,p.omega,d);
        p.x     = newcoord.x;
        p.omega = newcoord.omega;
        p.t     = p.t - d / c;
        return;
    } // stream_particle

    template <typename ParticleT, typename MeshT, typename OpacityT,
              typename VelocityT, typename fp_t>
    void
    apply_event(ParticleT & p,
                Event const & event,
                geom_t const distance,
                MeshT const & mesh,
                OpacityT const & opacity,
                VelocityT const & velocity,
                Tally<fp_t> & tally,
                fp_t const alpha = fp_t(2.0)
        )
    {
        cell_t const index = make_idx(p.cell,opacity.m_n_cells);

        stream_particle(p,distance,mesh);

        tally.accum_pl(distance);

        switch(event)
        {
        case Event::nucleon_abs:
            apply_nucleon_abs(p,tally);
            break;
        case Event::nucleon_elastic_scatter:
            apply_nucleon_elastic_scatter<ParticleT,Tally<fp_t>,VelocityT,MeshT>( p,tally,velocity);
            break;
        case Event::electron_scatter:
            {
                fp_t const ebar = opacity.m_T.T_e_minus[index];
                fp_t const e_e  = ebar;
                // fp_t const e_e  = gen_power_law_energy_alpha2(ebar,p.rng);
                apply_lepton_scatter<ParticleT,Tally<fp_t>,VelocityT,MeshT>(p,tally,e_e,velocity);
            }
            break;
        case Event::positron_scatter:
            {
                fp_t const ebar = opacity.m_T.T_e_plus[index];
                fp_t const e_e  = ebar;
                // fp_t const e_e  = gen_power_law_energy_alpha2(ebar,p.rng);
                apply_lepton_scatter<ParticleT,Tally<fp_t>,VelocityT,MeshT>(p,tally,e_e,velocity);
            }
            break;
        // case Event::nu_e_annhilation:
        //     break;
        // case Event::nu_x_annhilation:
        //     break;
        case Event::cell_low_x_boundary:
            apply_low_x_boundary(p,tally);
            break;
        case Event::cell_high_x_boundary:
            apply_hi_x_boundary(p,tally);
            break;
        case Event::escape:
            apply_escape(p,tally);
            break;
        case Event::reflect:
            apply_reflect(p,tally);
            break;
        case Event::step_end:
            apply_step_end(p,tally);
            break;
        // error if we get to an unresolved collision or boundary
        case Event::undetermined:
            p.kill(Fates::UNDETERMINED_EVENT);
            break;
        case Event::boundary:   // fall through to next
        case Event::collision:  //
            // underresolved_event(event);
            // break;
        default:
            p.kill(Fates::UNHANDLED_EVENT);
        } // switch
        if(event != Event::nucleon_elastic_scatter
           && event != Event::electron_scatter
           && event != Event::positron_scatter)
        {
            p.rng.random(); // compensate to maintain constant RNs/MC step
        }
        return;
    } // apply_event

    template <typename p_t, typename tally_t>
    void apply_nucleon_abs(p_t & p, tally_t & t)
    {
        t.deposit_energy(p.cell,p.weight,p.e);
        t.deposit_momentum_elastic(p.cell,p.omega,p.e,p.weight);
        t.count_nucleon_abs(p.cell,p.species,p.weight);
        p.kill(Fates::NUCLEON_ABS);
        return;
    } // apply_nucleon_abs


    template <typename p_t, typename tally_t,
              typename vel_t, typename Mesh_T>
    void apply_nucleon_elastic_scatter(p_t & p, tally_t & t, vel_t const & vel)
    {
        // typedef typename tally_t::FP_T fp_t;
        size_t const dim(tally_t::dim);
        cell_t const cell     = p.cell;
        vec_t<dim> const v    = vel.v(cell);
        geom_t const eli      = p.e;
        vec_t<dim> const oli  = p.omega;
        // LT to comoving frame (need init comoving e to compute
        // final lab e).

        EandOmega<dim> eno_cmi = Mesh_T::LT_to_comoving(v,eli,oli);
        geom_t const eci  = eno_cmi.first;
        // scatter: sample a new direction
        vec_t<dim> ocf = Mesh_T::sample_direction(p.rng);

        // LT comoving -> lab
        EandOmega<dim> const eno_lf = Mesh_T::LT_to_lab(v,eci,ocf);
        geom_t const elf = eno_lf.first;
        vec_t<dim> const olf = eno_lf.second;
        // tally
        t.count_nucleon_elastic_scatter(p.cell);
        t.deposit_inelastic_scat(cell, eli, elf, oli, olf,
                                 p.weight, p.species);
        // update particle
        p.omega = olf;
        p.e = elf;

        return;
    }

    template <typename p_t, typename tally_t,
              typename velocity_t, typename Mesh_T>
    void apply_lepton_scatter(p_t & p, tally_t & t, typename tally_t::FP_T const e_lep,
                              velocity_t const & vel)
    {
        typedef typename tally_t::FP_T fp_t;
        uint32_t const dim(p_t::dim);
        cell_t const cell     = p.cell;
        vec_t<dim> const v    = vel.v(cell);
        geom_t const eli      = p.e;
        vec_t<dim> const oli  = p.omega;

        // LT to comoving frame (need init comoving e to compute
        // final lab e).
        EandOmega<dim> eno_cmi = Mesh_T::LT_to_comoving(v,eli,oli);
        geom_t const eci  = eno_cmi.first;
        // scatter
        vec_t<dim> ocf = Mesh_T::sample_direction(p.rng);
        fp_t const de = (eci - e_lep)/4.;
        fp_t const ecf = eci - de;
        // LT comoving -> lab
        EandOmega<dim> const eno_lf = Mesh_T::LT_to_lab(v,ecf,ocf);
        geom_t const elf = eno_lf.first;
        vec_t<dim> const olf = eno_lf.second;
        // tally
        t.count_lepton_scatter(p.cell,p.species);
        t.deposit_inelastic_scat(cell,eli,elf,oli,olf,p.weight,
                                 p.species);
        // update
        p.omega = olf;
        p.e     = elf;
        return;
    } // apply_lepton_scatter


    template <typename p_t, typename tally_t>
    void apply_low_x_boundary(p_t & p, tally_t & t)
    {
        t.count_cell_bdy(p.cell);
        p.cell -= 1;
        return;
    }

    template <typename p_t, typename tally_t>
    void apply_hi_x_boundary(p_t & p, tally_t & t)
    {
        t.count_cell_bdy(p.cell);
        p.cell += 1;
        return;
    }

    template <typename p_t, typename tally_t>
    void apply_escape(p_t & p, tally_t & t)
    {
        t.count_escape(p.cell,p.weight,p.e);
        p.kill(Fates::ESCAPE);
        return;
    }

    template <typename p_t, typename tally_t>
    void apply_reflect(p_t & p, tally_t & t)
    {
        t.count_reflect(p.cell);
        p.omega = -p.omega;
        return;
    }

    template <typename p_t, typename tally_t>
    void apply_step_end(p_t & p, tally_t & t)
    {
        t.count_census(p.cell,p.weight,p.species);
        p.kill(Fates::STEP_END);
        return;
    }

} // nut::


#endif // include guard

// version
// $Id$

// End of file


/*
"collision","nucleon_abs","nucleon_elastic_scatter","electron_scatter","positron_scatter","nu_e_annhilation","nu_x_annhilation","boundary","cell_low_x_boundary","cell_high_x_boundary","escape","reflect","step_end"
*/
