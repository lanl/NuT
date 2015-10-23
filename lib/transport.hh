// transport.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef TRANSPORT_HH
#define TRANSPORT_HH

#include "apply_event.hh"
#include "decision.hh"
#include "Tally.hh"
#include "types.hh"


namespace nut
{

    template <typename ParticleT, typename MeshT, typename OpacityT,
              typename VelocityT, typename CensusT>
    ParticleT
    transport_particle(ParticleT const & in_p
                      ,MeshT const & mesh
                      ,OpacityT const & opacity
                      ,VelocityT const & velocity
                      ,Tally<typename OpacityT::fp_t> & tally
                      ,CensusT & census);

    /** transport a collection of particles, populating a collection of
     * new particles, accumulating statistics in a tally, banking particles
     * in a census, and optionally logging Monte Carlo steps. */
    template <typename PContainer, typename MeshT, typename OpacityT,
              typename VelocityT, typename fp_t = typename OpacityT::fp_t>
    void
    transport(PContainer const & p_source,
              uint32_t const n_particles,
              MeshT const & mesh,
              OpacityT const & opacity,
              VelocityT const & vel,
              Tally<fp_t> & tally,
              PContainer & p_sink,
              fp_t const alpha)
    {
        if(p_sink.size() != p_source.size())
        {
            p_sink.resize(p_source.size());
        }
        for(uint32_t ip = 0; ip < n_particles; ++ip)
        {
            p_sink[ip] =
                transport_particle(p_source[ip],mesh,opacity,vel,tally,alpha);
        }
        return;
    } // transport


    /** \brief transport a single particle until it's dead, gone, or its
        time is up. */
    template <typename ParticleT, typename MeshT, typename OpacityT,
              typename VelocityT, typename fp_t = typename OpacityT::fp_t>
    ParticleT
    transport_particle(ParticleT const & in_p,
                       MeshT const & mesh,
                       OpacityT const & opacity,
                       VelocityT const & vel,
                       Tally<fp_t> & tally,
                       fp_t const alpha)
    {
        ParticleT particle = in_p;
        while(particle.t >= 0.0 && particle.alive == true)
        {
            event_n_dist e_n_d = decide_event(particle,mesh,opacity,vel);
            Event const event = e_n_d.first;
            geom_t const dist         = e_n_d.second;
            apply_event(particle, event, dist, mesh, opacity, vel, tally, alpha);
        }
        return particle;
    } // transport_particle

} // nut::

#endif // TRANSPORT_HH

// End of file
