// transport.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef TRANSPORT_HH
#define TRANSPORT_HH

#include "Tally.hh"
#include "apply_event.hh"
#include "decision.hh"
#include "types.hh"
#include <sstream>
// #include <tr1/functional>  // bind & co
#include <functional>  // bind & co

namespace nut {

template <typename ParticleT,
          typename MeshT,
          typename OpacityT,
          typename VelocityT,
          typename CensusT,
          typename TallyT,
          typename LogT>
ParticleT
transport_particle(ParticleT const & in_p,
                   MeshT const & mesh,
                   OpacityT const & opacity,
                   VelocityT const & velocity,
                   TallyT & tally,
                   CensusT & census,
                   LogT & log,
                   typename TallyT::fp_t const alpha);

template <typename ParticleT,
          typename MeshT,
          typename OpacityT,
          typename VelocityT,
          typename CensusT,
          typename TallyT>
ParticleT
transport_particle_no_log(ParticleT const & in_p,
                          MeshT const & mesh,
                          OpacityT const & opacity,
                          VelocityT const & velocity,
                          TallyT & tally,
                          CensusT & census,
                          typename TallyT::fp_t const alpha);

/** transport a collection of particles, populating a collection of
 * new particles, accumulating statistics in a tally, banking particles
 * in a census, and optionally logging Monte Carlo steps. */
template <typename PContainer,
          typename MeshT,
          typename OpacityT,
          typename VelocityT,
          typename CensusT,
          typename TallyT,
          typename LogT>
void
transport(PContainer & p_source,
          MeshT const & mesh,
          OpacityT const & opacity,
          VelocityT const & vel,
          TallyT & tally,
          PContainer & p_sink,
          CensusT & census,
          LogT & log,
          typename TallyT::fp_t const alpha)
{
  typedef typename PContainer::value_type p_t;
  using namespace std::placeholders;

  if(p_sink.size() != p_source.size()) { p_sink.resize(p_source.size()); }

  // would be nice to make the log a functor that can be plugged
  // into the transport_particle fn.
  if(log.isNull()) {
    std::transform(
        p_source.begin(), p_source.end(), p_sink.begin(), [&](p_t & p) {
          return transport_particle_no_log<p_t, MeshT, OpacityT, VelocityT,
                                           CensusT, TallyT>(
              p, mesh, opacity, vel, tally, census, alpha);
        });
  }
  else {
    std::transform(p_source.begin(), p_source.end(), p_sink.begin(),
                   [&](p_t & p) {
                     return transport_particle<p_t, MeshT, OpacityT, VelocityT,
                                               CensusT, TallyT, LogT>(
                         p, mesh, opacity, vel, tally, census, log, alpha);
                   });
  }
  return;
}  // transport

template <typename ParticleT,
          typename MeshT,
          typename OpacityT,
          typename VelocityT,
          typename CensusT,
          typename TallyT,
          typename LogT>
ParticleT
transport_particle(ParticleT const & in_p,
                   MeshT const & mesh,
                   OpacityT const & opacity,
                   VelocityT const & vel,
                   TallyT & tally,
                   CensusT & census,
                   LogT & log,
                   typename TallyT::fp_t const alpha)
{
  ParticleT particle = in_p;
  cntr_t i = 1;
  while(particle.t >= 0.0 && particle.alive == true) {
    event_n_dist e_n_d = decide_event(particle, mesh, opacity, vel);
    events::Event const event = e_n_d.first;
    geom_t const dist = e_n_d.second;
    apply_event(particle, event, dist, mesh, opacity, vel, tally, census,
                alpha);
    i++;
  }
  return particle;
}  // transport_particle

template <typename ParticleT,
          typename MeshT,
          typename OpacityT,
          typename VelocityT,
          typename CensusT,
          typename TallyT>
ParticleT
transport_particle_no_log(ParticleT const & in_p,
                          MeshT const & mesh,
                          OpacityT const & opacity,
                          VelocityT const & vel,
                          TallyT & tally,
                          CensusT & census,
                          typename TallyT::fp_t const alpha)
{
  ParticleT particle = in_p;
  while(particle.t >= 0.0 && particle.alive == true) {
    event_n_dist e_n_d = decide_event(particle, mesh, opacity, vel);
    events::Event const event = e_n_d.first;
    geom_t const dist = e_n_d.second;
    apply_event(particle, event, dist, mesh, opacity, vel, tally, census,
                alpha);
  }
  return particle;
}  // transport_particle_no_log

}  // namespace nut

#endif  // TRANSPORT_HH

// End of file
