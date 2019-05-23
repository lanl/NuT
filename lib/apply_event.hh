// application.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef APPLY_EVENT_HH
#define APPLY_EVENT_HH

/**\file apply events to particles */
#include "Planck.hh"

#include "Tally.hh"
#include "Velocity.hh"
#include "constants.hh"
#include "lorentz.hh"
#include "types.hh"
#include <iomanip>
#include <sstream>
#include <stdexcept>

/* unknown_event and underresolved_event are not used in all translation units,
 * leading to unexciting warnings. Suppress those. */
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunneeded-internal-declaration"

namespace nut {
namespace {
void
underresolved_event(nut::events::Event const & event);
void
unknown_event(nut::events::Event const & event);
}  // namespace

template <typename p_t, typename tally_t, typename vel_t, typename Mesh_T>
void
apply_nucleon_elastic_scatter(p_t & p, tally_t & t, vel_t const & vel);

template <typename p_t, typename tally_t, typename velocity_t, typename Mesh_T>
void
apply_lepton_scatter(p_t & p,
                     tally_t & t,
                     typename tally_t::FP_T const e_lep,
                     velocity_t const & vel);

template <typename Mesh_T, typename Particle_T, typename fp_t>
void
stream_particle(Particle_T & p, fp_t const d, Mesh_T const & mesh)
{
  typename Mesh_T::coord_t newcoord = mesh.new_coordinate(p.x, p.omega, d);
  p.x = newcoord.x;
  p.omega = newcoord.omega;
  p.t = p.t - d / c;
  return;
}  // stream_particle

template <typename ParticleT,
          typename MeshT,
          typename OpacityT,
          typename VelocityT,
          typename CensusT,
          typename TallyT>
void
apply_event(ParticleT & p,
            events::Event const & event,
            geom_t const distance,
            MeshT const & mesh,
            OpacityT const & opacity,
            VelocityT const & velocity,
            TallyT & tally,
            CensusT & census,
            typename TallyT::fp_t const alpha = typename TallyT::fp_t(2.0))
{
  using namespace events;
  using fp_t = typename TallyT::fp_t;

  cell_t const index = make_idx(p.cell, opacity.m_n_cells);

  stream_particle(p, distance, mesh);

  tally.accum_pl(distance);

  switch(event) {
    case nucleon_abs: apply_nucleon_abs(p, tally); break;
    case nucleon_elastic_scatter:
      apply_nucleon_elastic_scatter<ParticleT, TallyT, VelocityT, MeshT>(
          p, tally, velocity);
      break;
    case electron_scatter: {
      fp_t const ebar = opacity.m_T.T_e_minus[index];
      fp_t const e_e = ebar;
      // fp_t const e_e  = gen_power_law_energy_alpha2(ebar,p.rng);
      apply_lepton_scatter<ParticleT, TallyT, VelocityT, MeshT>(p, tally, e_e,
                                                                velocity);
    } break;
    case positron_scatter: {
      fp_t const ebar = opacity.m_T.T_e_plus[index];
      fp_t const e_e = ebar;
      // fp_t const e_e  = gen_power_law_energy_alpha2(ebar,p.rng);
      apply_lepton_scatter<ParticleT, TallyT, VelocityT, MeshT>(p, tally, e_e,
                                                                velocity);
    } break;
    // case nu_e_annhilation:
    //     break;
    // case nu_x_annhilation:
    //     break;
    case cell_low_x_boundary: apply_low_x_boundary(p, tally); break;
    case cell_high_x_boundary: apply_hi_x_boundary(p, tally); break;
    case escape: apply_escape(p, tally); break;
    case reflect: apply_reflect(p, tally); break;
    case step_end: apply_step_end(p, tally, census); break;
    // error if we get to an unresolved collision or boundary
    case boundary:  // fall through to next
    case collision: underresolved_event(event); break;
    default: unknown_event(event);
  }  // switch
  if(event != nucleon_elastic_scatter && event != electron_scatter &&
     event != positron_scatter) {
    p.rng.random();  // compensate to maintain constant RNs/MC step
  }
  return;
}  // apply_event

template <typename p_t, typename tally_t>
void
apply_nucleon_abs(p_t & p, tally_t & t)
{
  t.deposit_energy(p.cell, p.weight, p.e);
  t.deposit_momentum_elastic(p.cell, p.omega, p.e, p.weight);
  t.count_nucleon_abs(p.cell, p.species, p.weight);
  p.kill();
  return;
}  // apply_nucleon_abs

template <typename p_t, typename tally_t, typename vel_t, typename Mesh_T>
void
apply_nucleon_elastic_scatter(p_t & p, tally_t & t, vel_t const & vel)
{
  // typedef typename tally_t::FP_T fp_t;
  size_t const dim(tally_t::dim);
  cell_t const cell = p.cell;
  vec_t<dim> const v = vel.v(cell);
  geom_t const eli = p.e;
  vec_t<dim> const oli = p.omega;
  // LT to comoving frame (need init comoving e to compute
  // final lab e).

  EandOmega<dim> eno_cmi = Mesh_T::LT_to_comoving(v, eli, oli);
  geom_t const eci = eno_cmi.first;
  // scatter: sample a new direction
  vec_t<dim> ocf = Mesh_T::sample_direction(p.rng);

  // LT comoving -> lab
  EandOmega<dim> const eno_lf = Mesh_T::LT_to_lab(v, eci, ocf);
  geom_t const elf = eno_lf.first;
  vec_t<dim> const olf = eno_lf.second;
  // tally
  t.count_nucleon_elastic_scatter(p.cell);
  t.deposit_inelastic_scat(cell, eli, elf, oli, olf, p.weight, p.species);
  // update particle
  p.omega = olf;
  p.e = elf;

  return;
}

template <typename p_t, typename tally_t, typename velocity_t, typename Mesh_T>
void
apply_lepton_scatter(p_t & p,
                     tally_t & t,
                     typename tally_t::FP_T const e_lep,
                     velocity_t const & vel)
{
  typedef typename tally_t::FP_T fp_t;
  uint32_t const dim(p_t::dim);
  cell_t const cell = p.cell;
  vec_t<dim> const v = vel.v(cell);
  geom_t const eli = p.e;
  vec_t<dim> const oli = p.omega;

  // LT to comoving frame (need init comoving e to compute
  // final lab e).
  EandOmega<dim> eno_cmi = Mesh_T::LT_to_comoving(v, eli, oli);
  geom_t const eci = eno_cmi.first;
  // scatter
  vec_t<dim> ocf = Mesh_T::sample_direction(p.rng);
  fp_t const de = (eci - e_lep) / 4.;
  fp_t const ecf = eci - de;
  // LT comoving -> lab
  EandOmega<dim> const eno_lf = Mesh_T::LT_to_lab(v, ecf, ocf);
  geom_t const elf = eno_lf.first;
  vec_t<dim> const olf = eno_lf.second;
  // tally
  t.count_lepton_scatter(p.cell, p.species);
  t.deposit_inelastic_scat(cell, eli, elf, oli, olf, p.weight, p.species);
  // update
  p.omega = olf;
  p.e = elf;
  return;
}  // apply_lepton_scatter

template <typename p_t, typename tally_t>
void
apply_low_x_boundary(p_t & p, tally_t & t)
{
  t.count_cell_bdy(p.cell);
  p.cell -= 1;
  return;
}

template <typename p_t, typename tally_t>
void
apply_hi_x_boundary(p_t & p, tally_t & t)
{
  t.count_cell_bdy(p.cell);
  p.cell += 1;
  return;
}

template <typename p_t, typename tally_t>
void
apply_escape(p_t & p, tally_t & t)
{
  t.count_escape(p.cell, p.weight, p.e);
  p.kill();
  return;
}

template <typename p_t, typename tally_t>
void
apply_reflect(p_t & p, tally_t & t)
{
  t.count_reflect(p.cell);
  p.omega = -p.omega;
  return;
}

template <typename p_t, typename tally_t, typename census_t>
void
apply_step_end(p_t & p, tally_t & t, census_t & c)
{
  t.count_census(p.cell, p.weight, p.species);
  p.alive = false;
  c.append(p);
  return;
}

namespace {
void
underresolved_event(nut::events::Event const & event)
{
  std::stringstream errstr;
  errstr << "Event not fully resolved: event type "
         << nut::events::event_name(event);
  throw std::runtime_error(errstr.str());
}

void
unknown_event(nut::events::Event const & event)
{
  std::stringstream errstr;
  errstr << "Unknown event: " << event;
  throw std::runtime_error(errstr.str());
}
}  // namespace

}  // namespace nut

#pragma clang diagnostic pop

#endif  // include guard

// End of file
