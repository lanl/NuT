// application.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef APPLY_EVENT_HH
#define APPLY_EVENT_HH

/**\file apply events to particles */
#include "Planck.hh"
#include "Tally.hh"
#include "constants.hh"
#include "events.hh"
#include "lorentz.hh"
#include "types.hh"
#include <iomanip>
#include <sstream>
#include <stdexcept>

/* unknown_event and underresolved_event are not used in all translation units,
 * leading to unexciting warnings. Suppress those. */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"

namespace nut {
namespace {
void
underresolved_event(nut::events::Event const & event);
void
unknown_event(nut::events::Event const & event);
}  // namespace

template <typename p_t,
          typename tally_t,
          typename Mesh_T,
          typename vector_t = typename Mesh_T::Vector>
void
apply_nucleon_elastic_scatter(p_t & p, tally_t & t, vector_t const & v);

template <typename p_t,
          typename tally_t,
          typename Mesh_T,
          typename vector_t = typename Mesh_T::Vector>
void
apply_lepton_scatter(p_t & p,
                     tally_t & t,
                     vector_t const & velocity,
                     typename tally_t::FP_T const e_lep);

template <typename Mesh_T, typename Particle_T, typename fp_t>
void
stream_particle(Particle_T & p, fp_t const d, Mesh_T const & mesh)
{
  typename Mesh_T::Ray newcoord = mesh.advance_point(p.x, p.omega, d);
  p.x = newcoord.position();
  p.omega = newcoord.direction();
  p.t = p.t - d / c;
  return;
}  // stream_particle

template <typename ParticleT,
          typename MeshT,
          typename OpacityT,
          typename CensusT,
          typename TallyT>
void
apply_event(ParticleT & p,
            event_data<typename MeshT::face_handle_t> const & event_data,
            MeshT const & mesh,
            OpacityT const & opacity,
            TallyT & tally,
            CensusT & census,
            typename TallyT::fp_t const alpha = typename TallyT::fp_t(2.0))
{
  using namespace events;
  using fp_t = typename TallyT::fp_t;
  using vector_t = typename MeshT::Vector;

  Event const & event = events::get_event(event_data);
  geom_t distance = events::get_distance(event_data);
  typename MeshT::face_handle_t face = events::get_face(event_data);
  // typename MeshT::Vector face_normal = mesh.get_normal(p.cell, face);

  stream_particle(p, distance, mesh);

  tally.accum_pl(distance);

  vector_t const & velocity = opacity.velocity(p.cell);
  typename OpacityT::Cell_Data_T const & cell_data = opacity.cell_data(p.cell);

  /*

  template <typename p_t, typename tally_t, typename vector_t, typename Mesh_T>
  void
  apply_nucleon_elastic_scatter(p_t & p, tally_t & t, vector_t const & v)

  */
  switch(event) {
    case nucleon_abs: apply_nucleon_abs(p, tally); break;
    case nucleon_elastic_scatter:
      apply_nucleon_elastic_scatter<ParticleT, TallyT, MeshT>(p, tally,
                                                              velocity);
      break;
    case electron_scatter: {
      fp_t const ebar = cell_data.T_e_minus;
      fp_t const e_e = ebar;
      // fp_t const e_e  = gen_power_law_energy_alpha2(ebar,p.rng);
      apply_lepton_scatter<ParticleT, TallyT, MeshT>(p, tally, velocity, e_e);
    } break;
    case positron_scatter: {
      fp_t const ebar = cell_data.T_e_plus;
      fp_t const e_e = ebar;
      // fp_t const e_e  = gen_power_law_energy_alpha2(ebar,p.rng);
      apply_lepton_scatter<ParticleT, TallyT, MeshT>(p, tally, velocity, e_e);
    } break;
    // case nu_e_annhilation:
    //     break;
    // case nu_x_annhilation:
    //     break;
    case cell_boundary: apply_cell_boundary(mesh, face, p, tally); break;
    // case cell_boundary: apply_hi_x_boundary(p, tally); break;
    case escape: apply_escape(p, tally); break;
    // case reflect: apply_reflect(p, tally, face_normal); break;
    case reflect: apply_reflect(p, tally, mesh, face); break;
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

template <typename p_t, typename tally_t, typename Mesh_T, typename vector_t>
void
apply_nucleon_elastic_scatter(p_t & p, tally_t & t, vector_t const & v)
{
  using vector4_t = typename Mesh_T::vector4_t;
  // size_t const dim(tally_t::dim);
  cell_t const cell = p.cell;
  // Vector const v = vel.v(cell);
  geom_t const eli = p.e;
  vector_t const oli = p.omega;
  // LT to comoving frame (need init comoving e to compute
  // final lab e).

  vector4_t eno_cmi = Mesh_T::LT_to_comoving(v, eli, oli);
  geom_t const eci = eno_cmi[0];
  // scatter: sample a new direction
  vector_t ocf = Mesh_T::sample_direction_isotropic(p.rng);

  // LT comoving -> lab
  geom_t elf;
  vector_t olf;
  try {
    vector4_t const eno_lf = Mesh_T::LT_to_lab(v, eci, ocf);
    elf = eno_lf[0];
    for(size_t i = 0; i < olf.size(); ++i) { olf[i] = eno_lf[i + 1]; }
  } catch(std::exception & e) {
    printf("%s:%i Caught exception \"%s\"\n", __FUNCTION__, __LINE__, e.what());
    printf("%s:%i ocf = %f\n", __FUNCTION__, __LINE__, ocf[0]);
    throw(e);
  }

  // tally
  t.count_nucleon_elastic_scatter(p.cell);
  t.deposit_inelastic_scat(cell, eli, elf, oli, olf, p.weight, p.species);
  // update particle
  p.omega = olf;
  p.e = elf;

  return;
}  // namespace nut

template <typename p_t, typename tally_t, typename Mesh_T, typename vector_t>
void
apply_lepton_scatter(p_t & p,
                     tally_t & t,
                     vector_t const & velocity,
                     typename tally_t::FP_T const e_lep)
{
  using vector4_t = typename Mesh_T::vector4_t;
  using fp_t = typename tally_t::FP_T;
  cell_t const cell = p.cell;
  // Vector const v = op.velocity(cell);
  geom_t const eli = p.e;
  vector_t const oli = p.omega;

  // LT to comoving frame (need init comoving e to compute
  // final lab e).
  vector4_t eno_cmi = Mesh_T::LT_to_comoving(velocity, eli, oli);
  geom_t const eci = eno_cmi[0];
  // scatter
  vector_t ocf = Mesh_T::sample_direction_isotropic(p.rng);
  fp_t const de = (eci - e_lep) / 4.;
  fp_t const ecf = eci - de;
  // LT comoving -> lab
  vector4_t const eno_lf = Mesh_T::LT_to_lab(velocity, ecf, ocf);
  geom_t const elf = eno_lf[0];
  vector_t olf{};
  for(size_t i = 0; i < olf.size(); ++i) { olf[i] = eno_lf[i + 1]; }
  // tally
  t.count_lepton_scatter(p.cell, p.species);
  t.deposit_inelastic_scat(cell, eli, elf, oli, olf, p.weight, p.species);
  // update
  p.omega = olf;
  p.e = elf;
  return;
}  // apply_lepton_scatter

template <typename p_t, typename tally_t, typename mesh_t, typename face_t>
void
apply_cell_boundary(mesh_t const & m, face_t const & face, p_t & p, tally_t & t)
{
  t.count_cell_bdy(p.cell);
  typename mesh_t::Cell mcell{p.cell};
  if constexpr(std::is_class<decltype(
                   m.cell_across(mcell, face, p.x))>::value) {
    auto new_mcell = m.cell_across(mcell, face, p.x);
    p.cell = new_mcell.as_id();
  }
  else {
    p.cell = m.cell_across(mcell, face, p.x);
  }
  // Require(new_mcell != m.null_cell(), "invalid cell from cell boundary");
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

template <typename p_t, typename tally_t, typename mesh_t>
void
apply_reflect(p_t & p,
              tally_t & t,
              mesh_t const & m, /*vector_t const & face_normal*/
              typename mesh_t::face_handle_t const face)
{
  using vector_t = typename mesh_t::Vector;
  t.count_reflect(p.cell);
  vector_t unused{0.0};
  vector_t reflected = m.reflect(face, p.omega, unused);
  p.omega = reflected;
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
