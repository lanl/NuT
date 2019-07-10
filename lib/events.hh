// Events.hh
// May 24, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#include "Assert.hh"
#include "types.hh"
#include <utility>  // make_pair

namespace nut {

namespace events {

enum Event : uint32_t {
  // physics events
  collision = 0,
  nucleon_abs = 1,
  nucleon_elastic_scatter = 2,
  electron_scatter = 3,
  positron_scatter = 4,
  nu_e_annhilation = 5,
  nu_x_annhilation = 6,
  // # compute events
  boundary = 100,
  cell_boundary = 101,
  escape = 110,
  reflect = 111,
  step_end = 115,
  weight_cutoff = 120,

  // testing only
  null = 0xFFFF'FFFF
};

std::string
event_name(Event const & e);

}  // namespace events

template <typename Face_T>
using event_data = std::tuple<events::Event, geom_t, Face_T>;

namespace events {
template <typename Face_T>
Event &
get_event(event_data<Face_T> & e)
{
  return std::get<0>(e);
}
template <typename Face_T>
Event const &
get_event(event_data<Face_T> const & e)
{
  return std::get<0>(e);
}
template <typename Face_T>
geom_t
get_distance(event_data<Face_T> const & e)
{
  return std::get<1>(e);
}
template <typename Face_T>
Face_T
get_face(event_data<Face_T> const & e)
{
  return std::get<2>(e);
}
}  // namespace events

}  // namespace nut

// End of file
