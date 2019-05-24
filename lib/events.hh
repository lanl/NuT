// Events.hh
// May 24, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#include "Assert.hh"
#include <utility>  // make_pair

namespace nut {

namespace events {

constexpr uint32_t face_bits{10};

constexpr uint32_t max_num_faces{0xFF};

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

  // Reserve bytes 5&6 for encoding face
  face_start = (1 << face_bits),
  face_mask = (max_num_faces << face_bits),

  // testing only
  null = 0xFFFF'FFFF
};

/**\brief Encode a face id into an Event.
 * \tparam Face_id: the face id type, must be integral.
 * \return New Event, with Face_id encode.
 */
template <typename Face_id>
Event
encode_face(Event const & e, Face_id const f)
{
  static_assert(std::is_integral<Face_id>::value,
                "Face_id must be integral type");
  Require(f >= 0 && f < max_num_faces, "encode_face: invalid face");
  uint32_t const f_in_place{f << face_bits};
  uint32_t const new_e = f_in_place | static_cast<uint32_t>(e);
  return static_cast<Event>(new_e);
}  // encode_face

/**\brief Decode a face id from an Event.
 * \tparam Face_id: the face id type, must be integral.
 * \return pair<Event,Face_id> (event w/out face_id, and the face_id)
 */
template <typename Face_id>
std::pair<Event, Face_id>
decode_face(Event const & e)
{
  static_assert(std::is_integral<Face_id>::value,
                "Face_id must be integral type");
  // project the face id out, and shift it back to normal
  uint32_t const f_in_place{e & face_mask};
  uint32_t f{f_in_place >> face_bits};
  LessThan(f, max_num_faces, "decode_face: face_id");
  uint32_t const new_eu = e & ~face_mask;
  Event new_e{static_cast<Event>(new_eu)};
  auto p = std::make_pair<Event, Face_id>(std::move(new_e), std::move(f));
  return p;
}  // encode_face

std::string
event_name(Event const & e);

}  // namespace events

using event_n_dist = std::pair<events::Event, geom_t>;

}  // namespace nut

// End of file
