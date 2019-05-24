// Events.cc
// May 24, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#include "events.hh"
#include <string>

namespace nut {

namespace events {

std::string
event_name(Event const & e)
{
  std::string name;
  switch(e) {
    case collision: name = "Collision"; break;
    case nucleon_abs: name = "Nucleon_abs"; break;
    case nucleon_elastic_scatter: name = "Nucleon_elastic_scatter"; break;
    case electron_scatter: name = "Electron_scatter"; break;
    case positron_scatter: name = "Positron_scatter"; break;
    case nu_e_annhilation: name = "Nu_e_annhilation"; break;
    case nu_x_annhilation:
      name = "Nu_x_annhilation";
      break;
      ;  // # compute events
    case boundary: name = "boundary"; break;
    case cell_boundary: name = "Cell_boundary"; break;
    case escape: name = "Escape"; break;
    case reflect: name = "Reflect"; break;
    case step_end: name = "Step_end"; break;
    case weight_cutoff:
      name = "Weight_cutoff";
      break;
      ;  // testing only
    case null: name = "Null"; break;
  };
  return name;
}  // event_name

}  // namespace events

}  // namespace nut

// End of file
