// Event.hh
// T. M. Kelley
// Oct 26, 2015
// (c) Copyright 2015 LANSLLC, all rights reserved


#ifndef EVENT_HH
#define EVENT_HH

#include "types.hh"

namespace nut
{
    enum struct Event
    {
        // physics events
        collision               = 0,
        nucleon_abs             = 1,
        nucleon_elastic_scatter = 2,
        electron_scatter        = 3,
        positron_scatter        = 4,
        nu_e_annhilation        = 5,
        nu_x_annhilation        = 6,
        // # compute events
        boundary             = 100,
        cell_low_x_boundary  = 101,
        cell_high_x_boundary = 102,
        escape        = 110,
        reflect       = 111,
        step_end      = 115,
        weight_cutoff = 120,
        // bad times
        undetermined  = 999,  // Could not make a decision--error!
        // testing only
        null          = 9999
    }; // events::

    typedef std::pair<Event, geom_t> event_n_dist;

    std::string event_name(Event const & e);

} // nut::

#endif // include guard


// End of file
