// Event.cc
// T. M. Kelley
// Oct 26, 2015
// (c) Copyright 2015 LANSLLC, all rights reserved

#include "Event.hh"
#include <string>

namespace nut
{

    std::string event_name(Event const & e)
    {
        std::string name;
        switch(e)
        {
        case Event::collision:
            name = "Collision";
            break;
        case Event::nucleon_abs:
            name = "Nucleon_abs";
            break;
        case Event::nucleon_elastic_scatter:
            name = "Nucleon_elastic_scatter";
            break;
        case Event::electron_scatter:
            name = "Electron_scatter";
            break;
        case Event::positron_scatter:
            name = "Positron_scatter";
            break;
        case Event::nu_e_annhilation:
            name = "Nu_e_annhilation";
            break;
        case Event::nu_x_annhilation:
            name = "Nu_x_annhilation";
            break;;            // # compute events
        case Event::boundary:
            name = "boundary";
            break;
        case Event::cell_low_x_boundary:
            name = "Cell_low_x_boundary";
            break;
        case Event::cell_high_x_boundary:
            name = "Cell_high_x_boundary";
            break;
        case Event::escape:
            name = "Escape";
            break;
        case Event::reflect:
            name = "Reflect";
            break;
        case Event::step_end:
            name = "Step_end";
            break;
        case Event::weight_cutoff:
            name = "Weight_cutoff";
            break;;            // testing only
        case Event::undetermined:
            name = "undetermined";
            break;
        case Event::null:
            name = "Null";
            break;
        };
        return name;
    } // event_name

} // nut::


// End of file
