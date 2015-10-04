// types.cc
// T. M. Kelley
// May 23, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#include "types.hh"
#include <string>

namespace nut
{

    namespace events
    {

        std::string event_name(Event const & e)
        {
            std::string name;
            switch(e)
            {
            case collision:
                name = "Collision";
                break;
            case nucleon_abs:
                name = "Nucleon_abs";
                break;
            case nucleon_elastic_scatter:
                name = "Nucleon_elastic_scatter";
                break;
            case electron_scatter:
                name = "Electron_scatter";
                break;
            case positron_scatter:
                name = "Positron_scatter";
                break;
            case nu_e_annhilation:
                name = "Nu_e_annhilation";
                break;
            case nu_x_annhilation:
                name = "Nu_x_annhilation";
                break;;            // # compute events
            case boundary:
                name = "boundary";
                break;
            case cell_low_x_boundary:
                name = "Cell_low_x_boundary";
                break;
            case cell_high_x_boundary:
                name = "Cell_high_x_boundary";
                break;
            case escape:
                name = "Escape";
                break;
            case reflect:
                name = "Reflect";
                break;
            case step_end:
                name = "Step_end";
                break;
            case weight_cutoff:
                name = "Weight_cutoff";
                break;;            // testing only
            case null:
                name = "Null";
                break;
            };
            return name;
        } // event_name

    } // events::


    seed_t species_seed( Species const s)
    {
        seed_t n(-1);
        switch(s)
        {
        case nu_e:       n = 0; break;
        case nu_e_bar:   n = 1; break;
        case nu_x:       n = 6; break;
        case nu_mu:      n = 2; break;
        case nu_mu_bar:  n = 3; break;
        case nu_tau:     n = 4; break;
        case nu_tau_bar: n = 5; break;
        }
        return n;
    }


    std::string species_name( Species const s)
    {
        std::string n;
        switch(s)
        {
        case nu_e:       n = "nu_e";       break;
        case nu_e_bar:   n = "nu_e_bar";   break;
        case nu_x:       n = "nu_x";       break;
        case nu_mu:      n = "nu_mu";      break;
        case nu_mu_bar:  n = "nu_mu_bar";  break;
        case nu_tau:     n = "nu_tau";     break;
        case nu_tau_bar: n = "nu_tau_bar"; break;
        }
        return n;
    }


    Species species_type( std::string const & n)
    {
        Species s;
        if(n == "nu_e")           s = nu_e;
        else if(n == "nu_e_bar")  s = nu_e_bar;
        else if(n == "nu_x")      s = nu_x;
        else if(n == "nu_mu")     s = nu_mu;
        else if(n == "nu_mu_bar") s = nu_mu_bar;
        else if(n == "nu_tau")    s = nu_tau;
        else                      s = nu_tau_bar;
        return s;
    }

} // nut::

// version
// $Id$

// End of file
