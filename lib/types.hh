// types.hh
// T. M. Kelley
// Jan 03, 2011
// Header for types
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef TYPES_H
#define TYPES_H


#include <stdint.h>
#include <stdexcept>
#include <utility>  // std::pair
#include <vector>
#include <string>
#include <iterator>
#include "RNG.hh"

namespace nut
{
    namespace bdy_types
    {
        enum descriptor
        {
            V, /** Vacuum */
            R, /** Reflective */
            T  /** Transmissive */
        };
    }

    namespace events
    {
        enum Event
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
            // testing only
            null          = 9999
        };
    } // events::

    typedef double geom_t;   /*^ type for geometry calculations */

    typedef std::vector<geom_t> vg;

    typedef std::pair<events::Event, geom_t> event_n_dist;
    
    typedef uint32_t cell_t; /*^ cell index. Use cell_t(-1) for null cell. */

    typedef uint32_t id_t;   /*^ particle index */

    typedef std::vector<id_t> vid;

    typedef uint32_t group_t;

    typedef uint32_t cntr_t;

    typedef std::invalid_argument arg_error;

    typedef std::ostream_iterator<size_t> sz_o_it;

    typedef std::ostream_iterator<geom_t> ge_o_it;
    
    typedef nut::Philox4x32_RNG         rng_t;

    typedef rng_t::ctr_t                ctr_t;

    typedef rng_t::key_t                key_t;

    typedef rng_t::seed_t               seed_t;

    /* Neutrino species */
    enum Species
    {
          nu_e
        , nu_e_bar
        , nu_x
        , nu_mu
        , nu_mu_bar
        , nu_tau
        , nu_tau_bar
    }; // Species; 
    
    std::string species_name( Species const s);
    
    Species species_type( std::string const & s);
    
    seed_t species_seed(Species const s);

    namespace events
    {
        namespace 
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

        } // anonymous::

    } // events::
    
} // nut::

#endif



// version
// $Id$

// End of file
