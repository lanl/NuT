// Census.hh
// T. M. Kelley
// Jan. 10, 2010
// Header for Census
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef CENSUS_H
#define CENSUS_H

#include "types.hh"
#include <vector>

namespace nut
{
    /*!\brief arrays of each type of particle, banked for next time step */
    template <typename p_t>
    struct Census
    {
        typedef p_t particle_t;
        typedef std::vector<particle_t> vp;
        vp nu_es;
        vp nu_e_bars;
        vp nu_mus;
        vp nu_mu_bars;
        vp nu_taus;
        vp nu_tau_bars;

        void append(p_t const & particle){
            switch(particle.species)
            {
            case nu_e:
                nu_es.push_back(particle);
                break;
            case nu_e_bar:
                nu_e_bars.push_back(particle);
                break;
            case nu_mu:
                nu_mus.push_back(particle);
                break;
            case nu_mu_bar:
                nu_mu_bars.push_back(particle);
                break;
            case nu_tau:
                nu_taus.push_back(particle);
                break;
            case nu_tau_bar:
                nu_tau_bars.push_back(particle);
                break;
            default: // should never get here
                std::cerr << "Census::append unknown particle species: " 
                          << particle.species << std::endl;
            } // switch
            return;
        } // append
        
    }; // Census

} // nut::

#endif



// version
// $Id$

// End of file
