// Particle.hh
// T. M. Kelley
// Dec 21, 2010
// Header for Particle
// (c) Copyright 2010 LANSLLC all rights reserved.

#ifndef PARTICLE_H
#define PARTICLE_H

#include <ostream>
#include "types.hh"

namespace nut
{

    template <typename fpt, typename rngt>
    struct Particle
    {
        typedef fpt  fp_t;
        typedef rngt rng_t;

        fp_t x;
        fp_t omega;
        fp_t e;
        fp_t t;
        fp_t weight;
        cell_t  cell;
        rng_t   rng;
        Species species;
        bool alive;

        // ctor
        Particle( fp_t x_, 
                  fp_t omega_,
                  fp_t e_,
                  fp_t t_,
                  fp_t weight_,
                  cell_t  cell_,
                  rng_t   rng_,
                  Species s_)
            : 
            x    (x_),
            omega(omega_),
            e    (e_),
            t    (t_),
            weight (weight_),
            cell   (cell_),
            rng    (rng_),
            species(s_),
            alive  (true)
            {}
            

        void kill(){
            alive = false;
        }

        bool operator==(Particle const & p) const {
            bool passed = 
                x       == p.x &&
                omega   == p.omega &&
                e       == p.e &&
                t       == p.t &&
                weight  == p.weight &&
                cell    == p.cell &&
                rng     == p.rng && 
                species == p.species &&
                alive   == p.alive;
            return passed;
        }

    }; // Particle

    // template <typename fp_t,typename rng_t>
    // std::ostream &
    // operator<<(std::ostream & os, Particle<fp_t,rng_t> const & p)
    // {
    //     os << "x: "       << p.x
    //        << ", omega: " << p.omega
    //        << ", t: " << p.t
    //        << ", e: " << p.e
    //        << ", weight: " << p.weight
    //        << ", cell: "   << (p.cell - 1)
    //        << ", rng: "    << p.rng
    //        << ", species: " << p.species
    //        << ", alive: "   << p.alive
    //         ;
    //     return os;
    // }

} // nut::

#endif



// version
// $Id$

// End of file
