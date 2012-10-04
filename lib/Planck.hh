// Planck.hh
// T. M. Kelley
// Jan 12, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef PLANCK_HH
#define PLANCK_HH

#include <cmath>

namespace nut
{
    /*\tparam rng_t: random number generator, must implement 
     *               fp_t random(). */

    template <typename fp_t>
    inline 
    fp_t pwr_law_reject(fp_t const alpha, fp_t const e_a,
                          fp_t const x)
    {
        return std::pow(x,alpha) * std::exp(-alpha*x)/e_a;
    } // pwr_law_reject


    template <typename rng_t, typename fp_t>
    inline
    fp_t gen_power_law_energy(fp_t const alpha,fp_t const ebar,
                              rng_t & rng)
    {
        fp_t energy = 0.0;

        fp_t const e_a = std::exp(-alpha);
        bool found = false;
        while(!found)
        {
            fp_t const xsd      = rng.random();
            fp_t const x        = -1.0 * std::log(xsd);
            fp_t const selector = rng.random();
            fp_t const rejector = pwr_law_reject(alpha,e_a,x);
            if(selector <= rejector)
            {
                energy = x * ebar;
                found  = true;
            }
        }
        return energy;
    } // gen_power_law_energy


    /* same functions, specialized to alpha = 2: */

    template <typename fp_t>
    inline 
    fp_t pwr_law_reject_alpha2(fp_t const x)
    {
        fp_t const exp2 = fp_t(0.1353352832366127); // e^{-2}
        return x*x*std::exp(-2*x)/exp2;
    }

    template <typename rng_t, typename fp_t>
    inline
    fp_t gen_power_law_energy_alpha2(fp_t const ebar, rng_t & rng)
    {
        fp_t energy = 0.0;

        bool found = false;
        while(!found)
        {
            fp_t const x        = -1.0 * std::log(rng.random());
            fp_t const selector = rng.random();
            fp_t const rejector = pwr_law_reject_alpha2<fp_t>(x);
            if(selector <= rejector)
            {
                energy = x * ebar;
                found  = true;
            }
        }
        return energy;
    } // gen_power_law_energy

} // nut::


#endif // PLANCK_HH

// version
// $Id$

// End of file
