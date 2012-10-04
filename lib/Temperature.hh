// Temperature.hh
// T. M. Kelley
// Jan. 10, 2010
// Header for Temperature
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "Assert.hh"
#include <vector>
#include <iostream>

namespace nut
{
    /*! \brief arrays to keep track of temperature in each mesh cell
     * \tparam <fp_t> {floating point type} 
     */
    template <typename fp_t>
    struct Temperature
    {
        typedef std::vector<fp_t> vf;

        vf T_p;
        vf T_e_minus;
        vf T_e_plus;
        
        /*!\brief structure of arrays to store temperature data.
         *
         * Temperatures are stored for the various material types. 
         * The temperature vectors must all be the same
         * size (one element per mesh cell).  
         */
        Temperature(vf const & T_p_,
                    vf const & T_e_minus_,
                    vf const & T_e_plus_
            )
            : T_p(T_p_),
              T_e_minus(T_e_minus_),
              T_e_plus(T_e_plus_)
            {
                this->check_sizes();
                return;
            }

        /*!\brief copy ctor */
        explicit Temperature ( Temperature const & t)
            : T_p(t.T_p),
              T_e_minus(t.T_e_minus),
              T_e_plus(t.T_e_plus)
            {}

        /*!\brief construct empty Temperature of given size */
        explicit Temperature(size_t const n)
            : T_p       (n, fp_t(0)),
              T_e_minus (n, fp_t(0)),
              T_e_plus  (n, fp_t(0))
            {}

        /*!\brief check sizes of arrays for consistency at construction */
        void check_sizes(){
            using nut::Equal;
            size_t n_cells = T_p.size();
            Equal(T_e_minus.size(), n_cells,"T_e_minus","n_cells");
            Equal(T_e_plus.size(),  n_cells,"T_e_plus","n_cells");
            return;
        } // check_sizes

        size_t size() const {return T_p.size();}

    }; // Temperature

    template <typename fp_t>
    struct temp_t
    {
        fp_t t_p, t_e_m, t_e_p;
    };


} // nut::

#endif



// version
// $Id$

// End of file
