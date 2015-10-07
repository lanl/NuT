// utilities.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "types.hh"
// #include <sstream>
#include <numeric>
#include <cmath>
// #include "Assert.hh"

/**!\file Useful operators. */

namespace nut
{
    template <typename fp_t>
    struct mult
    {
        fp_t operator()(fp_t const a) const {return a * t;}
        fp_t const t;
        explicit mult(fp_t const t_) : t(t_) {}
    }; // mult

    template <typename fp_t, typename int_t>
    struct int_mult
    {
        int_t operator()(fp_t const a) const
        {
            return int_t( std::floor(a * t + 0.5) );
        }
        fp_t const t;
        explicit int_mult(fp_t const t_) : t(t_) {}
    }; // mult

    template <typename fp_t>
    struct div_by
    {
        fp_t operator()(fp_t const a) const {return a / d;}
        fp_t const d;
        explicit div_by(fp_t const d_) : d(d_)
        {
            Require(d != fp_t(0),"div_by::ctor: divisor may not equal 0");
        }
    }; // div_by

    /*!\brief if divisor is not equal to 0, returns a / b, otherwise
     * returns 0. */
    template <typename fp_t>
    fp_t div_if_ne0(fp_t const a,size_t const b)
    {
        return b != 0 ? fp_t(a / b) : fp_t (0);
    }

    // sum contents of an STL container
    template <typename ContT>
    typename ContT::value_type
    sum(ContT const & v)
    {
        return std::accumulate(v.begin(),v.end(),
                               typename ContT::value_type(0));
    }

} // nut::


#endif // include_guard

// version
// $Id$

// End of file
