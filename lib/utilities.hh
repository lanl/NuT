// utilities.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "types.hh"
#include <numeric>
#include <algorithm>
#include <cmath>
#include "Assert.hh"

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

    /** Add each element of v1 to corr. element v2. */
    template <typename T>
    inline
    void merge_vectors(std::vector<T> const & v1, std::vector<T> & v2)
    {
        std::transform(v1.begin(),v1.end(),v2.begin(),v2.begin(),std::plus<T>());
    }

    /** Append the elements of v1 to v2. */
    template <typename T>
    inline
    void append_vector(std::vector<T> const & v1, std::vector<T> & v2)
    {
        v2.insert(v2.end(),v1.begin(),v1.end());
    }


    /** Compare arrays element-wise*/
    template <size_t sz>
    bool arrays_eq(vec_t<sz> const & a1, vec_t<sz> const & a2)
    {
        return std::equal(a1.begin(),a1.end(),a2.begin());
    }


} // nut::


#endif // include_guard

// version
// $Id$

// End of file
