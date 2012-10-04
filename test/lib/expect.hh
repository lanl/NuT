// expect.hh
// T. M. Kelley
// Jan 11, 2011
// Header for expect
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef EXPECT_H
#define EXPECT_H


#include "soft_equiv.hh"
#include <iostream>
#include <iomanip>


namespace test_aux
{
    // handy functor harness for running soft_equiv in std::equal
    template <typename fp_t>
    struct soft_eq_bound_tol
    {
        explicit soft_eq_bound_tol(fp_t const tol) : m_tol(tol){}
        bool operator()(fp_t const v, fp_t const r){
            return nut::soft_equiv(v,r,m_tol);
        }
        fp_t m_tol;
    };

    static 
    std::ostream & outstr(std::cout);

    template <typename t>
    void expect_fail(t const v, t const r, std::string const & name)
    {
        outstr << name << " = " << v << ", expected " << r << std::endl;
    }

    /*!\brief check equality within tolerance, print message if fails. */
    template <typename fp_t>
    bool soft_expect(fp_t const val, fp_t const ref, std::string const & name,
                     fp_t const tol = fp_t(1e-15))
    {
        bool passed = nut::soft_equiv(val,ref,tol);
        if(!passed)
        {
            // fp_t const prec_f = std::log10(tol);
            // std::streamsize const prec = 1+static_cast<std::streamsize>(prec_f);
            expect_fail(val, ref,name);
        }
        return passed;
    }

    template <typename int_t>
    bool
    expect(int_t const val, int_t const ref, std::string const & name)
    {
        bool passed = val == ref;
        if(!passed)
        {
            expect_fail(val,ref,name);
        }
        return passed;
    }

    // see expect.cc for explicit specializations of expect_fail

}// test_aux::

#endif



// version
// $Id$

// End of file
