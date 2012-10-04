// expect.cc
// T. M. Kelley
// Jan 11, 2011
// Implementation, expect
// (c) Copyright 2011 LANSLLC all rights reserved.

#include "expect.hh"

namespace test_aux
{
    // explicit specializations for floating point.
    template <>
    void expect_fail<double>(double const v, double const r, std::string const & name)
    {
        std::streamsize const old_prec = outstr.precision();
        outstr << std::setprecision(15);
        outstr << name << " = " << v << ", expected " << r << std::endl;
        outstr << std::setprecision(old_prec);
    }

    template <>
    void expect_fail<float>(float const v, float const r, std::string const & name)
    {
        std::streamsize const old_prec = outstr.precision();
        outstr << std::setprecision(7);
        outstr << name << " = " << v << ", expected " << r << std::endl;
        outstr << std::setprecision(old_prec);
    }

} // test_aux::


// version
// $Id$

// End of file
