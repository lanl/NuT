// lorentz.hh
// T. M. Kelley
// May 02, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef LORENTZ_HH
#define LORENTZ_HH

#include "constants.hh"
#include <cmath>

namespace nut
{
    typedef std::pair<geom_t,geom_t> EandOmega;

    // template <typename fp_t>
    geom_t 
    inline 
    gamma(geom_t const v)
    {
        return 1.0/std::sqrt(1 - v*v/(c*c));
    }

    EandOmega
    inline
    LT_to_comoving_sphere1D(geom_t const v_lab, 
                            geom_t const e_lab, 
                            geom_t const omega_lab)
    {
        geom_t const voc       = v_lab/c;
        geom_t const onemovoc  = (1 - omega_lab * voc);
        geom_t const e_com     = e_lab * gamma(v_lab) * onemovoc;
        geom_t const omega_com = (omega_lab - voc) / onemovoc;
        return EandOmega(e_com,omega_com);
    }
    
    // Yes, you should pass v_lab here--it's the material speed in the lab frame!
    EandOmega 
    inline
    LT_to_lab_sphere1D(geom_t const v_lab, 
                       geom_t const e_com, 
                       geom_t const omega_com)
    {
        geom_t const voc       = v_lab/c;
        geom_t const onepovoc  = 1 + omega_com * voc;
        geom_t const e_lab     = e_com * gamma(v_lab) * onepovoc;
        geom_t const omega_lab = (voc + omega_com) / onepovoc;
        return EandOmega(e_lab,omega_lab);
    }

} // nut::


#endif // LORENTZ_HH

// version
// $Id$

// End of file
