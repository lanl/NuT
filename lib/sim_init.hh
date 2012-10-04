// sim_init.hh
// T. M. Kelley
// Jun 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#include "Density.hh"
#include "Temperature.hh"
#include "Opacity.hh"
#include "Velocity.hh"
#include "Mesh.hh"

namespace nut
{
    /* Gather up opacity, luminosity, velocity into one convenient struct. */
    template <typename fp_t>
    struct Mat_State
    {
        typedef Density<fp_t> D_t;
        typedef Temperature<fp_t> T_t;
        typedef Opacity<fp_t> op_t;
        typedef Velocity<fp_t> vel_t;
        typedef Luminosity<fp_t> lum_t;
        typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;
        
        
    }; 

} // nut::


// version
// $Id$

// End of file
