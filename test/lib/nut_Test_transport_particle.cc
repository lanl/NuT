//T. M. Kelley (c) 2011 LANS LLC

#include "nut_Test_transport_particle.hh"
#include "test_aux.hh"
#include "Log.hh"
#include "transport.hh"
#include "RNG.hh"
#include "Census.hh"
#include "Particle.hh"
#include "Mesh.hh"
#include "Assert.hh"
#include <vector>


namespace Nut_Test
{
    namespace transport_particle_tests
    {
        // target describes the code being tested
        char target[] = "transport_particle";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "transport 1 particle";

    }

    bool test_transport_particle()
    {
        using namespace transport_particle_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        // call additional tests here.

        return passed1;
    }

    namespace transport_particle_tests
    {
        using nut::Null_Log;
        using nut::cell_t;
        using nut::Species;
        using test_aux::check_same;
        using test_aux::check_same_verb;
        using test_aux::comp_verb;
        using test_aux::check_one_changed;
        using test_aux::check_two_changed;
        using std::cerr;
        using std::endl;

        typedef double                    fp_t;
        typedef nut::Buffer_RNG<fp_t>     rng_t;
        typedef nut::Particle<fp_t,rng_t> p_t;
        typedef nut::Tally<fp_t>          t_t;
        typedef nut::Census<p_t>          c_t;
        typedef nut::Sphere_1D<cell_t,nut::geom_t,nut::bdy_types::descriptor> mesh_t;

        // bits for std particle
        fp_t const rns[] = {0.30897681610609407,
                            0.92436920169905545,
                            0.21932404923057958};
        rng_t rng(rns,3);
        fp_t const x = 0.5, omega = 1.0, e = 5.0, t = 1.0, wt = 1.0;
        cell_t const cell = 1;
        Species const s(nut::nu_e);

        p_t make_std_particle()
        {
            p_t p(x,omega,e,t,wt,cell,rng,s);
            return p;
        }

        // generate a mesh that's reflective at the inner boundary,
        // vacuum at the outer boundary, and transmissive at all others
        struct gen_bdy_types
        {
            nut::bdy_types::descriptor operator()(){
                if(ctr++ == 0) return nut::bdy_types::descriptor::R;
                if(ctr == nbdy) return nut::bdy_types::descriptor::V;
                nut::Insist(ctr <= nbdy,"called gen_bdy_types too often");
                return nut::bdy_types::descriptor::T;
            }
            explicit gen_bdy_types(cell_t const nbdy_) : ctr(0),nbdy(nbdy_) {}
            cell_t ctr;
            cell_t const nbdy;
        }; // gen_bdy_types


        // generate a uniformly-spaced mesh, with spacing dx
        struct gen_bounds
        {
            fp_t operator()(){return dx * ctr++;}
            explicit gen_bounds(fp_t const dx_) : dx(dx_),ctr(0) {}
            fp_t const dx;
            cell_t ctr;
        }; // gen_bdy_types


        bool test_1()
        {
            bool passed(true);

            p_t p(make_std_particle());

            // mesh
            using nut::events::Event;

            using nut::decide_boundary_event;
            typedef mesh_t::vb vb;
            typedef mesh_t::vbd vbd;
            // generate uniform mesh
            cell_t n_cells(10);
            vb bounds(n_cells+1);
            vbd b_types(n_cells+1);
            std::generate(b_types.begin(),b_types.end(),gen_bdy_types(n_cells+1));
            std::generate(bounds.begin(),bounds.end(),gen_bounds(1.0));

            // reflects from innermost face
            mesh_t mesh(bounds, b_types);

            // opacity

            // tally
            nut::Tally<fp_t> tally(n_cells), ref(n_cells);

            // census
            c_t c,c_ref;


            return passed;
        }

        // define additional tests here.

    } // transport_particle_tests::

} // Nut_Test::



// version
// $Id$

// End of file
