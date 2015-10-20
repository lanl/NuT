//T. M. Kelley (c) 2011 LANS LLC

#include "Particle.hh"
#include "RNG.hh"
#include "nut_Test_Particle.hh"
#include "test_aux.hh"


namespace Nut_Test
{
    namespace particle_tests
    {
        // target describes the code being tested
        char target[] = "Particle";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "instantiation/initialization";
        bool test_2();
        char aspect2[] = "assignment";

    }

    bool test_Particle()
    {
        using namespace particle_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        // bool passed2 = test( target, aspect2, test_2);
        // call additional tests here.

        return passed1
         // && passed2
            ;
    }

    namespace particle_tests
    {

        // bool test_2()
        // {
        //     typedef nut::Philox4x32_RNG rng_t;
        //     typedef double fp_t;
        //     typedef nut::Particle<fp_t, rng_t> p_t;

        //     rng_t r1({0,0,0,6},{7,7,8,9});
        //     rng_t r2({59,60,61,62},{801,802,803,804});

        //     // p_t p1();
        //     return false;
        // } // test_2


        bool test_1()
        {
            // for a very simple test, we can use a totally bogus RNG
            typedef uint32_t rng_t;
            typedef double fp_t;

            typedef nut::Particle<fp_t, rng_t, 1> part_t;

            bool passed(true);

            fp_t x(1.0);
            fp_t omega(1.0);
            fp_t e(1.0);
            fp_t t(1.0);
            fp_t wt(1.0);
            nut::cell_t cell(1);
            nut::Species s(nut::nu_e);

            rng_t rng(42);

            part_t particle({x},{omega},e,t,wt,cell,rng,s);

            return passed;
        }

        // define additional tests here.

    } // particle_tests::

} // Nut_Test::



// version
// $Id$

// End of file
