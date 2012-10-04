//T. M. Kelley (c) 2011 LANS LLC

#include "RNG.hh"
#include "nut_Test_RNG.hh"
#include "test_aux.hh"


namespace Nut_Test
{
    namespace rNG_tests
    {
        // target describes the code being tested
        char target[] = "random (and not-so-random) number generators";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "inst. & init. LCG_RNG";

    }

    bool test_RNG()
    {
        using namespace rNG_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        // call additional tests here.

        return passed1;
    }

    namespace rNG_tests
    {
        bool test_1()
        {
            return false;
        }

        // define additional tests here.

    } // rNG_tests::

} // Nut_Test::



// version
// $Id$

// End of file
