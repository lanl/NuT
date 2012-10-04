//T. M. Kelley (c) 2011 LANS LLC

#include "nut_Test_species_name.hh"
#include "test_aux.hh"
#include "types.hh"
#include "copyright.hh"

namespace Nut_Test
{
    namespace species_name_tests
    {
        // target describes the code being tested
        char target[] = "species_name";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "nu_e";

    }

    bool test_species_name()
    {
        using namespace species_name_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        // call additional tests here.

        return passed1;
    }

    namespace species_name_tests
    {
        bool test_1()
        {
            bool passed(true);

            std::cout << nut::copyright() << std::endl;

            std::string n = nut::species_name(nut::nu_e);

            passed = n == "nu_e";

            return passed;
        }

        // define additional tests here.

    } // species_name_tests::

} // Nut_Test::



// version
// $Id$

// End of file
