//T. M. Kelley (c) 2011 LANS LLC

#include "Temperature.hh"
#include "nut_Test_Temperature.hh"
#include "test_aux.hh"


namespace Nut_Test
{
    namespace Temperature_tests
    {
        // target describes the code being tested
        char target[] = "Temperature";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "inst. & init.";

    }

    bool test_Temperature()
    {
        using namespace Temperature_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        // call additional tests here.

        return passed1;
    }

    namespace Temperature_tests
    {
        bool test_1()
        {
            bool passed(true);

            typedef float fp_t;
            typedef std::vector<fp_t> v_t;

            size_t n_cells(100);

            v_t T_p(n_cells);
            v_t T_e_minus(n_cells);
            v_t T_e_plus(n_cells);
            
            nut::Temperature<fp_t> temperature(T_p,
                                               T_e_minus,
                                               T_e_plus
                );
            return passed;
        }

        // define additional tests here.

    } // Temperature_tests::

} // Nut_Test::



// version
// $Id$

// End of file
