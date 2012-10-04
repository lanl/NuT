//T. M. Kelley (c) 2011 LANS LLC

#include "Density.hh"
#include "nut_Test_Density.hh"
#include "test_aux.hh"


namespace Nut_Test
{
    namespace Density_tests
    {
        // target describes the code being tested
        char target[] = "Density";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "inst. & init.";

    }

    bool test_Density()
    {
        using namespace Density_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        // call additional tests here.

        return passed1;
    }

    namespace Density_tests
    {
        bool test_1()
        {
            bool passed(true);

            typedef float fp_t;
            typedef std::vector<fp_t> v_t;

            size_t n_cells(100);

            v_t rho_nu_e(n_cells);
            v_t rho_nu_e_bar(n_cells);
            v_t rho_nu_x(n_cells);
            v_t rho_nu_x_bar(n_cells);
            v_t rho_p(n_cells);
            v_t rho_n(n_cells);
            v_t rho_e_minus(n_cells);
            v_t rho_e_plus(n_cells);
            v_t rho_A(n_cells);
            v_t y_e(n_cells);
            v_t abar(n_cells);
            
            nut::Density<fp_t> density(rho_p,
                                       rho_e_minus,
                                       rho_e_plus,
                                       rho_A,
                                       y_e,
                                       abar);
            return passed;
        }

        // define additional tests here.

    } // Density_tests::

} // Nut_Test::



// version
// $Id$

// End of file
