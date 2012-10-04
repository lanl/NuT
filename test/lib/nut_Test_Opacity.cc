//T. M. Kelley (c) 2011 LANS LLC

#define USING_HFBC_SIGMAS
#include "Opacity.hh"
#undef  USING_HFBC_SIGMAS

#include "Density.hh"
#include "Temperature.hh"
#include "nut_Test_Opacity.hh"
#include "test_aux.hh"


namespace Nut_Test
{
    namespace Opacity_tests
    {
        // target describes the code being tested
        char target[] = "Opacity";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "inst. & init.";
        bool test_2();
        char aspect2[] = "sigma<nu_e_e_minus>";

    }

    bool test_Opacity()
    {
        using namespace Opacity_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        // bool passed2 = test( target, aspect2, test_2);

        // call additional tests here.

        return passed1
            // and passed2
            ;
    }

    namespace Opacity_tests
    {
        using nut::Opacity;

        bool test_1()
        {
            bool passed(true);
            typedef double fp_t;
            // bogus random number gen type: ok to test instantiation
            typedef uint32_t rng_t; 
            typedef Opacity<fp_t> op_t;
            typedef nut::Density<fp_t> rho_t;
            typedef nut::Temperature<fp_t> temp_t;
            
            std::vector<fp_t> nil(10,0);
            rho_t rho(nil,nil,nil,nil,nil,nil);
            temp_t T(nil,nil,nil);
            op_t op(rho,T);

            return passed;
        } // test_1


        // bool test_2()
        // {
        //     bool passed(false);
        //     typedef double fp_t;
        //     // bogus random number gen type: ok to test instantiation
        //     typedef uint32_t rng_t; 
        //     typedef Opacity<fp_t,rng_t> op_t;
        //     typedef nut::Density<fp_t> rho_t;
        //     typedef nut::Temperature<fp_t> temp_t;
            
        //     fp_t const e_nu = 1.0;
        //     fp_t const e_e_minus = 1.0;
        //     fp_t const rho = 1.0;
        //     fp_t sigma = op_t::sigma<nut::sigmas::nu_e_e_minus>(
        //         e_nu,e_e_minus,rho);

        //     std::cout << "sigma = " << sigma << std::endl;
                
        //     return passed;
        // } // test_2



        // define additional tests here.

    } // Opacity_tests::

} // Nut_Test::



// version
// $Id$

// End of file
