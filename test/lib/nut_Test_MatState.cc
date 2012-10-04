//T. M. Kelley (c) 2011 LANS LLC

#include "MatState.hh"

#include "nut_Test_MatState.hh"
#include "test_aux.hh"
#include "expect.hh"
#include "fileio.hh"
#include <string>
#include <sstream>

namespace Nut_Test
{
    namespace MatState_tests
    {
        // target describes the code being tested
        char target[] = "MatState";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "instantiation and initialization";

    }

    bool test_MatState()
    {
        using namespace MatState_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        // call additional tests here.

        return passed1;
    }

    namespace MatState_tests
    {
        std::string line(" 1  1.5000E-02  1.7348E+06  2.7438E+11  3.5650E+08"
                         "  8.1058E+00  0.0000E+00  5.1030E+02  1.7185E+02  "
                         "5.2054E+21  3.1954E+07  2.2382E+02  3.1954E+07  "
                         "2.2248E+02  1.1083E+08  4.7589E+02  1.1083E+08  "
                         "4.7589E+02  9.5198E+07  4.7589E+02");

        using nut::MatState;
        using test_aux::expect;

        bool test_1()
        {
            // cf nut_Test_fileio.cc, test_3:
            bool passed(true);
            typedef double fp_t;
            typedef nut::MatStateRowP<fp_t> row_t;
            typedef std::vector<row_t> vecrow;

            std::stringstream instr(line);
            vecrow rows( nut::read_mat_state_file<fp_t>(instr));

            MatState<fp_t> state(rows);

            // check mat state
            // sizes of everything:
            size_t const exp_sz(1u);
            passed = passed 
                and expect(state.density.size(), exp_sz, "density size")
                and expect(state.luminosity.size(), exp_sz, "luminosity size")
                and expect(state.temperature.size(), exp_sz, "temperature size")
                and expect(state.velocity.size(), exp_sz, "velocity size")
                ;
            // values

            // density
            passed = passed
                and expect(state.density.rho_p[0],1.6404204182659538e+35,"rho_p")
                and expect(state.density.rho_e_minus[0],1.329691982638017e+36,"rho_e_minus")
                and expect(state.density.rho_e_plus[0],0.0,"rho_e_plus")
                and expect(state.density.rho_A[0],0.0,"rho_A")
                and expect(state.density.y_e[0],0.0,"y_e")
                and expect(state.density.abar[0],0.0,"abar")
                ;

            // luminosity
            passed = passed
                and expect(state.luminosity.nue[0],1.42784e+59,"nue")
                and expect(state.luminosity.nueb[0],1.42784e+59,"nuebar")
                and expect(state.luminosity.nux[0],9.5198000000000001e+58,"nux")
                ;

            // temperature
            passed = passed 
                and expect(state.temperature.T_p[0],43.974286020000001,"T_p")
                and expect(state.temperature.T_e_minus[0],43.974286020000001,"T_e_minus")
                and expect(state.temperature.T_e_plus[0],43.974286020000001,"T_e_plus")
                ;
            // velocity
            passed = passed 
                and expect(state.velocity.vs[0],3.565e+08,"v");
            return passed;
        }

        // define additional tests here.

    } // MatState_tests::

} // Nut_Test::




// version
// $Id$

// End of file
