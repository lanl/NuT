//T. M. Kelley (c) 2011 LANS LLC

#include "fileio.hh"

#include "nut_Test_fileio.hh"
#include "test_aux.hh"
#include "expect.hh"
#include <iostream>
#include <sstream>

namespace Nut_Test
{
    namespace fileio_tests
    {
        // target describes the code being tested
        char target[] = "fileio";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "line_to_struct (float)";

        bool test_2();
        char aspect2[] = "line_to_struct (double)";

        bool test_3();
        char aspect3[] = "read_mat_state_file (double)";

    }

    bool test_fileio()
    {
        using namespace fileio_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        bool passed2 = test( target, aspect2, test_2);

        bool passed3 = test( target, aspect3, test_3);

        // call additional tests here.

        return passed1 and passed2 and passed3;
    }

    namespace fileio_tests
    {
        std::string line(" 1  1.5000E-02  1.7348E+06  2.7438E+11  3.5650E+08"
                         "  8.1058E+00  0.0000E+00  5.1030E+02  1.7185E+02  "
                         "5.2054E+21  3.1954E+07  2.2382E+02  3.1954E+07  "
                         "2.2248E+02  1.1083E+08  4.7589E+02  1.1083E+08  "
                         "4.7589E+02  9.5198E+07  4.7589E+02");

        using test_aux::expect;

        bool test_1()
        {
            bool passed(false);
            typedef float fp_t;
            typedef nut::MatStateRowP<fp_t> row_t;

            row_t row = nut::line_to_struct<fp_t>(line);

            row_t exp;
            exp.zone = 1;
            exp.m_encl = 1.5000e-2;
            exp.radius = 1.7348e+6;
            exp.density = 2.7438e+11;
            exp.velocity = 3.565e+08;
            exp.ye = 8.1058;
            
            exp.eta = 0.0;
            exp.temperature = 5.1030e+02 ;
            exp.entropy = 1.7185e+02;
            
            exp.u = 5.2054e+21;
            exp.lnue_capture = 3.1954e+07;
            exp.enue_capture = 2.2382e+02;
            
            exp.lnueb_capture = 3.1954e+7;
            exp.enueb_capture = 2.2248e+02;
            exp.lnue_pair = 1.1083e+08;
            
            exp.enue_pair = 4.7589e+02;
            exp.lnueb_pair = 1.1083e+08;
            exp.enueb_pair  = 4.7589e+02;
            
            exp.lnux_pair = 9.5198e+07;
            exp.enux_pair = 4.7589e+02;
            
            passed = row == exp;
            if(!passed)
            {
                std::cout << "row: ";
                std::cout << row.zone <<","<< row.m_encl <<","<< row.radius <<","<< 
                    row.density       <<","<< row.velocity <<","<< row.ye <<","<< 
                    row.eta <<","<< row.temperature   <<","<< row.entropy      <<","<< 
                    row.u   <<","<< row.lnue_capture  <<","<< row.enue_capture <<","<< 
                    row.lnueb_capture <<","<< row.enueb_capture <<","<< row.lnue_pair    <<","<<
                    row.enue_pair <<","<< row.lnueb_pair <<","<< row.enueb_pair <<","<< 
                    row.lnux_pair <<","<< row.enux_pair;
                std::cout << "\nexp: ";
                std::cout << exp.zone <<","<< exp.m_encl <<","<< exp.radius <<","<< 
                    exp.density       <<","<< exp.velocity <<","<< exp.ye <<","<< 
                    exp.eta <<","<< exp.temperature   <<","<< exp.entropy      <<","<< 
                    exp.u   <<","<< exp.lnue_capture  <<","<< exp.enue_capture <<","<< 
                    exp.lnueb_capture <<","<< exp.enueb_capture <<","<< exp.lnue_pair    <<","<<
                    exp.enue_pair <<","<< exp.lnueb_pair <<","<< exp.enueb_pair <<","<< 
                    exp.lnux_pair <<","<< exp.enux_pair;
                std::cout << std::endl;
            }
            
            return passed;
        } // test_1

        bool test_2()
        {
            bool passed(false);
            typedef double fp_t;
            typedef nut::MatStateRowP<fp_t> row_t;

            row_t row = nut::line_to_struct<fp_t>(line);

            row_t exp;
            exp.zone = 1;
            exp.m_encl = 1.5000e-2;
            exp.radius = 1.7348e+6;
            exp.density = 2.7438e+11;
            exp.velocity = 3.565e+08;
            exp.ye = 8.1058;
            
            exp.eta = 0.0;
            exp.temperature = 5.1030e+02 ;
            exp.entropy = 1.7185e+02;
            
            exp.u = 5.2054e+21;
            exp.lnue_capture = 3.1954e+07;
            exp.enue_capture = 2.2382e+02;
            
            exp.lnueb_capture = 3.1954e+7;
            exp.enueb_capture = 2.2248e+02;
            exp.lnue_pair = 1.1083e+08;
            
            exp.enue_pair = 4.7589e+02;
            exp.lnueb_pair = 1.1083e+08;
            exp.enueb_pair  = 4.7589e+02;
            
            exp.lnux_pair = 9.5198e+07;
            exp.enux_pair = 4.7589e+02;
            
            passed = row == exp;
            if(!passed)
            {
                std::cout << "row: ";
                std::cout << row.zone <<","<< row.m_encl <<","<< row.radius <<","<< 
                    row.density       <<","<< row.velocity <<","<< row.ye <<","<< 
                    row.eta <<","<< row.temperature   <<","<< row.entropy      <<","<< 
                    row.u   <<","<< row.lnue_capture  <<","<< row.enue_capture <<","<< 
                    row.lnueb_capture <<","<< row.enueb_capture <<","<< row.lnue_pair    <<","<<
                    row.enue_pair <<","<< row.lnueb_pair <<","<< row.enueb_pair <<","<< 
                    row.lnux_pair <<","<< row.enux_pair;
                std::cout << "\nexp: ";
                std::cout << exp.zone <<","<< exp.m_encl <<","<< exp.radius <<","<< 
                    exp.density       <<","<< exp.velocity <<","<< exp.ye <<","<< 
                    exp.eta <<","<< exp.temperature   <<","<< exp.entropy      <<","<< 
                    exp.u   <<","<< exp.lnue_capture  <<","<< exp.enue_capture <<","<< 
                    exp.lnueb_capture <<","<< exp.enueb_capture <<","<< exp.lnue_pair    <<","<<
                    exp.enue_pair <<","<< exp.lnueb_pair <<","<< exp.enueb_pair <<","<< 
                    exp.lnux_pair <<","<< exp.enux_pair;
                std::cout << std::endl;
            }
            
            return passed;
        } // test_2

        bool test_3()
        {
            bool passed(true);
            typedef double fp_t;
            typedef nut::MatStateRowP<fp_t> row_t;

            std::stringstream instr(line);
            std::vector<row_t> rows( nut::read_mat_state_file<fp_t>(instr));

            size_t const exp_sz(1u);

            passed = passed and expect(rows.size(),exp_sz,"vector size"); 
            row_t const & row(rows[0]);

            {
                row_t exp;
                exp.zone = 1;
                exp.m_encl = 1.5000e-2;
                exp.radius = 1.7348e+6;
                exp.density = 2.7438e+11;
                exp.velocity = 3.565e+08;
                exp.ye = 8.1058;
                
                exp.eta = 0.0;
                exp.temperature = 5.1030e+02 ;
                exp.entropy = 1.7185e+02;
                
                exp.u = 5.2054e+21;
                exp.lnue_capture = 3.1954e+07;
                exp.enue_capture = 2.2382e+02;
                
                exp.lnueb_capture = 3.1954e+7;
                exp.enueb_capture = 2.2248e+02;

                exp.lnue_pair = 1.1083e+08;                
                exp.enue_pair = 4.7589e+02;

                exp.lnueb_pair = 1.1083e+08;
                exp.enueb_pair  = 4.7589e+02;
                
                exp.lnux_pair = 9.5198e+07;
                exp.enux_pair = 4.7589e+02;
                
                passed = passed and exp == row;
            }

            return passed;
        } // test_3

        // define additional tests here.

    } // fileio_tests::

} // Nut_Test::



// version
// $Id$

// End of file
