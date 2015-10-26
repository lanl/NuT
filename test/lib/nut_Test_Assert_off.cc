//T. M. Kelley (c) 2011 LANS LLC

#include "nut_Test_Assert_off.hh"
#include "test_aux.hh"
#include <string>
#include <iostream>

#ifdef REQUIRE_ON
#undef REQUIRE_ON
#include "Assert.hh"
#define REQUIRE_ON
#else
#include "Assert.hh"
#endif

namespace Nut_Test
{
    namespace Assert_off_tests
    {
        // target describes the code being tested
        char target[] = "Assert (REQUIRE_ON turned off)";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "Require (OFF)";
        bool test_2();
        char aspect2[] = "Equal (OFF)";
        bool test_3();
        char aspect3[] = "InOpenRange<double> (OFF)";
        bool test_4();
        char aspect4[] = "GreaterThan<double> (OFF)";
        bool test_5();
        char aspect5[] = "InOpenRange<uint32_t> (OFF)";
        bool test_6();
        char aspect6[] = "GreaterThan<int64_t> (OFF)";

    }

    bool test_Assert_off()
    {
        using namespace Assert_off_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);


        bool passed2 = test( target, aspect2, test_2);

        bool passed3 = test( target, aspect3, test_3);

        bool passed4 = test( target, aspect4, test_4);

        bool passed5 = test( target, aspect5, test_5);

        bool passed6 = test( target, aspect6, test_6);

        // call additional tests here.

        return passed1 and passed2 and passed3 and passed4 and passed5
            and passed6;
    }

    namespace Assert_off_tests
    {

        bool test_1()
        {
            bool passed(true);
            std::string err("test 1");
            try
            {
                dbc::Require(false,err.c_str());
                passed = passed and true;
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_1: Threw assertion--not expected"
                          << std::endl;
            }
            try
            {
                dbc::Require(true,err.c_str());
                passed = passed and true;
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_1: Threw assertion--not expected"
                          << std::endl;
            }

            return passed;
        } // test_1




        bool test_2()
        {
            bool passed(true);
            std::string err("a != b");
            try
            {
                double a = 1.0, b = 2.0;
                dbc::Equal(a,b,"a","b");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_2:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }
            try
            {
                double a = 1.0, b = 1.0;
                dbc::Equal(a,b,"a","b");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_2:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }

            return passed;
        } // test_2


        bool test_3()
        {
            bool passed(true);
            std::string err("a (2) was not in range (1,2)");
            try
            {
                double a = 2.0, min=1.0, max=2.0;
                dbc::InOpenRange(a,min,max,"a");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_3:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }
            try
            {
                double a = 2.0, min=1.0, max=2.000000000001;
                dbc::InOpenRange(a,min,max,"a");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_3:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }

            return passed;
        } // test_3


        bool test_4()
        {
            bool passed(true);
            std::string err("a (2) <= 3");
            try
            {
                double a = 2.0, min=3.0;
                dbc::GreaterThan(a,min,"a");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_4:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }
            try
            {
                double a = 2.0, min=1.0;
                dbc::GreaterThan(a,min,"a");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_4:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }

            return passed;
        } // test_4


        bool test_5()
        {
            bool passed(true);
            std::string err("a (2) was not in range (1,2)");
            try
            {
                uint32_t a = 2, min=1, max=2;
                dbc::InOpenRange(a,min,max,"a");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_5:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }
            try
            {
                uint32_t a = 2, min=1, max=3;
                dbc::InOpenRange(a,min,max,"a");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_5:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }

            return passed;
        } // test_5


        bool test_6()
        {
            bool passed(true);
            std::string err("a (-2) <= 3");
            try
            {
                int64_t a = -2, min=3;
                dbc::GreaterThan(a,min,"a");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_6:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }
            try
            {
                int64_t a = 2, min=-1;
                dbc::GreaterThan(a,min,"a");
            }
            catch(std::exception & exc)
            {
                passed = false;
                std::cerr << "test_6:" << __LINE__
                          << " Caught assertion--not expected" << std::endl;
            }

            return passed;
        } // test_6


        // define additional tests here.

    } // Assert_off_tests::

} // Nut_Test::



// version
// $Id$

// End of file
