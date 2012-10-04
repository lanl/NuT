//T. M. Kelley (c) 2011 LANS LLC

#include "nut_Test_Assert.hh"
#include "test_aux.hh"
#include <string>
#include <iostream>


namespace Nut_Test
{
    namespace Assert_tests
    {
        // target describes the code being tested
        char target[] = "Assert";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "Require (ON)";
        bool test_2();
        char aspect2[] = "Equal (ON)";
        bool test_3();
        char aspect3[] = "InOpenRange<double> (ON)";
        bool test_4();
        char aspect4[] = "GreaterThan<double> (ON)";
        bool test_5();
        char aspect5[] = "InOpenRange<uint32_t> (ON)";
        bool test_6();
        char aspect6[] = "GreaterThan<int64_t> (ON)";
    }

    bool test_Assert()
    {
        using namespace Assert_tests;
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

} // Nut_Test::

#ifndef REQUIRE_ON
#define REQUIRE_ON
#include "Assert.hh"
#undef REQUIRE_ON
#else
#include "Assert.hh"
#endif

namespace Nut_Test
{
    namespace Assert_tests
    {

        bool test_1()
        {
            bool passed(true);
            std::string err("test 1");
            try
            {
                nut::Require(false,err.c_str());
                passed = false;
            }
            catch(std::exception & exc)
            {
                if(exc.what() != err)
                {
                    std::cerr << "test_1: Caught assertion as expected"
                              << " wrong what: " << exc.what() << std::endl;
                    passed = false;
                }
            }
            try
            {
                nut::Require(true,err.c_str());
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
                nut::Equal(a,b,"a","b");
                passed = false;
            }
            catch(std::exception & exc)
            {
                if(exc.what() != err)
                {
                    std::cerr << "test_2: Caught assertion as expected"
                              << " wrong what: " << exc.what() << std::endl;
                    passed = false;
                }
            }
            try
            {
                double a = 1.0, b = 1.0;
                nut::Equal(a,b,"a","b");
                passed = passed and true;
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
                nut::InOpenRange(a,min,max,"a");
                passed = false;
                std::cerr << "test_3:" << __LINE__ 
                          << " Didn't throw assertion as expected" << std::endl;
            }
            catch(std::exception & exc)
            {
                if(exc.what() != err)
                {
                    std::cerr << "test_3:" << __LINE__ 
                              << " Caught assertion as expected"
                              << " wrong what: " << exc.what() << std::endl;
                    passed = false;
                }
            }
            try
            {
                double a = 2.0, min=1.0, max=2.000000000001;
                nut::InOpenRange(a,min,max,"a");
                passed = passed and true;
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
                nut::GreaterThan(a,min,"a");
                passed = false;
                std::cerr << "test_4:" << __LINE__ 
                          << " Didn't throw assertion as expected" << std::endl;
            }
            catch(std::exception & exc)
            {
                if(exc.what() != err)
                {
                    std::cerr << "test_4:" << __LINE__ 
                              << " Caught assertion as expected wrong what: "
                              << exc.what() << std::endl;
                    passed = false;
                }
            }
            try
            {
                double a = 2.0, min=1.0;
                nut::GreaterThan(a,min,"a");
                passed = passed and true;
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
                nut::InOpenRange(a,min,max,"a");
                passed = false;
                std::cerr << "test_5:" << __LINE__ 
                          << " Didn't throw assertion as expected" << std::endl;
            }
            catch(std::exception & exc)
            {
                if(exc.what() != err)
                {
                    std::cerr << "test_5:" << __LINE__ 
                              << " Caught assertion as expected wrong what: " 
                              << exc.what() << std::endl;
                    passed = false;
                }
            }
            try
            {
                uint32_t a = 2, min=1, max=3;
                nut::InOpenRange(a,min,max,"a");
                passed = passed and true;
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
                nut::GreaterThan(a,min,"a");
                passed = false;
                std::cerr << "test_6:" << __LINE__ 
                          << " Didn't throw assertion as expected" << std::endl;
            }
            catch(std::exception & exc)
            {
                if(exc.what() != err)
                {
                    std::cerr << "test_6:" << __LINE__ 
                              << " Caught assertion as expected wrong what: "
                              << exc.what() << std::endl;
                    passed = false;
                }
            }
            try
            {
                int64_t a = 2, min=-1;
                nut::GreaterThan(a,min,"a");
                passed = passed and true;
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

    } // Assert_tests::

} // Nut_Test::


// version
// $Id$

// End of file
