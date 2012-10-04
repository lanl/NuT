// test_aux.cc
// T. M. Kelley
// Mar 23, 2009
// Implementation, test_aux
// (c) Copyright 2009 LANSLLC all rights reserved.

#include "test_aux.hh"
#include <stdio.h>

namespace test_aux
{

    bool test( const char *target, const char *aspect, bool (*test2run)())
    {
        pre_report( target, aspect);
        bool passed = test2run();
        post_report( target, aspect, passed);
        return passed;
    } // test

    void pre_report( const char* target, const char* aspect)
    {
        printf("----------------------------------------\n"
               "Testing %s, %s\n", target, aspect);
        return;
    }

    void post_report( const char* target, const char* aspect, bool passed)
    {
        if(passed)
        {
            printf("Test of %s, %s PASSED\n"
                   "----------------------------------------\n",
                   target,aspect);
        }
        else
        {
            printf("Test of %s, %s FAILED\n"
                   "----------------------------------------\n",
                   target,aspect);
        }
        return;
    } 

    bool test( std::string const &target, std::string const & aspect, bool (*test2run)())
    {
        pre_report( target, aspect);
        bool passed = test2run();
        post_report( target, aspect, passed);
        return passed;
    } // test

    void pre_report( std::string const & target, std::string const & aspect)
    {
        std::cout << "----------------------------------------\n"
                  << "Testing " << target << ", " << aspect << "\n";
        return;
    }

    void post_report( std::string const & target, std::string const & aspect, bool passed)
    {
        if(passed)
        {
            std::cout << "Test of " << target<< ", " << aspect << " PASSED\n"
                      << "----------------------------------------\n";
        }
        else
        {
            std::cout << "Test of " << target<< ", " << aspect << " FAILED\n"
                      << "----------------------------------------\n";
        }
        return;
    } 

} // test_aux::


// version
// $Id$

// End of file
