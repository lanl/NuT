// test-main.cc
// exa11 class
// Nov 24, 2010
// Implementation, test-main
// (c) Copyright 2010 LANSLLC all rights reserved.

#include <stdio.h>
#include <vector>
#include "nut_test.hh"
#include "test_aux.hh"

int main(int argc, char ** argv)
{
    using test_aux::test;
    using namespace NutTest_standalones;

    // get individual names of tests from command line
    std::vector<std::string> names;

    if(argc > 1)
    {
        for(int i = 1; i < argc; ++i)
        {
            names.push_back(argv[i]);
        }
    }

    if(names.size() != 0)
    {
        run_specific_tests(names);
    }
    else
    {
        // standalone tests
        test( NutTest_standalones::target, "all",
              NutTest_standalones::run_standalone_tests);
    }

    return 0;
}



// version
// $Id$

// End of file
