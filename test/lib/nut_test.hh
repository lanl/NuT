// nut_test.hh
// T. M. Kelley
// Jan 03, 2011
// Header for nut_test
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef NUT_TEST_H
#define NUT_TEST_H

#include <map>
#include <vector>
#include <string>

namespace NutTest_standalones
{
    extern char target[]; 

    extern char aspect1[]; 

    bool run_standalone_tests();

    // struct Test_Desc
    typedef std::pair<std::string, bool (*) ()> Test_Desc;
    // {
    //     std::string target;
    //     bool (*test) (); 
    //     Test_Desc() : target(""), test(NULL) {}
    // }; // Test_Desc

    typedef std::map<std::string,Test_Desc> test_map;
    typedef test_map::iterator test_it;
    typedef test_map::value_type test_t;

    void print_available_tests( test_map const & tests, std::ostream & o);

    test_map get_tests();

    bool run_specific_tests(std::vector<std::string> const & names);

} // NutTest


#endif



// version
// $Id$

// End of file
