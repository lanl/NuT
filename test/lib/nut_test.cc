// nut_test.cc
// T. M. Kelley
// Jan 03, 2011
// Implementation, nut_test
// (c) Copyright 2011 LANSLLC all rights reserved.

#include "nut_test.hh"
#include "nut_Test_apply_event.hh"
#include "nut_Test_Assert.hh"
#include "nut_Test_Assert_off.hh"
#include "nut_Test_decide_event.hh"
#include "nut_Test_Density.hh"
#include "nut_Test_Mesh.hh"
#include "nut_Test_Mesh3D.hh"
#include "nut_Test_Opacity.hh"
#include "nut_Test_Particle.hh"
#include "nut_Test_Planck.hh"
#include "nut_Test_RNG.hh"
#include "nut_Test_Tally.hh"
#include "nut_Test_Temperature.hh"
#include "nut_Test_transport_particle.hh"
#include "nut_Test_species_name.hh"
#include "nut_Test_fileio.hh"
#include "nut_Test_MatState.hh"
#include "nut_Test_partition.hh"
#include "test_aux.hh"

namespace NutTest_standalones
{
    char target[] = "neutrino transport";

    char aspect1[] = "tests of transport lib";

    namespace
    {
        bool check_requests(test_map & tests,
                            std::vector<std::string> const & names);
    }

    void print_available_tests( test_map & tests, std::ostream & o)
    {
        o << "available tests: " << std::endl;
        for(test_it curr = tests.begin(); curr != tests.end(); ++curr)
        {
            o << "'" << curr -> first <<"'" << std::endl;
        }
        o << "note: test names are case sensitive" << std::endl;
        return;
    }

    test_map get_tests()
    {
        using test_aux::test;
        using namespace Nut_Test;

        test_map tests;
        // register tests here
        // tests["name_used_on_command_line"] = Test_Desc(target, test);
        tests["assert"]     =  Test_Desc(Assert_tests::target, test_Assert);
        tests["density"]    =  Test_Desc(Density_tests::target, test_Density);
        tests["mesh"]       =  Test_Desc(Mesh_tests::target, test_Mesh);
        tests["mesh3D"]       =  Test_Desc(Mesh3D_tests::target, test_Mesh3D);
        tests["opacity"]    =  Test_Desc(Opacity_tests::target, test_Opacity);
        tests["planck"]     =  Test_Desc(Planck_tests::target, test_Planck);
        tests["ptcl"]       =  Test_Desc(particle_tests::target, test_Particle);
        tests["rng"]        =  Test_Desc(RNG_tests::target, test_RNG);
        tests["fileio"]     =  Test_Desc(fileio_tests::target, test_fileio);
        tests["tally"]      =  Test_Desc(Tally_tests::target, test_Tally);
        tests["assert_off"]   =  Test_Desc(Assert_off_tests::target,
                                           test_Assert_off);
        tests["spec_name"]    =  Test_Desc(species_name_tests::target,
                                           test_species_name);
        tests["temperature"]  =  Test_Desc(Temperature_tests::target,
                                           test_Temperature);
        tests["decide_event"] =  Test_Desc(decide_event_tests::target,
                                           test_decide_event);
        tests["apply_event"]  =  Test_Desc(apply_event_tests::target,
                                           test_apply_event);
        tests["mat_state"]    =  Test_Desc(MatState_tests::target,
                                           test_MatState);
        tests["transport_particle"] = Test_Desc(
            transport_particle_tests::target, test_transport_particle);
        tests["partition"]    =  Test_Desc(partition_tests::target, test_partition);

        return tests;
    } // get_tests

    bool run_specific_tests(std::vector<std::string> const & names)
    {
        bool all_passed(true);
        test_map tests = get_tests();

        bool reqsok = check_requests(tests,names);
        if(!reqsok)
        {
            print_available_tests(tests,std::cout);
            return -1;
        }
        for(size_t i = 0; i < names.size(); ++i)
        {
            test_it it = tests.find(names[i]);
            bool passed = test_aux::test( (it->second).first,"all",(it->second).second);
            all_passed = all_passed && passed;
        }
        return all_passed;
    }


    bool run_standalone_tests()
    {
        using test_aux::test;
        using namespace Nut_Test;

        bool assert_passed =  test(Assert_tests::target, "all",
                                   test_Assert);

        bool assert_off_passed =  test(Assert_off_tests::target,
                                       "all", test_Assert_off);

        bool density_passed =  test(Density_tests::target, "all",
                                       test_Density);

        bool mesh_passed  =  test(Mesh_tests::target, "all",
                                 test_Mesh);

        bool mesh3D_passed  =  test(Mesh3D_tests::target, "all",
                                 test_Mesh3D);

        bool opacity_passed  =  test(Opacity_tests::target, "all",
                                 test_Opacity);

        bool planck_passed  =  test(Planck_tests::target, "all",
                                 test_Planck);

        bool ptcl_passed =  test(particle_tests::target, "all",
                                 test_Particle);

        bool rng_passed  =  test(RNG_tests::target, "all",
                                 test_RNG);

        bool spec_name_passed =  test(species_name_tests::target, "all",
                                      test_species_name);

        bool tally_passed =  test(Tally_tests::target, "all",
                                       test_Tally);

        bool temperature_passed =  test(Temperature_tests::target,
                                        "all", test_Temperature);

        bool decide_event_passed  =  test(decide_event_tests::target,
                                          "all", test_decide_event);

        bool apply_event_passed =  test(apply_event_tests::target,
                                        "all", test_apply_event);

        bool transport_particle_passed = test(transport_particle_tests::target,
                                        "all", test_transport_particle);

        bool fileio_passed =  test(fileio_tests::target, "all", test_fileio);

        bool mat_state_passed =  test(MatState_tests::target, "all",
                                      test_MatState);

        return
            apply_event_passed
            and assert_passed
            and assert_off_passed
            and decide_event_passed
            and density_passed
            and mesh_passed
            and mesh3D_passed
            and opacity_passed
            and planck_passed
            and ptcl_passed
            and rng_passed
            and tally_passed
            and temperature_passed
            and transport_particle_passed
            and fileio_passed
            and mat_state_passed
            and spec_name_passed
            ;
    } // run_standalone_tests


    namespace
    {

        bool check_requests(NutTest_standalones::test_map & tests,
                            std::vector<std::string> const & names)
        {
            bool ok(true);

            if(names.size() > tests.size())
            {
                ok = false;
            }

            for(size_t i = 0; i < names.size(); ++i)
            {
                if(tests.count(names[i]) == 0)
                {
                    ok = false;
                    std::cerr << "could not find requested test '"
                              << names[i] << "'" << std::endl;
                }
            }

            return ok;
        } // check_requests

    } // anonymous::
} // NutTest::

// version
// $Id$

// End of file
