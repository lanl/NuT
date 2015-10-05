//T. M. Kelley (c) 2011 LANS LLC

#include "Tally.hh"
#include "nut_Test_Tally.hh"
#include "test_aux.hh"
#include "types.hh"
#include "soft_equiv.hh"
#include <iomanip>


namespace Nut_Test
{
    namespace Tally_tests
    {
        // target describes the code being tested
        char target[] = "Tally";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "inst. & init.";
        bool test_2();
        char aspect2[] = "deposit_inelastic_el_scat";
        bool test_3();
        char aspect3[] = "deposit_energy";
        bool test_4();
        char aspect4[] = "deposit_momentum_elastic";
        bool test_5();
        char aspect5[] = "count_electron_scatter (nu_e)";
        bool test_6();
        char aspect6[] = "count_electron_scatter (nu_e_bar)";
        bool test_7();
        char aspect7[] = "count_electron_scatter (nu_mu)";
        bool test_8();
        char aspect8[] = "count_electron_scatter (nu_mu_bar)";
        bool test_9();
        char aspect9[] = "count_electron_scatter (nu_tau)";
        bool test_10();
        char aspect10[] = "count_electron_scatter (nu_tau_bar)";
        bool test_11();
        char aspect11[] = "count_nucleon_abs (nu_e)";
        bool test_12();
        char aspect12[] = "count_nucleon_abs (nu_e_bar)";
        bool test_13();
        char aspect13[] = "count_nucleon_elastic_scatter";
        bool test_14();
        char aspect14[] = "count_escape";
        bool test_15();
        char aspect15[] = "count_reflect";
        bool test_16();
        char aspect16[] = "count_cell_bdy";
        bool test_17();
        char aspect17[] = "count_cutoff";
        bool test_18();
        char aspect18[] = "count_census (nu_e)";
        bool test_19();
        char aspect19[] = "count_census (nu_e_bar)";
        bool test_20();
        char aspect20[] = "count_census (nu_mu)";
        bool test_21();
        char aspect21[] = "count_census (nu_mu_bar)";
        bool test_22();
        char aspect22[] = "count_census (nu_tau)";
        bool test_23();
        char aspect23[] = "count_census (nu_tau_bar)";

    }

    bool test_Tally()
    {
        using namespace Tally_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        bool passed2 = test( target, aspect2, test_2);

        bool passed3 = test( target, aspect3, test_3);

        bool passed4 = test( target, aspect4, test_4);

        bool passed5 = test( target, aspect5, test_5);

        bool passed6 = test( target, aspect6, test_6);

        bool passed7 = test( target, aspect7, test_7);

        bool passed8 = test( target, aspect8, test_8);

        bool passed9 = test( target, aspect9, test_9);

        bool passed10 = test( target, aspect10, test_10);

        bool passed11 = test( target, aspect11, test_11);

        bool passed12 = test( target, aspect12, test_12);

        bool passed13 = test( target, aspect13, test_13);

        bool passed14 = test( target, aspect14, test_14);

        bool passed15 = test( target, aspect15, test_15);

        bool passed16 = test( target, aspect16, test_16);

        bool passed17 = test( target, aspect17, test_17);

        bool passed18 = test( target, aspect18, test_18);

        bool passed19 = test( target, aspect19, test_19);

        bool passed20 = test( target, aspect20, test_20);

        bool passed21 = test( target, aspect21, test_21);

        bool passed22 = test( target, aspect22, test_22);

        bool passed23 = test( target, aspect23, test_23);

        // call additional tests here.

        return passed1 and passed2 and passed3 and passed4 and passed5
            and passed6 and passed7 and passed8 and passed9 and passed10
            and passed11 and passed12 and passed13 and passed14 and passed15
            and passed16 and passed17 and passed18 and passed19 and passed20
            and passed21 and passed22 and passed23;
    }

    namespace Tally_tests
    {
        using test_aux::check_one_changed;
        using test_aux::check_two_changed;
        using test_aux::check_same;
        using test_aux::check_same_verb;
        using test_aux::comp_verb;


        bool test_1()
        {
            bool passed(true);

            typedef float fp_t;

            size_t n_cells(100);

            nut::Tally<fp_t> tally(n_cells);

            return passed;
        } // test_1


        bool test_2()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            nut::Tally<fp_t> tally(n_cells), ref(n_cells);

            fp_t const ei(2.0), ef(1.9);
            fp_t const omega_i(1.0), omega_f(-1.0);
            fp_t const wt(0.5);
            nut::Species const s(nut::nu_e);
            nut::cell_t const c(21);
            tally.deposit_inelastic_scat(c,ei,ef,omega_i,omega_f,wt,s);

            ref.momentum[c-1] = 0.5*(ei*omega_i - ef*omega_f);
            ref.energy[c-1] = 0.5*(ei-ef);

            passed = check_same(&tally.momentum,&ref.momentum) and passed;
            passed = check_same(&tally.energy,&ref.energy) and passed;

            return passed;
        } // test_2


        bool test_3()
        {
            bool passed(true);
            typedef float fp_t;
            typedef std::vector<fp_t> vf;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            fp_t const wt = 0.2;
            nut::cell_t const c(21);
            fp_t const e  = 3.0;
            tally.deposit_energy(c,wt,e);

            passed = check_one_changed<t_t,vf>(tally,ref,&tally.energy)
                and passed;

            ref.energy[20] = 0.6;

            passed = check_same(&tally.energy,&ref.energy) and passed;

            return passed;
        } // test_3


        bool test_4()
        {
            bool passed(true);
            typedef float fp_t;
            typedef std::vector<fp_t> vf;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            fp_t const o  = 0.3;
            fp_t const wt = 0.2;
            fp_t const e  = 4.0;
            nut::cell_t const c(21);
            tally.deposit_momentum_elastic(c,o,e,wt);

            passed = check_one_changed<t_t,vf>(tally,ref,&tally.momentum)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.momentum[c - 1] = 0.2 * 0.3 * 4.0;
            passed = check_same_verb(&tally.momentum,&ref.momentum,
                                     comp_verb<fp_t>()) and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_4


        bool test_5()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_e);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_lepton_scatter(c,s,n);

            passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_nu_e_el_scat)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_nu_e_el_scat[c - 1] = 1;
            passed = check_same(&tally.n_nu_e_el_scat, &ref.n_nu_e_el_scat)
                and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_5


        bool test_6()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_e_bar);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_lepton_scatter(c,s,n);

            passed = check_one_changed<t_t,vc>(tally,ref,
                                               &tally.n_nu_e_bar_pos_scat)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_nu_e_bar_pos_scat[c - 1] = 1;
            passed = check_same(&tally.n_nu_e_bar_pos_scat,
                                &ref.n_nu_e_bar_pos_scat) and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_6


        bool test_7()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_mu);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_lepton_scatter(c,s,n);

            passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_nu_x_el_scat)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_nu_x_el_scat[c - 1] = 1;
            passed = check_same(&tally.n_nu_x_el_scat, &ref.n_nu_x_el_scat)
                and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_7


        bool test_8()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_mu_bar);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_lepton_scatter(c,s,n);

            passed = check_one_changed<t_t,vc>(tally,ref,
                                               &tally.n_nu_x_bar_pos_scat)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_nu_x_bar_pos_scat[c - 1] = 1;
            passed = check_same(&tally.n_nu_x_bar_pos_scat,
                                &ref.n_nu_x_bar_pos_scat) and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_8


        bool test_9()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_tau);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_lepton_scatter(c,s,n);

            passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_nu_x_el_scat)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_nu_x_el_scat[c - 1] = 1;
            passed = check_same(&tally.n_nu_x_el_scat, &ref.n_nu_x_el_scat)
                and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_9


        bool test_10()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_tau_bar);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_lepton_scatter(c,s,n);

            passed = check_one_changed<t_t,vc>(tally,ref,
                                               &tally.n_nu_x_bar_pos_scat)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_nu_x_bar_pos_scat[c - 1] = 1;
            passed = check_same(&tally.n_nu_x_bar_pos_scat,
                                &ref.n_nu_x_bar_pos_scat) and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_10


        bool test_11()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_e);
            nut::cell_t const c(21);
            fp_t const wt(59.2);
            t_t::cntr_t const n(1);
            tally.count_nucleon_abs(c,s,wt,n);

            ref.n_nu_e_nucl_abs[c - 1] = 1;
            ref.ew_nu_e_nucl_abs[c - 1] = 59.2;

            passed = check_same(&tally.n_nu_e_nucl_abs,
                                  &ref.n_nu_e_nucl_abs) and passed;
            passed = check_same(&tally.ew_nu_e_nucl_abs,
                                  &ref.ew_nu_e_nucl_abs) and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_11


        bool test_12()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_e_bar);
            nut::cell_t const c(21);
            fp_t const wt(59.2);
            t_t::cntr_t const n(1);
            tally.count_nucleon_abs(c,s,wt,n);


            ref.n_nu_e_bar_nucl_abs[c - 1] = 1;
            ref.ew_nu_e_bar_nucl_abs[c - 1] = 59.2;

            passed = check_same(&tally.n_nu_e_bar_nucl_abs,
                                &ref.n_nu_e_bar_nucl_abs) and passed;
            passed = check_same(&tally.ew_nu_e_bar_nucl_abs,
                                &ref.ew_nu_e_bar_nucl_abs) and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_12


        bool test_13()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_nucleon_elastic_scatter(c,n);

            passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_nucl_el_scat)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_nucl_el_scat[c - 1] = 1;
            passed = check_same(&tally.n_nucl_el_scat, &ref.n_nucl_el_scat)
                and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_13


        bool test_14()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::cell_t const c(21);
            fp_t const ew(37.64);
            t_t::cntr_t const n(1);
            tally.count_escape(c,ew,n);

            // two change
            // passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_escape)
            //     and passed;
            // if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_escape[c - 1] = 1;
            ref.ew_escaped[c - 1] = ew;
            passed = check_same(&tally.n_escape, &ref.n_escape)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            passed = check_same(&tally.ew_escaped, &ref.ew_escaped)
                and passed;

            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_14


        bool test_15()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_reflect(c,n);

            // two change
            // passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_escape)
            //     and passed;
            // if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_reflect[c - 1] = 1;
            passed = check_same(&tally.n_reflect, &ref.n_reflect)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_15


        bool test_16()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_cell_bdy(c,n);

            passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_cell_bdy)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_cell_bdy[c - 1] = 1;
            passed = check_same(&tally.n_cell_bdy, &ref.n_cell_bdy)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_16


        bool test_17()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;
            typedef t_t::vc vc;

            t_t tally(n_cells), ref(n_cells);

            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            tally.count_cutoff(c,n);

            passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_cutoff)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_cutoff[c - 1] = 1;
            passed = check_same(&tally.n_cutoff, &ref.n_cutoff)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_17


        bool test_18()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_e);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            fp_t const ew(17.13);
            tally.count_census(c,ew,s,n);

            passed = check_two_changed(tally,ref,
                                       &tally.n_census_nu_e,
                                       &tally.ew_census_nu_e)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_census_nu_e[c - 1] = n;
            ref.ew_census_nu_e[c - 1] = ew;
            passed = check_same(&tally.n_census_nu_e, &ref.n_census_nu_e)
                and passed;
            passed = check_same(&tally.ew_census_nu_e, &ref.ew_census_nu_e)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_18


        bool test_19()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_e_bar);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            fp_t const ew(17.13);
            tally.count_census(c,ew,s,n);

            passed = check_two_changed(tally,ref,
                                       &tally.n_census_nu_e_bar,
                                       &tally.ew_census_nu_e_bar)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_census_nu_e_bar[c - 1] = n;
            ref.ew_census_nu_e_bar[c - 1] = ew;
            passed = check_same(&tally.n_census_nu_e_bar, &ref.n_census_nu_e_bar)
                and passed;
            passed = check_same(&tally.ew_census_nu_e_bar, &ref.ew_census_nu_e_bar)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_19


        bool test_20()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_mu);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            fp_t const ew(17.13);
            tally.count_census(c,ew,s,n);

            passed = check_two_changed(tally,ref,
                                       &tally.n_census_nu_x,
                                       &tally.ew_census_nu_x)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_census_nu_x[c - 1] = n;
            ref.ew_census_nu_x[c - 1] = ew;
            passed = check_same(&tally.n_census_nu_x, &ref.n_census_nu_x)
                and passed;
            passed = check_same(&tally.ew_census_nu_x, &ref.ew_census_nu_x)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_20


        bool test_21()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_mu_bar);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            fp_t const ew(17.13);
            tally.count_census(c,ew,s,n);

            passed = check_two_changed(tally,ref,
                                       &tally.n_census_nu_x_bar,
                                       &tally.ew_census_nu_x_bar)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_census_nu_x_bar[c - 1] = n;
            ref.ew_census_nu_x_bar[c - 1] = ew;
            passed = check_same(&tally.n_census_nu_x_bar, &ref.n_census_nu_x_bar)
                and passed;
            passed = check_same(&tally.ew_census_nu_x_bar, &ref.ew_census_nu_x_bar)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_21


        bool test_22()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_tau);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            fp_t const ew(17.13);
            tally.count_census(c,ew,s,n);

            passed = check_two_changed(tally,ref,
                                       &tally.n_census_nu_x,
                                       &tally.ew_census_nu_x)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_census_nu_x[c - 1] = n;
            ref.ew_census_nu_x[c - 1] = ew;
            passed = check_same(&tally.n_census_nu_x, &ref.n_census_nu_x)
                and passed;
            passed = check_same(&tally.ew_census_nu_x, &ref.ew_census_nu_x)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_22


        bool test_23()
        {
            bool passed(true);
            typedef float fp_t;
            size_t n_cells(100);

            typedef nut::Tally<fp_t> t_t;

            t_t tally(n_cells), ref(n_cells);

            nut::Species const s(nut::nu_tau_bar);
            nut::cell_t const c(21);
            t_t::cntr_t const n(1);
            fp_t const ew(17.13);
            tally.count_census(c,ew,s,n);

            passed = check_two_changed(tally,ref,
                                       &tally.n_census_nu_x_bar,
                                       &tally.ew_census_nu_x_bar)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            ref.n_census_nu_x_bar[c - 1] = n;
            ref.ew_census_nu_x_bar[c - 1] = ew;
            passed = check_same(&tally.n_census_nu_x_bar, &ref.n_census_nu_x_bar)
                and passed;
            passed = check_same(&tally.ew_census_nu_x_bar, &ref.ew_census_nu_x_bar)
                and passed;
            if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

            return passed;
        } // test_23


        // define additional tests here.


        // additional helpers

    } // Tally_tests::

} // Nut_Test::



// version
// $Id$

// End of file
