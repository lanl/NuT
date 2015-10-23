//T. M. Kelley (c) 2011 LANS LLC

#define USING_HFBC_SIGMAS

#include "nut_Test_decide_event.hh"
#include "test_aux.hh"
#include "Density.hh"
#include "Temperature.hh"
#include "Opacity.hh"
#include "decision.hh"
#include "RNG.hh"
#include "types.hh"
#include "Particle.hh"
#include "Mesh.hh"
#include "Velocity.hh"
#include "soft_equiv.hh"
#include <algorithm>
#include <iterator>

namespace Nut_Test
{
    namespace decide_event_tests
    {
        // target describes the code being tested
        char target[] = "event decision";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "decide_scatter_event: nucleon_abs";
        bool test_2();
        char aspect2[] = "decide_boundary_event";
        bool test_3();
        char aspect3[] = "decide_event: stream to cell boundary";
        bool test_4();
        char aspect4[] = "decide_event: stream through 10 steps, escapes on the last";
        bool test_5();
        char aspect5[] = "decide_scatter_event: nucleon_elastic_scatter";
        bool test_6();
        char aspect6[] = "decide_scatter_event: electron_scatter";

        bool test_7();
        char aspect7[] = "decide_event: nucleon_abs";
        bool test_8();
        char aspect8[] = "decide_event: nucleon_elastic_scatter";
        bool test_9();
        char aspect9[] = "decide_event: electron_scatter";
    }

    bool test_decide_event()
    {
        using namespace decide_event_tests;
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

        // call additional tests here.

        return passed1 and passed2 and passed3 and passed4 and passed5
            and passed6 and passed7 and passed8
            and passed9
            ;
    }

    namespace decide_event_tests
    {
        using nut::Species;
        using nut::cell_t;
        using nut::cell_t;
        using nut::geom_t;
        using nut::bdy_types::descriptor;
        using nut::vec_t;
        using nut::Event;

        typedef double                    fp_t;
        typedef std::vector<fp_t>         vf;
        typedef std::vector<vec_t<1>>         vec_vec;
        typedef nut::Density<fp_t>        rho_t;
        typedef nut::Temperature<fp_t>    T_t;
        typedef nut::Buffer_RNG<fp_t>     BRNG;
        typedef nut::Opacity<fp_t>  OpB;
        typedef nut::Particle<fp_t, BRNG> p_t;
        typedef nut::Sphere_1D<cell_t,geom_t,descriptor> mesh_t;
        typedef nut::Velocity<fp_t,1>       v_t;
        size_t ncells(10);
        vf nullv(ncells,fp_t(0));

        fp_t const pmg = 1.67262e-24;  // proton mass in gram

        // create empty Density and Temperature objects: tests can copy
        // and modify these.
        rho_t const rho( nullv, nullv, nullv, nullv, nullv, nullv);
        T_t const T( nullv, nullv, nullv);

        vec_vec zerovs(ncells);
        v_t const vel0s(zerovs);


        /* This test uses the setup of Python test_2 of decide_event to test
         * decide_scatter_event. We use the same random seeds, starting at the
         * third seed, since the test of decide_event burns the first one for
         * sampling the exponential to come up with d_collision, and the second
         * to choose the interaction channel. */
        bool test_1()
        {
            bool passed(true);

            using nut::decide_scatter_event;

            rho_t rho1(rho);
            T_t   T1(T);
            rho1.rho_p[0] = 1e14/pmg;

            OpB op(rho1,T1);
            // from Python test_2 of decide_event
            fp_t const rns[] = {// 0.30897681610609407, 1st seed burned for d_coll
                                // 0.92436920169905545,
                                0.21932404923057958};
            BRNG rng(rns,3);
            fp_t const e     = 5.0;
            cell_t const cell  = 1;
            Species const s(nut::nu_e);

            Event event = decide_scatter_event(rng,e,cell,op,s);

            Event event_exp = Event::nucleon_abs;

            passed = event == event_exp;
            if(!passed)
            {
                std::cout << "event was " << nut::event_name(event) << ", should have been "
                          << nut::event_name(Event::nucleon_abs) << std::endl;
            }
            return passed;
        } // test_1


        // generate a mesh that's reflective at the inner boundary,
        // vacuum at the outer boundary, and transmissive at all others
        struct gen_bdy_types
        {
            nut::bdy_types::descriptor operator()(){
                if(ctr++ == 0) return nut::bdy_types::descriptor::R;
                if(ctr == nbdy) return nut::bdy_types::descriptor::V;
                return nut::bdy_types::descriptor::T;
            }
            explicit gen_bdy_types(cell_t const nbdy_) : ctr(0),nbdy(nbdy_) {}
            cell_t ctr;
            cell_t const nbdy;
        }; // gen_bdy_types


        // generate a uniformly-spaced mesh, with spacing dx
        struct gen_bounds
        {
            fp_t operator()(){return dx * ctr++;}
            explicit gen_bounds(fp_t const dx_) : dx(dx_),ctr(0) {}
            fp_t const dx;
            cell_t ctr;
        }; // gen_bdy_types


        bool test_2()
        {
            bool passed(true);
            using nut::Event;

            using nut::decide_boundary_event;
            using nut::event_name;
            typedef mesh_t::vb vb;
            typedef mesh_t::vbd vbd;
            // generate uniform mesh
            cell_t n_cells(10);
            vb bounds(n_cells+1);
            vbd b_types(n_cells+1);
            std::generate(b_types.begin(),b_types.end(),gen_bdy_types(n_cells+1));
            std::generate(bounds.begin(),bounds.end(),gen_bounds(1.0));

            // reflect from innermost face
            mesh_t mesh(bounds, b_types);
            cell_t const c1(1);
            cell_t const f1(0);
            Event const exp_ev1(Event::reflect);
            Event ev1 = decide_boundary_event(mesh,c1,f1);
            passed = exp_ev1 == ev1 and passed;
            // transmit: high face of cell 1
            cell_t const c2(1);
            cell_t const f2(1);
            Event const exp_ev2(Event::cell_high_x_boundary);
            Event ev2 = decide_boundary_event(mesh,c2,f2);
            passed = exp_ev2 == ev2 and passed;
            if(!passed) std::cout << "FAILED: " << __LINE__ << std::endl;
            // transmit: low face of cell 2
            cell_t const c3(2);
            cell_t const f3(0);
            Event const exp_ev3(Event::cell_low_x_boundary);
            Event ev3 = decide_boundary_event(mesh,c3,f3);
            passed = exp_ev3 == ev3 and passed;
            if(!passed) std::cout << "FAILED: " << __LINE__ << std::endl;
            // transmit: low face of cell 10
            cell_t const c4(10);
            cell_t const f4(0);
            Event const exp_ev4(Event::cell_low_x_boundary);
            Event ev4 = decide_boundary_event(mesh,c4,f4);
            passed = exp_ev4 == ev4 and passed;
            if(!passed) std::cout << "FAILED: " << __LINE__ << std::endl;
            // escape: high face of cell 10
            cell_t const c5(10);
            cell_t const f5(1);
            Event const exp_ev5(Event::escape);
            Event ev5 = decide_boundary_event(mesh,c5,f5);
            passed = exp_ev5 == ev5 and passed;
            if(!passed) std::cout << "FAILED: " << event_name(ev5) << ":"
                                  << __LINE__ << std::endl;

            return passed;
        } // test_2


        /** Reworks the Python test_0 of decide_event. */
        bool test_3()
        {
            bool passed(true);

            typedef mesh_t::vb vb;
            typedef mesh_t::vbd vbd;

            using nut::decide_scatter_event;
            using nut::Event;
            using nut::event_name;

            rho_t rho1(rho);
            T_t   T1(T);
            rho1.rho_p[0] = 1.0;

            OpB op(rho1,T1);
            fp_t const rns[] = {0.30897681610609407,
                                0.92436920169905545,
                                0.21932404923057958};
            BRNG rng(rns,3);
            fp_t const x = 0.5, omega = 1.0, e = 5.0, t = 1.0, wt = 1.0;
            cell_t const cell = 1;
            Species const s(nut::nu_e);

            // generate uniform mesh
            cell_t n_cells(10);
            vb bounds(n_cells+1);
            vbd b_types(n_cells+1);
            std::generate(b_types.begin(),b_types.end(),gen_bdy_types(n_cells+1));
            std::generate(bounds.begin(),bounds.end(),gen_bounds(1.0));

            // reflect from innermost face
            mesh_t mesh(bounds, b_types);

            p_t p({x},{omega},e,t,wt,cell,rng,s);

            nut::event_n_dist e_n_d = decide_event(p,mesh,op,vel0s);

            Event const event_exp = Event::cell_high_x_boundary;
            Event const event = e_n_d.first;
            passed = event == event_exp;
            if(!passed)
            {
                std::cout << "event was " << event_name(event) << ", should have been "
                          << event_name(nut::Event::nucleon_abs) << std::endl;
            }
            return passed;
        } // test_3


        /** Reworks the Python test_1 of decide_event. */
        bool test_4()
        {
            bool passed(true);

            typedef mesh_t::vb vb;
            typedef mesh_t::vbd vbd;

            using nut::decide_scatter_event;
            using nut::Event;
            using nut::event_name;

            rho_t rho1(rho);
            T_t   T1(T);
            rho1.rho_p[0] = 1.0;

            OpB op(rho1,T1);
            fp_t const rns[] = {0,0,0,0,0,0,0,0,0,0};
            BRNG rng(rns,3);
            fp_t const x = 0.5, omega = 1.0, e = 5.0, t = 1.0, wt = 1.0;
            cell_t const cell = 1;
            Species const s(nut::nu_e);

            // generate uniform mesh
            cell_t n_cells(10);
            vb bounds(n_cells+1);
            vbd b_types(n_cells+1);
            std::generate(b_types.begin(),b_types.end(),gen_bdy_types(n_cells+1));
            std::generate(bounds.begin(),bounds.end(),gen_bounds(1.0));

            // reflect from innermost face
            mesh_t mesh(bounds, b_types);

            p_t p({x},{omega},e,t,wt,cell,rng,s);

            for(size_t i = 0; i < 9; ++i)
            {
                nut::event_n_dist e_n_d = decide_event(p,mesh,op,vel0s);
                std::cout << event_name(e_n_d.first)
                          << ", d = " << e_n_d.second << std::endl;
                p.x += e_n_d.second;
                p.cell += 1;
            }

            nut::event_n_dist e_n_d = decide_event(p,mesh,op,vel0s);

            Event const event_exp = Event::escape;
            Event const event = e_n_d.first;
            passed = event == event_exp;
            if(!passed)
            {
                std::cout << "event was " << event_name(event) << ", AKA "
                          << event_name(event)
                          << ", should have been "
                          << event_name( event_exp) << std::endl;
            }
            geom_t const d_exp = 1.0;
            passed = e_n_d.second == d_exp and passed;
            return passed;
        } // test_4


        /* This test uses the setup of Python test_3 of decide_event to test
         * decide_scatter_event. We use the same random seeds, starting at the
         * third seed, since the test of decide_event burns the first one for
         * sampling the exponential to come up with d_collision, and the second
         * is used to choose the particle interaction. This test
         * should result in nucleon elastic scatter. */
        bool test_5()
        {
            bool passed(true);

            using nut::decide_scatter_event;

            rho_t rho1(rho);
            T_t   T1(T);
            rho1.rho_p[0] = 1e14/pmg;

            OpB op(rho1,T1);
            // from Python test_2 of decide_event
            fp_t const rns[] = {// 0.21932404923057958, 1st used for d_coll
                                // 0.20867489035315723,
                                0.91525579001682567};
            BRNG rng(rns,3);
            fp_t const e     = 5.0;
            cell_t const cell  = 1;
            Species const s(nut::nu_e);

            Event event = decide_scatter_event(rng,e,cell,op,s);

            Event event_exp = Event::nucleon_elastic_scatter;

            passed = event == event_exp;
            if(!passed)
            {
                std::cout << "event was " << event_name(event)
                          << ", should have been "
                          << nut::event_name(Event::nucleon_elastic_scatter) << std::endl;
            }
            return passed;
        } // test_5


        /* Stack the deck to use electron scatter */
        bool test_6()
        {
            bool passed(true);

            using nut::Event;

            using nut::decide_scatter_event;

            rho_t rho1(rho);
            T_t   T1(T);
            rho1.rho_p[0] = nut::tiny;
            rho1.rho_e_minus[0] = 1e14/pmg;
            T1.T_e_minus[0] = 1.0;

            OpB op(rho1,T1);
            // from Python test_2 of decide_event
            fp_t const rns[] = {0.9,0.1};
            BRNG rng(rns,3);
            fp_t const e     = 5.0;
            cell_t const cell  = 1;
            Species const s(nut::nu_e);

            Event event = decide_scatter_event(rng,e,cell,op,s);

            Event event_exp = Event::electron_scatter;

            passed = event == event_exp;
            if(!passed)
            {
                std::cout << "event was " << event_name(event)
                          << ", should have been "
                          << event_name(event_exp) << std::endl;
            }
            return passed;
        } // test_6


        /* This test uses the setup of Python test_2 of decide_event to test
         * decide_event. We use the same random seeds, skipping the
         * second seed, since the test of decide_event burns the first one for
         * sampling the exponential to come up with d_collision. In the Python code,
         * the second seed is used to select which type of particle to scatter from;
         * NuT avoids this in keeping with McPhD. */
        bool test_7()
        {
            bool passed(true);

            using nut::Event;

            typedef mesh_t::vb vb;
            typedef mesh_t::vbd vbd;
            // generate uniform mesh
            cell_t n_cells(10);
            vb bounds(n_cells+1);
            vbd b_types(n_cells+1);
            std::generate(b_types.begin(),b_types.end(),gen_bdy_types(n_cells+1));
            std::generate(bounds.begin(),bounds.end(),gen_bounds(1.0e6));
            mesh_t mesh(bounds, b_types);

            // generate opacity
            rho_t rho1(rho);
            T_t   T1(T);
            rho1.rho_p[0] = 1e14/pmg;

            OpB op(rho1,T1);

            // create particle
            // from Python test_2 of decide_event
            fp_t const rns[] = {0.30897681610609407, // 1st seed burned for d_coll
                                // 0.92436920169905545,
                                0.21932404923057958};
            BRNG rng(rns,3);
            fp_t const e     = 5.0;
            cell_t const cell  = 1;
            Species const s(nut::nu_e);
            geom_t const x = 0.5;
            geom_t const omega = 1.0;
            fp_t const t = 100.0;
            fp_t const wt = 1.0;
            p_t p({x},{omega},e,t,wt,cell,rng,s);


            nut::event_n_dist e_n_d = decide_event(p,mesh,op,vel0s);

            Event const & event = e_n_d.first;
            Event event_exp = Event::nucleon_abs;

            passed = event == event_exp;
            if(!passed)
            {
                std::cout << "event was " << event_name(event)
                          << ", should have been "
                          << event_name(nut::Event::nucleon_abs) << std::endl;
            }
            geom_t const d_exp = 7343.827;
            geom_t const epsilon =  0.001;
            geom_t const d = e_n_d.second;
            passed = nut::soft_equiv(d,d_exp,epsilon) && passed;
            if(!passed)
            {
                std::cout << "distance to event was " << std::setprecision(15) <<  d
                          << "expected distance was " << std::setprecision(15) << d_exp
                          << std::endl;
            }

            return passed;
        } // test_7


        /* This test uses the setup of Python test_3 of decide_event to test
         * decide_scatter_event. We use the same random seeds, skipping the
         * second seed, since the test of decide_event burns the first one for
         * sampling the exponential to come up with d_collision. See comment to
         * test_7 for more on the rationale. This test
         * should result in nucleon elastic scatter. */
        bool test_8()
        {
            bool passed(true);

            using nut::Event;

            typedef mesh_t::vb vb;
            typedef mesh_t::vbd vbd;
            // generate uniform mesh
            cell_t n_cells(10);
            vb bounds(n_cells+1);
            vbd b_types(n_cells+1);
            std::generate(b_types.begin(),b_types.end(),gen_bdy_types(n_cells+1));
            std::generate(bounds.begin(),bounds.end(),gen_bounds(1.0e6));
            mesh_t mesh(bounds, b_types);

            rho_t rho1(rho);
            T_t   T1(T);
            rho1.rho_p[0] = 1e14/pmg;

            OpB op(rho1,T1);
            // from Python test_2 of decide_event
            fp_t const rns[] = {0.21932404923057958, // 1st used for d_coll
                                0.20867489035315723,
                                0.91525579001682567};
            BRNG rng(rns,3);
            fp_t const e     = 5.0;
            cell_t const cell  = 1;
            Species const s(nut::nu_e);
            geom_t const x = 0.5;
            geom_t const omega = 1.0;
            fp_t const t = 100.0;
            fp_t const wt = 1.0;
            p_t p({x},{omega},e,t,wt,cell,rng,s);

            nut::event_n_dist e_n_d  = nut::decide_event(p,mesh,op,vel0s);

            Event const event = e_n_d.first;
            Event const event_exp = Event::nucleon_elastic_scatter;

            passed = event == event_exp;
            if(!passed)
            {
                std::cout << "event was " << event_name(event)
                          << ", should have been "
                          << event_name(Event::nucleon_elastic_scatter) << std::endl;
            }
            geom_t const d_exp = 9486.756;
            geom_t const epsilon =  0.001;
            geom_t const d = e_n_d.second;
            passed = nut::soft_equiv(d,d_exp,epsilon) && passed;
            if(!passed)
            {
                std::cout << "distance to event was " << d
                          << "expected distance was " << d_exp << std::endl;
            }
            return passed;
        } // test_8


        /* Stack the deck to use electron scatter */
        bool test_9()
        {
            bool passed(true);

            using nut::Event;

            typedef mesh_t::vb vb;
            typedef mesh_t::vbd vbd;
            // generate uniform mesh
            cell_t n_cells(10);
            vb bounds(n_cells+1);
            vbd b_types(n_cells+1);
            std::generate(b_types.begin(),b_types.end(),gen_bdy_types(n_cells+1));
            std::generate(bounds.begin(),bounds.end(),gen_bounds(1.0e6));
            mesh_t mesh(bounds, b_types);

            rho_t rho1(rho);
            T_t   T1(T);
            rho1.rho_p[0] = nut::tiny;
            rho1.rho_e_minus[0] = 1e14/pmg;
            T1.T_e_minus[0] = 1.0;

            OpB op(rho1,T1);
            // from Python test_2 of decide_event
            fp_t const rns[] = {0.9,0.9,0.1};
            BRNG rng(rns,3);
            fp_t const e     = 5.0;
            cell_t const cell  = 1;
            Species const s(nut::nu_e);
            geom_t const x = 0.5;
            geom_t const omega = 1.0;
            fp_t const t = 1000.0;
            fp_t const wt = 1.0;
            p_t p({x},{omega},e,t,wt,cell,rng,s);

            nut::event_n_dist e_n_d  = nut::decide_event(p,mesh,op,vel0s);

            Event const event = e_n_d.first;

            Event const event_exp = Event::electron_scatter;

            passed = event == event_exp;
            if(!passed)
            {
                std::cout << "event was " << event_name(event)
                          << ", should have been "
                          << event_name(event_exp) << std::endl;
            }
            fp_t const d = e_n_d.second;
            fp_t const d_exp = 38310.5;
            geom_t const epsilon =  0.001;
            passed = nut::soft_equiv(d,d_exp,epsilon) && passed;
            if(!passed)
            {
                std::cout << "distance to event was " << d
                          << "expected distance was " << d_exp << std::endl;
            }
            return passed;
        } // test_9


        // define additional tests here.

    } // decide_event_tests::

} // Nut_Test::



// version
// $Id$

// End of file
