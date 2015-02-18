//T. M. Kelley (c) 2011 LANS LLC

#include "Planck.hh"
#include "RNG.hh"
#include "nut_Test_Planck.hh"
#include "test_aux.hh"
#include "expect.hh"
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>


namespace Nut_Test
{
    namespace Planck_tests
    {
        // target describes the code being tested
        char target[] = "Planck-ish functions";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "gen_power_law_energies: 1000 tests via Mathematica";
        bool test_2();
        char aspect2[] = "pwr_law_reject: 1000 tests via Mathematica";
        bool test_3();
        char aspect3[] = "gen_pwr_law specialized to alpha=2: 1000 tests via Mathematica";
        bool test_4();
        char aspect4[] = "prw_law_reject specialized to alpha=2: 1000 tests via Mathematica";

    }

    bool test_Planck()
    {
        using namespace Planck_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        bool passed2 = test( target, aspect2, test_2);

        bool passed3 = test( target, aspect3, test_3);

        bool passed4 = test( target, aspect4, test_4);

        // call additional tests here.

        return passed1 and passed2 and passed3 and passed4;
    }

    namespace Planck_tests
    {
        using test_aux::soft_eq_bound_tol;

        /* Read in 1000 test cases from a file, run them, check the results.
         * The test cases were generated using Mathematica and random sampling.
         * Here, we use a Buffer_RNG object loaded with the random values from
         * Mathematica as inputs and verify the outputs. */
        namespace
        {
            using std::bind;
            using namespace std::placeholders;  // for _1, _2, etc in bind

            /*\param rns: vector to hold random numbers
             *\param nrgs: vector to hold sampled energies
             *\param rejfs: vector to hold rejection function values */
            template <typename vt>
            bool read_input_t1(vt & all_rns, vt & x_rns, vt & nrgs, vt & rejfs)
            {
                typedef typename vt::value_type geom_t;
                std::string const fname("generate-power-law-energies-tests-cxx.txt");
                std::ifstream inf(fname.c_str());
                if(!inf.good())
                {
                    std::cerr << "Planck test_1 cannot read file " << fname
                              << std::endl;
                    return false;
                }
                std::string trash;
                std::getline(inf,trash);

                uint32_t count(0);
                while(1)
                {
                    geom_t r1,r2,nrg,rejf;
                    inf >> r1 >> r2 >> nrg >> rejf;
                    if(inf.eof()) break;
                    x_rns.push_back(r1);
                    all_rns.push_back(r1);
                    all_rns.push_back(r2);
                    nrgs.push_back(nrg);
                    rejfs.push_back(rejf);
                    count++;
                }
                std::cout << "read " << count << " input lines, will run that "
                          << " many problems."
                          << std::endl;
                inf.close();
                return true;
            } // read_input_t1

            /** - Log(x) */
            template <typename fp_t>
            fp_t neg_log(fp_t const x){ return -std::log(x);}

        } // anonymous::


        bool test_1()
        {
            bool passed(true);

            using nut::gen_power_law_energy;
            using namespace test_aux;

            typedef double fp_t;
            typedef std::vector<fp_t> vf;
            vf rns, x_rns, es_exp, rejs_exp;
            fp_t const ebar(2.0);
            fp_t const alpha(2.0);

            // load values
            read_input_t1(rns,x_rns,es_exp,rejs_exp);

            // construct Buffer_RNG from rns
            typedef nut::Buffer_RNG<fp_t> rng_t;
            rng_t rng(&rns[0],rns.size());

            vf es(es_exp.size());
            std::generate(es.begin(),es.end(),
                          std::bind( nut::gen_power_law_energy<rng_t,fp_t>,
                                          alpha,ebar,rng));

            bool es_passed = check_same_verb(&es,&es_exp,comp_verb<fp_t>("energies",1e-15));

            passed = es_passed and passed;

            return passed;
        } // test_1


        bool test_2()
        {
            bool passed(true);

            using nut::gen_power_law_energy;

            typedef double fp_t;
            typedef std::vector<fp_t> vf;
            vf rns, x_rns, es_exp, rejs_exp;

            // load values
            read_input_t1(rns,x_rns,es_exp,rejs_exp);

            fp_t const e_a(0.1353352832366127); // exp(-2)
            fp_t const alpha(2.0);

            // divide energies by ebar (= 2.0 MeV) to get xs.
            vf xs(es_exp.size());
            std::transform(x_rns.begin(),x_rns.end(),xs.begin(),
                           neg_log<fp_t>);

            // construct Buffer_RNG from rns
            typedef nut::Buffer_RNG<fp_t> rng_t;
            rng_t rng(&rns[0],rns.size());

            vf rejs(rejs_exp.size());
            // std::transform(xs.begin(),xs.end(),rejs.begin(),
            //                pl_rej<fp_t>(alpha,e_a));
            std::transform(xs.begin(),xs.end(),rejs.begin(),
                           bind(nut::pwr_law_reject<fp_t>,alpha,e_a,_1));

            // accept a relative error less than 10^(-15)
            // soft_eq_bound_tol<fp_t> seq(1e-15);
            // bool rejs_passed = std::equal(rejs.begin(),rejs.end(),
            //                               rejs_exp.begin(),seq);
            soft_eq_bound_tol<fp_t> seq(1e-15);
            bool rejs_passed = std::equal(rejs.begin(),rejs.end(),rejs_exp.begin(),
                                          bind(nut::soft_equiv<fp_t>,_1,_2,1e-15));

            passed = rejs_passed and passed;

            return passed;
        } // test_2


        bool test_3()
        {
            bool passed(true);

            using nut::gen_power_law_energy;
            using namespace test_aux;

            typedef double fp_t;
            typedef std::vector<fp_t> vf;
            vf rns, x_rns, es_exp, rejs_exp;
            fp_t const ebar(2.0);

            // load values
            read_input_t1(rns,x_rns,es_exp,rejs_exp);

            // construct Buffer_RNG from rns
            typedef nut::Buffer_RNG<fp_t> rng_t;
            rng_t rng(&rns[0],rns.size());

            vf es(es_exp.size());
            // std::generate(es.begin(),es.end(),pl_gen_a2<fp_t,rng_t>(ebar,rng));
            std::generate(es.begin(),es.end(),
                          bind(nut::gen_power_law_energy_alpha2<rng_t,fp_t>,ebar,rng));

            bool es_passed = check_same_verb(&es,&es_exp,comp_verb<fp_t>("energies",1e-15));

            passed = es_passed and passed;

            return passed;
        } // test_3


        bool test_4()
        {
            bool passed(true);

            using nut::gen_power_law_energy;

            typedef double fp_t;
            typedef std::vector<fp_t> vf;
            vf rns, x_rns, es_exp, rejs_exp;

            // load values
            read_input_t1(rns,x_rns,es_exp,rejs_exp);

            // divide energies by ebar (= 2.0 MeV) to get xs.
            vf xs(es_exp.size());
            std::transform(x_rns.begin(),x_rns.end(),xs.begin(),
                           neg_log<fp_t>);

            // construct Buffer_RNG from rns
            typedef nut::Buffer_RNG<fp_t> rng_t;
            rng_t rng(&rns[0],rns.size());

            vf rejs(rejs_exp.size());
            // std::transform(xs.begin(),xs.end(),rejs.begin(),
            //                pl_rej_a2<fp_t>(e_a));
            std::transform(xs.begin(),xs.end(),rejs.begin(),
                           bind(nut::pwr_law_reject_alpha2<fp_t>,_1));

            // accept a relative error less than 10^(-15)
            soft_eq_bound_tol<fp_t> seq(1e-15);
            // bool rejs_passed =
            //     std::equal(rejs.begin(),rejs.end(),rejs_exp.begin(),seq);
            bool rejs_passed =
                std::equal(rejs.begin(),rejs.end(),rejs_exp.begin(),
                           bind(nut::soft_equiv<fp_t>,_1,_2,1e-15));

            passed = rejs_passed and passed;

            return passed;
        } // test_4


        // define additional tests above here.


        // functors for use with older compilers without ::bind:
        namespace
        {

            // run the power law generator, wd be nice to replace with boost::bind
            template <typename fp_t, typename rng_t>
            struct pl_gen
            {
                pl_gen(fp_t const alpha, fp_t const ebar, rng_t & rng)
                    : m_a(alpha),m_e(ebar),m_rng(rng){}
                fp_t operator()(){return gen_power_law_energy(m_a,m_e,m_rng);}
                fp_t const m_a, m_e;    rng_t & m_rng;
            }; // pl_gen


            // run the power law rejection function,
            // wd be nice to replace with boost::bind
            template <typename fp_t>
            struct pl_rej
            {
                pl_rej(fp_t const alpha, fp_t const e_a)
                    : m_a(alpha),m_e(e_a){}
                fp_t operator()(fp_t const x){return nut::pwr_law_reject(m_a,m_e,x);}
                fp_t const m_a, m_e;
            }; // pl_rej


            // run the power law generator specialized to alpha = 2,
            // wd be nice to replace with boost::bind
            template <typename fp_t, typename rng_t>
            struct pl_gen_a2
            {
                pl_gen_a2(fp_t const ebar, rng_t & rng)
                    : m_e(ebar),m_rng(rng){}
                fp_t operator()(){return gen_power_law_energy_alpha2(m_e,m_rng);}
                fp_t const m_e;    rng_t & m_rng;
            }; // pl_gen


            // run the power law rejection function, specialized to alpha = 2
            // wd be nice to replace with boost::bind
            template <typename fp_t>
            struct pl_rej_a2
            {
                pl_rej_a2(fp_t const e_a)
                    : m_e(e_a){}
                fp_t operator()(fp_t const x){
                    return nut::pwr_law_reject_alpha2<fp_t>(x);}
                fp_t const m_e;
            }; // pl_rej_a2


        } // anonymous::



    } // Planck_tests::

} // Nut_Test::



// version
// $Id$

// End of file
