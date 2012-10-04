//T. M. Kelley (c) 2011 LANS LLC

#include "RNG.hh"
#include "nut_Test_RNG.hh"
#include "test_aux.hh"
#include <vector>
#include <algorithm>
#include <functional>

namespace Nut_Test
{
    namespace RNG_tests
    {
        // target describes the code being tested
        char target[] = "RNGs";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "init. & inst. Buffer_RNG<float>";
        bool test_2();
        char aspect2[] = "use Buffer_RNG<float>";
        bool test_3();
        char aspect3[] = "use Buffer_RNG<double>";
        bool test_4();
        char aspect4[] = "use Buffer_RNG<double> twice: look for \"Buffer_RNG rolled over from 3 to 0\" in output.";

        bool test_5();
        char aspect5[] = "init. & inst. LCG_RNG";
        bool test_6();
        char aspect6[] = "use LCG_RNG: print a few values";
        bool test_7();
        char aspect7[] = "use LCG_RNG: draw 1e6 values, check for 0 < x < 1";

        bool test_8();
        char aspect8[] = "init. & inst. MLCG";
        bool test_9();
        char aspect9[] = "use MLCG: print a few values (seed = 42)";
        bool test_10();
        char aspect10[] = "use MLCG: draw 1e6 values, check for 0 < x < 1";
        bool test_11();
        char aspect11[] = "use MLCG: test splitting";

        bool test_12();
        char aspect12[] = "init. & inst. Philox4x32";

        bool test_13();
        char aspect13[] = "print a few values";

        bool test_14();
        char aspect14[] = "compatible with Haskell/McPhD usage";
    }

    bool test_RNG()
    {
        using namespace RNG_tests;
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

        // call additional tests here.

        return passed1 and passed2 and passed3 and passed4 and passed5 and 
            passed6 and passed7 and passed8 and passed9 and passed11 and passed10
            and passed12 and passed13;
    }

    namespace RNG_tests
    {
        // make rng callable by STL algos
        template <typename fp_t, typename rng_t>
        class rng_adaptor
        {
        public:
            explicit rng_adaptor(rng_t & rng)
                : m_rng(rng){}
            fp_t operator()(){ return m_rng.random();}
            rng_t & m_rng;
            void dump_state(std::ostream & o){return m_rng.dump_state(o);}
        };

        bool test_1()
        {
            bool passed(true);

            typedef float fp_t;
            size_t const vals_sz = 4;
            fp_t const in_vals[vals_sz] = {0.,0.11,0.22,0.33};

            typedef nut::Buffer_RNG<fp_t> rng_t;
            rng_t rng(&in_vals[0],vals_sz);

            return passed;
        } // test_1


        bool test_2()
        {
            bool passed(true);

            typedef float fp_t;
            size_t const vals_sz = 4;
            fp_t const in_vals[vals_sz] = {0.,0.11,0.22,0.33};

            typedef nut::Buffer_RNG<fp_t> rng_t;
            rng_t rng_base(&in_vals[0],vals_sz);
            rng_adaptor<fp_t,rng_t> rng(rng_base);

            std::vector<fp_t> out_vals(vals_sz);

            std::generate(out_vals.begin(),out_vals.end(),rng);
            passed = passed and 
                std::equal(&in_vals[0],&in_vals[vals_sz],out_vals.begin());

            return passed;
        } // test_2


        bool test_3()
        {
            bool passed(true);

            typedef double fp_t;
            size_t const vals_sz = 4;
            fp_t const in_vals[vals_sz] = {0.,0.11,0.22,0.33};

            typedef nut::Buffer_RNG<fp_t> rng_t;
            rng_t rng_base(&in_vals[0],vals_sz);
            rng_adaptor<fp_t,rng_t> rng(rng_base);

            std::vector<fp_t> out_vals(vals_sz);

            std::generate(out_vals.begin(),out_vals.end(),rng);
            passed = passed and 
                std::equal(&in_vals[0],&in_vals[vals_sz],out_vals.begin());

            return passed;
        } // test_3


        bool test_4()
        {
            bool passed(true);

            typedef double fp_t;
            size_t const vals_sz = 4;
            fp_t const in_vals[vals_sz] = {0.,0.11,0.22,0.33};

            typedef nut::Buffer_RNG<fp_t> rng_t;
            rng_t rng_base(&in_vals[0],vals_sz);
            rng_adaptor<fp_t,rng_t> rng(rng_base);

            std::vector<fp_t> out_vals(2*vals_sz);

            std::generate(out_vals.begin(),out_vals.end(),rng);
            passed = passed and 
                std::equal(&in_vals[0],&in_vals[vals_sz],out_vals.begin());
            passed = passed and 
                std::equal(&in_vals[0],&in_vals[vals_sz],
                           out_vals.begin()+vals_sz);

            return passed;
        } // test_4


        bool test_5()
        {
            bool passed(true);

            typedef nut::LCG_RNG rng_t;
            int32_t const seed(-120211);
            rng_t rng(seed);

            return passed;
        } // test_5

        bool test_6()
        {
            bool passed(true);

            typedef nut::LCG_RNG rng_t;
            int32_t const seed(-120211);
            rng_t rng_base(seed);
            rng_adaptor<double,rng_t> rng(rng_base);

            size_t const vals_sz(4);

            std::vector<double> out_vals(vals_sz);
                        
            std::generate(out_vals.begin(),out_vals.end(),rng);

            std::copy(out_vals.begin(),out_vals.end(),
                      std::ostream_iterator<double>(std::cout,","));
            std::cout << std::endl;

            return passed;
        } // test_6

        namespace 
        {
            bool lessThanEq0(double const d){return d <= 0.0;}
            bool greaterThanEq1(double const d){return d >= 1.0;}
            std::string err_lte0("some values were < 0.0: FAIL");
            std::string err_gte1("some values were > 1.0: FAIL");

            /* check an iterator range against a unary predicate, print a 
             * string if predicate fails anywhere. */
            template <class it_t,class pred_t>
            bool checkPred(it_t first, it_t last, pred_t pred, std::string errstr)
            {
                it_t it = find_if(first,last,pred);
                bool passed = it == last;
                if(!passed)
                {
                    std::cerr << errstr << std::endl;
                }
                return passed;
            }

        } // anonymous

        bool test_7()
        {
            bool passed(true);

            typedef nut::LCG_RNG rng_t;
            int32_t const seed(-120211);
            rng_t rng_base(seed);
            rng_adaptor<double,rng_t> rng(rng_base);

            size_t const vals_sz(1000000);
            std::vector<double> out_vals(vals_sz);
            std::generate(out_vals.begin(),out_vals.end(),rng);

            // check each value greater than 0 and less than 1
            bool lt0_passed = checkPred(out_vals.begin(),out_vals.end(),
                                        lessThanEq0,err_lte0);
            passed = passed and lt0_passed;

            bool gt1_passed = checkPred(out_vals.begin(),out_vals.end(),
                                        greaterThanEq1,err_gte1);
            passed = passed and gt1_passed;

            return passed;
        } // test_7

       
        bool test_8()
        {
            bool passed(true);

            typedef nut::MLCG rng_t;
            int32_t const seed(42);
            rng_t rng(seed);

            return passed;
        } // test_8

        bool test_9()
        {
            bool passed(true);

            typedef nut::MLCG rng_t;
            int32_t const seed(42);
            rng_t rng_base(seed);
            rng_adaptor<double,rng_t> rng(rng_base);

            size_t const vals_sz(4);

            std::vector<double> out_vals(vals_sz);
                        
            std::generate(out_vals.begin(),out_vals.end(),rng);

            std::copy(out_vals.begin(),out_vals.end(),
                      std::ostream_iterator<double>(std::cout,","));
            std::cout << std::endl;

            return passed;
        } // test_9

        bool test_10()
        {
            bool passed(true);

            typedef nut::MLCG rng_t;
            int32_t const seed(42);
            rng_t rng_base(seed);
            rng_adaptor<double,rng_t> rng(rng_base);

            size_t const vals_sz(1000000);
            std::vector<double> out_vals(vals_sz);
            std::generate(out_vals.begin(),out_vals.end(),rng);

            // check each value greater than 0 and less than 1
            bool lt0_passed = checkPred(out_vals.begin(),out_vals.end(),
                                        lessThanEq0,err_lte0);
            passed = passed and lt0_passed;

            bool gt1_passed = checkPred(out_vals.begin(),out_vals.end(),
                                        greaterThanEq1,err_gte1);
            passed = passed and gt1_passed;

            return passed;
        } // test_10

        bool test_11()
        {
            bool passed(true);

            typedef nut::MLCG rng_t;
            int32_t const seed(42);
            rng_t rng_base1(seed);
            rng_adaptor<double,rng_t> rng1(rng_base1);
            
            nut::MLCG::new_gens ngens = rng_base1.split();
            rng_adaptor<double,rng_t> rng2(ngens.first);
            rng_adaptor<double,rng_t> rng3(ngens.second);
            
            // check old generator
            {
                std::cout << "after split, state of RNG 1: ";
                rng1.dump_state(std::cout);
                size_t const vals_sz(4);
                std::vector<double> out_vals(vals_sz);
                std::generate(out_vals.begin(),out_vals.end(),rng1);
                std::cout << "; some values: ";
                std::copy(out_vals.begin(),out_vals.end(),
                          std::ostream_iterator<double>(std::cout,","));
                std::cout << std::endl;
            }
            // check new generators
            {
                std::cout << "after split, state of RNG 2: ";
                rng2.dump_state(std::cout);
                size_t const vals_sz(4);
                std::vector<double> out_vals(vals_sz);
                std::generate(out_vals.begin(),out_vals.end(),rng2);
                std::cout << "; some values: ";
                std::copy(out_vals.begin(),out_vals.end(),
                          std::ostream_iterator<double>(std::cout,","));
                std::cout << std::endl;
            }

            {
                std::cout << "after split, state of RNG 3: ";
                rng3.dump_state(std::cout);
                size_t const vals_sz(4);
                std::vector<double> out_vals(vals_sz);
                std::generate(out_vals.begin(),out_vals.end(),rng3);
                std::cout << "; some values: ";
                std::copy(out_vals.begin(),out_vals.end(),
                          std::ostream_iterator<double>(std::cout,","));
                std::cout << std::endl;
            }

            return passed;
        } // test_11

        bool test_12(){
            bool passed(true);

            nut::Philox4x32_RNG::ctr_t c;
            nut::Philox4x32_RNG::key_t k;
            c[0] = c[1] = c[2] = 0;
            c[3] = 1;
            k[0] = 0xdeadbeef;
            k[1] = 0xcafecafe;

            nut::Philox4x32_RNG g(c,k);

            return passed;
        } // test_12

        bool test_13()
        {
            bool passed(true);

            typedef nut::Philox4x32_RNG rng_t;
            nut::Philox4x32_RNG::ctr_t c;
            nut::Philox4x32_RNG::key_t k;
            c[0] = c[1] = c[2] = 0;
            c[3] = 1;
            k[0] = 0xdeadbeef;
            k[1] = 0xcafecafe;
            rng_t rng_base(c,k);
            rng_adaptor<double,rng_t> rng(rng_base);

            size_t const vals_sz(10);

            std::vector<double> out_vals(vals_sz);
                        
            std::generate(out_vals.begin(),out_vals.end(),rng);

            std::copy(out_vals.begin(),out_vals.end(),
                      std::ostream_iterator<double>(std::cout,","));
            std::cout << std::endl;

            return passed;
        } // test_13


        bool test_14()
        {
            bool passed(true);

            typedef nut::Philox4x32_RNG rng_t;
            nut::Philox4x32_RNG::ctr_t c;
            nut::Philox4x32_RNG::key_t k;
            c[0] = c[1] = c[2] = 0;
            c[3] = 1;
            k[0] = 0xdeadbeef;
            k[1] = 0xcafecafe;
            rng_t rng_base(c,k);
            rng_adaptor<double,rng_t> rng(rng_base);

            size_t const vals_sz(12);

            std::vector<double> out_vals(vals_sz);
                        
            std::generate(out_vals.begin(),out_vals.end(),rng);

            // These values generated via the function randoms in the 
            // Philo2 module of McPhD/basic, using same counter & key
            // as this test.
            std::vector<double> exp_vals = { 
                0.9848038759115404,0.8800535828289092,0.8015360604392531,
                0.2025851484105382,0.27783947633403117,0.634825162264666,
                0.49563039019191557,0.769547123123007,0.7025144994130103,
                0.6585216467165459,0.677582395892678,0.8577909275075531
            };

            passed = passed and std::equal(out_vals.begin(),out_vals.end(),
                                           exp_vals.begin());

            // std::copy(out_vals.begin(),out_vals.end(),
            //           std::ostream_iterator<double>(std::cout,","));
            // std::cout << std::endl;

            return passed;
        } // test_14

        // define additional tests here.

    } // RNG_tests::

} // Nut_Test::



// version
// $Id$

// End of file
