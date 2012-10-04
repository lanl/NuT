//T. M. Kelley (c) 2011 LANS LLC

#include "nut_Test_partition.hh"
#include "test_aux.hh"
#include "partition.hh"
#include "mpi_helpers.hh"
#include "types.hh"
#include <numeric>     // for partial_sum
#include <iostream>
#include <iterator>
#include <random>


namespace Nut_Test
{
    namespace partition_tests
    {
        // target describes the code being tested
        char target[] = "partition";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "first particle id is 0";

        bool test_2();
        char aspect2[] = "last particle id is n_particles - 1";

        bool test_3();
        char aspect3[] = "each chunk end > chunk begin (except maybe last)";

        bool test_4();
        char aspect4[] = "each chunk start > previous chunk begin";        

        bool test_5();
        char aspect5[] = "chunks contiguous: c2.start == c1.end + 1";        

        bool test_6();
        char aspect6[] = "chunks disjoint";        

        bool test_7();
        char aspect7[] = "take_n_particles checks chkSz";        
    }

    bool test_partition()
    {
        using namespace partition_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        bool passed2 = test( target, aspect2, test_2);

        bool passed3 = test( target, aspect3, test_3);

        bool passed4 = test( target, aspect4, test_4);

        bool passed5 = test( target, aspect5, test_5);

        bool passed6 = test( target, aspect6, test_6);

        bool passed7 = test( target, aspect7, test_7);

        // call additional tests here.

        return 
            passed1 
            and passed2
            and passed3
            and passed4
            and passed5
            and passed6
            and passed7
            ;
    }

    namespace partition_tests
    {
        using nut::geom_t;
        using nut::id_t;
        using nut::Rank;
        using nut::CommSz;
        using nut::PtclId;
        // using nut::ps_per_rank;
        // using nut::gen_ptcl_ids;
        // using nut::ids_per_cell;

        typedef nut::src_stats_t<geom_t,id_t> stats_t;
        typedef std::shared_ptr<stats_t>      sp_stats;
        typedef stats_t::sz_t                 sz_t;
        using nut::Chunker;
        typedef Chunker<stats_t>              Chker;
        typedef Chker::ChunkVec               ChunkVec;
        typedef Chker::StatsVec               StatsVec;

        typedef std::minstd_rand gen_t;
        gen_t mkGen()
        {
            std::random_device rd;
            return gen_t(rd());
        }
        
        /** choose a number between min and max (inclusive-ish) */
        template <typename t>
        t choose(t min, t max, gen_t g)
        {
            double urd = std::generate_canonical<double,64>(g);
            return min + t(urd*(max-min));
        }

        template <typename t>
        struct Chooser
        {
            Chooser(t min, t max,gen_t gn): mn(min),mx(max),g(gn){}
            t operator()(){return choose(mn,mx,g);}
            t mn, mx;
            gen_t g;
        };

        struct Test_Props
        {
            uint32_t n;
            explicit Test_Props(uint32_t n_ts = 100) : n(n_ts) {}
        }; // Test_Props

        Test_Props t5(5u);
        Test_Props t100;

        struct StatsNSz
        {
            sp_stats ss;
            stats_t::sz_t chkSz;
            explicit StatsNSz(sp_stats ss_,
                              stats_t::sz_t sz_) 
                : ss(ss_),chkSz(sz_){}
            // for testing partitioning, we use only the
            // number of particles in cells
            static
            StatsNSz arbitrary() 
            {
                uint32_t const n_c_mx(4000u);
                uint32_t const n_p_mx(1000000u);

                gen_t g(mkGen());
                uint32_t n_clls(choose<uint32_t>(1u,n_c_mx,g));
                sp_stats ps(new stats_t(n_clls));
                std::generate(ps->ns.begin(),ps->ns.end(),
                              Chooser<sz_t>(1,n_p_mx,mkGen() ));
                sz_t n = ps->n_particles();
                sz_t chkSz(choose<sz_t>(1,n,g));
                return StatsNSz(ps,chkSz);
            }
        };

        struct prop_lastParticleNminus1
        {
            bool operator()(StatsNSz const & ssnsz)
            {
                ChunkVec chks;
                StatsVec stats;
                Chker::chunks(ssnsz.chkSz,*(ssnsz.ss),chks,stats);
                sz_t n_tot = ssnsz.ss->n_particles();
                bool passed = (n_tot-1 == chks.back().pend);
                if(!passed)
                {
                    std::cerr << "failed: chunk size: " << ssnsz.chkSz
                              << ", n_tot = " << n_tot
                              << ", n_chks = " << chks.size()
                              << ", last particle id: " 
                              << chks.back().pend
                              << std::endl;
                    std::cerr << "ss.ns: ";
                    std::copy(ssnsz.ss->ns.begin(),
                              ssnsz.ss->ns.end(),
                              std::ostream_iterator<sz_t>(
                                  std::cerr,","));
                    std::cerr << std::endl << "chunks:\n";
                    for(size_t i = 0; i < chks.size(); ++i)
                    {
                        std::cerr << chks[i].cid << ": start: "
                                  << chks[i].pstart << ", end: "
                                  << chks[i].pend << "\n";
                    }
                    std::cerr << std::endl << "chunked stats:\n";
                    for(size_t i = 0; i < stats.size(); ++i)
                    {
                        std::cerr << "stats chunk " << i << ":";
                        std::copy(stats[i]->ns.begin(),
                                  stats[i]->ns.end(),
                                  std::ostream_iterator<sz_t>(
                                      std::cerr,","));
                        std::cerr << "\n";
                    }
                }
                return passed;
            }
            std::string name()const{return "prop_lastParticleNMinus1";}
        };

        bool
        prop_firstParticle0(StatsNSz const & ssnsz)
        {
            ChunkVec chks;
            StatsVec stats;
            Chker::chunks(ssnsz.chkSz,*(ssnsz.ss),chks,stats);
            return (0 == chks[0].pstart);
        }

        struct prop_chunksIncrease
        {
            bool operator()(StatsNSz const & ssnsz)
            {
                ChunkVec chks;
                StatsVec stats;
                Chker::chunks(ssnsz.chkSz,*(ssnsz.ss),chks,stats);
                bool check(true);
                for(size_t i = 0; i < chks.size()-1; ++i)
                {
                    check = check and chks[i].pstart < chks[i].pend;
                }
                check = check and chks.back().pstart <= chks.back().pend;
                return check;
            }
            std::string name()const{return "prop_chunksIncrease";}
        };

        struct prop_chunksMonotonic
        {
            bool operator()(StatsNSz const & ssnsz)
            {
                ChunkVec chks;
                StatsVec stats;
                Chker::chunks(ssnsz.chkSz,*(ssnsz.ss),chks,stats);
                bool check(true);
                if(1 == chks.size())
                {
                    return true;
                }
                for(size_t i = 1; i < chks.size(); ++i)
                {
                    bool ok = chks[i].pstart > chks[i-1].pend;
                    check = check and ok;
                    if(!ok)
                    {
                        std::cerr << name() << ":" << i << ": failed, "
                                  << " start i: " << chks[i].pstart
                                  << " end i-1: " << chks[i-1].pend
                                  << "\n";
                    }
                }
                return check;
            }
            std::string name()const{return "prop_chunksMonotonic";}
        };

        struct prop_chunksContig
        {
            bool operator()(StatsNSz const & ssnsz)
            {
                ChunkVec chks;
                StatsVec stats;
                Chker::chunks(ssnsz.chkSz,*(ssnsz.ss),chks,stats);
                bool check(true);
                if(1 == chks.size())
                {
                    return true;
                }
                for(size_t i = 1; i < chks.size(); ++i)
                {
                    bool ok = chks[i].pstart == 1 + chks[i-1].pend
                        and chks[i].cid == 1 + chks[i-1].cid;
                    check = check and ok;
                    if(!ok)
                    {
                        std::cerr << name() << ":" << i << ": failed, "
                                  << " start i: " << chks[i].pstart
                                  << " end i-1: " << chks[i-1].pend
                                  << "\n";
                    }
                }
                return check;
            }
            std::string name()const{return "prop_chunksContiguous";}
        };

        struct prop_chunksDontOverlap
        {
            bool operator()(StatsNSz const & ssnsz)
            {
                ChunkVec chks;
                StatsVec stats;
                Chker::chunks(ssnsz.chkSz,*(ssnsz.ss),chks,stats);
                bool check(true);
                if(1 == chks.size())
                {
                    return true;
                }
                for(size_t i = 1; i < chks.size(); ++i)
                {
                    PtclId s1 = chks[i].pstart;
                    PtclId s2 = chks[i-1].pstart;
                    PtclId e1 = chks[i].pend;
                    PtclId e2 = chks[i-1].pend;
                    bool ok = s2 > s1 ? s2 > e1 : e2 < s1;
                    check = check and ok;
                    if(!ok)
                    {
                        std::cerr << name() << ":" << i << ": failed, "
                                  << " start i: " << chks[i].pstart
                                  << " end i: " << chks[i].pend
                                  << " start i-1: " << chks[i-1].pstart
                                  << " end i-1: " << chks[i-1].pend
                                  << "\n";
                    }
                }
                return check;
            }
            std::string name()const{return "prop_chunksDontOverlap";}
        };

        template <typename prop_t,typename test_props>
        bool
        test_n(test_props & tps)
        {
            bool passed(true);
            prop_t prop;
            std::cout << "running " << tps.n << " tests of "
                      << prop.name()
                      << " with random inputs\n";
            for(uint32_t i = 0; i < tps.n; ++i)
            {
                StatsNSz ssnsz(StatsNSz::arbitrary());
                bool pl = prop(ssnsz);
                if(!pl)
                {
                    std::cerr << "property " << prop.name() 
                              << " failed: ss size: "
                              << ssnsz.ss->size()
                              << "chunk size: " << ssnsz.chkSz
                              << std::endl;
                }
                passed = passed and pl;
            }
            if(passed)
            {
                std::cout << "PASSED " << tps.n << " tests.\n";
            }
            return passed;            
        }

        bool test_7()
        {
#ifndef REQUIRE_ON
#warning "Not really doing test_7: set -DREQUIRE_ON"
            bool passed(true);
            std::cerr << "Not really doing test_7: set -DREQUIRE_ON"
                      << std::endl;
#else
            bool passed(false);
            try
            {
                ChunkVec chks;
                StatsVec stats;
                StatsNSz ssnsz(StatsNSz::arbitrary());
                ssnsz.chkSz = 0;
                Chker::chunks(ssnsz.chkSz,*(ssnsz.ss),chks,stats);
                passed = false;
            }
            catch(std::exception & exc)
            {
                std::cerr << "test_7 caught expected exception:\n"
                          << exc.what()
                          << std::endl;
                passed = true;
            }
#endif
            return passed;
        } // test_7

        bool test_6()
        {
            return test_n<prop_chunksDontOverlap,Test_Props>(t100);
        } // test_6

        bool test_5()
        {
            return test_n<prop_chunksContig,Test_Props>(t100);
        } // test_5

        bool test_4()
        {
            return test_n<prop_chunksMonotonic,Test_Props>(t100);
        } // test_4

        bool test_3()
        {
            return test_n<prop_chunksIncrease,Test_Props>(t100);
        } // test_3


        bool test_2()
        {
            return test_n<prop_lastParticleNminus1,Test_Props>(t5);
        } // test_2

        bool test_1()
        {
            bool passed(true);
            {
                sp_stats ss(new stats_t(5));
                ss->ns[0] = 5;
                ss->ns[1] = 5;
                ss->ns[2] = 5;
                ss->ns[3] = 5;
                ss->ns[4] = 5;
                id_t chkSz(5);
                passed = passed and 
                    prop_firstParticle0(StatsNSz(ss,chkSz));
            }
            {
                uint32_t n_prop_tests(100);
                std::cout << "running " << n_prop_tests << " tests of "
                          << "prop_firstParticle0"
                          << " with random inputs\n";
                bool lpassed(true);
                for(uint32_t i = 0; i < n_prop_tests; ++i)
                {
                    StatsNSz ssnsz(StatsNSz::arbitrary());
                    bool pl = prop_firstParticle0(ssnsz);
                    if(!pl)
                    {
                        std::cerr << "property prop_firstParticle0 failed: "
                                  << ssnsz.ss->size()
                                  << "chunk size: " << ssnsz.chkSz
                                  << std::endl;
                    }
                    lpassed = lpassed and pl;
                }
                if(lpassed)
                {
                    std::cout << "PASSED " << n_prop_tests << " tests.\n";
                }
                passed = passed and lpassed;
            }
            return passed;
        }
        // define additional tests here.

    } // partition_tests::

} // Nut_Test::


// These are the tests for the old partitioning method:



        // struct RankNComm
        // {
        //     nut::Rank r;
        //     nut::CommSz csz;
        //     RankNComm(uint32_t rk,uint32_t c) : r(Rank(rk)),csz(CommSz(c)) {}
        // };

        // /** For given source statistics, the lists of particle ids, 
        //  * interleaved cell-by-cell, should be the same as the single list
        //  * generated by gen_ptcl_ids. This is what gives partition the 
        //  * property of keeping the simulation invariant even as
        //  * partitioning changes. Note that gen_ptcl_ids assembles
        //  * contiguously the particle ids for all cells for one rank--
        //  * we want contiguous for all ranks for each cell. The particle
        //  * ids in that order will be invariant wrt 
        //  * to the number of ranks. */
        // bool 
        // prop_genPtclIdsIndCommSz(RankNComm const & rnc,
        //                          stats_t const & stats)
        // {
        //     id_t const n_cells = stats.size();
        //     nut::vid cumNs(n_cells + 1,0);
        //     std::partial_sum(stats.ns.begin(),stats.ns.end(),&cumNs[1]);
        //     id_t totNPs = cumNs[n_cells];
        //     nut::vid idsManyRanks(totNPs,0);

        //     // for each cell
        //     for(uint32_t j = 0; j < n_cells; ++j)
        //     {
        //         // for each rank less than commSz, get list of particle ids 
        //         // for one cell
        //         id_t const cumN = cumNs[j];
        //         id_t const nc = stats.ns[j];
        //         id_t rank_off = 0;
        //         id_t idx_low  = 0;
        //         id_t idx_high = 0;

        //         for(uint32_t r = 0; r < fromCommSz(rnc.csz); ++r)
        //         {
        //             id_t const ppr = ps_per_rank(r,fromCommSz(rnc.csz),nc);
        //             idx_high += ppr;

        //             ids_per_cell(Rank(r),rnc.csz,cumN,nc,
        //                          &idsManyRanks[idx_low+cumN],
        //                          &idsManyRanks[idx_high+cumN]);
        //             // std::copy(&idsManyRanks[idx_low+cumN],
        //             //           &idsManyRanks[idx_high+cumN],
        //             //           std::ostream_iterator<id_t>(std::cout,","));

        //             idx_low = idx_high;
        //             rank_off += ppr;
        //         }
        //     }

        //     nut::vid idsOneRank(totNPs,0);
        //     gen_ptcl_ids(Rank(0),CommSz(1),stats,idsOneRank);
        //     bool passed = std::equal(idsManyRanks.begin(),idsManyRanks.end(),
        //                              idsOneRank.begin());
        //     if(!passed)
        //     {
        //         std::copy(idsManyRanks.begin(),idsManyRanks.end(),
        //                   std::ostream_iterator<id_t>(std::cout,","));
        //         std::cout << std::endl;
        //         std::copy(idsOneRank.begin(),idsOneRank.end(),
        //                   std::ostream_iterator<id_t>(std::cout,","));
        //         std::cout << std::endl;
        //     }
        //     return passed;
        // }

        // bool test_1()
        // {
        //     bool passed(true);

        //     {
        //         RankNComm rnc(0,1);
        //         stats_t stats(2);
        //         stats.ns[0] = 3;
        //         stats.ns[1] = 4;
        //         passed = passed and prop_genPtclIdsIndCommSz(rnc,stats);
        //     }

        //     {
        //         RankNComm rnc(0,1);
        //         stats_t stats(3);
        //         stats.ns[0] = 3;
        //         stats.ns[1] = 4;
        //         stats.ns[2] = 8;
        //         passed = passed and prop_genPtclIdsIndCommSz(rnc,stats);
        //     }

        //     {
        //         RankNComm rnc(0,2);
        //         stats_t stats(3);
        //         stats.ns[0] = 3;
        //         stats.ns[1] = 4;
        //         stats.ns[2] = 8;
        //         passed = passed and prop_genPtclIdsIndCommSz(rnc,stats);
        //     }

        //     {
        //         RankNComm rnc(0,3);
        //         stats_t stats(25);
        //         stats.ns[0] = 3;
        //         stats.ns[1] = 4;
        //         stats.ns[2] = 8;
        //         passed = passed and prop_genPtclIdsIndCommSz(rnc,stats);
        //     }

        //     {
        //         RankNComm rnc(0,20);
        //         stats_t stats(25);
        //         stats.ns[0] = 3;
        //         stats.ns[1] = 4;
        //         stats.ns[2] = 8;
        //         stats.ns[3] = 2;
        //         stats.ns[4] = 7;
        //         stats.ns[5] = 9;
        //         stats.ns[7] = 3;
        //         stats.ns[8] = 6;
        //         stats.ns[19] = 4;
        //         passed = passed and prop_genPtclIdsIndCommSz(rnc,stats);
        //     }

        //     return passed;
        // }



// version
// $Id$

// End of file
