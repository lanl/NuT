// T. M. Kelley (c) 2011 LANS LLC

#include "gtest/gtest.h"
#include "mpi_helpers.hh"
#include "partition.hh"
#include "test_aux.hh"
#include "types.hh"
#include <iostream>
#include <iterator>
#include <numeric>  // for partial_sum
#include <random>

using nut::CommSz;
using nut::geom_t;
using nut::id_t;
using nut::PtclId;
using nut::Rank;

typedef nut::src_stats_t<geom_t, id_t> stats_t;
typedef std::shared_ptr<stats_t> sp_stats;
typedef stats_t::sz_t sz_t;
using nut::Chunker;
typedef Chunker<stats_t> Chker;
typedef Chker::ChunkVec ChunkVec;
typedef Chker::StatsVec StatsVec;

typedef std::minstd_rand gen_t;
gen_t
mkGen()
{
  std::random_device rd;
  return gen_t(rd());
}

/** choose a number between min and max (inclusive-ish) */
template <typename t>
t
choose(t min, t max, gen_t g)
{
  double urd = std::generate_canonical<double, 64, gen_t>(g);
  return min + t(urd * (double)(max - min));
}

template <typename t>
struct Chooser {
  Chooser(t min, t max, gen_t gn) : mn(min), mx(max), g(gn) {}
  t operator()() { return choose(mn, mx, g); }
  t mn, mx;
  gen_t g;
};

struct Test_Props {
  uint32_t n;
  explicit Test_Props(uint32_t n_ts = 100) : n(n_ts) {}
};  // Test_Props

Test_Props t5(5u);
Test_Props t100;

struct StatsNSz {
  sp_stats ss;
  stats_t::sz_t chkSz;
  explicit StatsNSz(sp_stats ss_, stats_t::sz_t sz_) : ss(ss_), chkSz(sz_) {}
  // for testing partitioning, we use only the
  // number of particles in cells
  static StatsNSz arbitrary()
  {
    uint32_t const n_c_mx(4000u);
    uint32_t const n_p_mx(1000000u);

    gen_t g(mkGen());
    uint32_t n_clls(choose<uint32_t>(1u, n_c_mx, g));
    sp_stats ps(new stats_t(n_clls));
    std::generate(ps->ns.begin(), ps->ns.end(),
                  Chooser<sz_t>(1, n_p_mx, mkGen()));
    sz_t n = ps->n_particles();
    sz_t chkSz(choose<sz_t>(1, n, g));
    return StatsNSz(ps, chkSz);
  }
};

struct prop_lastParticleNminus1 {
  bool operator()(StatsNSz const & ssnsz)
  {
    ChunkVec chks;
    StatsVec stats;
    Chker::chunks(ssnsz.chkSz, *(ssnsz.ss), chks, stats);
    sz_t n_tot = ssnsz.ss->n_particles();
    bool passed = (n_tot - 1 == chks.back().pend);
    if(!passed) {
      std::cerr << "failed: chunk size: " << ssnsz.chkSz
                << ", n_tot = " << n_tot << ", n_chks = " << chks.size()
                << ", last particle id: " << chks.back().pend << std::endl;
      std::cerr << "ss.ns: ";
      std::copy(ssnsz.ss->ns.begin(), ssnsz.ss->ns.end(),
                std::ostream_iterator<sz_t>(std::cerr, ","));
      std::cerr << std::endl << "chunks:\n";
      for(size_t i = 0; i < chks.size(); ++i) {
        std::cerr << chks[i].cid << ": start: " << chks[i].pstart
                  << ", end: " << chks[i].pend << "\n";
      }
      std::cerr << std::endl << "chunked stats:\n";
      for(size_t i = 0; i < stats.size(); ++i) {
        std::cerr << "stats chunk " << i << ":";
        std::copy(stats[i]->ns.begin(), stats[i]->ns.end(),
                  std::ostream_iterator<sz_t>(std::cerr, ","));
        std::cerr << "\n";
      }
    }
    return passed;
  }
  std::string name() const { return "prop_lastParticleNMinus1"; }
};

bool
prop_firstParticle0(StatsNSz const & ssnsz)
{
  ChunkVec chks;
  StatsVec stats;
  Chker::chunks(ssnsz.chkSz, *(ssnsz.ss), chks, stats);
  return (0 == chks[0].pstart);
}

struct prop_chunksIncrease {
  bool operator()(StatsNSz const & ssnsz)
  {
    ChunkVec chks;
    StatsVec stats;
    Chker::chunks(ssnsz.chkSz, *(ssnsz.ss), chks, stats);
    bool check(true);
    for(size_t i = 0; i < chks.size() - 1; ++i) {
      check = check and chks[i].pstart < chks[i].pend;
    }
    check = check and chks.back().pstart <= chks.back().pend;
    return check;
  }
  std::string name() const { return "prop_chunksIncrease"; }
};

struct prop_chunksMonotonic {
  bool operator()(StatsNSz const & ssnsz)
  {
    ChunkVec chks;
    StatsVec stats;
    Chker::chunks(ssnsz.chkSz, *(ssnsz.ss), chks, stats);
    bool check(true);
    if(1 == chks.size()) { return true; }
    for(size_t i = 1; i < chks.size(); ++i) {
      bool ok = chks[i].pstart > chks[i - 1].pend;
      check = check and ok;
      if(!ok) {
        std::cerr << name() << ":" << i << ": failed, "
                  << " start i: " << chks[i].pstart
                  << " end i-1: " << chks[i - 1].pend << "\n";
      }
    }
    return check;
  }
  std::string name() const { return "prop_chunksMonotonic"; }
};

struct prop_chunksContig {
  bool operator()(StatsNSz const & ssnsz)
  {
    ChunkVec chks;
    StatsVec stats;
    Chker::chunks(ssnsz.chkSz, *(ssnsz.ss), chks, stats);
    bool check(true);
    if(1 == chks.size()) { return true; }
    for(size_t i = 1; i < chks.size(); ++i) {
      bool ok = chks[i].pstart == 1 + chks[i - 1].pend and
                chks[i].cid == 1 + chks[i - 1].cid;
      check = check and ok;
      if(!ok) {
        std::cerr << name() << ":" << i << ": failed, "
                  << " start i: " << chks[i].pstart
                  << " end i-1: " << chks[i - 1].pend << "\n";
      }
    }
    return check;
  }
  std::string name() const { return "prop_chunksContiguous"; }
};

struct prop_chunksDontOverlap {
  bool operator()(StatsNSz const & ssnsz)
  {
    ChunkVec chks;
    StatsVec stats;
    Chker::chunks(ssnsz.chkSz, *(ssnsz.ss), chks, stats);
    bool check(true);
    if(1 == chks.size()) { return true; }
    for(size_t i = 1; i < chks.size(); ++i) {
      PtclId s1 = chks[i].pstart;
      PtclId s2 = chks[i - 1].pstart;
      PtclId e1 = chks[i].pend;
      PtclId e2 = chks[i - 1].pend;
      bool ok = s2 > s1 ? s2 > e1 : e2 < s1;
      check = check and ok;
      if(!ok) {
        std::cerr << name() << ":" << i << ": failed, "
                  << " start i: " << chks[i].pstart
                  << " end i: " << chks[i].pend
                  << " start i-1: " << chks[i - 1].pstart
                  << " end i-1: " << chks[i - 1].pend << "\n";
      }
    }
    return check;
  }
  std::string name() const { return "prop_chunksDontOverlap"; }
};

template <typename prop_t, typename test_props>
bool
test_n(test_props & tps)
{
  bool passed(true);
  prop_t prop;
  std::cout << "running " << tps.n << " tests of " << prop.name()
            << " with random inputs\n";
  for(uint32_t i = 0; i < tps.n; ++i) {
    StatsNSz ssnsz(StatsNSz::arbitrary());
    bool pl = prop(ssnsz);
    if(!pl) {
      std::cerr << "property " << prop.name()
                << " failed: ss size: " << ssnsz.ss->size()
                << "chunk size: " << ssnsz.chkSz << std::endl;
    }
    passed = passed and pl;
  }
  if(passed) { std::cout << "PASSED " << tps.n << " tests.\n"; }
  return passed;
}

TEST(nut_partition, take_n_particles_checks_chkSz)
{
#ifndef REQUIRE_ON
#warning "Not really doing test_7: set -DREQUIRE_ON"
  bool passed(true);
  std::cerr << "Not really doing test_7: set -DREQUIRE_ON" << std::endl;
#else
  bool passed(false);
  try {
    ChunkVec chks;
    StatsVec stats;
    StatsNSz ssnsz(StatsNSz::arbitrary());
    ssnsz.chkSz = 0;
    Chker::chunks(ssnsz.chkSz, *(ssnsz.ss), chks, stats);
    passed = false;
  } catch(std::exception & exc) {
    std::cerr << "test_7 caught expected exception:\n"
              << exc.what() << std::endl;
    passed = true;
  }
#endif
  EXPECT_TRUE(passed);
  return;
}  // test_7

TEST(nut_partition, chunks_disjoint)
{
  bool passed = test_n<prop_chunksDontOverlap, Test_Props>(t100);
  EXPECT_TRUE(passed);
  return;
}  // test_6

TEST(nut_partition, chunks_contiguous)
{
  bool passed = test_n<prop_chunksContig, Test_Props>(t100);
  EXPECT_TRUE(passed);
  return;
}  // test_5

TEST(nut_partition, each_chunk_start_gt_previous_chunk_begin)
{
  bool passed = test_n<prop_chunksMonotonic, Test_Props>(t100);
  EXPECT_TRUE(passed);
  return;
}  // test_4

TEST(nut_partition, each_chunk_end_gt_chunk_begin)
{
  bool passed = test_n<prop_chunksIncrease, Test_Props>(t100);
  EXPECT_TRUE(passed);
  return;
}  // test_3

TEST(nut_partition, last_particle_id_correct)
{
  bool passed = test_n<prop_lastParticleNminus1, Test_Props>(t5);
  EXPECT_TRUE(passed);
  return;
}  // test_2

TEST(nut_partition, first_particle_id_is_0)
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
    passed = passed and prop_firstParticle0(StatsNSz(ss, chkSz));
  }
  {
    uint32_t n_prop_tests(100);
    std::cout << "running " << n_prop_tests << " tests of "
              << "prop_firstParticle0"
              << " with random inputs\n";
    bool lpassed(true);
    for(uint32_t i = 0; i < n_prop_tests; ++i) {
      StatsNSz ssnsz(StatsNSz::arbitrary());
      bool pl = prop_firstParticle0(ssnsz);
      if(!pl) {
        std::cerr << "property prop_firstParticle0 failed: " << ssnsz.ss->size()
                  << "chunk size: " << ssnsz.chkSz << std::endl;
      }
      lpassed = lpassed and pl;
    }
    if(lpassed) { std::cout << "PASSED " << n_prop_tests << " tests.\n"; }
    passed = passed and lpassed;
  }
  EXPECT_TRUE(passed);
  return;
}

// End of file
