// bh-3.cc
// T. M. Kelley
// June 19, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved

// C++ version of the bh-2.py program
// bh-2.cc Differs in that it reads the material state file directly, instead of
// using the MC state file.
// bh-3.cc further differs in attempting OpenMP & MPI execution, as well as
// changing to the Philox CBRNG (Salmon et al, SC 2011).

#include "Mesh.hh"
#include "RNG.hh"
#include "constants.hh"
#include "types.hh"
#define USING_HFBC_SIGMAS
#include "Assert.hh"
#include "Cell_Data.hh"
#include "Census.hh"
#include "Log.hh"
#include "MatState.hh"
#include "Opacity.hh"
#include "Particle.hh"
#include "Tally.hh"
#include "cl-args.hh"
#include "fileio.hh"
#include "mpi_helpers.hh"
#include "partition.hh"
#include "serialize.hh"
#include "sourcery.hh"
#include "transport.hh"
#include "utilities.hh"
#include <algorithm>  // transform
#include <fstream>
#include <iterator>
#include <numeric>  // accumulate
#include <stdexcept>
#include <execinfo.h>
// #include <omp.h>

using nut::geom_t;

#ifdef HAVE_MURMELN
using Mesh_T = murmeln_mesh::Spherical_1D_Mesh;
using Mesh_Interface_T = murmeln::Spherical_Mesh_Interface;
#else
using Mesh_Interface_T =
    nut::Sphere_1D<nut::cell_t, geom_t, nut::bdy_types::descriptor>;
using Mesh_T = Mesh_Interface_T;
#endif

using Boundary_Cond_T = nut::Boundary_Cond<Mesh_Interface_T::face_handle_t>;
using vector_t = Mesh_Interface_T::Vector;
using op_t = nut::Opacity<geom_t, vector_t>;
using p_t = nut::Particle<geom_t, nut::rng_t, Mesh_Interface_T::Vector>;
using tally_t = nut::Tally<geom_t, 1>;
using census_t = nut::Census<p_t>;

using vec_geom = std::vector<geom_t>;
using vec_vec = tally_t::vv;
using vsz = std::vector<size_t>;
using log_t = nut::Null_Log;
// using log_t = nut::Std_Log        ;
// using MatState_t = nut::MatState<geom_t, vector_t>;
// Used to construct a mesh and material state
using cell_d_t = nut::Cell_Data<geom_t, vector_t>;
using cell_data_t = std::vector<cell_d_t>;
using cons_state_t = std::pair<cell_data_t, Mesh_T>;
// Used to share a material state and mesh interface
using state_t = std::pair<cell_data_t, Mesh_Interface_T>;
using src_stat_t = nut::src_stats_t<geom_t, nut::id_t>;
using Chnker = nut::Chunker<src_stat_t>;

void
run_cycle(src_stat_t const & stats,
          Mesh_Interface_T const & mesh,
          Boundary_Cond_T const & bcs,
          op_t const & op,
          geom_t const alpha,
          nut::Species const s,
          uint32_t const seed,
          tally_t & tally,
          census_t & census,
          nut::key_t const & key,
          uint32_t const chkSz,
          uint32_t rank,
          uint32_t commSz)
{
  using nut::cell_t;
  using nut::Chunk;
  using nut::ChunkId;
  using nut::ChunkIdVec;
  using nut::cntr_t;
  using nut::LessThan;
  using nut::PtclId;
  using nut::Require;

  cell_t const n_cells = mesh.num_cells();

  Require(stats.ns.size() == n_cells && stats.es.size() == n_cells &&
              stats.ews.size() == n_cells,
          "correct array lenghts");
  // very large for particles in postprocess--there's really no time
  // CAUTION: this can make sims take a very long time!
  geom_t const particle_dt = 1e12;

  nut::cntr_t ctr(0);  // for counting events
  nut::cntr_t n_tot = std::accumulate(stats.ns.begin(), stats.ns.end(), 0);
  cntr_t const rep_frac = n_tot > 5 ? n_tot / 5 : 1;  // report frequency

  log_t log;

  Chnker::StatsVec schunks;
  Chnker::ChunkVec chunks;
  Chnker::chunks(chkSz, stats, chunks, schunks);
  uint32_t const n_chks_tot(chunks.size());
  nut::ChunkIdVec chkIds;
  nut::getChunks(rank, commSz, n_chks_tot, chkIds);

  std::cout << "We have " << n_chks_tot << " total chunks." << std::endl;

  std::cout << "This PE has " << chkIds.size() << " chunks." << std::endl;

  /* This PE does the chunks scheduled by getChunks. Those chunks are listed
   * in chkIds; chkIds are indices into chunks and schunks arrays. For each
   * chunk, we look at the corresponding vector of source numbers in
   * schunks[chkId]: for each entry (cell) we run the corresponding number of
   * particles.
   */

  for(size_t ci = 0; ci < n_chks_tot; ++ci) {
    // for each chunk that this PE does
    nut::ChunkId const chkId = chkIds[ci];
    Chunk const & chunk = chunks[chkId];
    src_stat_t const & ss = *schunks[chkId];
    PtclId curr(chunk.pstart);

    // std::cout << "Processing chunk " << chkId
    //           << ", pstart: " << chunk.pstart
    //           << ", pend: "  << chunk.pend
    //           // << ", taken by thread " << tid
    //           // << " of " << n_threads
    //           << std::endl;

    // for each cell in the chunk, process the particles that
    // originate in that cell.
    for(uint32_t celli = 0; celli < ss.size(); ++celli) {
      Chnker::sz_t const n_ps = ss.ns[celli];
      cell_t const cidx(ss.cidxs[celli]);  // cell idx for this entry
      geom_t const ew(ss.ews[celli]);
      // std::cout << "*** Processing cell " << celli
      //           << ", n = " << n_ps << ", cell = " << cidx
      //           << ", ew = " << ew << std::endl;
      for(uint32_t pi = 0; pi < n_ps; ++pi) {
        LessThan(curr, chunk.pend + 1, "current ptcl id", "last ptcl id");
        // std::cout << "\nProcessing particle " << curr << std::endl;
        id_t const ptcl_id(curr);
        nut::ctr_t ptcl_ctr(nut::rng_t::make_ctr(ptcl_id, 0u, 0u, 0u));
        nut::rng_t ptcl_rng(ptcl_ctr, key);  // for generating the particle
        auto velocity = op.velocity(cidx);
        p_t p_in =
            nut::gen_init_particle<Mesh_Interface_T, geom_t, nut::rng_t, p_t>(
                mesh, cidx, particle_dt, alpha, s, ew, op.T_p(cidx), velocity,
                ptcl_rng);
        nut::ctr_t evt_ctr(nut::rng_t::make_ctr(0u, 0u, ptcl_id, 0u));
        nut::rng_t evt_rng(evt_ctr, key);  // for generating events
        p_in.rng = evt_rng;
        nut::transport_particle(p_in, mesh, op, tally, census, log, bcs, alpha);
        // std::cout << "final state: " << p_out << std::endl;
        ctr++;
        curr++;
        if(ctr % rep_frac == 0) {
          std::cout << ctr << "/" << n_tot << " " << species_name(s)
                    << "'s complete" << std::endl;
        }
      }  // loop over particles
    }    // loop over cells
  }      // chunk loop

  // std::cout << "run_cycle: Transport complete\n";

  // fix up momenta
  // vec_vec new_momenta(n_cells);
  geom_t const one_over_c = 1.0 / nut::c;
  std::transform(tally.momentum.begin(), tally.momentum.end(),
                 tally.momentum.begin(),
                 [&](auto const & v) { return v * one_over_c; });

  // std::transform(tally.momentum.begin(),tally.momentum.end(),tally.momentum.begin(),
  //                nut::div_by<geom_t>(nut::c) );
  // std::copy(tally.momentum
  return;
}  // run_cycle

void
run_cycle_buffer(src_stat_t const & stats,
                 Mesh_Interface_T const & mesh,
                 Boundary_Cond_T const & bcs,
                 op_t const & op,
                 geom_t const alpha,
                 nut::Species const s,
                 uint32_t const seed,
                 tally_t & tally,
                 census_t & census,
                 nut::key_t const & key,
                 uint32_t const chkSz,
                 uint32_t rank,
                 uint32_t commSz)
{
  using nut::cell_t;
  using nut::Chunk;
  using nut::ChunkId;
  using nut::ChunkIdVec;
  using nut::cntr_t;
  using nut::LessThan;
  using nut::PtclId;
  using nut::Require;

  cell_t const n_cells = mesh.num_cells();

  Require(stats.ns.size() == n_cells && stats.es.size() == n_cells &&
              stats.ews.size() == n_cells,
          "correct array lenghts");
  // very large for particles in postprocess--there's really no time
  // CAUTION: this can make sims take a very long time!
  geom_t const particle_dt = 1e12;

  nut::cntr_t ctr(0);  // for counting events
  nut::cntr_t n_tot = std::accumulate(stats.ns.begin(), stats.ns.end(), 0);
  cntr_t const rep_frac = n_tot > 5 ? n_tot / 5 : 1;  // report frequency

  log_t log;

  Chnker::StatsVec schunks;
  Chnker::ChunkVec chunks;
  Chnker::chunks(chkSz, stats, chunks, schunks);
  uint32_t const n_chks_tot(chunks.size());
  nut::ChunkIdVec chkIds;
  nut::getChunks(rank, commSz, n_chks_tot, chkIds);

  std::cout << "We have " << n_chks_tot << " total chunks." << std::endl;

  std::cout << "This PE has " << chkIds.size() << " chunks." << std::endl;

  /* This PE does the chunks scheduled by getChunks. Those chunks are listed
   * in chkIds; chkIds are indices into chunks and schunks arrays. For each
   * chunk, we look at the corresponding vector of source numbers in
   * schunks[chkId]: for each entry (cell) we run the corresponding number of
   * particles.
   */

  std::vector<p_t> p_buff_in(chkSz);
  std::vector<p_t> p_buff_out(chkSz);

  for(size_t ci = 0; ci < n_chks_tot; ++ci) {
    // for each chunk that this PE does
    nut::ChunkId const chkId = chkIds[ci];
    Chunk const & chunk = chunks[chkId];
    src_stat_t const & ss = *schunks[chkId];
    PtclId curr(chunk.pstart);

    // std::cout << "Processing chunk " << chkId
    //           << ", pstart: " << chunk.pstart
    //           << ", pend: "  << chunk.pend
    //           // << ", taken by thread " << tid
    //           // << " of " << n_threads
    //           << std::endl;

    uint32_t p_buff_idx(0);

    // for each cell in the chunk, process the particles that
    // originate in that cell.
    for(uint32_t celli = 0; celli < ss.size(); ++celli) {
      Chnker::sz_t const n_ps = ss.ns[celli];
      cell_t const cidx(ss.cidxs[celli]);  // cell idx for this entry
      geom_t const ew(ss.ews[celli]);
      // std::cout << "*** Processing cell " << celli
      //           << ", n = " << n_ps << ", cell = " << cidx
      //           << ", ew = " << ew << std::endl;
      for(uint32_t pi = 0; pi < n_ps; ++pi) {
        LessThan(curr, chunk.pend + 1, "current ptcl id", "last ptcl id");
        // std::cout << "\nProcessing particle " << curr << std::endl;
        id_t const ptcl_id(curr);
        nut::ctr_t ptcl_ctr(nut::rng_t::make_ctr(ptcl_id, 0u, 0u, 0u));
        nut::rng_t ptcl_rng(ptcl_ctr, key);  // for generating the particle
        p_buff_in[p_buff_idx] =
            nut::gen_init_particle<Mesh_Interface_T, geom_t, nut::rng_t, p_t>(
                mesh, cidx, particle_dt, alpha, s, ew, op.T_p(cidx),
                op.velocity(cidx), ptcl_rng);
        nut::ctr_t evt_ctr(nut::rng_t::make_ctr(0u, 0u, ptcl_id, 0u));
        nut::rng_t evt_rng(evt_ctr, key);  // for generating events
        p_buff_in[p_buff_idx].rng = evt_rng;
        // nut::transport_particle(p_in,mesh,op,
        //                         vel,tally,census,
        //                         log,alpha);
        // std::cout << "final state: " << p_out << std::endl;
        p_buff_idx++;
        ctr++;
        curr++;
        if(ctr % rep_frac == 0) {
          std::cout << ctr << "/" << n_tot << " " << species_name(s)
                    << "'s complete" << std::endl;
        }
      }  // loop over particles
    }    // loop over cells
    // transport particles
    transport(p_buff_in, mesh, bcs, op, tally, p_buff_out, census, log, alpha);

    // dispose of particles, tally escape spectrum

  }  // chunk loop

  // std::cout << "run_cycle: Transport complete\n";
  geom_t const one_over_c{1.0 / nut::c};
  // fix up momenta--divide all by c to derive momentum.
  std::transform(tally.momentum.begin(), tally.momentum.end(),
                 tally.momentum.begin(),
                 [&](auto const & v) { return v * one_over_c; });
  return;
}  // run_cycle_buffered

/*!\brief generate a mesh & material state info by reading a material
 * state file and parsing it into MatState and Mesh objects. */
cons_state_t
get_mat_state(std::string const filename,
              geom_t const llimit,
              geom_t const ulimit)
{
  using nut::Require;
  using vecrows = std::vector<nut::MatStateRowP<geom_t, vector_t>>;
  std::ifstream infile(filename.c_str());
  if(!infile) {
    std::stringstream errstr;
    errstr << "Unable to open file \"" << filename << "\"";
    throw(std::runtime_error(errstr.str()));
  }
  vecrows rows(nut::read_mat_state_file<geom_t, vector_t>(infile));
  infile.close();

  // get a mesh that includes only those cells within the
  // specified limits; also, get back the indices corresponding
  // to the limits.
  size_t llimitIdx(0), ulimitIdx(0);
#ifdef HAVE_MURMELN
  Mesh_T mesh = nut::rows_to_murmeln_mesh<geom_t>(rows, llimit, ulimit,
                                                  llimitIdx, ulimitIdx);
#else
  Mesh_Interface_T mesh =
      nut::rows_to_mesh<geom_t>(rows, llimit, ulimit, llimitIdx, ulimitIdx);
#endif
  Require(ulimitIdx >= llimitIdx, "invalid limits");
  size_t const nrows(ulimitIdx - llimitIdx);
  Require(mesh.num_cells() == nrows,
          "get_mat_state: mesh size and nrows disagree");
  // get the subset of rows within the limits
  vecrows limitedRows(nrows);
  std::copy(&rows[llimitIdx], &rows[ulimitIdx], limitedRows.begin());
  // return the mesh & mat state within the limits
  cell_data_t cell_data{make_cell_data(limitedRows)};
  return cons_state_t(cell_data, mesh);
}  // get_mat_state

void
run_one_species(nut::Species const spec, args_t const & args, state_t & state)
{
  using nut::Check;
  Mesh_Interface_T const & mesh = state.second;
  Boundary_Cond_T bcs{nut::make_vacuum_boundary_1D(mesh)};

  op_t const op(std::move(state.first));

  size_t const ncells(mesh.num_cells());
  size_t const ncen(0);

  std::vector<nut::cell_t> cidxs(mesh.num_cells());
  for(size_t i = 0; i < cidxs.size(); ++i) { cidxs[i] = i + 1; }
  Check(cidxs.size() == op.num_cells(), "Cell indexes size != velocity size");

  src_stat_t stats(ncells);
  tally_t tally(ncells);
  census_t census;
  uint32_t const cycle(1);

  auto const & cell_data = state.first;
  auto nu_es = nut::make_nu_e_iterator<cell_d_t>(cell_data.cbegin());
  auto nu_e_end = nut::make_nu_e_iterator<cell_d_t>(cell_data.end());
  auto nu_ebars = nut::make_nu_ebar_iterator<cell_d_t>(cell_data.begin());
  auto nu_ebar_end = nut::make_nu_ebar_iterator<cell_d_t>(cell_data.end());
  auto nu_xs = nut::make_nu_x_iterator<cell_d_t>(cell_data.begin());
  auto nu_x_end = nut::make_nu_x_iterator<cell_d_t>(cell_data.end());
  switch(spec) {
    case nut::nu_e:
      nut::calc_src_stats_lum(nu_es, nu_e_end, cidxs, stats, args.dt,
                              args.n_particles, ncen);
      break;
    case nut::nu_e_bar:
      nut::calc_src_stats_lum(nu_ebars, nu_ebar_end, cidxs, stats, args.dt,
                              args.n_particles, ncen);
      break;
    case nut::nu_x:
      nut::calc_src_stats_lum(nu_xs, nu_x_end, cidxs, stats, args.dt,
                              args.n_particles, ncen);
      break;
    default:
      std::cerr << "run_one_species: unhandled species"
                << nut::species_name(spec) << std::endl;
      throw std::runtime_error("run_one_species: unhandled species");
  }

  summarize_stats(stats, std::cout, spec, cycle);

  nut::seed_t const lo = species_seed(spec);
  nut::key_t key = nut::rng_t::make_key(args.seed, lo);

  uint32_t const rank(0);
  uint32_t const commSz(1);

  run_cycle(stats, mesh, bcs, op, args.alpha, spec, args.seed, tally, census,
            key, args.chunkSz, rank, commSz);

  std::string outfname(args.outputF + "_" + nut::species_name(spec));
  std::ofstream outf(outfname.c_str());
  if(!outf) {
    std::stringstream errstr;
    errstr << "Unable to open output file \"" << outfname << "\"" << std::endl;
    throw std::runtime_error(errstr.str());
  }
  write_tally_mcphd(outf, tally, spec, cycle);
  summarize_tally(tally, std::cout, spec, cycle);

  return;
}  // run_species

// exception stack trace generator
void
handler()
{
  void * trace_elems[20];
  int trace_elem_count(backtrace(trace_elems, 20));
  char ** stack_syms(backtrace_symbols(trace_elems, trace_elem_count));
  for(int i = 0; i < trace_elem_count; ++i) {
    std::cout << stack_syms[i] << "\n";
  }
  free(stack_syms);

  exit(1);
}

int
main(int argc, char ** argv)
{
  // std::set_terminate(handler);

  args_t args = parseCL(argc, argv);

  if(args.halp){
    std::cout << help();
    return -1;
  }

  cons_state_t c_state = get_mat_state(args.inputF, args.llimit, args.ulimit);
  Mesh_Interface_T mesh{std::get<1>(c_state)};
  cell_data_t mat_state{std::get<0>(c_state)};
  state_t state{mat_state, mesh};

  // for each species, compute source stats, run particles, write tally
  // nu_e

  run_one_species(nut::nu_e, args, state);

  run_one_species(nut::nu_e_bar, args, state);

  run_one_species(nut::nu_x, args, state);

  return 0;
}

// version
// $Id$

// End of file
