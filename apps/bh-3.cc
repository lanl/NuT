// bh-3.cc
// T. M. Kelley
// June 19, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved

// C++ version of the bh-2.py program
// bh-2.cc Differs in that it reads the material state file directly, instead of
// using the MC state file.
// bh-3.cc further differs in attempting OpenMP & MPI execution, as well as
// changing to the Philox CBRNG (Salmon et al, SC 2011).

#include "Assert.hh"
#include "Census.hh"
#include "Fates.hh"
#include "Log.hh"
#include "MatState.hh"
#include "Mesh.hh"
#define USING_HFBC_SIGMAS
#include "Opacity.hh"
#include "Particle.hh"
#include "RNG.hh"
#include "Tally.hh"
#include "Velocity.hh"
#include "cl-args.hh"
#include "constants.hh"
#include "fileio.hh"
#include "mpi_helpers.hh"
#include "partition.hh"
#include "serialize.hh"
#include "sourcery.hh"
#include "transport.hh"
#include "types.hh"
#include "utilities.hh"

#include <fstream>
#include <execinfo.h>
#include <omp.h>


typedef nut::Opacity<nut::geom_t>        op_t;
typedef nut::Velocity<nut::geom_t,1>     Velocity_t;

typedef nut::Particle<nut::geom_t,nut::rng_t> p_t;
typedef std::vector<p_t>                      p_buff_t;
typedef std::vector<p_buff_t>                 vec_p_buff_t;

typedef nut::Tally<nut::geom_t>     tally_t;
typedef nut::Census<p_t>            census_t;

typedef nut::MatState<nut::geom_t> MatState_t;

typedef nut::Sphere_1D<nut::cell_t,nut::geom_t,nut::bdy_types::descriptor> Mesh_t;
typedef std::pair<MatState_t,Mesh_t> state_t;
typedef nut::src_stats_t<nut::geom_t,nut::id_t> src_stat_t;
typedef nut::Chunker<src_stat_t> Chnker;

typedef std::function<p_t(nut::cell_t const,nut::rng_t &,nut::geom_t const ew)> PGen_T;

/** \brief Generate the particles in the chunk.
    \return the number of particles generated for the buffer. */
uint32_t gen_particle_buffer(p_buff_t & ps
  , src_stat_t const & ss
  , nut::Chunk const & chunk
  , nut::key_t const & key
  , PGen_T gen_p
    )
{
    nut::PtclId curr(chunk.pstart);
    uint32_t p_buff_idx(0u);
    for(uint32_t celli = 0; celli < ss.size(); ++celli)
    {
        Chnker::sz_t const n_ps = ss.ns[celli];
        nut::cell_t const cidx(ss.cidxs[celli]);     // cell idx for this entry
        nut::geom_t const ew(ss.ews[celli]);
        for(uint32_t pi = 0; pi < n_ps; ++pi)
        {
            nut::LessThan(curr,chunk.pend+1,"current ptcl id","last ptcl id");
            id_t const ptcl_id(curr);
            nut::ctr_t ptcl_ctr(nut::rng_t::make_ctr(ptcl_id,0u,0u,0u));
            nut::rng_t ptcl_rng(ptcl_ctr,key); // for generating the particle
            ps[p_buff_idx] = gen_p(cidx,ptcl_rng,ew);
            nut::ctr_t evt_ctr(nut::rng_t::make_ctr(0u,0u,ptcl_id,0u));
            nut::rng_t evt_rng(evt_ctr,key);   // for generating events
            ps[p_buff_idx].rng = evt_rng;
            p_buff_idx++;
            curr++;
        } // loop over particles
    } // loop over cells
    return p_buff_idx;
} // gen_particle_buffer


/** \brief dispose of the particles in the buffer, pushing census particles
    onto the census, and recording escaped particles in the escape spectrum.
    \param p_buff_out: buffer of particles to dispose of.
    \param n_ps: the number of particles in the buffer
    \param census: a census instance
    \param tally: a tally instance */
void dispose_particle_buffer(p_buff_t const & p_buff_out, uint32_t const n_ps,
    census_t & census, tally_t & tally)
{
    using nut::Fates;
    for(uint32_t ip = 0; ip < n_ps; ++ip)
    {
        p_t const & p(p_buff_out[ip]);
        nut::Require(p.fate != Fates::NOT_DEAD_YET,"Output particle did not meet fate");
        switch(p.fate)
        {
        case Fates::STEP_END:
            census.append(p);
            break;
        case Fates::ESCAPE:
            tally.count_escape_spectrum(p.e,p.weight);
            break;
        case Fates::NUCLEON_ABS:
            break; // no action required
        default:
            nut::Require(false,"Unknown particle fate!");
        }
    } // loop over particles
    return;
}


void run_cycle_buffered(src_stat_t const & stats
               , Mesh_t const & mesh
               , op_t const & op
               , Velocity_t const & vel
               , nut::geom_t const alpha
               , nut::Species const s
               , uint32_t const seed
               , tally_t & tally
               , census_t & census
               , nut::key_t const & key
               , uint32_t const chkSz
               , uint32_t rank
               , uint32_t commSz
    )
{
    using nut::cell_t;
    using nut::geom_t;
    using nut::vec_t;
    using nut::Require;

    static size_t const dim(p_t::dim);

    cell_t const n_cells = mesh.n_cells();

    Require(stats.ns.size()  == n_cells &&
            stats.es.size()  == n_cells &&
            stats.ews.size() == n_cells, "correct array lenghts");
    // very large for particles in postprocess--there's really no time
    // CAUTION: this can make sims take a very long time!
    geom_t const particle_dt = 1e12;

    Chnker::StatsVec schunks;
    Chnker::ChunkVec chunks;
    Chnker::chunks(chkSz,stats,chunks,schunks);
    uint32_t const n_chks_tot(chunks.size());
    nut::ChunkIdVec chkIds;
    nut::getChunks(rank,commSz,n_chks_tot,chkIds);

    // create pools of tallies, particle buffers, and censuseses, one per thread
    uint32_t const n_threads(omp_get_max_threads());
    std::vector<census_t> censuses(n_threads);
    std::vector<tally_t> tallies(n_threads);
    vec_p_buff_t p_buffs_in(n_threads);
    vec_p_buff_t p_buffs_out(n_threads);
    for(uint32_t tid = 0; tid < n_threads; ++tid)
    {
        tallies[tid].resize(n_cells);
        p_buffs_in[tid].resize(chkSz);
        p_buffs_out[tid].resize(chkSz);
    }

    std::cout << "We have " << n_chks_tot << " total chunks." << std::endl;

    std::cout << "This PE has " << chkIds.size() << " chunks for "
        << n_threads << " threads."
        << std::endl;

    /* This PE does the chunks scheduled by getChunks. Those chunks are listed
     * in chkIds; chkIds are indices into chunks and schunks arrays. For each
     * chunk, we look at the corresponding vector of source numbers in
     * schunks[chkId]: for each entry (cell) we run the corresponding number of
     * particles.
     */

#pragma omp parallel for
    for(size_t ci = 0; ci < n_chks_tot; ++ci)
    {
        // for each chunk that this PE does
        nut::ChunkId const chkId = chkIds[ci];
        nut::Chunk const & chunk = chunks[chkId];
        src_stat_t const & ss = *schunks[chkId];

        uint32_t const tid(omp_get_thread_num());
        tally_t & thr_tally(tallies[tid]);
        p_buff_t & p_buff_in(p_buffs_in[tid]);
        p_buff_t & p_buff_out(p_buffs_out[tid]);
        census_t & thr_census(censuses[tid]);

        uint32_t const p_buff_idx = gen_particle_buffer(p_buff_in,ss,chunk,key,
            [&](cell_t const cidx,nut::rng_t & rng, geom_t const ew){
                return nut::gen_init_particle<Mesh_t,geom_t,nut::rng_t,p_t>(
                    mesh,cidx,particle_dt,alpha,s,ew,
                    op.temp(cidx),vel.v(cidx),rng);}
            );

        transport(p_buff_in,p_buff_idx,mesh,op,vel,thr_tally,p_buff_out,alpha);
        dispose_particle_buffer(p_buff_out, p_buff_idx, thr_census, thr_tally);
    } // chunk loop

    // std::cout << "run_cycle: Transport complete\n";
    for(uint32_t tid = 0; tid < n_threads; ++tid)
    {
        tally.merge(tallies[tid]);
        census.merge(censuses[tid]);
    }

    // std::cout << "run_cycle: Merge complete \n";

    // fix up momenta: divide by c
    std::transform(tally.momentum.begin(),tally.momentum.end(),
                   tally.momentum.begin(),
                   [&](vec_t<dim> & v){return v.div_by(nut::c);});
    return;
} // run_cycle_buffered


/*!\brief generate a mesh & material state info by reading a material
 * state file and parsing it into MatState and Mesh objects. */
state_t get_mat_state(std::string const filename,
                      nut::geom_t const llimit,
                      nut::geom_t const ulimit)
{
    using nut::Require;
    typedef std::vector<nut::MatStateRowP<nut::geom_t> > vecrows;
    std::ifstream infile(filename.c_str());
    if(!infile)
    {
        std::stringstream errstr;
        errstr << "Unable to open file \"" << filename << "\"";
        throw(std::runtime_error(errstr.str()));
    }
    vecrows rows(nut::read_mat_state_file<nut::geom_t>(infile));
    infile.close();

    // get a mesh that includes only those cells within the
    // specified limits; also, get back the indices corresponding
    // to the limits.
    size_t llimitIdx(0),ulimitIdx(0);
    Mesh_t mesh = nut::rows_to_mesh<Mesh_t,nut::geom_t>(
        rows, llimit, ulimit, llimitIdx, ulimitIdx);
    Require( ulimitIdx >= llimitIdx,"invalid limits");
    size_t const nrows(ulimitIdx - llimitIdx);
    Require(mesh.n_cells() == nrows,
            "get_mat_state: mesh size and nrows disagree");
    // get the subset of rows within the limits
    vecrows limitedRows(nrows);
    std::copy(&rows[llimitIdx],&rows[ulimitIdx],limitedRows.begin());
    // return the mesh & mat state within the limits
    return state_t(MatState_t(limitedRows),mesh);
} // get_mat_state


void run_one_species(nut::Species const spec,
                     args_t const & args,
                     state_t const & state)
{
    using nut::Check;
    MatState_t const & mat = state.first;
    Mesh_t const & mesh    = state.second;

    MatState_t::Density_T     const & d(mat.density);
    MatState_t::Luminosity_T  const & l(mat.luminosity);
    MatState_t::Temperature_T const & t(mat.temperature);
    MatState_t::Velocity_T    const & v(mat.velocity);
    op_t const op(d,t);

    size_t const ncells(mesh.n_cells());
    size_t const ncen(0);

    std::vector<nut::cell_t> cidxs(mesh.n_cells());
    for(size_t i = 0; i < cidxs.size(); ++i){ cidxs[i] = i+1; }
    Check(cidxs.size() == v.size(),"Cell indexes size != velocity size");

    src_stat_t stats(ncells);
    tally_t    tally(ncells);
    census_t   census;
    uint32_t const cycle(1);

    switch(spec)
    {
    case nut::nu_e:
        nut::calc_src_stats_lum(l.nue,cidxs,stats,args.dt,args.n_particles,ncen);
        break;
    case nut::nu_e_bar:
        nut::calc_src_stats_lum(l.nueb,cidxs,stats,args.dt,args.n_particles,ncen);
        break;
    case nut::nu_x:
        nut::calc_src_stats_lum(l.nux,cidxs,stats,args.dt,args.n_particles,ncen);
        break;
    default:
        std::cerr << "run_one_species: unhandled species"
                  << nut::species_name(spec) << std::endl;
        throw std::runtime_error("run_one_species: unhandled species");
    }

    summarize_stats(stats,std::cout,spec,cycle);

    nut::seed_t const lo = species_seed(spec);
    nut::key_t key = nut::rng_t::make_key(args.seed,lo);

    uint32_t const rank(0);
    uint32_t const commSz(1);

    run_cycle_buffered( stats,
               mesh,
               op,
               v,
               args.alpha,
               spec,
               args.seed,
               tally,
               census,
               key
               ,args.chunkSz
               ,rank
               ,commSz
               );

    std::string outfname(args.outputF + "_" + nut::species_name(spec));
    std::ofstream outf(outfname.c_str());
    if(!outf)
    {
        std::stringstream errstr;
        errstr << "Unable to open output file \"" << outfname
               << "\"" << std::endl;
        throw std::runtime_error(errstr.str());
    }
    write_tally_mcphd(outf,tally,spec,cycle);
    summarize_tally(tally,std::cout,spec,cycle);

    return;
} // run_species


// exception stack trace generator
void handler()
{
    void *trace_elems[20];
    int trace_elem_count(backtrace( trace_elems, 20 ));
    char **stack_syms(backtrace_symbols( trace_elems, trace_elem_count ));
    for ( int i = 0 ; i < trace_elem_count ; ++i )
    {
        std::cout << stack_syms[i] << "\n";
    }
    free( stack_syms );

    exit(1);
}


int main(int argc, char ** argv)
{
    // std::set_terminate(handler);

    args_t args = parseCL(argc,argv);

    state_t state = get_mat_state(args.inputF,args.llimit,args.ulimit);

    // for each species, compute source stats, run particles, write tally
    // nu_e

    run_one_species(nut::nu_e, args, state);

    run_one_species(nut::nu_e_bar, args, state);

    run_one_species(nut::nu_x, args, state);

    return 0;
}

// End of file
