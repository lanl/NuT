// bh-3.cc
// T. M. Kelley
// June 19, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved

// C++ version of the bh-2.py program
// bh-2.cc Differs in that it reads the material state file directly, instead of
// using the MC state file.
// bh-3.cc further differs in attempting OpenMP & MPI execution, as well as
// changing to the Philox CBRNG (Salmon et al, SC 2011).

#include "types.hh"
#include "constants.hh"
#include "Mesh.hh"
#include "RNG.hh"
#define USING_HFBC_SIGMAS
#include "Opacity.hh"
#include "Velocity.hh"
#include "Particle.hh"
#include "Assert.hh"
#include "Tally.hh"
#include "Census.hh"
#include "sourcery.hh"
#include "transport.hh"
#include "Log.hh"
#include "MatState.hh"
#include "fileio.hh"
#include "cl-args.hh"
#include "serialize.hh"
#include "utilities.hh"
#include "mpi_helpers.hh"
#include "partition.hh"
#include <numeric>    // accumulate
#include <algorithm>  // transform
#include <fstream>
#include <stdexcept>
#include <execinfo.h>
#include <iterator>
// #include <omp.h>

typedef nut::Opacity<nut::geom_t>        op_t;
typedef nut::Velocity<nut::geom_t>       Velocity_t;

typedef nut::Particle<nut::geom_t,nut::rng_t> p_t;
typedef nut::Tally<nut::geom_t>          tally_t;
typedef nut::Census<p_t>            census_t;

typedef std::vector<nut::geom_t>   vg;
typedef std::vector<size_t>   vsz;
typedef nut::Null_Log         log_t;
// typedef nut::Std_Log         log_t;
typedef nut::MatState<nut::geom_t> MatState_t;

typedef nut::Sphere_1D<nut::cell_t,nut::geom_t,nut::bdy_types::descriptor> Mesh_t;
typedef std::pair<MatState_t,Mesh_t> state_t;
typedef nut::src_stats_t<nut::geom_t,nut::id_t> src_stat_t;
typedef nut::Chunker<src_stat_t> Chnker;


void run_cycle(src_stat_t const & stats
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
    using nut::Require;
    using nut::cntr_t;
    using nut::Chunk;
    using nut::ChunkId;
    using nut::ChunkIdVec;
    using nut::PtclId;
    using nut::LessThan;

    cell_t const n_cells = mesh.n_cells();

    Require(stats.ns.size()  == n_cells &&
            stats.es.size()  == n_cells &&
            stats.ews.size() == n_cells, "correct array lenghts");
    // very large for particles in postprocess--there's really no time
    // CAUTION: this can make sims take a very long time!
    geom_t const particle_dt = 1e12;

    nut::cntr_t ctr(0);  // for counting events
    nut::cntr_t n_tot = std::accumulate(stats.ns.begin(),stats.ns.end(),0);
    cntr_t const rep_frac = n_tot > 5 ? n_tot/5 : 1; // report frequency

    log_t log;

    Chnker::StatsVec schunks;
    Chnker::ChunkVec chunks;
    Chnker::chunks(chkSz,stats,chunks,schunks);
    uint32_t const n_chks_tot(chunks.size());
    nut::ChunkIdVec chkIds;
    nut::getChunks(rank,commSz,n_chks_tot,chkIds);

    std::cout << "We have " << n_chks_tot << " total chunks." << std::endl;

    std::cout << "This PE has " << chkIds.size() << " chunks."
        << std::endl;

    /* This PE does the chunks scheduled by getChunks. Those chunks are listed
     * in chkIds; chkIds are indices into chunks and schunks arrays. For each
     * chunk, we look at the corresponding vector of source numbers in
     * schunks[chkId]: for each entry (cell) we run the corresponding number of
     * particles.
     */


    for(size_t ci = 0; ci < n_chks_tot; ++ci)
    {
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
        for(uint32_t celli = 0; celli < ss.size(); ++celli)
        {
            Chnker::sz_t const n_ps = ss.ns[celli];
            cell_t const cidx(ss.cidxs[celli]);     // cell idx for this entry
            geom_t const ew(ss.ews[celli]);
            // std::cout << "*** Processing cell " << celli
            //           << ", n = " << n_ps << ", cell = " << cidx
            //           << ", ew = " << ew << std::endl;
            for(uint32_t pi = 0; pi < n_ps; ++pi)
            {
                LessThan(curr,chunk.pend+1,"current ptcl id","last ptcl id");
                // std::cout << "\nProcessing particle " << curr << std::endl;
                id_t const ptcl_id(curr);
                nut::ctr_t ptcl_ctr(nut::rng_t::make_ctr(ptcl_id,0u,0u,0u));
                nut::rng_t ptcl_rng(ptcl_ctr,key); // for generating the particle
                p_t p_in = nut::gen_init_particle<Mesh_t,geom_t,nut::rng_t,p_t>(
                    mesh,cidx,particle_dt,alpha,s,ew,
                    op.temp(cidx),vel.v(cidx),
                    ptcl_rng);
                nut::ctr_t evt_ctr(nut::rng_t::make_ctr(0u,0u,ptcl_id,0u));
                nut::rng_t evt_rng(evt_ctr,key);   // for generating events
                p_in.rng = evt_rng;
                nut::transport_particle(p_in,mesh,op,
                                        vel,tally,census,
                                        log,alpha);
                // std::cout << "final state: " << p_out << std::endl;
                ctr++;
                curr++;
                if(ctr % rep_frac == 0)
                {
                    std::cout << ctr << "/" << n_tot << " " << species_name(s)
                              << "'s complete" << std::endl;
                }
            } // loop over particles
        } // loop over cells
    } // chunk loop

    // std::cout << "run_cycle: Transport complete\n";

    // fix up momenta
    vg new_momenta(n_cells,0);
    std::transform(tally.momentum.begin(),tally.momentum.end(),
                   new_momenta.begin(),nut::div_by<geom_t>(nut::c) );
    tally.momentum = new_momenta;
    return;
} // run_cycle


void run_cycle_buffer(
                 src_stat_t const & stats
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
    using nut::Require;
    using nut::cntr_t;
    using nut::Chunk;
    using nut::ChunkId;
    using nut::ChunkIdVec;
    using nut::PtclId;
    using nut::LessThan;

    cell_t const n_cells = mesh.n_cells();

    Require(stats.ns.size()  == n_cells &&
            stats.es.size()  == n_cells &&
            stats.ews.size() == n_cells, "correct array lenghts");
    // very large for particles in postprocess--there's really no time
    // CAUTION: this can make sims take a very long time!
    geom_t const particle_dt = 1e12;

    nut::cntr_t ctr(0);  // for counting events
    nut::cntr_t n_tot = std::accumulate(stats.ns.begin(),stats.ns.end(),0);
    cntr_t const rep_frac = n_tot > 5 ? n_tot/5 : 1; // report frequency

    log_t log;

    Chnker::StatsVec schunks;
    Chnker::ChunkVec chunks;
    Chnker::chunks(chkSz,stats,chunks,schunks);
    uint32_t const n_chks_tot(chunks.size());
    nut::ChunkIdVec chkIds;
    nut::getChunks(rank,commSz,n_chks_tot,chkIds);

    std::cout << "We have " << n_chks_tot << " total chunks." << std::endl;

    std::cout << "This PE has " << chkIds.size() << " chunks."
        << std::endl;

    /* This PE does the chunks scheduled by getChunks. Those chunks are listed
     * in chkIds; chkIds are indices into chunks and schunks arrays. For each
     * chunk, we look at the corresponding vector of source numbers in
     * schunks[chkId]: for each entry (cell) we run the corresponding number of
     * particles.
     */


    std::vector<p_t> p_buff_in(chkSz);
    std::vector<p_t> p_buff_out(chkSz);

    for(size_t ci = 0; ci < n_chks_tot; ++ci)
    {
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
        for(uint32_t celli = 0; celli < ss.size(); ++celli)
        {
            Chnker::sz_t const n_ps = ss.ns[celli];
            cell_t const cidx(ss.cidxs[celli]);     // cell idx for this entry
            geom_t const ew(ss.ews[celli]);
            // std::cout << "*** Processing cell " << celli
            //           << ", n = " << n_ps << ", cell = " << cidx
            //           << ", ew = " << ew << std::endl;
            for(uint32_t pi = 0; pi < n_ps; ++pi)
            {
                LessThan(curr,chunk.pend+1,"current ptcl id","last ptcl id");
                // std::cout << "\nProcessing particle " << curr << std::endl;
                id_t const ptcl_id(curr);
                nut::ctr_t ptcl_ctr(nut::rng_t::make_ctr(ptcl_id,0u,0u,0u));
                nut::rng_t ptcl_rng(ptcl_ctr,key); // for generating the particle
                p_buff_in[p_buff_idx] = nut::gen_init_particle<Mesh_t,geom_t,nut::rng_t,p_t>(
                    mesh,cidx,particle_dt,alpha,s,ew,
                    op.temp(cidx),vel.v(cidx),
                    ptcl_rng);
                nut::ctr_t evt_ctr(nut::rng_t::make_ctr(0u,0u,ptcl_id,0u));
                nut::rng_t evt_rng(evt_ctr,key);   // for generating events
                p_buff_in[p_buff_idx].rng = evt_rng;
                // nut::transport_particle(p_in,mesh,op,
                //                         vel,tally,census,
                //                         log,alpha);
                // std::cout << "final state: " << p_out << std::endl;
                p_buff_idx++;
                ctr++;
                curr++;
                if(ctr % rep_frac == 0)
                {
                    std::cout << ctr << "/" << n_tot << " " << species_name(s)
                              << "'s complete" << std::endl;
                }
            } // loop over particles
        } // loop over cells
        // transport particles
        transport(p_buff_in,mesh,op,vel,tally,p_buff_out,census,log,alpha);

        // dispose of particles, tally escape spectrum

    } // chunk loop

    // std::cout << "run_cycle: Transport complete\n";

    // fix up momenta
    vg new_momenta(n_cells,0);
    std::transform(tally.momentum.begin(),tally.momentum.end(),
                   new_momenta.begin(),nut::div_by<geom_t>(nut::c) );
    tally.momentum = new_momenta;
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

    run_cycle( stats,
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


// version
// $Id$

// End of file
