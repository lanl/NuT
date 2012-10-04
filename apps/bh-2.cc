// bh-2.cc
// T. M. Kelley
// May 23, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

// C++ version of the bh-2.py program
// Differs in that it reads the material state file directly, instead of using 
// the MC state file.

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

#include <numeric>    // accumulate
#include <algorithm>  // transform
#include <fstream>
#include <stdexcept>
#include <execinfo.h>


using nut::cell_t;
using nut::geom_t;
using nut::Require;

// typedef nut::LCG_RNG                rng_t;
typedef nut::MLCG                   rng_t;
// typedef nut::Buffer_RNG<geom_t>             rng_t;
typedef nut::Opacity<geom_t>        op_t;
typedef nut::Velocity<geom_t>       vel_t;

typedef nut::Particle<geom_t,rng_t> p_t;
typedef nut::Tally<geom_t>          tally_t;
typedef nut::Census<p_t>            census_t;

typedef std::vector<geom_t>   vg;
typedef nut::Null_Log         log_t;
// typedef nut::Std_Log         log_t;
typedef nut::MatState<geom_t> MatState_t;

typedef nut::Sphere_1D<cell_t,geom_t,nut::bdy_types::descriptor> Mesh_t;
typedef std::pair<MatState_t,Mesh_t> state_t;
typedef nut::src_stats_t<geom_t,uint32_t> src_stats_t;
typedef src_stats_t::vsz   vsz;


void run_cycle(vsz const & ns, 
               vg const & es,
               vg const & ews,
               Mesh_t const & mesh,
               op_t const & op,
               vel_t const & vel,
               geom_t const alpha,
               nut::Species const s,
               uint32_t const seed,
               tally_t & tally,
               census_t & census)
{
    using nut::cntr_t;
    cell_t const n_cells = mesh.n_cells();
    Require(ns.size()  == n_cells && 
            es.size()  == n_cells && 
            ews.size() == n_cells, "correct array lenghts");
    // very large for particles in postprocess--there's really no time
    geom_t const particle_dt = 1e12;
    cntr_t ctr(0);
    cntr_t n_tot = std::accumulate(ns.begin(),ns.end(),0);
    cntr_t const rep_frac = n_tot > 5 ? n_tot/5 : 1; // report frequency

    log_t log;

    // geom_t rns[] = {0.1};
    // bool silent(true);
    // rng_t r0(rns,1,silent);
    rng_t r0(seed);

    rng_t::new_gens r3_r4 = r0.split();
    rng_t::new_gens r1_r2 = r3_r4.first.split();

    /* managing RNG streams: we take a stream and split it into two parts.    *
     * Within the loop over cells, the first part will be used for the        *
     * particles of the present cell. The second part is saved to be split    *
     * again in the next cell. This pattern is repeated in the loop over      *
     * particles: the inital RNG stream, rsdp, is split into two pieces.      *
     * The first is assigned the current particle, while the second is saved. *
     * On the next trip (particle i+1), that RNG will be split, with          *
     * the first RNG used for particle i+1, and the second saved for          *
     * trip i+2. This scheme is used to maintain compatibility with McPhD,    *
     * the Haskell version of this code.                                      */
    rng_t rsdc = r1_r2.first;
    for(size_t ic = 0; ic < n_cells; ++ic)
    {
        cell_t const cidx = ic + 1;
        geom_t const ew = ews[ic];
        rng_t::new_gens rp_rsdc = rsdc.split();
        rng_t rsdp = rp_rsdc.first;
        rsdc       = rp_rsdc.second;
        for(size_t ip = 0; ip < ns[ic]; ++ip)
        {
            rng_t::new_gens r_rsd = rsdp.split();
            rng_t rng = r_rsd.first;
            rsdp      = r_rsd.second;
            p_t p_in = nut::gen_init_particle<Mesh_t,geom_t,rng_t,p_t>(
                mesh,cidx,particle_dt,alpha,s,
                ew,
                op.m_T.T_p[ic],vel.v(cidx),
                rng);
            // p_t p_out = 
            nut::transport_particle(p_in,mesh,op,
                                    vel,tally,census,
                                    log,alpha);
            ctr++;
            if(ctr % rep_frac == 0)
            {
                std::cout << ctr << "/" << n_tot << " " << species_name(s)
                          << "'s complete" << std::endl;
            }
        }
    }
    // fix up momenta
    vg new_momenta(n_cells,0);
    std::transform(tally.momentum.begin(),tally.momentum.end(),
                   new_momenta.begin(),nut::div_by<geom_t>(nut::c) );
    tally.momentum = new_momenta;
    return;
} // run_cycle


/*!\brief generate a mesh & material state info by reading a material
 * state file and parsing it into MatState and Mesh objects. */
state_t get_mat_state(std::string const filename,
                      geom_t const llimit,
                      geom_t const ulimit)
{
    typedef std::vector<nut::MatStateRowP<geom_t> > vecrows;
    std::ifstream infile(filename.c_str());
    if(!infile)
    {
        std::stringstream errstr;
        errstr << "Unable to open file \"" << filename << "\"";
        throw(std::runtime_error(errstr.str()));
    }
    vecrows rows(nut::read_mat_state_file<geom_t>(infile));
    infile.close();
    
    // get a mesh that includes only those cells within the 
    // specified limits; also, get back the indices corresponding
    // to the limits. 
    size_t llimitIdx(0),ulimitIdx(0);
    Mesh_t mesh = nut::rows_to_mesh<Mesh_t,geom_t>(
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
    MatState_t const & mat = state.first;
    Mesh_t const & mesh    = state.second;

    MatState_t::Density_T     const & d(mat.density);
    MatState_t::Luminosity_T  const & l(mat.luminosity);
    MatState_t::Temperature_T const & t(mat.temperature);
    MatState_t::Velocity_T    const & v(mat.velocity);
    op_t const op(d,t);

    size_t const ncells(mesh.n_cells());
    size_t const ncen(0);

    // std::cout << "mesh boundaries: ";
    // std::copy(mesh.m_bdys.begin(),mesh.m_bdys.end(),
    //           nut::ge_o_it(std::cout,","));
    // std::cout << "mesh boundary descriptors: ";
    // std::copy(mesh.m_descs.begin(),mesh.m_descs.end(),
    //           nut::sz_o_it(std::cout,","));
    // std::cout << std::endl;

    src_stats_t stats(ncells);
    tally_t     tally(ncells);
    census_t    census;
    uint32_t const cycle(1);

    switch(spec)
    {
    case nut::nu_e:
        nut::calc_src_stats_lum(l.nue,stats,args.dt,args.n_particles,ncen);
        break;
    case nut::nu_e_bar:
        nut::calc_src_stats_lum(l.nueb,stats,args.dt,args.n_particles,ncen);
        break;
    case nut::nu_x:
        nut::calc_src_stats_lum(l.nux,stats,args.dt,args.n_particles,ncen);
        break;
    default:
        std::cerr << "run_one_species: unhandled species" 
                  << nut::species_name(spec) << std::endl;
        throw std::runtime_error("run_one_species: unhandled species");
    }
    summarize_stats(stats,std::cout,spec,cycle);

    run_cycle( stats.ns, stats.es, stats.ews,
               mesh, 
               op,
               v,
               args.alpha,
               spec,
               args.seed,
               tally,
               census);
    std::string outfname(args.outputF + "_" + nut::species_name(spec));
    std::ofstream outf(outfname.c_str());
    if(!outf)
    {
        std::stringstream errstr;
        errstr << "Unable to open output file \"" << outfname
               << "\"" << std::endl;
        throw std::runtime_error(errstr.str());
    }
    write_tally(outf,tally,spec,cycle);
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

    // run_one_species(nut::nu_e_bar, args, state);

    // run_one_species(nut::nu_x, args, state);

    return 0;
}


// version
// $Id$

// End of file
