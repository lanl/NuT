// partition.hh
// T. M. Kelley
// Jun 19, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved


#ifndef PARTITION_HH
#define PARTITION_HH

#include "Assert.hh"
#include "sourcery.hh"  // for src_stats_t
#include "mpi_helpers.hh"
#include "types.hh"     // for id_t
#include <iterator>     // for distance
#include <memory>
#include <algorithm>
#include <iostream>

// functions for partitioning NuT across processes and threads.
namespace nut
{
    /** The new way: match up with McPhD's new method for partitioning. Here,
     * chunks are explicitly generated; chunks are explicitly scheduled onto
     * different ranks via getChunks.
     */
    typedef uint32_t ChunkId;

    typedef uint32_t PtclId;

    typedef std::vector<ChunkId> ChunkIdVec;

    template <typename Tally_t>
    struct IdxTally
    {
        ChunkId cid;
        Tally_t & t;
    }; // IdxTally

    struct Chunk
    {
        ChunkId cid;
        PtclId pstart;
        PtclId pend;
        Chunk(ChunkId const c,PtclId const s, PtclId const e)
            : cid(c),pstart(s),pend(e) {}
        Chunk()
            : cid(uint32_t(-1)),pstart(uint32_t(-1)),pend(uint32_t(-1)) {}
    };

    /** Compute the chunks for which a rank is responsible, given
     * the size of the communicator dividing work and the
     * total number of chunks.*/
    void getChunks(uint32_t const rank, uint32_t const commSz,
                   uint32_t const nc, ChunkIdVec & chunks);

    /** Split source statistics into chunks. Provides types via its *
     * template parameter.                                          */
    template <typename src_stats_t>
    struct Chunker
    {
        typedef std::shared_ptr<src_stats_t> sp_stats;
        typedef std::vector<sp_stats> StatsVec;
        typedef std::vector<Chunk>   ChunkVec;
        typedef src_stats_t src_stats;
        typedef typename src_stats::fp_t fp_t;
        typedef typename src_stats::sz_t sz_t;
        typedef typename src_stats::vsz  vsz;

        /** Make a vector of chunk objects, and a parallel vector of *
         * chunked source statistics.                                */
        static
        void
        chunks(uint32_t const chkSz,
               src_stats const & statsIn,
               ChunkVec & chunks,
               StatsVec & statsVOut)
        {
            take_n_particles(chkSz, statsIn,statsVOut);
            size_t const n_chks(statsVOut.size());
            Check(n_chks > 0,"Chunker::chunks() Got 0 chunks??");
            vsz ns(n_chks);
            std::transform(statsVOut.begin(),statsVOut.end(),ns.begin(),
                           &get_n_ps);
            vsz cums(n_chks,0);
            // if n_chks == 1, next line ok (&cums[1])
            // b/c then ns.begin() == --ns.end()
            std::partial_sum(ns.begin(),--(ns.end()),&cums[1]);
            if(chunks.size() != n_chks)
            {
                chunks.resize(n_chks);
            }
            for(sz_t i = 0; i < chunks.size(); ++i)
            {
                mkChunk(i,cums[i],ns[i],chunks[i]);
            }
            return;
        } // chunks

        static
        sz_t get_n_ps(sp_stats const ss){return ss->n_particles();}

        static
        void
        mkChunk( uint32_t const i,
                 sz_t const pid,
                 sz_t const chkSz,
                 Chunk & chunk)
        {
            chunk.cid = i;
            chunk.pstart = pid;
            chunk.pend   = pid + chkSz - 1;
        }


        /** Given source statistics and chunk size, define a set of chunks.*
         * The last chunk will have fewer particles (unless n_particles is *
         * a multiple of n_chunks).                                        */
        static
        void
        take_n_particles(uint32_t const chkSz,
                         src_stats const & statsIn,
                         StatsVec & statsVOut)
        {
            GreaterThan(chkSz,0u,"take_n_particles:chunk size");
            sz_t const n_ps_tot(statsIn.n_particles());
            uint32_t const nchks(n_ps_tot / chkSz +
                                 ((n_ps_tot % chkSz != 0) ? 1 : 0));
            statsVOut.reserve(nchks);
            statsVOut.resize(0);
            // anticipate that chunk will be size chkSz
            sp_stats chk(new src_stats_t(0));
            chk -> reserve(chkSz);
            sz_t n_curr(0); // # of particles in chunk
            sz_t curr_chunk(0);
            for(sz_t i = 0; i < statsIn.size(); ++i)
            {
                sz_t const n  = statsIn.ns[i];
                sz_t const c  = statsIn.cidxs[i];
                fp_t const e  = statsIn.es[i];
                fp_t const ew = statsIn.ews[i];
                sz_t n_rem(n),n_used(0);
                while(n_rem > 0)
                {
                    sz_t const n_take(chkSz - n_curr);
                    if(n_rem < n_take)
                    {
                        chk -> push_back(n_rem,c,e,ew);
                        n_curr += n_rem;
                        n_used += n_rem;
                        n_rem = 0;
                    }
                    else
                    {
                        // finalize this chunk, start new one
                        chk -> push_back(n_take,c,e,ew);
                        Check(chk->size() <= chkSz,"chunk took too many");
                        n_used += n_take;
                        n_rem -= n_take;
                        statsVOut.push_back(chk);
                        chk = sp_stats(new src_stats_t(0));
                        chk -> reserve(chkSz);
                        n_curr = 0;
                        curr_chunk++;
                    }
                } // while(n_rem)
                Check(n_used == n,"Failed loop inv n_used == n");
            } // loop over statistics / cells
            if(n_curr != 0) // particles in last chunk...
            {
                statsVOut.push_back(chk);  // ...so commit chunk
            }
            return;
        } // takeNParticle



        /* Alternate implementation of take_n_particles.
           Meh. Looks like this won't really be much savings. Could be
           a bit better if we had better notation for dealing with
           tuples and lists/vectors.
         */
        // // the information in the "head" of the input list. This changes as
        // // we "transfer" particles from input to output.
        // struct head
        // {
        //     sz_t n;
        //     fp_t e,ew;
        // };
        // // pointers/indices into various data structures
        // struct comm_ptr
        // {
        //     std::shared_ptr<head> ph;
        //     size_t issIn; // index into src_stats_t arrays
        //     src_stats const & ssIn;
        //     StatsVec & ssOut;
        //     comm_ptr(std::shared_ptr<head> phi,size_t issIni,
        //              src_stats const & ssIni,StatsVec & ssOuti)
        //         : ph(phi),issIn(issIni),ssIn(ssIni),ssOut(ssOuti)
        //     {}
        // };


        // static
        // void take_n_ps(uint32_t const chkSz,
        //                src_stats const & statsIn,
        //                StatsVec & statsVOut)
        // {
        //     std::shared_ptr<head> ph(new head);
        //     ph->n=sz_t(0); ph->e=fp_t(0); ph->ew=fp_t(0);
        //     comm_ptr cp(

        // }

        // static
        // comm_ptr &
        // take_n(n,comm_ptr & curr)
        // {
        //     if(0 == curr.ph -> n) // then move to next cell
        //     {
        //         curr.issIn++;
        //         if(curr.ssIn.size() == curr.issIn)
        //         {
        //             return curr;
        //         }
        //         curr.ph->n = curr.ssIn.ns[issIn];
        //         curr.ph->e = curr.ssIn.es[issIn];
        //         curr.ph->ew = curr.ssIn.ews[issIn];
        //     }
        //     sz_t & nc(curr.ph->n);
        //     if(n < nc)
        //     {
        //         ssOut.back()->push_back(n,e,ew);
        //         nc -= n;
        //         return curr;
        //     }
        //     else
        //     {
        //         ssOut.back()->push_back(nc,e,ew);
        //         nc = 0;
        //         return take_n(n - nc,curr);
        //     }
        // } // take_n

    }; // Chunker

} // :: nut

// /** Parallel decomposition: divide up particles in each cell between
//  * processors. Every rank gets nc / commSz particles, plus one extra iff
//  *
//  *    r < nc % comm_sz.
//  *
//  * This potentially puts more work onto processors with lower rank: each
//  * time, the lowest ranks get an extra particle. Alternatives are to
//  * sacrifice reproducibility (simply replicate), or move toward a domain-
//  * decomposed approach. The latter has its own load imbalances. Or could
//  * get a little more clever about assigning the extra patricles in each cell.
//  *
//  * To get the same result regardless of the number of MPI ranks, we
//  * associate a counter with every particle. The mapping from rank and
//  * communicator size to particle ids is performed by genPtclIds.
//  */
// namespace nut
// {

//     /** Number of particles per cell per rank, based on comm size */
//     template <typename sz_t>
//     inline
//     sz_t ps_per_rank(sz_t const rank, sz_t const csz, sz_t const n_cell)
//     {
//         sz_t const nbase = n_cell / csz;
//         sz_t const nextra = rank < (n_cell % csz) ? 1 : 0;
//         return nbase + nextra;
//     }


//     /** renumber_stats: partitions particle numbers among ranks in *
//      * domain replicated MC. Straightforward transcription of the  *
//      * McPhD function renumberStats. The statsout array must be    *
//      * sized correctly (same size as statsin).                     */
//     template <typename fp_t,typename sz_t>
//     void renumber_stats(Rank const & r , CommSz const & c,
//                         src_stats_t<fp_t,sz_t> const & statsin,
//                         src_stats_t<fp_t,sz_t> & statsout
//         )
//     {
//         Equal(statsout.size(),statsin.size(),"statsout size","statsin size");
//         // copy the energy and energy weight arrays verbatim
//         statsout.ews.assign(statsin.ews.begin(),statsin.ews.end());
//         statsout.es.assign(statsin.es.begin(),statsin.es.end());
//         std::transform(statsin.ns.begin(),statsin.ns.end(),statsout.ns.begin(),
//                        [&](sz_t const n){return ps_per_rank<sz_t>(r,c,n);});
//         return;
//     } // renumber_stats


//     /** Compute a rank's offset in the global set of particle ids (per cell).
//      *    offset = \sum_{i=0}^{i<r} ps_per_rank(i,c,n_cell) */
//     template <typename sz_t>
//     inline
//     sz_t
//     rank_off( Rank const & rank, CommSz const & csz, sz_t const n_cll)
//     {
//         sz_t const c = sz_t(fromCommSz(csz));
//         sz_t const r = sz_t(fromRank(rank));
//         sz_t acc = 0;
//         for(uint32_t i = 0; i < r; ++i)
//         {
//             acc += ps_per_rank(i,c,n_cll);
//         }
//         return acc;
//     }


//     /** Generate particle ids for one cell, implementing id of kth particle
//      * in cell j for rank r given by
//      *
//      *     cumN + rank_offset + k
//      *
//      * where cumN = number of all particles in cells < j
//      *       rank_offset = number of all particles in cell j with rank < p
//      *       k = all k in [0, ps_per_rank r )
//      */

//     template <class It_t>
//     inline
//     void
//     ids_per_cell(Rank const & rk, CommSz const & csz, size_t const n_cum,
//                  id_t const nc,
//                  It_t const & id_begin,
//                  It_t const & id_end)
//     {
//         id_t const ppr = ps_per_rank<id_t>(fromRank(rk),fromCommSz(csz),nc);
//         GreaterEqual(ppr,id_t(0),"particles per cell");
//         auto dist = std::distance(id_begin,id_end);
//         Equal(ppr,id_t(dist),"particles per cell","number of ids");
//         It_t ids = id_begin;
//         id_t const r_offset = rank_off(rk,csz,nc);

//         for( uint32_t i = 0; i < ppr; ++i)
//         {
//             *ids = i + n_cum + r_offset;
//             ++ids;
//         }
//         return;
//     }

//     /** Compute a set of unique particle ids for a given rank, total
//      * number of ranks, and global statistics. Correctly arranged, the
//      * output for all ranks for a given communicator size should be the
//      * the same as the output for all ranks for a different communicator
//      * size. This is the core property of the invariant partition. */
//     template <typename fp_t, typename sz_t>
//     void
//     gen_ptcl_ids(Rank const & rk, CommSz const & csz,
//                  src_stats_t<fp_t,sz_t> const & stats,
//                  vid & ids)
//     {
//         size_t n_cells = stats.size();
//         typename src_stats_t<fp_t,sz_t>::vsz cumNs(n_cells+1);
//         cumNs[0] = 0;
//         std::partial_sum(stats.ns.begin(),stats.ns.end(),&cumNs[1]);
//         sz_t const totNPs = cumNs[n_cells];
//         if(totNPs != ids.size()){ ids.resize(totNPs);}
//         id_t idx_high = cumNs[0];
//         for(uint32_t i = 0; i < n_cells; ++i)
//         {
//             id_t const idx_low  = idx_high;
//             idx_high = cumNs[i+1];
//             ids_per_cell(rk,csz,idx_low,stats.ns[i],
//                          &ids[idx_low],&ids[idx_high]);
//         }

//     } // gen_ptcl_ids

// } // nut::

#endif // include guard


// End of file
