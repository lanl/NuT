// serialize.hh
// T. M. Kelley
// May 23, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#include "types.hh"
#include "Assert.hh"
#include "utilities_io.hh"
#include <iterator>
#include <string>
#include <ostream>
#include <istream>
#include <sstream>
#include <iomanip>

namespace nut
{
    template <typename stats_t>
    void summarize_stats(stats_t const & stats,
                         std::ostream & o,
                         Species const s,
                         uint32_t const cycle_no)
    {
        typedef typename stats_t::fp_t fp_t;
        fp_t const e_tot   = sum(stats.es);
        fp_t const n_p     = sum(stats.ns);
        fp_t const mn_e_wt = sum(stats.ews)/stats.size();
        std::string sname( species_name(s) );
        o << "For species " << sname << ":\n" << std::setprecision(15)
          << "\t" << n_p << " particles\n"
          << "\t" << e_tot << " ergs total energy\n"
          << "\t" << mn_e_wt << " mean energy weight\n"
          << "--------------------------------------------------"
          << std::endl;
    } // summarize_stats

    template <typename Tally_T>
    void summarize_tally(Tally_T const & tally,
                         std::ostream & o,
                         Species const s,
                         uint32_t const cycle_no)
    {
        typedef typename Tally_T::FP_T fp_t;
        typedef typename Tally_T::cntr_t cntr_t;

        static const size_t dim = Tally_T::dim;

        std::string sname = species_name(s);
        fp_t const e_tot   = sum(tally.energy);
        vec_t<dim> const mom_tot = sum(tally.momentum);
        cntr_t const n_nucl_el     = sum(tally.n_nucl_el_scat);
        cntr_t const n_cell_bdy    = sum(tally.n_cell_bdy);
        cntr_t const n_reflections = sum(tally.n_reflect);
        cntr_t const n_escapes     = sum(tally.n_escape);
        cntr_t n_lepton(0), n_nucl_abs(0);
        switch(s)
        {
        case nu_e:
            n_lepton   = sum(tally.n_nu_e_el_scat);
            n_nucl_abs = sum(tally.n_nu_e_nucl_abs);
            break;
        case nu_e_bar:
            n_lepton = sum(tally.n_nu_e_bar_pos_scat);
            n_nucl_abs = sum(tally.n_nu_e_bar_nucl_abs);
            break;
        case nu_x:
            n_lepton = sum(tally.n_nu_x_el_scat);
            n_nucl_abs = sum(tally.n_nu_x_nucl_abs);
            break;
        default:
            std::cerr << "summarize_tally: unhandled species " << sname
                      << std::endl;
        }
        o << "==================================================\n"
          << "Results for cycle " << cycle_no << ", " << sname << ":\n" << std::setprecision(15)
          << "Total energy deposited: " << e_tot << "\n"
          << "Net radial momentum deposited: " << mom_tot << "\n"
          << "Total path length: " << tally.path_length << "\n"
          << "Scatters:\n"
          << "\tnucleon elastic: " << n_nucl_el << "\n"
          << "\tlepton scatters: " << n_lepton << "\n"
          << "Absorptions:\n"
          << "\tnucleon absorptions: " << n_nucl_abs << "\n"
          << "Mesh:\n"
          << "\tcell boundary crossings: " << n_cell_bdy << "\n"
          << "\treflections: " << n_reflections << "\n"
          << "\tescapes: " << n_escapes << "\n"
          << "Timeouts:\n"
          << "\ttimeouts: " << tally.total_n_census() << "\n"
          << "Total number of MC steps: " << tally.total_mc_steps()
          << "\n"
            ;
    } // summarize_tally

    template <typename Tally_T>
    void write_deposition(std::ostream & o,
                          Tally_T const & tally
                          )
    {
        // header
        o << "\nCell        Momentum                Energy    \n";
        // body of deposition
        for(uint32_t cidx = 0; cidx < tally.n_cells(); ++cidx)
        {
            o << std::setw(4) << std::setfill('0') << cidx << ": "
              << std::setprecision(15) << tally.energy[cidx] << ", "
              << tally.momentum[cidx] << "\n";
        }
        return;
    } // write_deposition

    /** write tally in style of mcphd*/
    template <typename Tally_T>
    void write_tally_mcphd(std::ostream & o,
                           Tally_T & tally,
                           Species const s,
                         uint32_t const cycle_no)
    {
        summarize_tally(tally,o,s,cycle_no);
        write_deposition(o,tally);
        return;
    }

    template <typename Tally_T>
    void write_tally(std::ostream & o,
                     Tally_T & tally,
                     Species const s,
                     uint32_t const cycle_no)
    {
        typedef typename Tally_T::FP_T fp_t;
        typedef typename Tally_T::cntr_t cntr_t;
        typedef typename Tally_T::vf::iterator fpItT;
        typedef typename Tally_T::vc::iterator cItT;

        std::string const n = species_name(s);
        o << "# Tally for cycle " << cycle_no
          << "\n# " << tally.m_n_cells << " cells\n";

        print_per_cell_field<fp_t,fpItT>(
            "energy",tally.energy.begin(),
            tally.energy.end(),o);
        print_per_cell_field<fp_t,fpItT>(
            "momentum",tally.momentum.begin(),
            tally.momentum.end(),o);

        // number and weight created in other fields
        print_per_cell_field<cntr_t,cItT>(
            "n_n",tally.n_n.begin(),
            tally.n_n.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_p",tally.n_p.begin(),
            tally.n_p.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_e_minus",tally.n_e_minus.begin(),
            tally.n_e_minus.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_e_plus",tally.n_e_plus.begin(),
            tally.n_e_plus.end(),o);

        print_per_cell_field<fp_t,fpItT>(
            "ew_n",tally.ew_n.begin(),
            tally.ew_n.end(),o);
        print_per_cell_field<fp_t,fpItT>(
            "ew_p",tally.ew_p.begin(),
            tally.ew_p.end(),o);
        print_per_cell_field<fp_t,fpItT>(
            "ew_e_minus",tally.ew_e_minus.begin(),
            tally.ew_e_minus.end(),o);
        print_per_cell_field<fp_t,fpItT>(
            "ew_e_plus",tally.ew_e_plus.begin(),
            tally.ew_e_plus.end(),o);

        // event counters
        print_per_cell_field<cntr_t,cItT>(
            "n_escape",tally.n_escape.begin(),
            tally.n_escape.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_reflect",tally.n_reflect.begin(),
            tally.n_reflect.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_cell_bdy",tally.n_cell_bdy.begin(),
            tally.n_cell_bdy.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_cutoff",tally.n_cutoff.begin(),
            tally.n_cutoff.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_nucl_el_scat",tally.n_nucl_el_scat.begin(),
            tally.n_nucl_el_scat.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_nu_e_el_scat",tally.n_nu_e_el_scat.begin(),
            tally.n_nu_e_el_scat.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_nu_e_bar_pos_scat",tally.n_nu_e_bar_pos_scat.begin(),
            tally.n_nu_e_bar_pos_scat.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_nu_x_el_scat",tally.n_nu_x_el_scat.begin(),
            tally.n_nu_x_el_scat.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_nu_x_bar_pos_scat",tally.n_nu_x_bar_pos_scat.begin(),
            tally.n_nu_x_bar_pos_scat.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_nu_e_nucl_abs",tally.n_nu_e_nucl_abs.begin(),
            tally.n_nu_e_nucl_abs.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_nu_e_bar_nucl_abs",tally.n_nu_e_bar_nucl_abs.begin(),
            tally.n_nu_e_bar_nucl_abs.end(),o);

        // some associated energy weights (???)
        print_per_cell_field<fp_t,fpItT>(
            "ew_nu_e_nucl_abs",tally.ew_nu_e_nucl_abs.begin(),
            tally.ew_nu_e_nucl_abs.end(),o);
        print_per_cell_field<fp_t,fpItT>(
            "ew_nu_e_bar_nucl_abs",tally.ew_nu_e_bar_nucl_abs.begin(),
            tally.ew_nu_e_bar_nucl_abs.end(),o);

        // census counts & weights
        print_per_cell_field<cntr_t,cItT>(
            "n_census_nu_e",tally.n_census_nu_e.begin(),
            tally.n_census_nu_e.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_census_nu_e_bar",tally.n_census_nu_e_bar.begin(),
            tally.n_census_nu_e_bar.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_census_nu_x",tally.n_census_nu_x.begin(),
            tally.n_census_nu_x.end(),o);
        print_per_cell_field<cntr_t,cItT>(
            "n_census_nu_x_bar",tally.n_census_nu_x_bar.begin(),
            tally.n_census_nu_x_bar.end(),o);

        print_per_cell_field<fp_t,fpItT>(
            "ew_census_nu_e",tally.ew_census_nu_e.begin(),
            tally.ew_census_nu_e.end(),o);
        print_per_cell_field<fp_t,fpItT>(
            "ew_census_nu_e_bar",tally.ew_census_nu_e_bar.begin(),
            tally.ew_census_nu_e_bar.end(),o);
        print_per_cell_field<fp_t,fpItT>(
            "ew_census_nu_x",tally.ew_census_nu_x.begin(),
            tally.ew_census_nu_x.end(),o);
        print_per_cell_field<fp_t,fpItT>(
            "ew_census_nu_x_bar",tally.ew_census_nu_x_bar.begin(),
            tally.ew_census_nu_x_bar.end(),o);

        // escape spectrum
        print_esc_spectrum<typename Tally_T::esc,
                           typename Tally_T::vesc::iterator>(
            "escape_spectrum", tally.escape_spectrum.begin(),
            tally.escape_spectrum.end(),o);
        o << "# End of tally for cycle " << cycle_no << std::endl;
        return;
    }


    template <typename Particle_T>
    void write_particles(std::ostream & o,
                         std::vector<Particle_T> const & ps,
                         Species const s)
    {
        std::string const n = species_name(s);
        o << "## " << n << "particle data\nn_particles = "
          << ps.size() << "\n";
        o << "## x \t omega \t weight \t energy \t cell \t species\n";
        o <<  n << "_pdata = [ \n";
        for(size_t i = 0; i < ps.size(); i++)
        {
            Particle_T const & p = ps[i];
            o << "    (" << p.x << ", " << p.omega << ", " << p.weight << ", "
              << p.e << ", " << p.cell << ", " << n << "),\n";
        }
        o << "    \n";
        return;
    } // write_particles

    template <typename fp_t>
    struct P_arrays
    {
        explicit P_arrays( size_t n)
            : xs(n,fp_t(0)),
              os(n,fp_t(0)),
              ws(n,fp_t(0)),
              es(n,fp_t(0)),
              cells(n,cell_t(0)),
              specs(n)
            {}
        std::vector<fp_t> xs,os,ws,es;
        std::vector<cell_t> cells;
        std::vector<Species> specs;
    }; // P_arrays

    void advance_to_string(std::istream & i, const char * s)
    {
        size_t const linesz(256);
        char line[linesz];
        bool advance(true);
        std::string sline;
        sline.reserve(linesz);
        while(advance)
        {
            i.getline(line,linesz);
            sline.assign(line);
            if(sline.find(s) != std::string::npos)
            {
                advance = false;
            }
        }
        return;
    }

    // template <typename Particle_T>
    // P_arrays
    // read_particles(std::istream & i)
    // {
    //     typedef typename Particle_T::fp_t fpt;

    //     advance_to_string(i,"particle data");

    //     std::stringstream s = o.getline();
    //     size_t n_particles(0u);
    //     s >> n_particles;
    //     Insist(n_particles > 0, "read_particles: cannot read 0 particles!");

    //     P_arrays pas(n_particles);



    //     // skip line with column headers
    //     o.getline(); o.getline();
    //     fpt x,o,w,e;
    //     cell_t c;
    //     std::string sname;
    //     for(size_t k = 0; k < n_particles; ++k)
    //     {
    //         i >> x >> o >> w >> e >> c >> sname;
    //         pas.xs[k] = x;
    //         pas.os[k] = o;
    //         pas.ws[k] = w;
    //         pas.es[k] = e;
    //         pas.cells[k] = c;
    //         pas.specs[k] = spec_type(sname);
    //     }
    //     o.getline(); // skip termination of Python list

    //     return pas;

    // } // read_particles

    // template <typename fp_t>
    // Density
    // read_density(std::istream & i)
    // {
    //     advance_to_string(i,"density class");


    // }

} // nut::


// version
// $Id$

// End of file
