// utilities.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "types.hh"
#include <sstream>
#include <numeric>
#include <cmath>
#include "Assert.hh"

/**!\file commonly used assertions, checks, exceptions, etc. */


namespace nut
{
    // print out arrays Python style
    template <typename T, typename ItT>
    void print_per_cell_field(std::string const & name,
                              ItT begin,
                              ItT end,
                              std::ostream & o)
    {
        o << name << " = [";
        std::copy(begin,end,std::ostream_iterator<T>(o,","));
        o << "]\n";
    }

    template <typename esc, typename vescIt>
    void print_esc_spectrum(std::string const & name,
                            vescIt const begin,
                            vescIt const end,
                            std::ostream & o)
    {
        o << name << " = [";
        for(vescIt curr = begin; curr != end; ++curr)
        {
            o << "(" << (*curr).first << "," 
              << (*curr).second << "),";
        }
        o << "]\n";
    }
    
    template <typename fp_t>
    struct 
    mult{
        fp_t operator()(fp_t const a) const {return a * t;}
        fp_t const t;
        explicit mult(fp_t const t_) : t(t_) {}
    }; // mult

    template <typename fp_t, typename int_t>
    struct 
    int_mult{
        int_t operator()(fp_t const a) const {
            return int_t( std::floor(a * t + 0.5) );
        }
        fp_t const t;
        explicit int_mult(fp_t const t_) : t(t_) {}
    }; // mult

    template <typename fp_t>
    struct 
    div_by{
        fp_t operator()(fp_t const a) const {return a / d;}
        fp_t const d;
        explicit div_by(fp_t const d_) : d(d_) {
            Require(d != fp_t(0),"div_by::ctor: divisor may not equal 0");
        }
    }; // div_by

    /*!\brief if divisor is not equal to 0, returns a / b, otherwise
     * returns 0. */
    template <typename fp_t>
    fp_t div_if_ne0(fp_t const a,size_t const b) 
    {
        return b != 0 ? fp_t(a / b) : fp_t (0);
    }

    // sum contents of an STL container
    template <typename ContT>
    typename ContT::value_type
    sum(ContT const & v)
    {
        return std::accumulate(v.begin(),v.end(),
                               typename ContT::value_type(0));
    }

    // template <typename Num_T,typename ContT>
    // Num_T
    // sum2(ContT const & v)
    // {
    //     return std::accumulate(v.begin(),v.end(),Num_T(0));
    // }


    /**! Assert that the cell c > 0 and c <= n_cells */
    inline void 
    cellOK(cell_t const c, cell_t const n_cells) 
    {
        InOpenRange(c,cell_t(0),n_cells+1,"cell");
    }

    /**! Assert that cell c is ok (above) and return index into arrays */
    inline cell_t
    make_idx(cell_t const c, cell_t const n_cells) 
    {
        cellOK(c,n_cells); 
        return c - 1;
    }


    // throw an exception when an unknown species is encountered
    inline void 
    unhandled_species_exc(Species const s, std::string const & method)
    {
        std::stringstream errstr;
        errstr << method << ": Unhandled neutrino species " << s
               << ", better known as " << species_name(s);
        throw arg_error(errstr.str());            
    }

    template <typename fp_t>
    inline 
    void nrgOK(fp_t const e,char const * const s) 
    {
        GreaterThan(e,fp_t(0),s);
    }
    
    template <typename fp_t>
    inline 
    void nrgOK(fp_t const e) 
    {
        char const * const s = "neutrino energy";
        GreaterThan(e,fp_t(0),s);
    }
    



} // nut::


#endif // include_guard

// version
// $Id$

// End of file
