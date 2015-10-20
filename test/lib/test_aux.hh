// test_aux.hh
// T. M. Kelley
// Mar 23, 2009
// Header for test_aux
// (c) Copyright 2009 LANSLLC all rights reserved.

#ifndef TEST_AUX_H
#define TEST_AUX_H

#include "soft_equiv.hh"
#include "expect.hh"
#include <algorithm>        // std::equal
#include <iostream>
#include <iomanip>
#include <string>

namespace test_aux
{
    // the test framework:

    /* run a test, and report on it. For SPE, dump everything on stdout. */
    bool test( const char *target, const char *aspect, bool (*test2run)());

    void pre_report( const char* target, const char* aspect);

    void post_report( const char* target, const char* aspect, bool passed);

    /* run test and report on it. Use modern, non-SPE facilities. */
    bool test( std::string const & target, std::string const & aspect,
               bool (*test2run)());

    void pre_report( std::string const & target, std::string const & aspect);

    void post_report( std::string const & target, std::string const & aspect, bool passed);



    // some helper functions

    /** Compare two tallies, making sure that only one field has changed
     * from the reference. */
    template <typename tally_t, typename v_t>
    bool
    check_one_changed(tally_t const & tally, tally_t const & ref,
                      v_t const * const changed_v);

    /** Compare two tallies, making sure that only two fields have changed
     * from the reference. */
    template <typename tally_t, typename v_t1, typename v_t2>
    bool
    check_two_changed(tally_t const & tally, tally_t const & ref,
                      v_t1 const * const chgd_v1,
                      v_t2 const * const chgd_v2);


    /** compare vector v with r*/
    template <typename v_t>
    bool
    check_same(v_t const * const pv, v_t const * const pr)
    {
        v_t const & v(*pv);
        v_t const  & r(*pr);
        bool same(true);
        if(v.size() != r.size())
        {
            same = false;
        }
        same = same and std::equal(v.begin(),v.end(),r.begin());
        return same;
    }

    /** compare vector v with r*/
    template <typename v_t>
    bool
    check_same_v(v_t const * const pv, v_t const * const pr,
               std::string const & name)
    {
        std::cerr << "comparing " << name << std::endl;
        v_t const & v(*pv);
        v_t const  & r(*pr);
        bool same(true);
        if(v.size() != r.size())
        {
            same = false;
            std::cerr << "vector " << name << " has size " << v.size()
                      << "; expected size "
                      << r.size() << std::endl;
        }
        same = same and std::equal(v.begin(),v.end(),r.begin());
        if(!same)
        {
            std::cerr << "vector " << name << " does not agree with "
                      << "expected value using strict equality." << std::endl;
        }
        return same;
    }

    /** if vector v is not the one that changed (c), then compare with r*/
    template <typename v_t1, typename v_t2>
    bool
    check_same_unless_changed(v_t1 const * const v, v_t1 const * const r,
                              v_t2 const * const c)
    {
        return ( (void *)v != (void *)c) ?
            check_same(v,r) : true;
    }

    /** if vector v is not one of those that changed (c1 or c2), then
     * compare with r*/
    template <typename v_t1, typename v_t2, typename v_t3>
    bool
    check_same_unless_2changed(v_t1 const * const v, v_t1 const * const r,
                               v_t2 const * const c1, v_t3 const * const c2)
    {
        return ( (void *)v != (void *)c1 and (void *)v != (void *)c2) ?
            check_same(v,r) : true;
    }

    /** check_same, with user-supplied comparator function */
    template <typename v_t, typename comp>
    bool
    check_same_verb(v_t const * const pv, v_t const * const pr,
                    comp c)
    {
        v_t const & v(*pv);
        v_t const  & r(*pr);
        bool same(true);
        if(v.size() != r.size())
        {
            same = false;
        }
        same = same and std::equal(v.begin(),v.end(),r.begin(),c);
        return same;
    }

    // a verbose comparator for POD's
    template <typename fp_t>
    struct comp_verb
    {
        bool operator()(fp_t const & val, fp_t const & ref){
            bool passed = nut::soft_equiv(val,ref,m_tol);
            if(!passed)
            {
                std::cout << this -> name_
                          << std::setprecision(15)
                          << ": val = " << val << ", ref = " << ref
                          << std::endl;
            }
            return passed;
        }

        explicit comp_verb(std::string const & name, fp_t tol = 1e-7)
            : name_(name),m_tol(tol){}

        comp_verb() : name_(""),m_tol(1e-7) {}

        std::string name_;

        fp_t const m_tol;
    };

    template <typename iter_t>
    struct comp_verb_iter
    {
        typedef typename iter_t::value_type v_t;

        bool operator()(iter_t const & val, iter_t const & ref){
            return check_same_verb(&val,&ref,comp_verb<v_t>(name_,m_tol));
        }

        explicit comp_verb_iter(std::string const & name, v_t tol = 1e-7)
            : name_(name),m_tol(tol){}

        comp_verb_iter() : name_(""),m_tol(1e-7) {}

        std::string name_;

        v_t const m_tol;
    }; // comp_verb_iter


    template <typename tally_t>
    bool tallies_same(tally_t const & t1, tally_t const & t2)
    {
        bool same(true);
        same = same && check_same(&t1.energy,&t2.energy);
        same = same && check_same(&t1.momentum,&t2.momentum);
        same = same && check_same(&t1.n_n,&t2.n_n);
        same = same && check_same(&t1.n_p,&t2.n_p);
        same = same && check_same(&t1.n_e_minus,&t2.n_e_minus);
        same = same && check_same(&t1.n_e_plus,&t2.n_e_plus);
        same = same && check_same(&t1.ew_n,&t2.ew_n);
        same = same && check_same(&t1.ew_p,&t2.ew_p);
        same = same && check_same(&t1.ew_e_minus,&t2.ew_e_minus);
        same = same && check_same(&t1.ew_e_plus,&t2.ew_e_plus);
        same = same && check_same(&t1.n_escape,&t2.n_escape);
        same = same && check_same(&t1.n_reflect,&t2.n_reflect);
        same = same && check_same(&t1.n_cell_bdy,&t2.n_cell_bdy);
        same = same && check_same(&t1.n_cutoff,&t2.n_cutoff);
        same = same && check_same(&t1.n_nucl_el_scat,&t2.n_nucl_el_scat);
        same = same && check_same(&t1.n_nu_e_el_scat,&t2.n_nu_e_el_scat);
        same = same && check_same(&t1.n_nu_e_bar_pos_scat,&t2.n_nu_e_bar_pos_scat);
        same = same && check_same(&t1.n_nu_x_el_scat,&t2.n_nu_x_el_scat);
        same = same && check_same(&t1.n_nu_x_bar_pos_scat,&t2.n_nu_x_bar_pos_scat);
        same = same && check_same(&t1.ew_escaped,&t2.ew_escaped);
        same = same && check_same(&t1.n_nu_e_nucl_abs,&t2.n_nu_e_nucl_abs);
        same = same && check_same(&t1.n_nu_e_bar_nucl_abs,&t2.n_nu_e_bar_nucl_abs);
        same = same && check_same(&t1.n_nu_x_nucl_abs,&t2.n_nu_x_nucl_abs);
        same = same && check_same(&t1.ew_nu_e_nucl_abs,&t2.ew_nu_e_nucl_abs);
        same = same && check_same(&t1.ew_nu_e_bar_nucl_abs,&t2.ew_nu_e_bar_nucl_abs);
        same = same && check_same(&t1.ew_nu_x_nucl_abs,&t2.ew_nu_x_nucl_abs);
        same = same && check_same(&t1.n_census_nu_e,&t2.n_census_nu_e);
        same = same && check_same(&t1.n_census_nu_e_bar,&t2.n_census_nu_e_bar);
        same = same && check_same(&t1.n_census_nu_x,&t2.n_census_nu_x);
        same = same && check_same(&t1.n_census_nu_x_bar,&t2.n_census_nu_x_bar);
        same = same && check_same(&t1.ew_census_nu_e,&t2.ew_census_nu_e);
        same = same && check_same(&t1.ew_census_nu_e_bar,&t2.ew_census_nu_e_bar);
        same = same && check_same(&t1.ew_census_nu_x,&t2.ew_census_nu_x);
        same = same && check_same(&t1.ew_census_nu_x_bar,&t2.ew_census_nu_x_bar);
        same = same && soft_expect(t1.path_length,t2.path_length,"path length");
        return same;
    }


    template <typename tally_t, typename v_t>
    bool
    check_one_changed(tally_t const & tally, tally_t const & ref,
                      v_t const * const changed_v)
    {
        bool same(true);
        // for each tally field, if not the one in question,
        // make sure that the tally and ref fields are the same
        same = same && check_same_unless_changed( &tally.energy,
                                           & ref.energy, changed_v);
        same = same && check_same_unless_changed( &tally.momentum,
                                           & ref.momentum, changed_v);
        same = same && check_same_unless_changed( &tally.ew_escaped,
                                           & ref.ew_escaped, changed_v);
        same = same && check_same_unless_changed( &tally.n_n,
                                           & ref.n_n, changed_v);
        same = same && check_same_unless_changed( &tally.n_p,
                                           & ref.n_p, changed_v);
        same = same && check_same_unless_changed( &tally.n_e_minus,
                                           & ref.n_e_minus, changed_v);
        same = same && check_same_unless_changed( &tally.n_e_plus,
                                           & ref.n_e_plus, changed_v);
        same = same && check_same_unless_changed( &tally.ew_n,
                                           & ref.ew_n, changed_v);
        same = same && check_same_unless_changed( &tally.ew_p,
                                           & ref.ew_p, changed_v);
        same = same && check_same_unless_changed( &tally.ew_e_minus,
                                           & ref.ew_e_minus, changed_v);
        same = same && check_same_unless_changed( &tally.ew_e_plus,
                                           & ref.ew_e_plus, changed_v);
        same = same && check_same_unless_changed( &tally.n_escape,
                                           & ref.n_escape, changed_v);
        same = same && check_same_unless_changed( &tally.n_reflect,
                                           & ref.n_reflect, changed_v);
        same = same && check_same_unless_changed( &tally.n_cell_bdy,
                                           & ref.n_cell_bdy, changed_v);
        same = same && check_same_unless_changed( &tally.n_cutoff,
                                           & ref.n_cutoff, changed_v);
        same = same && check_same_unless_changed( &tally.n_nucl_el_scat,
                                           & ref.n_nucl_el_scat, changed_v);
        same = same && check_same_unless_changed( &tally.n_nu_e_el_scat,
                                           & ref.n_nu_e_el_scat, changed_v);
        same = same && check_same_unless_changed( &tally.n_nu_e_bar_pos_scat,
                                           & ref.n_nu_e_bar_pos_scat, changed_v);
        same = same && check_same_unless_changed( &tally.n_nu_x_el_scat,
                                           & ref.n_nu_x_el_scat, changed_v);
        same = same && check_same_unless_changed( &tally.n_nu_x_bar_pos_scat,
                                           & ref.n_nu_x_bar_pos_scat, changed_v);
        same = same && check_same_unless_changed( &tally.n_nu_e_nucl_abs,
                                           & ref.n_nu_e_nucl_abs, changed_v);
        same = same && check_same_unless_changed( &tally.n_nu_e_bar_nucl_abs,
                                           & ref.n_nu_e_bar_nucl_abs, changed_v);
        same = same && check_same_unless_changed( &tally.ew_nu_e_nucl_abs,
                                           & ref.ew_nu_e_nucl_abs, changed_v);
        same = same && check_same_unless_changed( &tally.ew_nu_e_bar_nucl_abs,
                                           & ref.ew_nu_e_bar_nucl_abs, changed_v);
        same = same && check_same_unless_changed( &tally.n_census_nu_e,
                                           & ref.n_census_nu_e, changed_v);
        same = same && check_same_unless_changed( &tally.n_census_nu_e_bar,
                                           & ref.n_census_nu_e_bar, changed_v);
        same = same && check_same_unless_changed( &tally.n_census_nu_x,
                                           & ref.n_census_nu_x, changed_v);
        same = same && check_same_unless_changed( &tally.n_census_nu_x_bar,
                                           & ref.n_census_nu_x_bar, changed_v);
        same = same && check_same_unless_changed( &tally.ew_census_nu_e,
                                           & ref.ew_census_nu_e, changed_v);
        same = same && check_same_unless_changed( &tally.ew_census_nu_e_bar,
                                           & ref.ew_census_nu_e_bar, changed_v);
        same = same && check_same_unless_changed( &tally.ew_census_nu_x,
                                           & ref.ew_census_nu_x, changed_v);
        same = same && check_same_unless_changed( &tally.ew_census_nu_x_bar,
                                           & ref.ew_census_nu_x_bar, changed_v);

        return same;
    } // check_one_changed


    template <typename tally_t, typename v_t1, typename v_t2>
    bool
    check_two_changed(tally_t const & tally, tally_t const & ref,
                      v_t1 const * const chgd_v1,
                      v_t2 const * const chgd_v2)
    {
        bool same(true);
        // for each tally field, if not the one in question,
        // make sure that the tally and ref fields are the same
        same = same && check_same_unless_2changed( &tally.energy,
                                            & ref.energy, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.momentum,
                                            & ref.momentum, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_escaped,
                                            & ref.ew_escaped, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_n,
                                            & ref.n_n, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_p,
                                            & ref.n_p, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_e_minus,
                                            & ref.n_e_minus, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_e_plus,
                                            & ref.n_e_plus, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_n,
                                            & ref.ew_n, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_p,
                                            & ref.ew_p, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_e_minus,
                                            & ref.ew_e_minus, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_e_plus,
                                            & ref.ew_e_plus, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_escape,
                                            & ref.n_escape, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_reflect,
                                            & ref.n_reflect, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_cell_bdy,
                                            & ref.n_cell_bdy, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_cutoff,
                                            & ref.n_cutoff, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_nucl_el_scat,
                                            & ref.n_nucl_el_scat, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_nu_e_el_scat,
                                            & ref.n_nu_e_el_scat, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_nu_e_bar_pos_scat,
                                            & ref.n_nu_e_bar_pos_scat, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_nu_x_el_scat,
                                            & ref.n_nu_x_el_scat, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_nu_x_bar_pos_scat,
                                            & ref.n_nu_x_bar_pos_scat, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_nu_e_nucl_abs,
                                            & ref.n_nu_e_nucl_abs, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_nu_e_bar_nucl_abs,
                                            & ref.n_nu_e_bar_nucl_abs, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_nu_e_nucl_abs,
                                            & ref.ew_nu_e_nucl_abs, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_nu_e_bar_nucl_abs,
                                            & ref.ew_nu_e_bar_nucl_abs, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_census_nu_e,
                                            & ref.n_census_nu_e, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_census_nu_e_bar,
                                            & ref.n_census_nu_e_bar, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_census_nu_x,
                                            & ref.n_census_nu_x, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.n_census_nu_x_bar,
                                            & ref.n_census_nu_x_bar, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_census_nu_e,
                                            & ref.ew_census_nu_e, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_census_nu_e_bar,
                                            & ref.ew_census_nu_e_bar, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_census_nu_x,
                                            & ref.ew_census_nu_x, chgd_v1, chgd_v2);
        same = same && check_same_unless_2changed( &tally.ew_census_nu_x_bar,
                                            & ref.ew_census_nu_x_bar, chgd_v1, chgd_v2);

        return same;

    } // check_two_changed




} // test_aux




#endif



// version
// $Id$

// End of file
