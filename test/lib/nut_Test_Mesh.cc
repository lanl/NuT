//T. M. Kelley (c) 2011 LANS LLC

#include "Mesh.hh"
#include "types.hh"
#include "nut_Test_Mesh.hh"
#include "test_aux.hh"
#include "expect.hh"
#include <vector>
#include <deque>
#include <algorithm>
#include <numeric>   // accumulate
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>


namespace Nut_Test
{
    namespace Mesh_tests
    {
        // target describes the code being tested
        char target[] = "Mesh";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "Sphere_1D: inst. & init.";
        bool test_2();
        char aspect2[] = "Sphere_1D: num_cells";
        bool test_3();
        char aspect3[] = "Sphere_1D: volume";
        bool test_4();
        char aspect4[] = "Sphere_1D: cell_across_face";
        bool test_5();
        char aspect5[] = "Sphere_1D: distance_to_boundary, theta = 0, cell 1";
        bool test_6();
        char aspect6[] = "Sphere_1D: distance_to_boundary, theta = pi, cell 1";
        bool test_7();
        char aspect7[] = "Sphere_1D: distance_to_boundary, theta = pi, cell 2";
        bool test_8();
        char aspect8[] = "Sphere_1D: distance_to_boundary cell 2, grazing inner sphere";
        bool test_9();
        char aspect9[] = "Sphere_1D: distance_to_boundary_impl: 100 random tests created with Mathematica";

        bool test_10();
        char aspect10[] = "Sphere_1D: sample position in cell";

    }

    bool test_Mesh()
    {
        using namespace Mesh_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        bool passed2 = test( target, aspect2, test_2);

        bool passed3 = test( target, aspect3, test_3);

        bool passed4 = test( target, aspect4, test_4);

        bool passed5 = test( target, aspect5, test_5);

        bool passed6 = test( target, aspect6, test_6);

        bool passed7 = test( target, aspect7, test_7);

        bool passed8 = test( target, aspect8, test_8);

        bool passed9 = test( target, aspect9, test_9);

        bool passed10 = test( target, aspect10, test_10);

        // call additional tests here.

        return passed1 and passed2 and passed3 and passed4 and passed5
            and passed6 and passed7 and passed8 and passed9 and passed10;
    }

    namespace Mesh_tests
    {
        using test_aux::soft_expect;
        using test_aux::expect;
        using nut::soft_equiv;

        bool test_1()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,35.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            return passed;
        }


        bool test_2()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,35.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            passed = passed and
                mesh.n_cells() == n_cells;

            return passed;
        } // test_2


        bool test_3()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,35.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            geom_t const vols_exp[]={4.1887902047863905, 29.321531433504735,
                                     179560.86970857819};
            Sp1D::vb vols(n_cells);
            vols[0] = mesh.volume(1);
            vols[1] = mesh.volume(2);
            vols[2] = mesh.volume(3);

            passed = passed and
                std::equal(vols.begin(),vols.end(),&vols_exp[0]);
            if(!passed)
            {
                std::cout << std::setprecision(16) << std::scientific
                          << "volumes: ";
                std::copy(vols.begin(),vols.end(),
                          std::ostream_iterator<geom_t>(std::cout,","));
                std::cout << std::endl << "expected: ";
                std::copy(&vols_exp[0],&vols_exp[n_cells],
                          std::ostream_iterator<geom_t>(std::cout,","));
                std::cout << std::endl;
            }
            return passed;
        } // test_3


        bool test_4()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,35.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            cell_t const cells_a_exp[]={1,2, 1,3, 2,0};
            std::vector<cell_t> cells_across(n_cells*2);
            cells_across[0] = mesh.cell_across_face(1,0);
            cells_across[1] = mesh.cell_across_face(1,1);
            cells_across[2] = mesh.cell_across_face(2,0);
            cells_across[3] = mesh.cell_across_face(2,1);
            cells_across[4] = mesh.cell_across_face(3,0);
            cells_across[5] = mesh.cell_across_face(3,1);

            passed = passed and
                std::equal(cells_across.begin(),cells_across.end(),&cells_a_exp[0]);
            if(!passed)
            {
                std::cout << std::setprecision(16) << std::scientific
                          << "cells across: ";
                std::copy(cells_across.begin(),cells_across.end(),
                          std::ostream_iterator<cell_t>(std::cout,","));
                std::cout << std::endl << "expected: ";
                std::copy(&cells_a_exp[0],&cells_a_exp[n_cells*2],
                          std::ostream_iterator<cell_t>(std::cout,","));
                std::cout << std::endl;
            }
            return passed;
        } // test_4


        bool test_5()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,35.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            Sp1D::d_to_b_t d_n_face = mesh.distance_to_bdy(0.5,1.0,1);

            passed = passed and soft_expect(d_n_face.d,0.5,"distance");
            passed = passed and expect(d_n_face.face,1u,"face");

            return passed;
        } // test_5


        bool test_6()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,35.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            Sp1D::d_to_b_t d_n_face = mesh.distance_to_bdy(0.5,-1.0,1);

            passed = passed and soft_expect(d_n_face.d,1.5,"distance");
            passed = passed and expect(d_n_face.face,1u,"face");

            return passed;
        } // test_6


        bool test_7()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,35.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            Sp1D::d_to_b_t d_n_face = mesh.distance_to_bdy(1.5,-1.0,2);

            passed = soft_expect(d_n_face.d,0.5,"distance") and passed;
            passed = expect(d_n_face.face,0u,"face") and passed;

            return passed;
        } // test_7


        bool test_8()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,35.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            cell_t const cell = 2;

            geom_t const x = 2;
            geom_t const omega = -std::sqrt(3)/2;
            geom_t const d_exp = std::sqrt(3);
            // This test has limited tolerance because of cancellation in
            // the code when computing the determinant.  The determinant is
            // 1.0*1.33333333333 - 4.0*0.33333333333. The exact value is
            // 1*4/3-4*1/3 = 0.  The floating point version fails pretty badly.
            geom_t const tol   = 3.0e-8;

            Sp1D::d_to_b_t d_n_face = mesh.distance_to_bdy(x,omega,cell);

            passed = soft_expect(d_n_face.d,d_exp,"distance",tol) and passed;
            passed = expect(d_n_face.face,0u,"face") and passed;

            return passed;
        } // test_8

        namespace
        {
            template <typename vb_t, typename vc_t>
            bool read_input_t9(vb_t &xs, vb_t &omegas, vb_t &r_los, vb_t &r_his,
                               vb_t & d_exps, vc_t &f_exps)
            {
                typedef typename vb_t::value_type geom_t;
                typedef typename vc_t::value_type cell_t;

                std::string const fname("spherical-d-to-b-tests-cxx.txt");
                std::ifstream inf(fname.c_str());
                if(!inf.good())
                {
                    std::cerr << "mesh test_9 cannot read file " << fname
                              << std::endl;
                    std::cerr << "Are you running in the directory in which "
                              << "the tests were built?" << std::endl;
                    return false;
                }

                uint32_t count(0);
                while(1)
                {
                    geom_t x, omega, r_lo, r_hi, d_exp;
                    cell_t f_exp;
                    inf >> x >> omega >> r_lo >> r_hi >> d_exp >> f_exp;
                    if(inf.eof()) break;
                    xs.push_back(x);
                    omegas.push_back(omega);
                    r_los.push_back(r_lo);
                    r_his.push_back(r_hi);
                    d_exps.push_back(d_exp);
                    f_exps.push_back(f_exp - 1);
                    count++;
                }
                std::cout << "read " << count << " input lines, will run that "
                          << " many problems."
                          << std::endl;
                inf.close();
                return true;
            } // read_input
        } //anonymous

        bool test_9()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            Sp1D::vb xs, omegas, r_los, r_his, d_exps;
            std::vector<cell_t> f_exps;

            if(!read_input_t9(xs, omegas, r_los, r_his, d_exps, f_exps))
            {
                std::cerr << "Failed to read input, this test SKIPPED" << std::endl;
                return true;
            }
            // This test has limited tolerance because of cancellation in
            // the code when computing the determinant.  The determinant is
            // 1.0*1.33333333333 - 4.0*0.33333333333. The exact value is
            // 1*4/3-4*1/3 = 0.  The floating point version fails pretty badly.
            geom_t const tol   = 3.0e-8;
            test_aux::soft_eq_bound_tol<geom_t> s_eq(tol);
            Sp1D::vb ds(xs.size());
            std::vector<cell_t> fs(xs.size());
            for(size_t i = 0; i < xs.size();++i)
            {
                geom_t x = xs[i], o=omegas[i], rl=r_los[i], rh=r_his[i];
                Sp1D::d_to_b_t dnf = Sp1D::dist_to_bdy_impl(x,o,rl,rh);
                ds[i] = (dnf.d);
                fs[i] = (dnf.face);
            }

            bool ds_passed = std::equal(ds.begin(),ds.end(),d_exps.begin(),s_eq);
            bool fs_passed = std::equal(fs.begin(),fs.end(),f_exps.begin());
            passed = ds_passed and passed;
            passed = fs_passed and passed;

            return passed;
        } // test_9

        bool test_10()
        {
            bool passed(true);

            using nut::cell_t;
            using nut::geom_t;
            typedef nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor> Sp1D;

            size_t const n_cells(3);
            size_t const n_bdys(n_cells+1);

            geom_t const bdys_in[n_bdys] = {0.0,1.0,2.0,3.0};
            nut::bdy_types::descriptor const bdy_ts[n_bdys] =
                {nut::bdy_types::R, nut::bdy_types::T,
                 nut::bdy_types::T, nut::bdy_types::V};

            Sp1D::vb bdys(&bdys_in[0],&bdys_in[n_bdys]);
            Sp1D::vbd bdy_types(&bdy_ts[0],&bdy_ts[n_bdys]);

            Sp1D mesh(bdys, bdy_types);

            geom_t position = mesh.sample_position(0.1,2);

            // Sp1D::extents_t xs = mesh.cell_extents(2);

            passed = soft_expect(position,1.193483191927337,"sampled position")
                and passed;
            // if(!passed)
            // {
            //     std::cerr << ""
            // }
            return passed;
        } // test_10



        // define additional tests here.

    } // Mesh_tests::

} // Nut_Test::



// version
// $Id$

// End of file
