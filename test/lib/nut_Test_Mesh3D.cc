//T. M. Kelley (c) 2011 LANS LLC

#include "Mesh3DCar.hh"
#include "nut_Test_Mesh3D.hh"
#include "test_aux.hh"
#include "types.hh"
#include "expect.hh"

namespace Nut_Test
{
    namespace Mesh3D_tests
    {
        // target describes the code being tested
        char target[] = "3D Cartesian mesh";

        // for each test of target, add a declaration and
        // a description of the aspect being tested.
        bool test_1();
        char aspect1[] = "instantiate";
        bool test_2();
        char aspect2[] = "instantiate correctly (n_cells)";
        bool test_3();
        char aspect3[] = "instantiate correctly (volume)";
        bool test_4();
        char aspect4[] = "isBoundary: 3 cell line";
        bool test_5();
        char aspect5[] = "isBoundary: 6 cell cluster";
        bool test_6();
        char aspect6[] = "cell_across_face: 6 cell cluster";
        bool test_7();
        char aspect7[] = "sample_position: buffer rng";
        bool test_8();
        char aspect8[] = "distance_to_bdy 1";
    }

    bool test_Mesh3D()
    {
        using namespace Mesh3D_tests;
        using test_aux::test;

        bool passed1 = test( target, aspect1, test_1);

        bool passed2 = test( target, aspect2, test_2);

        bool passed3 = test( target, aspect3, test_3);

        bool passed4 = test( target, aspect4, test_4);

        bool passed5 = test( target, aspect5, test_5);

        bool passed6 = test( target, aspect6, test_6);

        bool passed7 = test( target, aspect7, test_7);

        bool passed8 = test( target, aspect8, test_8);

        // call additional tests here.

        return passed1 and  passed2 and  passed3 and  passed4 and  passed5
            and passed6 and passed7 and passed8;
    }

    namespace Mesh3D_tests
    {
        using nut::geom_t;
        using nut::Equal;
        typedef uint64_t cell_t;
        typedef nut::Cartesian_3D<cell_t,double> mesh_t;
        typedef mesh_t::vbd vbd;
        using test_aux::expect;
        using test_aux::soft_expect;

        bool test_1()
        {
            cell_t const nx = 1;
            cell_t const ny = 1;
            cell_t const nz = 1;
            vbd bds;
            nut::mkReflectBCs<mesh_t,cell_t>(bds,nx,ny,nz);
            mesh_t mesh(1.0,nx,1.0,ny,1.0,nz,bds);
            return true;
        }

        template <typename cell_t>
        bool test_2_core(cell_t const nx,cell_t const ny,cell_t const nz)
        {
            bool passed(true);
            vbd bds;
            nut::mkReflectBCs<mesh_t,cell_t>(bds,nx,ny,nz);
            mesh_t mesh(1.0,nx,1.0,ny,1.0,nz,bds);
            cell_t const n_cs(nx*ny*nz);
            passed = passed && n_cs == mesh.n_cells();
            return passed;
        }

        bool test_2()
        {
            bool passed(true);
            cell_t const nx = 3;
            cell_t const ny = 5;
            cell_t const nz = 7;
            passed = passed && test_2_core(nx,ny,nz);
            return passed;
        } //test_2

        // for any 3D cartesian mesh, the volume should be constant
        // (and simple!)
        template <typename cell_t>
        bool test_3_core(geom_t const dx, cell_t const nx,
                         geom_t const dy, cell_t const ny,
                         geom_t const dz, cell_t const nz)
        {
            bool passed(true);
            vbd bds;
            nut::mkReflectBCs<mesh_t,cell_t>(bds,nx,ny,nz);
            mesh_t mesh(dx,nx,dy,ny,dz,nz,bds);
            geom_t const exp_vol(dx*dy*dz);
            for(cell_t c = 0; c < mesh.n_cells(); ++c)
            {
                geom_t const volume = mesh.volume(c);
                bool cok = volume == exp_vol;
                passed = passed && cok;
            }
            return passed;
        }

        bool test_3()
        {
            bool passed(true);
            geom_t const dx(2.0);
            geom_t const dy(3.0);
            geom_t const dz(5.0);
            cell_t const nx = 3;
            cell_t const ny = 5;
            cell_t const nz = 7;

            passed = passed && test_3_core(dx,nx,dy,ny,dz,nz);

            return passed;
        } //test_3


        bool test_4()
        {
            using test_aux::expect;
            bool passed(true);
            geom_t const dx(2.0);
            geom_t const dy(3.0);
            geom_t const dz(5.0);

            cell_t const nx = 3;
            cell_t const ny = 1;
            cell_t const nz = 1;

            vbd bds;
            nut::mkReflectBCs<mesh_t,cell_t>(bds,nx,ny,nz);
            mesh_t mesh(dx,nx,dy,ny,dz,nz,bds);

            {
                cell_t c(0);
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x), true,"c0 low_x" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_x),false,"c0 high_x")&&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),true,"c0 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),true,"c0 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),true,"c0 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),true,"c0 high_z");
            }
            {
                cell_t c(1);
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x ),false,"c1 low_x" )&&
                    expect(mesh.isBoundary(c,mesh_t::high_x),false,"c1 high_x")&&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),true,"c1 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),true,"c1 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),true,"c1 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),true,"c1 high_z");
            }
            {
                cell_t c(2);
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x ),false,"c2 low_x" )&&
                    expect(mesh.isBoundary(c,mesh_t::high_x),true,"c2 high_x") &&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),true,"c2 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),true,"c2 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),true,"c2 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),true,"c2 high_z");
            }
            return passed;
        } // test_4

        bool test_5()
        {
            using test_aux::expect;
            using nut::ijk_t;
            bool passed(true);
            geom_t const dx(2.0);
            geom_t const dy(3.0);
            geom_t const dz(5.0);

            cell_t const nx = 1;
            cell_t const ny = 3;
            cell_t const nz = 2;

            vbd bds;
            nut::mkReflectBCs<mesh_t,cell_t>(bds,nx,ny,nz);
            mesh_t mesh(dx,nx,dy,ny,dz,nz,bds);

            {
                cell_t c(0); // (0,0,0)
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x), true,"c0 low_x" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_x),true,"c0 high_x")&&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),true,"c0 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),false,"c0 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),true,"c0 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),false,"c0 high_z");
            }
            {
                cell_t c(1); // (0,1,0)
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x ),true,"c1 low_x" )&&
                    expect(mesh.isBoundary(c,mesh_t::high_x),true,"c1 high_x")&&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),false,"c1 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),false,"c1 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),true,"c1 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),false,"c1 high_z");
            }
            {
                cell_t c(2); // (0,2,0)
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x ),true,"c2 low_x" )&&
                    expect(mesh.isBoundary(c,mesh_t::high_x),true,"c2 high_x") &&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),false,"c2 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),true,"c2 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),true,"c2 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),false,"c2 high_z");
            }
            {
                cell_t c(3); // (0,0,1)
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x), true,"c3 low_x" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_x),true,"c3 high_x")&&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),true,"c3 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),false,"c3 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),false,"c3 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),true,"c3 high_z");
            }
            {
                cell_t c(4); // (0,1,1)
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x ),true,"c4 low_x" )&&
                    expect(mesh.isBoundary(c,mesh_t::high_x),true,"c4 high_x")&&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),false,"c4 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),false,"c4 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),false,"c4 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),true,"c4 high_z");
            }
            {
                cell_t c(5); // (0,2,1)
                passed = passed &&
                    expect(mesh.isBoundary(c,mesh_t::low_x ),true,"c5 low_x" )&&
                    expect(mesh.isBoundary(c,mesh_t::high_x),true,"c5 high_x") &&
                    expect(mesh.isBoundary(c,mesh_t::low_y ),false,"c5 low_y" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_y),true,"c5 high_y") &&
                    expect(mesh.isBoundary(c,mesh_t::low_z ),false,"c5 low_z" ) &&
                    expect(mesh.isBoundary(c,mesh_t::high_z),true,"c5 high_z");
            }
            return passed;
        } //test_5

        bool test_6()
        {
            using nut::ijk_t;
            bool passed(true);
            geom_t const dx(2.0);
            geom_t const dy(3.0);
            geom_t const dz(5.0);

            cell_t const nx = 1;
            cell_t const ny = 3;
            cell_t const nz = 2;

            vbd bds;
            nut::mkReflectBCs<mesh_t,cell_t>(bds,nx,ny,nz);
            mesh_t mesh(dx,nx,dy,ny,dz,nz,bds);
            {
                cell_t c(0); // (0,0,0)
                auto f = [c,mesh](mesh_t::face_t const f,cell_t const e, const char *s)
                    {return expect(mesh.cell_across_face(c,f),e,s);};
                passed = passed &&
                    f(mesh_t::low_x, cell_t(0),"c0 low x") &&
                    f(mesh_t::high_x,cell_t(0),"c0 high x") &&
                    f(mesh_t::low_y, cell_t(0),"c0 low y") &&
                    f(mesh_t::high_y,cell_t(1),"c0 high y") &&
                    f(mesh_t::low_z, cell_t(0),"c0 low z") &&
                    f(mesh_t::high_z,cell_t(3),"c0 high z");
            }
            {
                cell_t c(1);
                auto f = [c,mesh](mesh_t::face_t const f,cell_t const e, const char *s)
                    {return expect(mesh.cell_across_face(c,f),e,s);};
                passed = passed &&
                    f(mesh_t::low_x, cell_t(1),"c1 low x") &&
                    f(mesh_t::high_x,cell_t(1),"c1 high x") &&
                    f(mesh_t::low_y, cell_t(0),"c1 low y") &&
                    f(mesh_t::high_y,cell_t(2),"c1 high y") &&
                    f(mesh_t::low_z, cell_t(1),"c1 low z") &&
                    f(mesh_t::high_z,cell_t(4),"c1 high z");
            }
            {
                cell_t c(2);
                auto f = [c,mesh](mesh_t::face_t const f,cell_t const e, const char *s)
                    {return expect(mesh.cell_across_face(c,f),e,s);};
                passed = passed &&
                    f(mesh_t::low_x, cell_t(2),"c1 low x") &&
                    f(mesh_t::high_x,cell_t(2),"c1 high x") &&
                    f(mesh_t::low_y, cell_t(1),"c1 low y") &&
                    f(mesh_t::high_y,cell_t(2),"c1 high y") &&
                    f(mesh_t::low_z, cell_t(2),"c1 low z") &&
                    f(mesh_t::high_z,cell_t(5),"c1 high z");
            }
            {
                cell_t c(3);
                auto f = [c,mesh](mesh_t::face_t const f,cell_t const e, const char *s)
                    {return expect(mesh.cell_across_face(c,f),e,s);};
                passed = passed &&
                    f(mesh_t::low_x, cell_t(3),"c1 low x") &&
                    f(mesh_t::high_x,cell_t(3),"c1 high x") &&
                    f(mesh_t::low_y, cell_t(3),"c1 low y") &&
                    f(mesh_t::high_y,cell_t(4),"c1 high y") &&
                    f(mesh_t::low_z, cell_t(0),"c1 low z") &&
                    f(mesh_t::high_z,cell_t(3),"c1 high z");
            }
            {
                cell_t c(4);
                auto f = [c,mesh](mesh_t::face_t const f,cell_t const e, const char *s)
                    {return expect(mesh.cell_across_face(c,f),e,s);};
                passed = passed &&
                    f(mesh_t::low_x, cell_t(4),"c1 low x") &&
                    f(mesh_t::high_x,cell_t(4),"c1 high x") &&
                    f(mesh_t::low_y, cell_t(3),"c1 low y") &&
                    f(mesh_t::high_y,cell_t(5),"c1 high y") &&
                    f(mesh_t::low_z, cell_t(1),"c1 low z") &&
                    f(mesh_t::high_z,cell_t(4),"c1 high z");
            }
            {
                cell_t c(5);
                auto f = [c,mesh](mesh_t::face_t const f,cell_t const e, const char *s)
                    {return expect(mesh.cell_across_face(c,f),e,s);};
                passed = passed &&
                    f(mesh_t::low_x, cell_t(5),"c1 low x") &&
                    f(mesh_t::high_x,cell_t(5),"c1 high x") &&
                    f(mesh_t::low_y, cell_t(4),"c1 low y") &&
                    f(mesh_t::high_y,cell_t(5),"c1 high y") &&
                    f(mesh_t::low_z, cell_t(2),"c1 low z") &&
                    f(mesh_t::high_z,cell_t(5),"c1 high z");
            }
            return passed;
        } // test_6

        bool test_7()
        {
            bool passed(true);
            typedef nut::Buffer_RNG<double> rng_t;
            size_t const szRns(5);
            double rns[szRns] = {0.5,0.5,0.5,0.5,0.5};
            rng_t rng(rns,szRns);

            geom_t const dx(2.0);
            geom_t const dy(3.0);
            geom_t const dz(5.0);
            cell_t const nx = 1;
            cell_t const ny = 3;
            cell_t const nz = 2;
            vbd bds;
            nut::mkReflectBCs<mesh_t,cell_t>(bds,nx,ny,nz);
            mesh_t mesh(dx,nx,dy,ny,dz,nz,bds);

            mesh_t::coord_t newcrd = mesh.sample_position(rng,0);
            mesh_t::vec_t x = newcrd.x;
            mesh_t::vec_t o = newcrd.omega;

            passed = passed &&
                expect(x.v[0],1.0,"cell 0, x") &&
                expect(x.v[1],1.5,"cell 0, y") &&
                expect(x.v[2],2.5,"cell 0, z") &&
                soft_expect(o.v[0],-1.0,"cell 0, omega_x") &&
                soft_expect(o.v[1], 0.0,"cell 0, omega_y") &&
                soft_expect(o.v[2], 0.0,"cell 0, omega_z");
            return passed;
        } // test_7

        bool test_8()
        {
            bool passed(true);
            typedef nut::Buffer_RNG<double> rng_t;
            size_t const szRns(5);
            double rns[szRns] = {0.5,0.5,0.5,0.5,0.5};
            rng_t rng(rns,szRns);

            geom_t const dx(2.0);
            geom_t const dy(3.0);
            geom_t const dz(5.0);
            cell_t const nx = 1;
            cell_t const ny = 3;
            cell_t const nz = 2;
            vbd bds;
            nut::mkReflectBCs<mesh_t,cell_t>(bds,nx,ny,nz);
            mesh_t mesh(dx,nx,dy,ny,dz,nz,bds);

            mesh_t::coord_t crd = {{0.8,2.9,4.9},{1.0,0.0,0.0}};
            mesh_t::d_to_b_t d2b = mesh.distance_to_bdy(crd,cell_t(0));

            passed = passed
                && expect(d2b.d,1.2,"distance")
                && expect(d2b.face,mesh_t::high_x,"face")
                ;

            return passed;
        } // test_8

        // define additional tests here.

    } // Mesh3D_tests::

} // Nut_Test::


// End of file
