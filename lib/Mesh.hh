// Mesh.hh
// T. M. Kelley
// Jan 11, 2011
// Header for Mesh
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef MESH_H
#define MESH_H

#include "lorentz.hh"
#include "soft_equiv.hh"
#include "utilities_io.hh"
#include "Assert.hh"
#include "constants.hh"
#include "types.hh"
#include <vector>
#include <string>
#include <cmath>


namespace nut
{
    /*!\brief Mesh functions for 1D spherical geometry.
     * \tparam <cell_t> {cell index type}
     * \tparam <boundary_t> {geometry (numerical) type}
     * \tparam <bdy_descriptor_t> {boundary descriptor type}
     */
    template <typename cell_t,typename geometry_t, typename bdy_descriptor_t>
    struct Sphere_1D
    {
    public:
        static const size_t dim = 1;

        typedef geometry_t geom_t;
        typedef bdy_descriptor_t bdy_desc_t;
        typedef std::vector<geom_t> vb;
        typedef std::vector<bdy_desc_t> vbd;
        typedef std::pair<geom_t,geom_t> extents_t;

        struct coord_t
        {
            vec_t<dim> x;
            vec_t<dim> omega;
        };
        struct d_to_b_t
        {
            geom_t d;    /** distance to boundary */
            cell_t face; /** which face will be intersected: 0 (low) or 1 (high) */
        };


        //ctor
        Sphere_1D(vb const & bdys_, vbd const & descs_)
            : m_bdys(bdys_), m_descs(descs_), m_ncells(m_bdys.size() - 1){
            nut::Require(m_bdys.size() >= 2,"must have at least two boundaries");
            nut::Equal(m_bdys.size(),m_descs.size(),"bdys size","descs size");
        }


        cell_t
        n_cells() const {return m_ncells;}


        /*!\brief volume of spherical shell 'cell'. 0 < cell <= n_cells. */
        geom_t
        volume(cell_t const cell) const {
            cellOK(cell);
            cell_t const index(cell-1);
            geom_t const lo  = m_bdys.at(index);
            geom_t const hi  = m_bdys.at(index+1);
            geom_t const vol = 4./3.*pi*(hi*hi*hi - lo*lo*lo);
            nut::GreaterThan(vol,geom_t(0),"volume");
            return vol;
        } // volume


        /*!\brief which cell is across a face. Computed.
         *
         * \param cell: 0 < cell <= n_cells
         * \param face: 0 (left) or 1 (right)
         */
        cell_t
        cell_across_face(cell_t const cell,cell_t face) const
        {
            using namespace bdy_types;
            cellOK(cell);
            LessThan(face,cell_t(2),"face");
            // compute index into descriptors
            cell_t idx = cell - 1;
            idx += face;
            bdy_desc_t const & btype = m_descs[idx];
            cell_t result(-1);
            switch(btype)
            {
            case descriptor::T:
                result = cell - 1 + 2*face;
                break;
            case descriptor::R:
                result = cell;
                break;
            case descriptor::V:
                result = cell_t(0);
                break;
            default:
                std::stringstream errstr;
                errstr << "Sphere_1D::cell_across_face" << __LINE__
                       << " unknown boundary type: "<< btype;
                Require(false,errstr.str().c_str());
            }
            return result;
        } // cell_across_face

        geom_t
        sample_position( geom_t const urd,cell_t const cell) const {
            extents_t extents = this->cell_extents(cell);
            geom_t const lo = extents.first;
            geom_t const hi = extents.second;
            geom_t const lo3 = lo * lo * lo;
            geom_t const hi3 = hi * hi * hi;
            geom_t const r = std::pow( lo3 + (hi3 - lo3)*urd, 1.0/3.0);
            return r;
        }

        template <typename RNG_T>
        vec_t<dim>
        static
        sample_direction( RNG_T & rng)
        {
            geom_t const ctheta = geom_t(2)*rng.random()-geom_t(1);
            return ctheta;
        }

        /*!\brief type of cell boundary */
        bdy_desc_t
        bdy_type(cell_t const c, cell_t const face) const {
            cell_t const idx = make_idx(c,m_ncells);
            cell_t const fidx = face_to_index(face);
            return m_descs[idx + fidx];
        }

        /*!\brief get the lower and upper bounds of a cell. *
         * 0 < cell <= n_cells                              */
        extents_t
        cell_extents(cell_t const cell) const {
            cellOK(cell);
            return extents_t(m_bdys[cell-1],m_bdys[cell]);
        }


        /*!\brief compute the face part of the index into m_descs. */
        cell_t face_to_index(cell_t const f) const {
            LessThan(f,2u,"Mesh1D::face_to_index: face");
            return f;
        }


        /*!\brief calculate new coordinate and new direction cosine at a given
         *        distance along direction cosine omega.  */
        coord_t
        new_coordinate(vec_t<dim> const x, vec_t<dim> const omega,
                       geom_t const distance) const {
            geom_t const theta = std::acos(omega.v[0]);
            geom_t const s     = std::sin(theta);
            geom_t const new_x = x.v[0] + distance * omega.v[0];
            geom_t const new_y = distance * s;
            geom_t const new_r = std::sqrt(new_x*new_x + new_y*new_y);
            coord_t coord;
            coord.x.v[0] = new_r;
            coord.omega.v[0] = std::cos(theta - std::asin(distance/new_r*s));
            return std::move(coord);
        }


        d_to_b_t
        distance_to_bdy(vec_t<dim> const x, vec_t<dim> const omega,
                        cell_t const cell) const {
            cellOK(cell);
            extents_t extents = this->cell_extents(cell);
            return dist_to_bdy_impl(x.v[0],omega.v[0],extents.first,extents.second);
        } // distance_to_bdy


        static
        d_to_b_t
        dist_to_bdy_impl(geom_t const x, geom_t const omega,
                         geom_t const rlo, geom_t const rhi) {
            geom_t const rhisq = rhi * rhi;
            geom_t const rlosq = rlo * rlo;
            geom_t const xsq   = x * x;

            // need to (1-2) check for intersection with the two spheres that
            // bound the current cell, and (3) choose the first intersection
            // along the direction of travel.

            // 1. Compute intersections with outer sphere

            // get Tan(theta) from inverting omega=Cos(theta). We work with a
            // reduced polar angle, "abs_theta", in the range 0 <= theta <= pi/2,
            // then select the correct intercept between the line that the
            // particle travels and the spherical shell using the full value of
            // theta.
            geom_t const theta = std::acos(omega);
            geom_t const abs_theta = (theta > pi/2) ?  pi - theta : theta;
            geom_t const t = std::tan(abs_theta);
            geom_t const tsq = t*t;
            geom_t const tcub = tsq * t;
            geom_t const one_plus_tsq = 1 + tsq;
            geom_t const one_on_one_plus_tsq = 1.0/one_plus_tsq;
            // the determinant is r^2 + r^2 * t^2 - x^2 * t^2. Include the
            // square root and divide by 1+Tan(theta)^2 for convenience here.
            geom_t const dethi =
                std::sqrt( rhisq * one_plus_tsq - xsq * tsq) * one_on_one_plus_tsq;
            // These parts don't depend on the radius of the sphere.
            geom_t const xterm1 = x*tsq*one_on_one_plus_tsq;
            geom_t const yterm1 = -x * t + x * tcub * one_on_one_plus_tsq;
            // These are the two intersection with the outer sphere
            geom_t const xhip = xterm1 + dethi;
            geom_t const yhip = yterm1 + t * dethi;
            geom_t const xhim = xterm1 - dethi;
            geom_t const yhim = yterm1 - t * dethi;
            // if the polar angle is less than pi/2, we want the solution
            // that adds the determinant.
            geom_t const xhi = (theta < pi/2) ? xhip : xhim;
            geom_t const yhi = (theta < pi/2) ? yhip : yhim;
            geom_t const d_hi = std::sqrt( (xhi-x)*(xhi-x) + (yhi*yhi));

            // 2. Look for intersection with inner sphere
            // There is if the absolute value of the polar angle is
            // less than atan(1/(b^2 -1)), where b = x/r_lo. In this case, we can
            // just use the "plus" solution, since that's always the closest.
            // Note that you may want to know a real solution exists before
            // computing it--complicates branch removal.
            geom_t d_lo = huge;
            if(rlo > 0.0)
            {
                geom_t const b = x/rlo;
                // if the particle is on the inner sphere, and headed inward...
                if(soft_equiv(b,1.0,1e-11))
                {
                    if(theta > pi/2)
                    {
                        d_lo = 0.0;
                    }
                }
                else
                {
                    geom_t const tan_theta_lim = std::sqrt(1/(b*b - 1));
                    geom_t const theta_lim = std::atan(tan_theta_lim);
                    if(abs_theta <= theta_lim and theta > pi/2)
                    {
                        geom_t const detlo =
                            std::sqrt( rlosq * one_plus_tsq - xsq * tsq)
                            * one_on_one_plus_tsq;
                        geom_t const xlo = xterm1 + detlo;
                        geom_t const ylo = yterm1 + t * detlo;
                        d_lo = std::sqrt( (xlo-x)*(xlo-x) + ylo*ylo);
                    }
                }
            }
            // now select the shorter distance and the corresponding face
            geom_t d_bdy = huge;
            cell_t face  = -1;
            if( d_hi < d_lo)
            {
                d_bdy = d_hi;
                face  = 1;     // will intersect outer sphere
            }
            else
            {
                d_bdy = d_lo;
                face  = 0;     // will intersect inner sphere
            }
            d_to_b_t d2b;
            d2b.d = d_bdy;
            d2b.face = face;
            return d2b;
        } // dist_to_bdy_impl


        static
        inline
        EandOmega<1>
        LT_to_comoving(vec_t<1> const v_lab,
            geom_t const & e_lab,
            vec_t<1> const & omega_lab)
        {
            return spec_1D::LT_to_comoving_sphere1D(v_lab.v[0],e_lab,omega_lab);
        } // LT_to_comoving

        static
        EandOmega<1>
        inline
        LT_to_lab(vec_t<1> const v_lab,
            geom_t const & e_com,
            vec_t<1> const & omega_com)
        {
            return spec_1D::LT_to_lab_sphere1D(v_lab.v[0],e_com,omega_com);
        }

        vb const m_bdys;

        vbd const m_descs;

        cell_t const m_ncells;

        void cellOK(cell_t const cell_idx) const {
            nut::InOpenRange(cell_idx,cell_t(0),m_ncells+1,"cell id");
        }

    }; // Sphere_1D

} // nut::

#endif



// version
// $Id$

// End of file
