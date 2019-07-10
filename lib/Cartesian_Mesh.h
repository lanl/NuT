// Cartesian_Mesh.h
// Jan 17, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#pragma once

#include "constants.hh"
#include "detail/Cartesian_Mesh_detail.h"
#include "detail/Cell.h"
#include "detail/Cell_Face_Descriptor.h"
#include "detail/Face.h"
#include "detail/Ray.h"
#include "detail/Vertex.h"
#include "lorentz.hh"
#include "types.hh"
#include <cmath>
#include <iterator>
#include <map>
#include <unordered_set>
#include <vector>

/* Planes are ordered to match directions. */
enum Plane : uint32_t {
  YZ = 0, /*X*/
  XZ = 1, /*Y*/
  XY = 2, /*Z*/
  NUM_PLANES = 3
};

namespace nut {

class Cartesian_Mesh {
public:
  // TYPES
  using Ray = nut::Ray;
  using Vector = nut::Vector;
  using Point = nut::Vector;
  using Cell = nut::Cell;
  using intersect_t = std::pair<Cartesian_Face, geom_t>;
  using face_handle_t = Cartesian_Face;
  using mesh_t = Cartesian_Mesh;

  using Index = index_t;
  using Index_T = Index;
  using Geom_T = double;
  using vector4_t = nut::Vector4<Geom_T>;

  // useful containers
  // types
  using cell_list_t = std::unordered_set<Cell>;
  using vertex_list_t = std::vector<Vertex>;
  using face_list_t = std::vector<face_handle_t>;

  // element to element maps
  using cell_to_vertex_map = std::map<Cell, vertex_list_t>;
  using vertex_to_cell_map = std::map<Vertex, cell_list_t>;
  using cell_to_face_map = std::map<Cell, face_list_t>;

  enum Face_Name : uint32_t { LOW_X = 0, LOW_Y, LOW_Z, HIGH_X, HIGH_Y, HIGH_Z };

  enum Dir_Name : uint32_t { X = 0, Y, Z };

  static bool is_high(Face_Name const f) { return f > 2; }

  static Dir_Name to_dir(Face_Name const f) { return Dir_Name(f % 3); }

  static Plane to_plane(Dir_Name const d) { return static_cast<Plane>(d); }

  // convenient alias for the low- and high-extents of a cell in any direction
  using Extents_T = std::pair<geom_t, geom_t>;

  static constexpr index_t void_cell_idx{max_index_t};

  static const Cell null_cell_;

  static Cell null_cell() { return null_cell_; }

  static const face_handle_t null_face_;

  static face_handle_t null_face() { return null_face_; }

  static Geom_T get_distance(intersect_t const & i) { return i.second; }

  static face_handle_t get_face(intersect_t const & i) { return i.first; }

  /**\brief Compute vector v reflected in face f */
  Vector reflect(face_handle_t const & f, Vector const & v) const;

  /**\brief Compute vector v reflected in face f */
  Vector reflect(face_handle_t const & f,
                 Vector const & v,
                 Point const & /*unused*/) const
  {
    return this->reflect(f, v);
  }

  /**\brief Get inward-facing normal for f */
  Vector get_normal(Cell const & c, face_handle_t const & f) const;

  static Vector get_normal(Face_Name const);

  /**\brief Lorentz transform energy and direction from co-moving to lab
   * frame.
   * \param*/
  static vector4_t LT_to_lab(Vector const & v,
                             Geom_T const & ec,
                             Vector const & oc)
  {
    return nut::spec_3D_Cartesian::LT_to_lab(v, ec, oc);
  }

  static vector4_t LT_to_comoving(Vector const & v,
                                  Geom_T const el,
                                  Vector const & ol)
  {
    return nut::spec_3D_Cartesian::LT_to_comoving(v, el, ol);
  }

public:
  // INTERFACE
  /**\brief Identify cell across the face at given point. */
  Cell cell_across(Cell const & c,
                   face_handle_t const & f,
                   Point const & /*p*/) const;

  /**\brief What is the volume of the indicated cell? */
  Geom_T volume(Cell const & /*c*/) const;

  /**\brief Is this face a boundary of the problem domain? */
  bool is_boundary(face_handle_t f) const;

  /**\brief Find the cell in which the Point lies.
   *
   * Note that this will always return *something*. It will return junk if you
   * provide a point that is not in the mesh. You might want to check the
   * returned cell value with in_cell() (though you may get consistent junk).
   * So you may really want to check valid_cell.
   */
  Cell find_cell(Point const & p) const;

  bool valid_cell(Cell const & c) const;

  /**\brief Compute the x-extents {low, high} for a cell. This is where the
   * actual encoding of the cell index is important--we treat the linear index
   * as encoding position information. */
  Extents_T get_x_extents(Cell const & c) const;

  /**\brief Compute the y-extents {low, high} for a cell. */
  Extents_T get_y_extents(Cell const & c) const;

  /**\brief Compute the z-extents {low, high} for a cell. */
  Extents_T get_z_extents(Cell const & c) const;

  /**\brief Given a direction cosine and current coordinate, compute distance
   * to the coordinate of the face. */
  static geom_t compute_distance(geom_t d_cos,
                                 geom_t coord,
                                 geom_t face_coord_lo,
                                 geom_t face_coord_hi);

  /**\brief Compute intersection face and distance of ray r in Cell c.
   * \param r: current position & direction
   * \param c: current cell, must be consistent with r
   * \note The cell is a convenience, the position in the Ray must lie within
   * the Cell.
   *
   * Note that if a ray is pointed directly at the meeting of two or three
   * faces, then this algorithm always prefers the lesser of {X,Y,Z}. */
  intersect_t intersection(Ray const & r, Cell const & c) const;

  /**\brief A point is in a cell if each of its coordinates is greater than or
   * equal to the cell's lower extents and less than the cell's upper
   * extents.
   * \param p: a point
   * \param c: a cell
   * \return is p in c?
   */
  bool in_cell(Vector const & p, Cell const & c) const;

  /**\brief Compute a point a given distance along a given direction */
  Ray advance_point(Point const & init_point,
                    Vector const & direction,
                    double const distance) const
  {
    Point new_point = init_point + direction * distance;
    return {new_point, direction};
  }

  /**\brief Sample a random point in the Cell.
   * \tparam Any class that supports geom_t random(), returning a URD on [0,1)
   * \param: c: the cell in which to sample a point
   */
  template <typename RNG_T>
  Vector sample_position(RNG_T & r, Cell const & c) const;

  /**\brief Sample a random direction isotropically.
   * \tparam Any class that supports geom_t random(), returning a URD on [0,1)
   */
  template <typename RNG_T>
  static Vector sample_direction_isotropic(RNG_T & r);

  /**\brief Constructor
   * \param nx: number cells in x-dimension
   * \param ny: number cells in y-direction
   * \param nz: number cells in z-direction
   * \note Default cell size is 1 unit in each dimension
   */
  Cartesian_Mesh(index_t nx, index_t ny, index_t nz)
      : nx_(nx),
        ny_(ny),
        nz_(nz),
        dx_(1.0),
        dy_(1.0),
        dz_(1.0),
        xmin_{0.0},
        ymin_{0.0},
        zmin_{0.0},
        n_faces{num_yz_faces(nx_, ny_, nz_), num_xz_faces(nx_, ny_, nz_),
                num_xy_faces(nx_, ny_, nz_)}
  {
    this->make_cartesian_mesh();
  }

  /**\brief Constructor
   * \param nx: number cells in x-dimension
   * \param ny: number cells in y-direction
   * \param nz: number cells in z-direction
   * \param dx: extent of cells in x-dimension
   * \param dy: extent of cells in y-direction
   * \param dz: extent of cells in z-direction
   * \note Default cell size is 1 unit in each dimension
   */
  Cartesian_Mesh(index_t nx,
                 index_t ny,
                 index_t nz,
                 geom_t dx,
                 geom_t dy,
                 geom_t dz)
      : nx_(nx),
        ny_(ny),
        nz_(nz),
        dx_(dx),
        dy_(dy),
        dz_(dz),
        xmin_{0.0},
        ymin_{0.0},
        zmin_{0.0},
        n_faces{num_yz_faces(nx_, ny_, nz_), num_xz_faces(nx_, ny_, nz_),
                num_xy_faces(nx_, ny_, nz_)}
  {
    this->make_cartesian_mesh();
  }

  /**\brief Constructor
   * \param nx: number cells in x-dimension
   * \param ny: number cells in y-direction
   * \param nz: number cells in z-direction
   * \param dx: extent of cells in x-dimension
   * \param dy: extent of cells in y-direction
   * \param dz: extent of cells in z-direction
   * \param xmin: minimum of mesh in x-dimension
   * \param ymin: minimum of mesh in y-direction
   * \param zmin: minimum of mesh in z-direction
   * \note Default cell size is 1 unit in each dimension
   */
  Cartesian_Mesh(index_t nx,
                 index_t ny,
                 index_t nz,
                 geom_t dx,
                 geom_t dy,
                 geom_t dz,
                 geom_t xmin,
                 geom_t ymin,
                 geom_t zmin)
      : nx_(nx),
        ny_(ny),
        nz_(nz),
        dx_(dx),
        dy_(dy),
        dz_(dz),
        xmin_{xmin},
        ymin_{ymin},
        zmin_{zmin},
        n_faces{num_yz_faces(nx_, ny_, nz_), num_xz_faces(nx_, ny_, nz_),
                num_xy_faces(nx_, ny_, nz_)}
  {
    this->make_cartesian_mesh();
  }

  /**\brief Number of cells in this Cartesian mesh */
  index_t num_cells() const;

  /**\brief Convert a tuple of Cartesian indices to a Cell */
  Cell make_cell(index_t ix, index_t iy, index_t iz) const;

  virtual ~Cartesian_Mesh() {}

  /**\brief Generate a Cartesian_Face parallel to indicated plane for indicated
   * cell. \param ix, iy, iz: Cartesian indices of cell \param nx, ny: Number of
   * cells in X- and Y-dimension \param p: plane parallel to desired face.
   *
   * Note: These are the least faces for the indicated cell (lower, front, or
   * left-most). To get the upper face index, go up a cell. See
   * Cartesian_Mesh::make_faces for an example.
   *
   * The face index is the cell's Cartesian index, mapped into an array
   * index. The array index is the offset for the set of faces parallel to
   * that plane, plus whatever the linear cell index would be. Thus all faces
   * are indexed uniquely in the same
   * set. Think of this as returning the floor of the cell.
   */
  Cartesian_Face cell_to_face(index_t ix,
                              index_t iy,
                              index_t iz,
                              Plane const p) const;

  /**\brief Compute a Cartesian_Face index from a cell and Face_Name. The
   * numbering of faces in this function is consistent with the numbering in
   * cell_to_face(ix, iy, iz, plane); see description of that function.
   */
  Cartesian_Face cell_to_face(Cell const & c, Face_Name const fn) const;

  /**\brief Get a cell index corresponding to this face.
   *
   * Note: the cell index will be such that this face is the lower bound of
   * the cell. Indexed that way, the Cell index returned may not be a valid
   * cell in the mesh. But this is still useful. I think. It's necessary
   * because there is not a one-to-one mapping of faces to cells.
   */
  Cell face_to_cell(Cartesian_Face const & f) const;

  /**\brief Get the plane for a face (YZ, XZ, or XY)
   *
   * The plane is the plane with the greatest offset less than the index */
  Plane face_to_plane(index_t face_index) const;

  /**\brief Get the plane offset for a face
   *
   * Note: the plane offset is the largest offset less than the index. */
  index_t face_to_plane_offset(index_t const face_idx) const;

  /**\brief Compute the Face_Name for the face with respect to the cell. */
  Face_Name face_to_face_name(Cell const & c, face_handle_t const & f) const;

  /**\brief Get number cells in X-dimension */
  index_t get_num_x() const { return nx_; }

  /**\brief Get number cells in Y-dimension */
  index_t get_num_y() const { return ny_; }

  /**\brief Get number cells in Z-dimension */
  index_t get_num_z() const { return nz_; }

  /**\brief Iterate over one of the sets of boundary faces */
  struct boundary_face_iterator;

  boundary_face_iterator boundary_faces_begin(Face_Name const f) const;

  boundary_face_iterator boundary_faces_end(Face_Name const f) const;

  /**\brief Identify all the face indices on a given Face_Name */
  std::vector<face_handle_t> boundary_faces(Face_Name const f) const;

  void add_cell(Cell & c, face_list_t & faces, vertex_list_t & vtcs)
  {
    cells_.insert(c);
    cell_to_vertex_[c] = vtcs;
    cell_to_face_[c] = faces;
  }

  /**\brief get the vertices connected to a cell*/
  vertex_list_t const & get_vertices(Cell const & c)
  {
    // Require(valid_cell(c));
    auto it = cell_to_vertex_.find(c);
    return (*it).second;
  }  // get_vertices

  /**\brief Get the faces that bound a cell */
  face_list_t const & get_faces(Cell const & c) const
  {
    // Require(valid_cell(c));
    auto it = cell_to_face_.find(c);
    return (*it).second;
  }  // get_faces

  cell_list_t const & cells() const { return cells_; }

private:
  // state
  index_t nx_, ny_, nz_;

  geom_t dx_, dy_, dz_;

  geom_t xmin_, ymin_, zmin_;

  index_t n_faces[Plane::NUM_PLANES];

  // implementation methods
  Vertex cell_to_vertex(index_t ix, index_t iy, index_t iz) const;

  /** \brief Construct the set of faces adjecent to the cell with Cartesian
   * indices {ix, iy, iz} */
  std::vector<Cartesian_Face> make_faces(index_t ix,
                                         index_t iy,
                                         index_t iz) const;

  /** \brief Construct the set of vertices adjecent to the cell with Cartesian
   * indices {ix, iy, iz} */
  vertex_list_t make_vertices(index_t ix, index_t iy, index_t iz) const;

  bool valid_cell(Cell const & c) { return cells_.find(c) != cells_.end(); }

  cell_list_t cells_;

  cell_to_vertex_map cell_to_vertex_;
  cell_to_face_map cell_to_face_;

  /**\brief For each cell, construct its adjacent faces and vertices; add cell
   * to Mesh.
   */
  void make_cartesian_mesh();

public:
};  // class Cartesian_Mesh

}  // namespace nut

#define IM_ALLOWED_TO_INCLUDE_CARTESIAN_MESH_I_H
#include "Cartesian_Mesh.i.h"
#undef IM_ALLOWED_TO_INCLUDE_CARTESIAN_MESH_I_H

// End of file
