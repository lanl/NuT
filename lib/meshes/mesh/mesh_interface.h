/**
 * Interface definition for data access classes
 */

#pragma once

#include <utility>

namespace nut {

template <class MeshAdaptor>
class MeshInterface : public MeshAdaptor {
public:
  // type definitions
  /////////////////////////////////////////////////////////
  /// the actual mesh object
  using mesh_t = typename MeshAdaptor::mesh_t;
  /// index type
  using index_t = typename MeshAdaptor::index_t;
  /// geometron type
  using geometron_t = typename MeshAdaptor::geometron_t;
  /// data type to describe an intersection
  using intersect_t = typename MeshAdaptor::intersect_t;
  /// opaque data used for particle communication
  using ticket_t = typename MeshAdaptor::ticket_t;

  /// distance type
  using distance_t = typename MeshAdaptor::distance_t;
  /// area type
  using area_t = typename MeshAdaptor::area_t;
  /// volume type
  using volume_t = typename MeshAdaptor::volume_t;

  /// vector type
  using vector_t = typename MeshAdaptor::vector_t;
  /// 4 vector type
  using vector4_t = typename MeshAdaptor::vector4_t;
  /// point type
  using point_t = typename MeshAdaptor::point_t;

  /// handle to mesh vertex type
  using vertex_handle_t = typename MeshAdaptor::vertex_handle_t;
  /// handle to mesh face type
  using face_handle_t = typename MeshAdaptor::face_handle_t;
  /// cell handle type
  using cell_handle_t = typename MeshAdaptor::cell_handle_t;

  /// index for vertex type
  using vertex_index_t = typename MeshAdaptor::vertex_index_t;
  /// index for face type
  using face_index_t = typename MeshAdaptor::face_index_t;
  /// index for cell type
  using cell_index_t = typename MeshAdaptor::cell_index_t;
  /// index for boundary types
  using boundary_index_t = typename MeshAdaptor::boundary_index_t;

  /// index for vertex type
  using global_vertex_index_t = typename MeshAdaptor::global_vertex_index_t;
  /// index for face type
  using global_face_index_t = typename MeshAdaptor::global_face_index_t;
  /// index for cell type
  using global_cell_index_t = typename MeshAdaptor::global_cell_index_t;
  /// index for boundary types
  using global_boundary_index_t = typename MeshAdaptor::global_boundary_index_t;

  // constructor
  /////////////////////////////////////////////////////////
  /// default constructor
  MeshInterface(mesh_t const & _mesh) : MeshAdaptor(_mesh) {}
  /// copy constructor
  MeshInterface(MeshInterface const & _rhs) = default;
  /// move constructor
  MeshInterface(MeshInterface && _rhs) noexcept = default;
  /// destructor
  ~MeshInterface() = default;

  /// copy assignment operator
  MeshInterface & operator=(MeshInterface const & _rhs) = default;
  /// move assignment operator
  MeshInterface & operator=(MeshInterface && _rhs) noexcept = default;

  // connectivity
  /////////////////////////////////////////////////////////

  /// vertices connected to cell
  decltype(auto) cell_vertices(cell_handle_t _cell) const
  {
    return MeshAdaptor::cell_vertices(_cell);
  }

  /// vertices connected to cell
  decltype(auto) cell_faces(cell_handle_t _cell) const
  {
    return MeshAdaptor::cell_faces(_cell);
  }

  /////////////////////////////////////////////////////////
  /// get the number of vertices in local mesh
  size_t num_vertices() const { return MeshAdaptor::num_vertices(); }

  /// get number of faces in local mesh
  size_t num_faces() const { return MeshAdaptor::num_faces(); }

  /// number of cells in the local mesh
  size_t num_cells() const { return MeshAdaptor::num_cells(); }

  /// number of boundaries in the local mesh
  size_t num_boundaries() const { return MeshAdaptor::num_boundaries(); }

  /// get the number of vertices in global mesh
  size_t num_global_vertices() const { return MeshAdaptor::num_vertices(); }

  /// get number of faces in global mesh
  size_t num_global_faces() const { return MeshAdaptor::num_faces(); }

  /// number of cells in the global mesh
  size_t num_global_cells() const { return MeshAdaptor::num_cells(); }

  /// number of boundaries in the global mesh
  size_t num_global_boundaries() const { return MeshAdaptor::num_boundaries(); }

  // index function defintion
  /////////////////////////////////////////////////////////

  /// get index for a cell handle
  vertex_index_t vertex_idx(vertex_handle_t const & _vertex) const
  {
    return MeshAdaptor::vertex_idx(_vertex);
  }

  /// get index for a cell handle
  face_index_t face_idx(face_handle_t const & _face) const
  {
    return MeshAdaptor::face_idx(_face);
  }

  /// get index for a cell handle
  cell_index_t cell_idx(cell_handle_t const & _cell) const
  {
    return MeshAdaptor::cell_idx(_cell);
  }

  /// get global index for a vertex handle
  global_vertex_index_t global_vertex_idx(const vertex_handle_t & _vertex) const
  {
    return MeshAdaptor::global_vertex_idx(_vertex);
  }

  /// get global index for a face handle
  global_face_index_t global_face_idx(const face_handle_t & _face) const
  {
    return MeshAdaptor::global_face_idx(_face);
  }

  /// get global index for a cell handle
  global_cell_index_t global_cell_idx(const cell_handle_t & _cell) const
  {
    return MeshAdaptor::global_cell_idx(_cell);
  }

  // vertex function defintion
  /////////////////////////////////////////////////////////
  /// Get coordinated of vertex; point_t or point_t &
  decltype(auto) coordinates(vertex_handle_t const & _vertex) const
  {
    return MeshAdaptor::get_coordinates(_vertex);
  }

  /// get a handle that can be interpreted as the null vertex
  vertex_handle_t null_vertex() const { return MeshAdaptor::null_vertex(); }

  // face function defintion
  /////////////////////////////////////////////////////////

  /// the area of said face; area_t or area_t &
  decltype(auto) area(face_handle_t const & _face) const
  {
    return MeshAdaptor::area(_face);
  }

  /// the center of said face; point_t or point_t &
  decltype(auto) face_center(face_handle_t const & _face) const
  {
    return MeshAdaptor::face_center(_face);
  }

  /// shortest distance from point to said face; distance_t or distance_t &
  decltype(auto) shortest_distance(face_handle_t const & face,
                                   const point_t & p) const
  {
    return MeshAdaptor::shortest_distance(face, p);
  }

  /// Obtain a handle to a null face
  static face_handle_t null_face() { return MeshAdaptor::null_face(); }

  // cell function defintion
  /////////////////////////////////////////////////////////

  /// volume of said cell; volume_t or volume_t &
  decltype(auto) volume(cell_handle_t const & _cell) const
  {
    return MeshAdaptor::volume(_cell);
  }

  /// center of said cell; point_t or point_t &
  decltype(auto) cell_center(cell_handle_t const & _cell) const
  {
    return MeshAdaptor::cell_center(_cell);
  }

  /** \brief Whether the point is in the cell
   * \remark True if geometron's point is within a cell; true if on the surface,
   * as long as internally consistent with geometron state. In particular, the
   * geometron may define being in one cell and not the other cells that share
   * that point.
   */
  bool in_cell(geometron_t const & _geometron,
               cell_handle_t const & _cell) const
  {
    return MeshAdaptor::in_cell(_geometron, _cell);
  }

  /** \brief Find the cell that encompasses the geometron's point.
   * \remark Undefined if point lies on a cell boundary.
   */
  cell_handle_t find_cell(point_t const & _point) const
  {
    return MeshAdaptor::find_cell(_point);
  }

  /// Obtain a handle to a null cell
  cell_handle_t null_cell() const { return MeshAdaptor::null_cell(); }

  /// Is said cell null?
  bool is_null(cell_handle_t const & cell) const
  {
    return MeshAdaptor::is_null(cell);
  }

  // boundary function defintion
  /////////////////////////////////////////////////////////

  /// Is said face a boundary of the mesh?
  bool is_boundary(face_handle_t const & _face) const
  {
    return MeshAdaptor::is_boundary(_face);
  }

  /** \brier Obtain index of boundary set that face is associated with.
   * \remark If face is not on a boundary, return null bdy_index.
   */
  boundary_index_t boundary(face_handle_t const & _face) const
  {
    return MeshAdaptor::boundary(_face);
  }

  /// Obtain an index of a null boundary
  boundary_index_t null_bdy_index() const
  {
    return MeshAdaptor::null_bdy_index();
  }

  /// Is said boundary index null?
  bool is_null(boundary_index_t const & bdy) const
  {
    return MeshAdaptor::is_null(bdy);
  }

  // internal boundary function defintion
  /////////////////////////////////////////////////////////

  /// Is said face a subdomain boundary?
  bool is_subdomain_boundary(face_handle_t const & _face) const
  {
    return MeshAdaptor::is_subdomain_boundary(_face);
  }

  // Geometric Operations
  /////////////////////////////////////////////////////////

  /**\brief Reflect at a surface
   *
   * \param direction: vector to be reflected
   * \param face: face from which to reflect
   * \param point: point at which to reflect
   * \remark Undefined if point is not interpretable as lying on face.
   */
  vector_t reflect(vector_t const & _direction,
                   face_handle_t const & _face,
                   point_t const & _point) const
  {
    return MeshAdaptor::reflect(_direction, _face);
  }

  /**\brief Compute the point resulting from moving translation away from
   * starting point.
   *
   * \param point: starting point
   * \param translation: displacement vector
   * \remark This is a purely geometric operation: for example, the resultant
   * point may not lie in the same cell.
   */
  static point_t displace(point_t const & _point, vector_t const & _translation)
  {
    return MeshAdaptor::displace(_point, _translation);
  }

  /** \brief Calculate the next intersection with a mesh face.
   *
   * \param Geometron: current geometric state
   */
  intersect_t intersect(geometron_t const & _geometron) const
  {
    return MeshAdaptor::intersect(_geometron);
  }

  /** \brief Advance geometron by distance, no face crossings performed
   *
   * \remark As with \c displace, care must be taken that this operation does
   * not produce an inconsistent state. For example, this method does not update
   * the \c geometron's cell or surface information. This method is typically
   * used combination with \c find_cell or \c cross_face to implement a
   * complete, consistent update of the geometron.
   */
  void advance(geometron_t & _geometron, distance_t const & _distance) const
  {
    MeshAdaptor::advance(_geometron, _distance);
  }

  /** \briefReflect the geometron's direction about the intersection point.
   *
   * \remark Does not update the geometron's position; this method is typically
   * used in conjunction with \c advance to move the geometron to the face.
   * Undefined if geometron has not been moved to the intersection point. */
  void reflect(geometron_t & _geometron, intersect_t const & _intersect) const
  {
    MeshAdaptor::reflect(_geometron, _intersect);
  }

  /** \brief Move geometron across a mesh face.
   *
   * \remark Typically used after \c advance is called to move geometron to the
   * face. Undefined is geometron has not been moved to the intersection point.
   */
  void cross_face(geometron_t & _geometron,
                  intersect_t const & _intersect) const
  {
    return MeshAdaptor::cross_face(_geometron, _intersect);
  }

  /** \brief Cross a face that is also a subdomain boundary.
   *
   * Moves geometron across the face given by intersect, updating any
   * internal cell and surface information. Returns a communication token that
   * is both carryies all mesh information needed to communicate the particle as
   * well as being opaque to the user. For example, in MPI, the \c ticket_t
   * might hold the destination rank and communicator. Undefined if geometron is
   * not on the face being crossed.
   */
  ticket_t cross_subdomain_boundary(geometron_t & _geometron,
                                    intersect_t const & _intersect) const
  {
    return MeshAdaptor::cross_subdomain_boundary(_geometron, _intersect);
  }

  // intersection access functions
  /////////////////////////////////////////////////////////
  /// Get point of intersection; point_t or point_t &
  decltype(auto) intersection_point(intersect_t const & _intersect) const
  {
    return MeshAdaptor::intesection_point(_intersect);
  }

  /// Get face intersected; face_handle_t or face_handle_t &
  decltype(auto) intersection_face(intersect_t const & _intersect) const
  {
    return MeshAdaptor::intersection_face(_intersect);
  }

  /// Get distance to intersection; distance_t or distance_t &
  decltype(auto) intersection_distance(intersect_t const & _intersect) const
  {
    return MeshAdaptor::intersection_distance(_intersect);
  }

  // Random sampling functions
  /////////////////////////////////////////////////////////
  /**\brief Sample a random point in the Cell.
   * \tparam Any class that supports geom_t random(), returning a URD on [0,1)
   * \param: c: the cell in which to sample a point
   */
  template <typename RNG_T>
  point_t sample_position(RNG_T & _random_number_generator,
                          cell_handle_t const & _cell) const
  {
    return MeshAdaptor::sample_position(_random_number_generator, _cell);
  }

  /**\brief Sample a random direction in the Cell.
   * \tparam Any class that supports geom_t random(), returning a URD on [0,1)
   */
  template <typename RNG_T>
  static vector_t sample_direction_isotropic(RNG_T & _random_number_generator)
  {
    return MeshAdaptor::sample_direction_isotropic(_random_number_generator);
  }

  // Covariant geometric functions
  /////////////////////////////////////////////////////////

  /**\brief Transform input to lab frame, given ordinary (3) velocity measured
   * in lab frame.
   *
   * \remark Velocity is typically the material velocity measured in lab frame.
   */
  vector4_t transform_to_lab_frame(vector_t const & _velocity,
                                   const vector4_t _input) const
  {
    return MeshAdaptor::transform_to_lab_frame(_velocity, _input);
  }

  /**\brief Transform input to comoving frame, given ordinary (3) velocity
   * measured in lab frame.
   *
   * \remark Velocity is typically the material velocity measured in lab frame.
   */
  vector4_t transform_to_comoving_frame(vector_t const & _velocity,
                                        const vector4_t _input)
  {
    return MeshAdaptor::transform_to_comoving_frame(_velocity, _input);
  }

};  // Mesh Interface

}  // end namespace nut
