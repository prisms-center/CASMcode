#ifndef CASM_SymInfo
#define CASM_SymInfo

#include <iostream>
#include <string>
#include <iomanip>

#include "casm/external/Eigen/Dense"
#include "casm/crystallography/Coordinate.hh"
#include "casm/casm_io/EnumIO.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  /** \ingroup Symmetry
   *  @{
   */

  enum class symmetry_type {
    identity_op,
    mirror_op,
    glide_op,
    rotation_op,
    screw_op,
    inversion_op,
    rotoinversion_op,
    invalid_op
  };

  ENUM_TRAITS(symmetry_type)

  ENUM_IO(symmetry_type)

  /// \brief Simple struct to be used as return type for SymOp::info().
  struct SymInfo {

    SymInfo(const SymOp &op, const Lattice &lat);

    /// One of: identity_op, mirror_op, glide_op, rotation_op, screw_op,
    ///         inversion_op, rotoinversion_op, or invalid_op
    symmetry_type op_type;

    /// Rotation axis if operation S is rotation/screw operation
    /// If improper operation, rotation axis of inversion*S
    /// (implying that axis is normal vector for a mirror plane)
    /// normalized to length 1
    /// axis is zero if operation is identity or inversion
    Coordinate axis;

    /// Rotation angle, if operation S is rotation/screw operation
    /// If improper operation, rotation angle of inversion*S
    double angle;

    /// Component of tau parallel to 'axis' (for rotation)
    /// or perpendicular to 'axis', for mirror operation
    Coordinate screw_glide_shift;

    /// A Cartesian coordinate that is invariant to the operation (if one exists)
    Coordinate location;


  private:

    typedef SymOp::vector_type vector_type;
    typedef SymOp::matrix_type matrix_type;

    void _set(const vector_type &_axis,
              const vector_type &_screw_glide_shift,
              const vector_type &_location,
              const Lattice &lat);
  };


  /// \brief Print SymInfo to string
  std::string to_string(const SymInfo &info, COORD_TYPE mode);

  /// \brief Print SymInfo to string
  std::string description(const SymOp &op, const Lattice &lat, COORD_TYPE mode);

  /// \brief Add to existing JSON object
  void add_sym_info(const SymInfo &info, jsonParser &j);

  /** @} */
}

#endif
