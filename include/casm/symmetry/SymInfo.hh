#ifndef CASM_SymInfo
#define CASM_SymInfo

#include "casm/crystallography/Coordinate.hh"
#include "casm/external/Eigen/Dense"
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

/// \brief Simple struct to be used as return type for SymOp::info().
struct SymInfo {
  SymInfo(const SymOp &op, const xtal::Lattice &lat);

  /// One of: identity_op, mirror_op, glide_op, rotation_op, screw_op,
  ///         inversion_op, rotoinversion_op, or invalid_op
  symmetry_type op_type;

  /// Rotation axis if operation S is rotation/screw operation
  /// If improper operation, rotation axis of inversion*S
  /// (implying that axis is normal vector for a mirror plane)
  /// normalized to length 1
  /// axis is zero if operation is identity or inversion
  xtal::Coordinate axis;

  /// Rotation angle, if operation S is rotation/screw operation
  /// If improper operation, rotation angle of inversion*S
  double angle;

  /// Component of tau parallel to 'axis' (for rotation)
  /// or perpendicular to 'axis', for mirror operation
  xtal::Coordinate screw_glide_shift;

  /// A Cartesian coordinate that is invariant to the operation (if one exists)
  xtal::Coordinate location;

  /// If time reversal symmetry
  bool time_reversal;

 private:
  typedef SymOp::vector_type vector_type;
  typedef SymOp::matrix_type matrix_type;

  void _set(const vector_type &_axis, const vector_type &_screw_glide_shift,
            const vector_type &_location, const xtal::Lattice &lat);
};

/** @} */
}  // namespace CASM

#endif
