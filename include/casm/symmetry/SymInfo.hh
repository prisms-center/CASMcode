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
  namespace xtal {
    class Lattice;
    class Coordinate;
  }
  using xtal::Lattice;
  using xtal::Coordinate;


  class SymGroup;

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

  ENUM_IO_DECL(symmetry_type)

  struct SymInfoOptions {
    SymInfoOptions(COORD_TYPE _coord_type = FRAC, double _tol = TOL, Index _prec = 7, bool _print_matrix_tau = false) :
      coord_type(_coord_type), tol(_tol), prec(_prec), print_matrix_tau(_print_matrix_tau) {}
    COORD_TYPE coord_type;
    double tol;
    Index prec;
    bool print_matrix_tau;
  };

  jsonParser &to_json(const SymInfoOptions &opt, jsonParser &json);

  /// \brief Read from JSON
  void from_json(SymInfoOptions &opt, const jsonParser &json);

  template<>
  struct jsonConstructor<SymInfoOptions> {
    static SymInfoOptions from_json(const jsonParser &json);
  };

  /// \brief Simple struct to be used as return type for SymOp::info().
  struct SymInfo {

    SymInfo(const SymOp &op, const Lattice &lat, SymInfoOptions opt = SymInfoOptions());

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

  /// \brief Print SymInfo
  void print_sym_info(Log &log, const SymInfo &info, SymInfoOptions opt = SymInfoOptions());

  /// \brief Print SymInfo to string
  std::string to_string(const SymInfo &info, SymInfoOptions opt = SymInfoOptions());

  /// \brief Print symmetry symbol to string
  std::string to_brief_unicode(const SymInfo &info, SymInfoOptions opt = SymInfoOptions());

  /// \brief Print SymInfo to string
  std::string description(const SymOp &op, const Lattice &lat, SymInfoOptions opt = SymInfoOptions());

  /// \brief Print SymGroup with matrix / tau
  void description(Log &log, const SymGroup &g, const Lattice &lat, SymInfoOptions opt = SymInfoOptions());

  /// \brief Print SymInfo to brief string
  std::string brief_description(const SymOp &op, const Lattice &lat, SymInfoOptions opt = SymInfoOptions());

  /// \brief Print SymGroup with brief string
  void brief_description(Log &log, const SymGroup &g, const Lattice &lat, SymInfoOptions opt = SymInfoOptions());

  /// \brief Add to existing JSON object
  void add_sym_info(const SymInfo &info, jsonParser &j);

  /** @} */
}

#endif
