#ifndef SYMOP_HH
#define SYMOP_HH

#include <cmath>
#include <iostream>

#include "casm/symmetry/SymOpRepresentation.hh"

namespace CASM {

class Log;

/** \defgroup Symmetry
 *  \brief Relates to symmetry operations and groups
 */

/** \defgroup SymOp
 *  \ingroup Symmetry
 *  \brief Relates to symmetry operations
 *  @{
 */

class MasterSymGroup;

///\brief SymOp is the Coordinate representation of a symmetry operation
/// it keeps fraction (FRAC) and Cartesian (CART) information about how
/// a symetry operation transforms 3D spatial coordinates
class SymOp : public SymOpRepresentation {
 public:
  typedef Eigen::Matrix3d matrix_type;
  typedef Eigen::Vector3d vector_type;

  /// static method to create operation that describes pure translation
  static SymOp translation(const Eigen::Ref<const vector_type> &_tau) {
    return SymOp(matrix_type::Identity(), _tau, false, TOL, 0, nullptr);
  }

  /// static method to create time_reversal operation
  static SymOp time_reversal_op() {
    return SymOp(matrix_type::Identity(), vector_type::Zero(), true, TOL, -1,
                 nullptr);
  }

  /// static method to create point operation (no translation)
  static SymOp point_op(const Eigen::Ref<const matrix_type> &_mat,
                        double _map_error = TOL) {
    return SymOp(_mat, vector_type::Zero(), false, _map_error, -1, nullptr);
  }

  SymOp() : SymOp(matrix_type::Identity(), vector_type::Zero(), false, TOL) {}

  /// Create new SymOp from Matrix3 and tau translation
  /// by default, assume no translation
  SymOp(const Eigen::Ref<const matrix_type> &_mat,
        const Eigen::Ref<const vector_type> &_tau, bool _time_reversal,
        double _map_error)
      : SymOp(_mat, _tau, _time_reversal, _map_error, -1, nullptr) {}

  /// Const access of entire cartesian symmetry matrix
  inline const matrix_type &matrix() const { return m_mat; }

  /// Const access of the cartesian translation vector, 'tau'
  inline const vector_type &tau() const { return m_tau; }

  /// Const access of the time-reversal flag (true if operation reverses time)
  inline bool time_reversal() const { return m_time_reversal; }

  /// Const access of the sym op's cartesian shift from its MasterSymGroup
  inline const vector_type &integral_tau() const {
    if (m_valid_integral_tau) {
      return m_integral_tau;
    }

    throw std::runtime_error("Error in SymOp::integral_tau(), not valid");
  }

  ///\brief returns true if matrix part of operation is identity
  inline bool is_identity() const {
    return matrix().isIdentity(TOL) && !time_reversal();
  }

  /// set master_group and op_index for the SymOp
  /// Performs simple checks to ensure that (*this) is compatible with
  /// new_group[new_index]
  void set_index(const MasterSymGroup &new_group, Index new_index);

  /// Allows access to the map_error.
  const double &map_error() const;
  void set_map_error(const double &value);

  // Routines that relate symmetry options:
  /// get a new SymOp that is equivalent to subsequent application of both
  /// SymOps
  SymOp operator*(const SymOp &RHS) const;

  /// Cartesian translation of the origin of the symmetry operation by
  /// Coordinate RHS
  SymOp &operator+=(const Eigen::Ref<const vector_type> &RHS);
  SymOp &operator-=(const Eigen::Ref<const vector_type> &RHS);

  /// get the inverse of this SymOp
  SymOp inverse() const;

  /// Get copy of the SymOp without translation
  SymOp no_trans() const;

  /// Check equality of SymOps, (matrix and translation). Does not necessarily
  /// return true for translational equivalents
  bool operator==(const SymOp &RHS) const;

  /// Performs conjugation of this SymOp with SymOp 'op'
  /// In other words, transforms coordinates of this SymOp by
  /// transformation specified by 'op':  op.inverse()*(*this)*op
  SymOp &apply_sym(const SymOp &op);

  /// calculate and return character of this SymOp (neglecting tau)
  double character() const override { return matrix().trace(); }

  /// Print this SymOP as header, followed by symmetry_mat and tau_vec (aligned
  /// vertically) \param c2fmat is a transformation matrix that can be passed to
  /// convert from Cartesian to fractional coordinates
  void print(std::ostream &stream, const Eigen::Ref<const matrix_type> &c2fmat =
                                       matrix_type::Identity()) const;

  /// Prints abridged description of SymOp, including its type, angle, and
  /// eigenvector does NOT print out the entire symmetry operation matrix
  /// \param c2fmat is a transformation matrix that can be passed to convert
  /// from Cartesian to fractional coordinates
  void print_short(std::ostream &stream,
                   const Eigen::Ref<const matrix_type> &c2fmat =
                       matrix_type::Identity()) const;

  /// Return pointer to a new copy of this SymOp
  SymOpRepresentation *copy() const override { return new SymOp(*this); }

  SymOp unregistered_copy() const {
    return SymOp(matrix(), tau(), time_reversal(), map_error());
  }

 private:
  /// Create new SymOp from Matrix3 and tau translation
  /// by default, assume no translation
  SymOp(const Eigen::Ref<const matrix_type> &_mat,
        const Eigen::Ref<const vector_type> &_tau, bool _time_reversal,
        double _map_error, Index _op_index, MasterSymGroup const *_master_ptr)
      : SymOpRepresentation(_master_ptr, SymGroupRepID(), _op_index),
        m_mat(_mat),
        m_tau(_tau),
        m_time_reversal(_time_reversal),
        m_map_error(_map_error) {
    _set_integral_tau();
  }

  /*** Inherited from SymOpRepresentation ***/
  /*
  ///Index into MasterSymGroup that specifies the operation
  Index m_op_index;

  SymGroupRepID m_rep_ID;

  /// Pointer to the MasterSymGroup where prototype of this SymOp lives
  MasterSymGroup const *m_master_group;
  */

  /// \brief Set the difference between the translation of this compared to
  // its MasterSymGroup correspondent
  void _set_integral_tau() override;

  /// matrix representation of symettry operation in
  /// Cartesian (symmetry_mat[CART])
  /// Fractional (symmetry_mat[FRAC])
  matrix_type m_mat;

  /// translation vector that is applied to
  /// a point after matrix transformation
  vector_type m_tau;

  /// time-reversal bit
  /// true if operation includes time-reversal,
  /// false otherwise
  bool m_time_reversal;

  /// translation vector that maps this op's correspondent in its Master Group
  /// to this
  mutable vector_type m_integral_tau;

  bool m_valid_integral_tau;

  /// This stores the mapping error associated with this SymOp, which will
  /// depend on the tolerances
  ///  you choose when you attempt to generate symmetry groups.
  double m_map_error;
};

/// Accessor to allow conversion to xtal::SymOpType. Returns Cartesian
/// transformation matrix.
const SymOp::matrix_type &get_matrix(const SymOp &op);
/// Accessor to allow conversion to xtal::SymOpType. Returns translation vector
/// (tau).
const SymOp::vector_type &get_translation(const SymOp &op);
/// Accessor to allow conversion to xtal::SymOpType. Returns whether the
/// symmetry operation is time reversal active.
bool get_time_reversal(const SymOp &op);

/// \brief Print formatted SymOp matrix and tau
void print_matrix_tau_col(Log &log, const SymOp &op, Index prec);

namespace adapter {

template <typename ToType, typename FromType>
struct Adapter;

/// Convert symmetry types to CASM::SymOp
///
/// Works for any symmetry type that has the get_matrix, get_translation, and
/// get_time_reversal accessors defined.
///
/// If there is a different symmetry type you would like to adapt for
/// CASM::SymOp, simply declare and define accessors get_matrix,
/// get_translation, and get_time_reversal
///
template <typename FromType>
struct Adapter<SymOp, FromType> {
  SymOp operator()(const FromType &adaptable) {
    // TODO: How to get mapping error tol?
    return SymOp(get_matrix(adaptable), get_translation(adaptable),
                 get_time_reversal(adaptable), CASM::TOL);
  }
};

/// Convert CASM::SymOp -> CASM::SymOp
///
/// Note:
/// - This just forwards a reference
template <>
struct Adapter<SymOp, SymOp> {
  SymOp const &operator()(SymOp const &adaptable) const { return adaptable; }

  SymOp &operator()(SymOp &adaptable) const { return adaptable; }
};
}  // namespace adapter

/** @}*/

}  // namespace CASM
#endif
