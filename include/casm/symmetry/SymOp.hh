#ifndef SYMOP_HH
#define SYMOP_HH

#include <iostream>
#include <cmath>

#include "casm/symmetry/SymOpRepresentation.hh"


namespace CASM {
  class MasterSymGroup;

  ///\brief SymOp is the Coordinate representation of a symmetry operation
  /// it keeps fraction (FRAC) and Cartesian (CART) information about how
  /// a symetry operation transforms 3D spatial coordinates
  class SymOp : public SymOpRepresentation {
  public:
    struct SymInfo {
      symmetry_type op_type;
      Eigen::Vector3d axis;
      double angle;
      Eigen::Vector3d screw_glide_shift;
      Eigen::Vector3d location;
    };

    typedef Eigen::Matrix3d matrix_type;
    typedef Eigen::Vector3d vector_type;

    /// static method to create operation that describes pure translation
    static SymOp translation(const Eigen::Ref<const vector_type> &_tau) {
      return SymOp(matrix_type::Identity(), _tau);
    }

    ///Create new SymOp from Matrix3 and tau translation
    /// by default, assume no translation
    SymOp(const Eigen::Ref<const matrix_type> &_mat = matrix_type::Identity(),
          const Eigen::Ref<const vector_type> &_tau = vector_type::Zero(),
          double _map_error = TOL) :
      m_mat(_mat),
      m_tau(_tau),
      m_map_error(_map_error) {
    }

    SymOp(const Eigen::Ref<const matrix_type> &_mat,
          double _map_error) :
      SymOp(_mat, vector_type::Zero(), _map_error) {

    }

    /// Const access of entire cartesian symmetry matrix
    inline
    const matrix_type &matrix() const {
      return m_mat;
    }

    /// Const access of the cartesian translation vector, 'tau'
    inline
    const vector_type &tau() const {
      return m_tau;
    }

    ///\brief returns true if matrix part of operation is identity
    inline
    bool is_identity() const {
      return matrix().isIdentity(TOL);
    }

    /// set master_group and op_index for the SymOp
    /// Performs simple checks to ensure that (*this) is compatible with new_group[new_index]
    void set_index(const MasterSymGroup &new_group, Index new_index);

    /// Allows access to the map_error.
    const double &map_error() const;
    void set_map_error(const double &value);

    //Routines that relate symmetry options:
    ///get a new SymOp that is equivalent to subsequent application of both SymOps
    SymOp operator*(const SymOp &RHS) const;

    /// Cartesian translation of the origin of the symmetry operation by Coordinate RHS
    SymOp &operator+=(const Eigen::Ref<const vector_type> &RHS);
    SymOp &operator-=(const Eigen::Ref<const vector_type> &RHS);

    /// get the inverse of this SymOp
    SymOp inverse() const;

    /// Get copy of the SymOp without translation
    SymOp no_trans() const;

    /// Check equality of SymOps, (matrix and translation). Does not necessarily return true for translational equivalents
    bool operator==(const SymOp &RHS) const;

    /// Performs conjugation of this SymOp with SymOp 'op'
    /// In other words, transforms coordinates of this SymOp by
    /// transformation specified by 'op':  op.inverse()*(*this)*op
    SymOp &apply_sym(const SymOp &op);

    /// Use SymOp::matrix() and SymOp::tau() to populate a SymInfo object with type of symmetry,
    /// axis (i.e., rotation axis or mirror normal), screw/glide shift, and location
    /// because SymOp doesn't know any lattice info, it can't tell between rotation and screw or mirror and glide
    /// it assumes that a rotation is a screw and that a mirror is a glide when populating SymInfo::screw_glide_shift
    SymInfo info() const;

    /// calculate and return character of this SymOp (neglecting tau)
    double character() const override {
      return matrix().trace();
    }

    /// Print this SymOP as header, followed by symmetry_mat and tau_vec (aligned vertically)
    /// \param c2fmat is a transformation matrix that can be passed to convert from Cartesian to fractional coordinates
    void print(std::ostream &stream, const Eigen::Ref<const Eigen::Matrix3d> &c2fmat = matrix_type::Identity()) const;

    ///Prints abridged description of SymOp, including its type, angle, and eigenvector
    ///does NOT print out the entire symmetry operation matrix
    /// \param c2fmat is a transformation matrix that can be passed to convert from Cartesian to fractional coordinates
    void print_short(std::ostream &stream, const Eigen::Ref<const Eigen::Matrix3d> &c2fmat = matrix_type::Identity()) const;

    /// Return pointer to a new copy of this SymOp
    SymOpRepresentation *copy() const override {
      return new SymOp(*this);
    }

    SymOp unregistered_copy() const {
      return SymOp(matrix(), tau(), map_error());
    }

    jsonParser &to_json(jsonParser &json) const override;

    void from_json(const jsonParser &json) override;

  private:

    /*** Inherited from SymOpRepresentation ***/
    /*
    ///Index into MasterSymGroup that specifies the operation
    Index m_op_index;

    SymGroupRepID m_rep_ID;

    /// Pointer to the MasterSymGroup where prototype of this SymOp lives
    MasterSymGroup const *m_master_group;
    */

    ///matrix representation of symettry operation in
    ///Cartesian (symmetry_mat[CART])
    ///Fractional (symmetry_mat[FRAC])
    matrix_type m_mat;

    ///translation vector that is applied to
    ///a point after matrix transformation
    vector_type m_tau;

    /// This stores the mapping error associated with this SymOp, which will depend on the tolerances
    ///  you choose when you attempt to generate symmetry groups.
    double m_map_error;

  };


  jsonParser &to_json(const SymOp &sym, jsonParser &json);
  void from_json(SymOp &sym, const jsonParser &json);

  //std::istream& operator>> (std::istream &stream, SymOp &op){op.read(stream); return stream;};
  //std::ostream& operator<< (std::ostream &stream, SymOp &op){op.print(stream); return stream;};

}
#endif
