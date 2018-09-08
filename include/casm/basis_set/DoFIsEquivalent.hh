#ifndef CASM_DoFIsEquivalent
#define CASM_DoFIsEquivalent

#include "casm/external/Eigen/Dense"
#include "casm/basis_set/DoFSet.hh"

namespace CASM {
  class DoFSet;
  class SymOp;


  /// \brief Class for checking equivalence of two DoFSet objects, with respect to symmetry transformations
  ///
  /// DoFSets dofA and dofB are considered equivalent if
  /// - dofA.type_name()==dofB.type_name() AND
  /// - dofA.size() == dofB.size() AND
  /// - dofA[i].var_name()==dofB[i].var_name() for all i<dofA.size() AND
  /// - the matrix equation dofB.basis() * U = dofA.basis() has a unique solution U
  ///
  /// \ingroup IsEquivalent
  ///
  class DoFIsEquivalent {
  public:

    DoFIsEquivalent(DoFSet const &_dof) : m_dof(_dof) {}

    /// returns true if m_dof and _other have matching labels, and m_dof.basis() = _other.basis()*U
    bool operator()(DoFSet const &_other) const;

    /// returns true and copy_apply(_op,m_dof.basis()) = m_dof.basis()*U
    bool operator()(SymOp const &_op) const;

    /// returns true and copy_apply(_op,m_dof.basis()) = _other.basis()*U
    bool operator()(SymOp const &_op, DoFSet const &_other) const;

    /// return transformation matrix U calculated during last successful comparison
    Eigen::MatrixXd const &U() const {
      return m_U;
    }

  private:

    /// returns true if m_dof and _other are same type, same size, and have same variable names
    bool _label_equiv(DoFSet const &_other) const;

    /// returns true if the matrix equation _other_basis * U = m_dof.basis() has a unique solution U
    bool _vector_equiv(Eigen::Ref<const Eigen::MatrixXd> const &_other_basis) const;

    DoFSet m_dof;

    mutable Eigen::MatrixXd m_U;
  };
}

#endif
