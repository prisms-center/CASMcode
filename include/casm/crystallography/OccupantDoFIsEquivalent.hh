#ifndef CASM_OccupantDoFIsEquivalent
#define CASM_OccupantDoFIsEquivalent

#include "casm/container/Permutation.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/external/Eigen/Dense"

namespace CASM {
namespace xtal {

class Molecule;

/// \brief Class for checking equivalence of two OccupantDoF objects, with
/// respect to symmetry transformations
///
/// OccupantDoFs dofA and dofB are considered equivalent if
/// - dofA and dofB have the same number of allowed occupants AND
/// - each allowed occupant of dofA is also an allowed occupant of dofB
///
/// \ingroup IsEquivalent
///
class OccupantDoFIsEquivalent {
 public:
  OccupantDoFIsEquivalent(std::vector<xtal::Molecule> const &_dof,
                          double tol = TOL)
      : m_dof(_dof), m_tol(tol), m_P(_dof.size()) {}

  /// returns true if m_dof and _other have matching labels, and m_dof =
  /// P.permute(_other)
  bool operator()(std::vector<xtal::Molecule> const &_other) const;

  /// returns true if copy_apply(_op,m_dof) = P.permute(m_dof)
  bool operator()(xtal::SymOp const &_op) const;

  /// returns true if copy_apply(_op,m_dof) =  P.permute(_other)
  bool operator()(xtal::SymOp const &_op,
                  std::vector<xtal::Molecule> const &_other) const;

  /// return transformation permutation P calculated during last successful
  /// comparison
  Permutation const &perm() const { return m_P; }

 private:
  std::vector<xtal::Molecule> m_dof;

  double m_tol;

  mutable Permutation m_P;
};

}  // namespace xtal
}  // namespace CASM

#endif
