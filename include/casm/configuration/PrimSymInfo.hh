#ifndef CASM_config_PrimSymInfo
#define CASM_config_PrimSymInfo

#include "casm/configuration/Group.hh"
#include "casm/configuration/UnitCellCoordRep.hh"
#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

/// \brief Data structure describing how prim DoF values transform under
/// application of symmetry
struct PrimSymInfo {
  PrimSymInfo(BasicStructure const &basicstructure);

  /// \brief Structure factor group
  std::shared_ptr<Group const> factor_group;

  /// \brief Factor group elements without translations
  std::shared_ptr<Group const> point_group;

  /// \brief Describes how sublattices permute under symmetry
  BasisPermutationSymGroupRep basis_permutation_symgroup_rep;

  /// \brief True if any permutation in occ_symgroup_rep is non-trivial
  bool has_aniso_occs;

  /// \brief Permutations describe occupant index transformation under symmetry
  ///
  /// Usage:
  /// - occ_symgroup_rep[master_group_index][sublattice_index] -> Permutation
  ///
  /// Note:
  /// - This describes cases such as discrete molecular orientations or
  /// occupants with spin where a symmetry operation may transform one discrete
  /// occupant into another *before* permutating among sites
  OccSymGroupRep occ_symgroup_rep;

  /// \brief Matrices describe local DoF value transformation under symmetry
  ///
  /// Usage:
  /// - local_dof_symgroup_rep[dof_type][master_group_index][sublattice_index]
  /// -> Eigen::MatrixXd
  /// - For each group element there is one matrix representation per sublattice
  std::map<DoFKey, LocalDoFSymGroupRep> local_dof_symgroup_rep;

  /// \brief Matrices describe local DoF value transformation under symmetry
  ///
  /// Usage:
  /// - local_dof_symgroup_rep[dof_type][master_group_index][sublattice_index]
  /// -> Eigen::MatrixXd
  /// - There is one matrix representation per sublattice
  std::map<DoFKey, GlobalDoFSymGroupRep> global_dof_symgroup_rep;
};

}  // namespace config
}  // namespace CASM

#endif
