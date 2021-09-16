#ifndef CASM_config_SupercellSymInfo
#define CASM_config_SupercellSymInfo

#include "casm/symmetry2/StructureSymInfo.hh"

namespace CASM {
namespace config {

/// \brief Data structure describing application of symmetry in a supercell
struct SupercellSymInfo {
  SupercellSymInfo(std::shared_ptr<Prim const> const &prim,
                   Superlattice const &_superlattice);

  /// \brief The subgroup of the prim factor group that leaves
  /// the supercell lattice vectors invariant
  std::shared_ptr<Group const> factor_group;

  /// \brief Describes how sites permute due to translations within the
  /// supercell
  ///
  /// The number of translations is equal the supercell volume (as an integer
  /// multiple of the prim unit cell)
  std::vector<Permutation> translation_permutations;

  /// \brief Describes how sites permute due to supercell factor group
  /// operations.
  ///
  /// There is one element for each element in the supercell factor group.
  std::vector<Permutation> factor_group_permutations;
};

}  // namespace config
}  // namespace CASM

#endif
