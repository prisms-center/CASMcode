#ifndef CASM_ClexBasisInfo_impl
#define CASM_ClexBasisInfo_impl

#include "casm/clex/ClexBasis.hh"
#include "casm/clex/ClexBasisInfo.hh"

namespace CASM {

/// Make ClexBasisFunctionInfo for all ClexBasis functions
///
/// \param clex_basis The generated ClexBasis
/// \param begin, end The range of orbits used to generate clex_basis
///
/// \returns A vector of ClexBasisFunctionInfo of size equal to
/// `clex_basis.n_functions()`
///
template <typename ClusterOrbitIterator>
std::vector<ClexBasisFunctionInfo> make_clex_basis_function_info(
    ClexBasis const &clex_basis, ClusterOrbitIterator begin,
    ClusterOrbitIterator end) {
  std::vector<ClexBasisFunctionInfo> info;

  // collect basis function info
  Index orbit_index = 0;
  Index function_index = 0;

  for (auto orbit_it = begin; orbit_it != end; ++orbit_it) {
    std::vector<Index> invariant_group_indices;
    for (SymOp const &op : orbit_it->equivalence_map()[0]) {
      invariant_group_indices.push_back(op.index());
    }

    BasisSet const &bset_prototype = clex_basis.bset_orbit(orbit_index)[0];
    for (Index i = 0; i < bset_prototype.size(); ++i) {
      info.emplace_back(orbit_index, function_index, orbit_it->prototype(),
                        orbit_it->size(), invariant_group_indices);
      function_index++;
    }
    orbit_index++;
  }

  return info;
}

}  // namespace CASM

#endif
