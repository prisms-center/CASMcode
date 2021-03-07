#include "casm/symmetry/VectorSpaceSymReport.hh"

namespace CASM {

namespace SymRepTools_v2 {

/// Construct VectorSpaceSymReport
///
/// \param irrep_decomposition An IrrepDecomposition
/// \param calc_wedges If true, 'irreducible_wedge' is constructed. If false,
///     'irreducible_wedge' is empty.
VectorSpaceSymReport vector_space_sym_report(
    IrrepDecomposition const &irrep_decomposition, bool calc_wedges) {
  VectorSpaceSymReport result;
  result.irreps = irrep_decomposition.irreps;
  result.symmetry_adapted_subspace =
      irrep_decomposition.symmetry_adapted_subspace;
  for (Index element_index : irrep_decomposition.head_group) {
    auto matrix_rep = irrep_decomposition.fullspace_rep[element_index];
    result.symgroup_rep.push_back(matrix_rep);
  }
  result.axis_glossary =
      std::vector<std::string>(result.symmetry_adapted_subspace.rows(), "x");
  Index i = 0;
  for (std::string &x : result.axis_glossary) {
    x += std::to_string(++i);
  }
  if (calc_wedges) {
    result.irreducible_wedge = make_symrep_subwedges(irrep_decomposition);
  }
  return result;
}

}  // namespace SymRepTools_v2

/// Construct the VectorSpaceSymReport
///
/// \param _rep Matrix representation of head_group, this defines group action
///     on the underlying vector space
/// \param head_group group for which the sym report is to be generated
/// \param _subspace matrix such that _subspace.rows()==_rep.dim() and whose
///     columns specify a subspace of underlying vector space
/// \param calc_wedges if true, 'irreducible_wedge' of returned object is
///     initialized, if false, 'irreducible_wedge' is empty
SymRepTools_v2::VectorSpaceSymReport vector_space_sym_report_v2(
    SymGroupRep const &rep, SymGroup const &head_group,
    Eigen::MatrixXd const &subspace, bool calc_wedges) {
  bool allow_complex = true;
  SymRepTools_v2::IrrepDecomposition irrep_decomposition =
      make_irrep_decomposition(rep, head_group, subspace, allow_complex);
  return SymRepTools_v2::vector_space_sym_report(irrep_decomposition,
                                                 calc_wedges);
}

}  // namespace CASM
