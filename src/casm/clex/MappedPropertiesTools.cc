#include "casm/clex/MappedPropertiesTools.hh"

#include "casm/clex/MappedProperties.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {
MappedProperties copy_apply(PermuteIterator const &op,
                            MappedProperties const &props) {
  SymOp symop = op.sym_op();

  MappedProperties result;
  for (auto it = props.global.begin(); it != props.global.end(); ++it) {
    AnisoValTraits traits(AnisoValTraits::name_suffix(it->first));
    result.global[it->first] =
        traits.symop_to_matrix(symop.matrix(), symop.tau(),
                               symop.time_reversal()) *
        it->second;
  }

  Permutation tperm(op.combined_permute());
  for (auto it = props.site.begin(); it != props.site.end(); ++it) {
    AnisoValTraits traits(AnisoValTraits::name_suffix(it->first));
    Eigen::MatrixXd new_matrix =
        traits.symop_to_matrix(symop.matrix(), symop.tau(),
                               symop.time_reversal()) *
        it->second;
    auto it2 = result.site
                   .emplace(std::make_pair(
                       it->first,
                       Eigen::MatrixXd(new_matrix.rows(), new_matrix.cols())))
                   .first;
    for (Index i = 0; i < (it->second).cols(); i++) {
      (it2->second).col(i) = new_matrix.col(tperm[i]);
    }
  }
  return result;
}
}  // namespace CASM
