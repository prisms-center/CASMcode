#ifndef SYMMETRYADAPTER_HH
#define SYMMETRYADAPTER_HH

#include "casm/external/Eigen/Dense"
#include <tuple>
#include <vector>

/* namespace */
/* { */
/*     SymOpMatrixType symop_to_matrix(const CASM::SymOp &op); */

/*     template<typename SymGroupType> */
/*     std::vector<SymOpMatrixType> symop_to_matrix(const SymGroupType &group) { */
/*       std::vector<SymOpMatrixType> casted_group; */
/*       casted_group.reserve(group.size()); */
/*       for(const auto &op : group) { */
/*         casted_group.push_back(symop_to_matrix(op)); */
/*       } */

/*       return casted_group; */
/*     } */
/* } */

namespace CASM {
  namespace Adapter {
    typedef Eigen::Matrix3d SymOpMatrixType;
    typedef Eigen::Vector3d SymOpTranslationType;
    typedef bool SymOpTimeReversalType;
    typedef std::tuple<SymOpMatrixType, SymOpTranslationType, SymOpTimeReversalType> SymOpType;

    SymOpMatrixType &get_matrix(SymOpType &op);
    const SymOpMatrixType &get_matrix(const SymOpType &op);

    SymOpTranslationType &get_translation(SymOpType &op);
    const SymOpTranslationType &get_translation(const SymOpType &op);

    SymOpTimeReversalType get_time_reversal(const SymOpType &op);

    SymOpType construct_sym_op(const SymOpMatrixType &sym_matrix,
                               const SymOpTranslationType &sym_translation,
                               SymOpTimeReversalType &sym_time_reversal);

  } // namespace Adapter

} // namespace CASM

#endif
