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
  namespace xtal {
    namespace Adapter {
      typedef Eigen::Matrix3d SymOpMatrixType;
      typedef Eigen::Vector3d SymOpTranslationType;
      typedef bool SymOpTimeReversalType;

      /// Within the scope of crystallography, this struct will serve as the symmetry
      /// object, which holds a transformation matrix, translation vector, and time
      /// reversal boolean, whithout any other overhead.
      struct SymOpData {
        SymOpData(const SymOpMatrixType &mat, const SymOpTranslationType &translation, SymOpTimeReversalType time_reversal)
          : matrix(mat), translation(translation), is_time_reversal_active(time_reversal) {
        }

        SymOpMatrixType matrix;
        SymOpTranslationType translation;
        SymOpTimeReversalType is_time_reversal_active;
      };

      /// This defines the type of the object representing symmetry operations within the crystallography
      /// classes. Any symmetry related operations within the crystallography module must be in terms
      /// of this type.
      /* typedef std::tuple<SymOpMatrixType, SymOpTranslationType, SymOpTimeReversalType> SymOpType; */
      typedef SymOpData SymOpType;
      typedef std::vector<SymOpType> SymGroupType;

      // TODO:
      // How do we want to go about this exactly? The initial idea was to have these accessors defined
      // outside of the SymOpType, so that if we decide to swap the SymOpType (e.g. to a tuple or something)
      // we just have to change the accessors here.
      // I made these templates because there's another option instead. For situations where a function
      // accepts a SymOp, we can either change the function signature so that it accepts a SymOpType,
      // or we can template the function so that it takes any type, and will work as long as accessors
      // like the ones below are defined. To an extent, this can eliminate the need for a concrete definition
      // of SymOpType.

      template <typename ExternSymOpType>
      SymOpMatrixType &get_matrix(ExternSymOpType &op);
      template <typename ExternSymOpType>
      const SymOpMatrixType &get_matrix(const ExternSymOpType &op);

      template <typename ExternSymOpType>
      SymOpTranslationType &get_translation(ExternSymOpType &op);
      template <typename ExternSymOpType>
      const SymOpTranslationType &get_translation(const ExternSymOpType &op);

      template <typename ExternSymOpType>
      SymOpTimeReversalType get_time_reversal(const SymOpType &op);

      /// For cases where you have a symmetry operation type that is not SymOpType, you should
      /// have definitions of these template functions, so that
      template <typename ExternSymOpType>
      SymOpType to_symop_type(const ExternSymOpType &op);

      template <typename ExternSymGroupType>
      SymGroupType to_symgroup_type(const ExternSymGroupType &group);
      template <typename ExternSymGroupTypeIt>
      SymGroupType to_symgroup_type(ExternSymGroupTypeIt begin, ExternSymGroupTypeIt end) {
        SymGroupType casted_group;
        for(auto it = begin; it != end; ++it) {
          casted_group.emplace_back(to_symop_type(*it));
        }
        return casted_group;
      }
    } // namespace Adapter

  } // namespace xtal
} // namespace CASM

#endif
