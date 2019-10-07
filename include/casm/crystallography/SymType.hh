#ifndef XTALSYMTYPE_HH
#define XTALSYMTYPE_HH

#include "casm/external/Eigen/Dense"
#include <tuple>
#include <vector>

namespace CASM {
  namespace xtal {
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

    /// Accessor for SymOpType. Returns transformation matrix (Cartesian).
    const SymOpMatrixType &get_matrix(const SymOpType &op);
    /// Accessor for SymOpType. Returns translation vector (tau).
    const SymOpTranslationType &get_translation(const SymOpType &op);
    /// Accessor for SymOpType. Returns whether the symmetry operation is time reversal active.
    SymOpTimeReversalType get_time_reversal(const SymOpType &op);

  } // namespace xtal
} // namespace CASM

#endif
