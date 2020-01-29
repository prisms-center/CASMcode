#ifndef XTALSYMTYPE_HH
#define XTALSYMTYPE_HH

#include "casm/external/Eigen/Dense"
#include <functional>
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
    struct SymOp {
      SymOp(const SymOpMatrixType &mat, const SymOpTranslationType &translation, SymOpTimeReversalType time_reversal)
        : matrix(mat), translation(translation), is_time_reversal_active(time_reversal) {
      }

      static SymOp identity() {
        return SymOp(SymOpMatrixType::Identity(), SymOpTranslationType::Zero(), false);
      }

      static SymOp time_reversal() {
        return SymOp(SymOpMatrixType::Identity(), SymOpTranslationType::Zero(), true);
      }

      static SymOp translation_operation(const SymOpTranslationType &translation) {
        return SymOp(SymOpMatrixType::Identity(), translation, false);
      }

      static SymOp point_operation(const SymOpMatrixType &mat) {
        return SymOp(mat, SymOpTranslationType::Zero(), false);
      }

      SymOpMatrixType matrix;
      SymOpTranslationType translation;
      SymOpTimeReversalType is_time_reversal_active;
    };

    /// Get a new SymOp that is equivalent to subsequent application of both SymOps
    SymOp operator*(const SymOp &LHS, const SymOp &RHS);

    /// This defines the type of the object representing symmetry operations within the crystallography
    /// classes. Any symmetry related operations within the crystallography module must be in terms
    /// of this type.
    /* typedef std::tuple<SymOpMatrixType, SymOpTranslationType, SymOpTimeReversalType> SymOpType; */
    typedef std::vector<SymOp> SymOpVector;

    /// Accessor for SymOpType. Returns transformation matrix (Cartesian).
    const SymOpMatrixType &get_matrix(const SymOp &op);
    /// Accessor for SymOpType. Returns translation vector (tau).
    const SymOpTranslationType &get_translation(const SymOp &op);
    /// Accessor for SymOpType. Returns whether the symmetry operation is time reversal active.
    SymOpTimeReversalType get_time_reversal(const SymOp &op);

    //*********************************************************************************************************************//

    class Lattice;

    /**
     * Unary predicate for finding SymOp in a container such as std::vector.
     * Compares the transformation matrix and translation to within the CASM
     * tolerance, and also checks for time reversal match.
     *
     * Members provided at construction must still exist externally when
     * using this functor.
     */

    struct SymOpCompare_f : public std::unary_function<SymOp, bool> {
      explicit SymOpCompare_f(const SymOp &target_operation, double tolerance) : m_target_operation(&target_operation), m_tolerance(tolerance) {}
      bool operator()(const SymOp &possible_match);

    private:
      const SymOp *m_target_operation;
      double m_tolerance;
    };

    /**
     * Unary predicate for finding SymOp in a container such as std::vector.
     * Compares the transformation matrix and translation to within the CASM
     * tolerance, and also checks for time reversal match.
     * A lattice is given at construction so that comparisons return true if
     * the operation is equivalent by lattice translations.
     *
     * Members provided at construction must still exist externally when
     * using this functor.
     */

    struct SymOpPeriodicCompare_f : public std::unary_function<SymOp, bool> {
      explicit SymOpPeriodicCompare_f(const SymOp &target_operation, const Lattice &periodicity_lattice, double tolerance)
        : m_target_operation(&target_operation), m_periodicity_lattice(&periodicity_lattice), m_tolerance(tolerance) {
      }
      bool operator()(const SymOp &possible_match);

    private:
      const SymOp *m_target_operation;
      const Lattice *m_periodicity_lattice;
      double m_tolerance;
    };

    /**
     * Unary predicate for finding SymOp in a container by considering
     * only the transformation matrix. When using this comparator
     * both translation and time reversal are ignored. Returns true
     * if the transformation matrix matches.
     *
     * Members provided at construction must still exist externally when
     * using this functor.
     */

    struct SymOpMatrixCompare_f : public std::unary_function<SymOp, bool> {
      explicit SymOpMatrixCompare_f(const SymOp &target_operation, double tolerance) : m_target_operation(&target_operation), m_tolerance(tolerance) {}
      bool operator()(const SymOp &possible_match);

    private:
      const SymOp *m_target_operation;
      double m_tolerance;
    };

    //*********************************************************************************************************************//

    /// Combines every pair of symmetry operations and adds any missing resulting operations to the group.
    /// Comparisons do not take any sort of periodicity into account.
    void close_group(SymOpVector *partial_group);

    /// Combines every pair of symmetry operations and adds any missing resulting operations to the group.
    /// Operations are considered equivalent if the translations are equivalent by lattice tranlsations.
    void close_group(SymOpVector *partial_group, const Lattice &periodicity_lattice);
  } // namespace xtal
} // namespace CASM

#endif
