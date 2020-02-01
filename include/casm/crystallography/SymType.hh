#ifndef XTALSYMTYPE_HH
#define XTALSYMTYPE_HH

#include "casm/external/Eigen/Dense"
#include <algorithm>
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
      explicit SymOpCompare_f(const SymOp &target_operation, double tolerance) : m_target_operation(target_operation), m_tolerance(tolerance) {
      }
      bool operator()(const SymOp &possible_match) const;

    private:
      const SymOp &m_target_operation;
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
        : m_target_operation(target_operation), m_periodicity_lattice(periodicity_lattice), m_tolerance(tolerance) {
      }
      bool operator()(const SymOp &possible_match) const;

    private:
      const SymOp &m_target_operation;
      const Lattice &m_periodicity_lattice;
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
      explicit SymOpMatrixCompare_f(const SymOp &target_operation, double tolerance)
        : m_target_operation(target_operation), m_tolerance(tolerance) {
      }
      bool operator()(const SymOp &possible_match) const;

    private:
      const SymOp &m_target_operation;
      double m_tolerance;
    };

    //*********************************************************************************************************************//

    // TODO: Sort out how to handle a maximum group size? A bad closure could result in an infinite group.
    /// Combines every pair of symmetry operations and adds any missing resulting operations to the group.
    /// In order to determine what a missing operation is, a comparator type must be provided, and any
    /// arguments needed to construct it must be passed as final arguments.
    ///
    /// For example, when closing a factor group, we consider operations to be the same when the transformation
    /// matrix matches to within a tolerance, and the translation components are equivalent under translational symmetry.
    /// We therefore want to use the SymOpPeriodicCompare_f comparator, which requires a lattice and a tolerance.
    /// The group closure call for this situation might look like:
    ///
    /// close_group<SymOpPeriodicCompare_f> close_group(&my_partial_group, my_lattice, CASM::TOL);
    template <typename SymOpCompareType, typename... CompareArgs>
    void close_group(SymOpVector *partial_group, const CompareArgs &... args) {
      bool is_closed = false;
      while(!is_closed) {
        is_closed = true;
        int current_size = partial_group->size();

        for(int i = 0; i < current_size; ++i) {
          for(int j = 0; j < current_size; ++j) {
            const SymOp &l_op = partial_group->at(i);
            const SymOp &r_op = partial_group->at(j);
            SymOp candidate = l_op * r_op;
            SymOpCompareType compare_candidate(candidate, args...);

            //If you can't find find the operation then the group wasn't closed.
            //Add it and make sure to start over to continue closing with the new operations.
            if(std::find_if(partial_group->begin(), partial_group->end(), compare_candidate) == partial_group->end()) {
              partial_group->push_back(candidate);
              is_closed = false;
            }
          }
        }
      }

      return;
    }
  } // namespace xtal
} // namespace CASM

#endif
