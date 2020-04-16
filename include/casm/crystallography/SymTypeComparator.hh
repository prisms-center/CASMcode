#ifndef SYMTYPECOMPARATOR_HH
#define SYMTYPECOMPARATOR_HH

#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {
  namespace xtal {
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
      const SymOp m_target_operation;
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
      const SymOp m_target_operation;
      const Lattice m_periodicity_lattice;
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
      const SymOp m_target_operation;
      double m_tolerance;
    };

    //*********************************************************************************************************************//

  }
}

#endif
