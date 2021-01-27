#ifndef SYMTYPECOMPARATOR_HH
#define SYMTYPECOMPARATOR_HH

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymType.hh"

namespace CASM {
namespace xtal {
class Lattice;

/**
 * Binary functor for checking if two SymOp are equal.
 * Compares the transformation matrix and translation to within the CASM
 * tolerance, and also checks for time reversal match.
 */

struct SymOpCompare_f {
  explicit SymOpCompare_f(double tolerance) : m_tolerance(tolerance) {}
  bool operator()(const SymOp &lhs, const SymOp &rhs) const;

 private:
  double m_tolerance;
};

/**
 * Binary functor for checking if two SymOp are equal.
 * Compares the transformation matrix and translation to within the CASM
 * tolerance, and also checks for time reversal match.
 * A lattice is given at construction so that comparisons return true if
 * the operation is equivalent by lattice translations.
 */

struct SymOpPeriodicCompare_f {
  explicit SymOpPeriodicCompare_f(const Lattice &periodicity_lattice,
                                  double tolerance)
      : m_periodicity_lattice(periodicity_lattice), m_tolerance(tolerance) {}
  bool operator()(const SymOp &lhs, const SymOp &rhs) const;

 private:
  const Lattice m_periodicity_lattice;
  double m_tolerance;
};

/**
 * Binary functor for checking if two SymOp are equal.
 * only the transformation matrix. When using this comparator
 * both translation and time reversal are ignored. Returns true
 * if the transformation matrix matches.
 */

struct SymOpMatrixCompare_f {
  explicit SymOpMatrixCompare_f(double tolerance) : m_tolerance(tolerance) {}
  bool operator()(const SymOp &lhs, const SymOp &rhs) const;

 private:
  double m_tolerance;
};

//*********************************************************************************************************************//

}  // namespace xtal
}  // namespace CASM

#endif
