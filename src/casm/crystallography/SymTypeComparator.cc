#include "casm/crystallography/SymTypeComparator.hh"

#include "casm/crystallography/Coordinate.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {

bool SymOpCompare_f::operator()(const SymOp &lhs, const SymOp &rhs) const {
  if (lhs.is_time_reversal_active != rhs.is_time_reversal_active) {
    return false;
  }

  if (!almost_equal(lhs.translation, rhs.translation, this->m_tolerance)) {
    return false;
  }

  if (!almost_equal(lhs.matrix, rhs.matrix, this->m_tolerance)) {
    return false;
  }

  return true;
}

bool SymOpMatrixCompare_f::operator()(const SymOp &lhs,
                                      const SymOp &rhs) const {
  return almost_equal(lhs.matrix, rhs.matrix, this->m_tolerance);
}

bool SymOpPeriodicCompare_f::operator()(const SymOp &lhs,
                                        const SymOp &rhs) const {
  if (lhs.is_time_reversal_active != rhs.is_time_reversal_active) {
    return false;
  }

  if (!almost_equal(lhs.matrix, rhs.matrix, this->m_tolerance)) {
    return false;
  }

  xtal::Coordinate this_translation_coord(lhs.translation,
                                          this->m_periodicity_lattice, CART);
  xtal::Coordinate match_translation_coord(rhs.translation,
                                           this->m_periodicity_lattice, CART);

  return this_translation_coord.min_dist(match_translation_coord) <
         this->m_tolerance;
}

}  // namespace xtal
}  // namespace CASM
