#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
  namespace xtal {
    SymOp operator*(const SymOp &LHS, const SymOp &RHS) {
      return SymOp(LHS.matrix * RHS.matrix,
                   LHS.translation + LHS.matrix * RHS.translation,
                   LHS.is_time_reversal_active != RHS.is_time_reversal_active); //This is an XOR operation
    }

    //**********************************************************************************//

    const SymOpMatrixType &get_matrix(const SymOp &op) {
      return op.matrix;
    }

    const SymOpTranslationType &get_translation(const SymOp &op) {
      return op.translation;
    }

    SymOpTimeReversalType get_time_reversal(const SymOp &op) {
      return op.is_time_reversal_active;
    }

    //**********************************************************************************//

    bool SymOpCompare_f::operator()(const SymOp &possible_match) const {
      if(this->m_target_operation.is_time_reversal_active != possible_match.is_time_reversal_active) {
        return false;
      }

      if(!almost_equal(this->m_target_operation.translation, possible_match.translation, this->m_tolerance)) {
        return false;
      }

      if(!almost_equal(this->m_target_operation.matrix, possible_match.matrix, this->m_tolerance)) {
        return false;
      }

      return true;
    }

    bool SymOpMatrixCompare_f::operator()(const SymOp &possible_match) const {
      return almost_equal(this->m_target_operation.matrix, possible_match.matrix, this->m_tolerance);
    }

    bool SymOpPeriodicCompare_f::operator()(const SymOp &possible_match) const {
      if(this->m_target_operation.is_time_reversal_active != possible_match.is_time_reversal_active) {
        return false;
      }

      if(!almost_equal(this->m_target_operation.matrix, possible_match.matrix, this->m_tolerance)) {
        return false;
      }

      xtal::Coordinate this_translation_coord(this->m_target_operation.translation, this->m_periodicity_lattice, CART);
      xtal::Coordinate match_translation_coord(possible_match.translation, this->m_periodicity_lattice, CART);

      return this_translation_coord.min_dist(match_translation_coord) < this->m_tolerance;
    }

  } // namespace xtal
} // namespace CASM
