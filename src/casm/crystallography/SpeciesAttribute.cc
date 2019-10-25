#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/SpeciesAttribute.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  SpeciesAttribute &SpeciesAttribute::apply_sym(SymOp const &_op) {
    m_value = traits().symop_to_matrix(_op.matrix(), _op.tau(), _op.time_reversal()) * m_value;
    return *this;
  }

  //*******************************************************************

  bool SpeciesAttribute::identical(SpeciesAttribute const &other, double _tol) const {
    return name() == other.name() && almost_equal(value(), other.value(), _tol);
  }
}
