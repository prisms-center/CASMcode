#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/MoleculeAttribute.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {
  namespace MoleculeAttribute_impl {
  }

  template<>
  ParsingDictionary<MoleculeAttribute::BasicTraits>  make_parsing_dictionary<MoleculeAttribute::BasicTraits>() {
    ParsingDictionary<MoleculeAttribute::BasicTraits> dict;
    return dict;
  }

  //*******************************************************************

  bool MoleculeAttribute::identical(MoleculeAttribute const &other, double _tol) const {
    return name() == other.name() && almost_equal(value(), other.value(), _tol);
  }

  //*******************************************************************

  MoleculeAttribute &MoleculeAttribute::apply_sym(SymOp const &_op) {
    if(m_rep_ID.empty())
      _generate_symrep(_op.master_group());
    m_value = (*(_op.get_matrix_rep(m_rep_ID))) * m_value;
    return *this;
  }

  //*******************************************************************

  void MoleculeAttribute::_generate_symrep(MasterSymGroup const &_group) {
    m_rep_ID = _traits().generate_symrep(_group);
  }

}
