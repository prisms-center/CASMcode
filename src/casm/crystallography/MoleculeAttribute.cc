#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/MoleculeAttribute.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {
  namespace MoleculeAttribute_impl {
  }

  //*******************************************************************

  MoleculeAttribute_impl::TraitsDictionary &MoleculeAttribute::_traits_dict() {
    static MoleculeAttribute_impl::TraitsDictionary _static_dict;

    return _static_dict;
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

  //*******************************************************************

  void MoleculeAttribute::_load_traits() const {
    m_traits_ptr = _traits_dict().lookup(name());
  }
}
