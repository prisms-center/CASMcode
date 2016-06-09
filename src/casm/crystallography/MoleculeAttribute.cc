#include "casm/misc/CASM_Math.hh"
#include "casm/crystallography/MoleculeAttribute.hh"

namespace CASM {
  namespace MoleculeAttribute_impl {
    /// \brief Equivalent to find, but throw error with suggestion if _name not found
    notstd::cloneable_ptr<BasicTraits> TraitsDictionary::lookup(const key_type &_name) const {

      auto res = this->find(_name);
      if(res != this->end()) {
        return notstd::cloneable_ptr<BasicTraits>(*res);
      }
      else {

        // If no match, try to use demerescau-levenshtein distance to make a helpful suggestion
        auto it = this->begin();
        int min_dist(-1);
        for(; it != this->end(); ++it) {
          int dist = dl_string_dist(_name, it->name());
          if(min_dist < 0 || dist < min_dist) {
            min_dist = dist;

            res = it;
          }
        }

        throw std::runtime_error("CRITICAL ERROR: Invalid Molecule attribute \"" + _name + "\" specified.\n"
                                 + "                Did you mean \"" + res->name() + "\"?\n");

      }

    }
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
    m_traits_ptr = MoleculeAttribute_impl::traits_dict().lookup(name());
  }
}
