#include "casm/misc/CASM_Math.hh"
#include "casm/crystallography/MoleculeAttribute.hh"

namespace CASM {
  namespace MoleculeAttribute_impl {
    /// \brief Equivalent to find, but throw error with suggestion if _name not found
    notstd::cloneable_ptr<BasicTraits> TraitsDictionary::lookup(const key_type &_name) const {

      auto res = find(_name);
      if(res != end()) {
        return res->clone();
      }
      else {

        // If no match, try to use demerescau-levenshtein distance to make a helpful suggestion
        auto it = begin();
        int min_dist(-1);
        for(; it != end(); ++it) {
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

  MoleculeAttribute_impl::TraitsDictionary &MoleculeAttribute::_traits_dict() {
    static MoleculeAttribute_impl::TraitsDictionary _static_dict;

    return _static_dict;
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
