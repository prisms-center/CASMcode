#include "casm/symmetry/SymPermutation.hh"

#include "casm/container/Permutation.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  double SymPermutation::character() const {
    int n_fix(0);
    for(Index i = 0; i < m_permute.size(); i++) {
      n_fix += int(i == m_permute[i]);
    }
    return n_fix;
  }

  //*******************************************************************************************

  //Makes permutation matrix from m_permute
  //Permute matrix P reorders elements of a column vector v via
  // v' = P*v
  //Therefore, we loop over m_permute and assign ones to the
  //appropriate row/column

  void SymPermutation::_calc_mat() const {
    m_mat.resize(m_permute.size(), m_permute.size());
    m_mat.setZero();
    for(Index i = 0; i < m_permute.size(); i++)
      m_mat(i, m_permute[i]) = 1;
    return;
  }

  //*******************************************************************************************

  jsonParser &SymPermutation::to_json(jsonParser &json) const {
    json.put_obj();

    // Members not included:
    //
    // From SymOpRepresentation:
    //   MasterSymGroup const *head_group;

    json["SymOpRep_type"] = "SymPermutation";

    json["op_index"] = index();
    json["rep_ID"] = rep_ID();

    json["m_permute"] = m_permute;
    if(m_has_mat)
      json["m_mat"] = m_mat;
    return json;
  }

  //*******************************************************************************************

  void SymPermutation::from_json(const jsonParser &json) {
    try {
      CASM::from_json(m_op_index, json["op_index"]);
      CASM::from_json(m_rep_ID, json["rep_ID"]);

      CASM::from_json(m_permute, json["m_permute"]);
      CASM::from_json(m_mat, json["m_mat"]);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //*******************************************************************************************

  jsonParser &to_json(const SymPermutation &sym, jsonParser &json) {
    return sym.to_json(json);
  }

  //*******************************************************************************************

  void from_json(SymPermutation &sym, const jsonParser &json) {
    try {
      sym.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

}
