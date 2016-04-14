#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {


  jsonParser &SymMatrixXd::to_json(jsonParser &json) const {
    json.put_obj();

    // Members not included:
    //
    // From SymOpRepresentation:
    //   MasterSymGroup const *head_group;

    json["SymOpRep_type"] = "SymMatrixXd";

    json["op_index"] = index();
    json["rep_ID"] = rep_ID();

    json["mat"] = mat;
    return json;
  }

  //*******************************************************************************************

  void SymMatrixXd::from_json(const jsonParser &json) {
    try {
      CASM::from_json(m_op_index, json["op_index"]);
      CASM::from_json(m_rep_ID, json["rep_ID"]);

      CASM::from_json(mat, json["mat"]);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //*******************************************************************************************

  jsonParser &to_json(const SymMatrixXd &sym, jsonParser &json) {
    return sym.to_json(json);
  }

  //*******************************************************************************************

  void from_json(SymMatrixXd &sym, const jsonParser &json) {
    try {
      sym.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

}
