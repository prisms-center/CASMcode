#include "casm/symmetry/SymMatrixXd.hh"

namespace CASM {


  jsonParser &SymMatrixXd::to_json(jsonParser &json) const {
    json.put_obj();

    // Members not included:
    //
    // From SymOpRepresentation:
    //   MasterSymGroup const *head_group;

    json["SymOpRep_type"] = "SymMatrixXd";

    ///type of symmetry, given by one of the allowed values of symmetry_type
    json["symmetry"] = symmetry;
    json["op_index"] = op_index;
    json["rep_ID"] = rep_ID;

    json["mat"] = mat;
    return json;
  }

  //*******************************************************************************************

  void SymMatrixXd::from_json(const jsonParser &json) {
    try {
      CASM::from_json(symmetry, json["symmetry"]);
      CASM::from_json(op_index, json["op_index"]);
      CASM::from_json(rep_ID, json["rep_ID"]);

      Eigen::from_json(mat, json["mat"]);
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
