#include "casm/symmetry/io/json/SymGroup_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymInfo.hh"

namespace CASM {

// --------- SymmetryIO Declarations
// --------------------------------------------------

void write_symop(const SymGroup &grp, Index i, jsonParser &j) {
  j = jsonParser::object();

  const SymOp &op = grp[i];

  to_json(op.matrix(), j["CART"]["matrix"]);
  to_json_array(op.tau(), j["CART"]["tau"]);
  to_json(op.time_reversal(), j["CART"]["time_reversal"]);

  to_json(grp.lattice().inv_lat_column_mat() * op.matrix() *
              grp.lattice().lat_column_mat(),
          j["FRAC"]["matrix"]);
  to_json_array(grp.lattice().inv_lat_column_mat() * op.tau(),
                j["FRAC"]["tau"]);
  to_json(op.time_reversal(), j["FRAC"]["time_reversal"]);

  to_json(grp.class_of_op(i) + 1, j["info"]["conjugacy_class"]);
  to_json(grp.ind_inverse(i) + 1, j["info"]["inverse_operation"]);

  add_sym_info(grp.info(i), j["info"]);
}

void write_symgroup(const SymGroup &grp, jsonParser &json) {
  json = jsonParser::object();
  {
    jsonParser &json_ops = json["group_operations"];
    for (int i = 0; i < grp.size(); i++) {
      std::string op_name = "op_" + to_sequential_string(i + 1, grp.size());
      write_symop(grp, i, json_ops[op_name]);
    }
  }
  {
    jsonParser &json_info = json["group_classification"];
    json_info["name"] = grp.get_name();
    json_info["latex_name"] = grp.get_latex_name();
    json_info["periodicity"] = grp.periodicity();
    if (grp.periodicity() == PERIODIC)
      json_info["possible_space_groups"] = grp.possible_space_groups();
  }

  {
    jsonParser &json_struc = json["group_structure"];
    auto const &classes = grp.get_conjugacy_classes();

    for (Index c = 0; c < classes.size(); ++c) {
      std::string class_name =
          "class_" + to_sequential_string(c + 1, classes.size());
      jsonParser &json_class = json_struc["conjugacy_classes"][class_name];

      SymInfo info = grp.info(classes[c][0]);
      json_class["operation_type"] = info.op_type;
      if (info.op_type == symmetry_type::rotation_op ||
          info.op_type == symmetry_type::screw_op ||
          info.op_type == symmetry_type::rotoinversion_op) {
        json_class["rotation_angle"] = info.angle;
      }
      json_class["operations"].put_array();
      for (Index o : classes[c]) {
        json_class["operations"].push_back(o + 1);
      }
    }
    auto multi_table = grp.get_multi_table();
    for (auto &v : multi_table) {
      for (auto &i : v) {
        ++i;
      }
    }
    json_struc["multiplication_table"] = multi_table;
  }

  // Character table excluded for now. Will revisit after SymGroup is refactored
  // json["character_table"] = grp.character_table();
}

}  // namespace CASM
