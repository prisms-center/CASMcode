#include "casm/symmetry/io/json/SymGroup_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/container/io/PermutationIO.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/io/json/SymInfo_json_io.hh"

namespace CASM {

// --------- SymmetryIO Declarations
// --------------------------------------------------

void write_symop(SymGroup const &grp, Index i, jsonParser &j) {
  j = jsonParser::object();

  const SymOp &op = grp[i];

  j["master_group_index"] = op.master_group_index();

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

  to_json(grp.info(i), j["info"]);

  to_json(to_brief_unicode(grp.info(i), SymInfoOptions(CART)),
          j["info"]["brief"]["CART"]);
  to_json(to_brief_unicode(grp.info(i), SymInfoOptions(FRAC)),
          j["info"]["brief"]["FRAC"]);
}

void write_symgroup(SymGroup const &grp, jsonParser &json) {
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

/// \brief Describes how integral site coordinates transform under application
/// of symmetry.
///
/// Writes an array, with one element for each group element, containing:
/// - matrix
/// - sublattice_permute
/// - sublattice_shift
///
/// The ith factor group operation transforms a basis site:
///
///     (b, r_frac) -> (b', r_frac')
///
/// according to:
///
///     b' = sublattice_permute[b]
///     r_frac' = matrix * r_frac + sublattice_shift[b]`
///
/// where `b` is the basis index and `r_frac` is the integer unit cell
/// coordinate of a site.
void write_basis_permutation_rep(SymGroup const &grp,
                                 jsonParser &group_rep_json,
                                 SymGroupRepID symgrouprep_id) {
  group_rep_json = jsonParser::array();
  for (auto const &op : grp) {
    jsonParser op_rep_json = jsonParser::object();
    SymBasisPermute const &op_rep = *op.get_basis_permute_rep(symgrouprep_id);
    op_rep_json["master_group_index"] = op_rep.master_group_index();
    op_rep_json["matrix"] = op_rep.matrix();
    op_rep_json["sublattice_permute"] = jsonParser::array();
    op_rep_json["sublattice_shift"] = jsonParser::array();
    for (auto const integral_coordinate : op_rep.data()) {
      op_rep_json["sublattice_permute"].push_back(
          integral_coordinate.sublattice());

      jsonParser tjson;
      to_json(integral_coordinate.unitcell(), tjson, jsonParser::as_flattest());
      op_rep_json["sublattice_shift"].push_back(tjson);
    }
    group_rep_json.push_back(op_rep_json);
  }
}

/// Writes a 3d array, where occ_permutation_rep[b][i] is the permutation array
/// for occupant values on sublattice 'b', due to group operation grp[i],
/// according to: occ(l) = occ_permutation_rep[b][i][occ(l)]. Occupation values
/// are transformed on each site in this way before begin permuted among sites.
void write_occ_permutation_rep(SymGroup const &grp, jsonParser &json,
                               std::vector<SymGroupRepID> occupant_symrep_IDs) {
  json = jsonParser::array();
  for (auto symrep_ID : occupant_symrep_IDs) {
    jsonParser sublattice_json = jsonParser::array();
    for (auto const &op : grp) {
      Permutation const &permute = *op.get_permutation_rep(symrep_ID);
      jsonParser permute_json;
      to_json(permute, permute_json);
      sublattice_json.push_back(permute_json);
    }
    json.push_back(sublattice_json);
  }
}

void write_matrix_rep(SymGroupRepHandle const &grp, jsonParser &json) {
  json = jsonParser::array();
  for (Index i = 0; i < grp.size(); ++i) {
    json.push_back(*grp[i]->MatrixXd());
  }
}

}  // namespace CASM
