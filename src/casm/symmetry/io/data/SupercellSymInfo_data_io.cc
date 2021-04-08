#include "casm/symmetry/io/data/SupercellSymInfo_data_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/container/io/PermutationIO.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/io/json/SymGroup_json_io.hh"

namespace CASM {

namespace SupercellSymInfo_dataformatter_impl {

// use for ValueType=bool, int, double, std::string, jsonParser
template <typename ValueType>
using SupercellSymInfoFormatter =
    GenericDatumFormatter<ValueType, SupercellSymInfo>;

typedef Generic1DDatumFormatter<Eigen::VectorXd, SupercellSymInfo>
    VectorXdSupercellSymInfoFormatter;

typedef Generic1DDatumFormatter<Eigen::VectorXi, SupercellSymInfo>
    VectorXiSupercellSymInfoFormatter;

typedef Generic2DDatumFormatter<Eigen::MatrixXd, SupercellSymInfo>
    MatrixXdSupercellSymInfoFormatter;

typedef Generic2DDatumFormatter<Eigen::MatrixXi, SupercellSymInfo>
    MatrixXiSupercellSymInfoFormatter;

MatrixXiSupercellSymInfoFormatter transformation_matrix_to_super() {
  return MatrixXiSupercellSymInfoFormatter(
      "transformation_matrix_to_super",
      "Transformation matrix T, defining the supercell lattice vectors, S, in "
      "terms of the prim lattice vectors, P: `S = T * P`, where S and P are "
      "column vector matrices.",
      [](SupercellSymInfo const &supercell_sym_info) -> Eigen::MatrixXi {
        return supercell_sym_info.transformation_matrix_to_super().cast<int>();
      });
}

VectorXdSupercellSymInfoFormatter supercell_lattice() {
  return VectorXdSupercellSymInfoFormatter(
      "supercell_lattice",
      "Supercell lattice vectors, unrolled: (a0, a1, a2, b0, ...)",
      [](SupercellSymInfo const &supercell_sym_info) -> Eigen::VectorXd {
        Eigen::Matrix<double, 3, 3, Eigen::ColMajor> L =
            supercell_sym_info.supercell_lattice().lat_column_mat();
        return Eigen::Map<Eigen::VectorXd>(L.data(), L.size());
      });
}

MatrixXdSupercellSymInfoFormatter supercell_lattice_column_matrix() {
  return MatrixXdSupercellSymInfoFormatter(
      "supercell_lattice_column_matrix",
      "Supercell lattice vectors, as column vector matrix",
      [](SupercellSymInfo const &supercell_sym_info) -> Eigen::MatrixXd {
        return supercell_sym_info.supercell_lattice().lat_column_mat();
      });
}

MatrixXdSupercellSymInfoFormatter supercell_lattice_vectors() {
  return MatrixXdSupercellSymInfoFormatter(
      "supercell_lattice_vectors",
      "Supercell lattice vectors, as row vector matrix",
      [](SupercellSymInfo const &supercell_sym_info) -> Eigen::MatrixXd {
        Lattice const &supercell_lattice =
            supercell_sym_info.supercell_lattice();
        return supercell_lattice.lat_column_mat().transpose();
      });
}

VectorXdSupercellSymInfoFormatter supercell_lattice_params() {
  return VectorXdSupercellSymInfoFormatter(
      "supercell_lattice_params",
      "Supercell lattice parameters, as: (a, b, c, alpha, beta, gamma)",
      [](SupercellSymInfo const &supercell_sym_info) -> Eigen::VectorXd {
        Lattice const &lat = supercell_sym_info.supercell_lattice();
        Eigen::VectorXd res(6);
        res << lat.length(0), lat.length(1), lat.length(2), lat.angle(0),
            lat.angle(1), lat.angle(2);
        return res;
      });
}

SupercellSymInfoFormatter<Index> supercell_size() {
  return SupercellSymInfoFormatter<Index>(
      "supercell_size",
      "Supercell size, given as the integer number of primitive cells",
      [](SupercellSymInfo const &supercell_sym_info) -> Index {
        return supercell_sym_info.unitcell_index_converter().total_sites();
      });
}

SupercellSymInfoFormatter<double> supercell_volume() {
  return SupercellSymInfoFormatter<double>(
      "supercell_volume", "Supercell volume (length^3)",
      [](SupercellSymInfo const &supercell_sym_info) -> double {
        return supercell_sym_info.supercell_lattice().volume();
      });
}

SupercellSymInfoFormatter<jsonParser> unitcells() {
  return SupercellSymInfoFormatter<jsonParser>(
      "unitcells",
      "Integer coordinates of unit cells in the supercell. The order "
      "corresponds to the both the linear order of site degree of freedom "
      "values within one sublattice of a configuration, and the order in which "
      "translational symmetry operations are applied.",
      [](SupercellSymInfo const &supercell_sym_info) -> jsonParser {
        jsonParser json;
        to_json(xtal::make_lattice_points(
                    supercell_sym_info.transformation_matrix_to_super()),
                json, jsonParser::as_flattest());
        return json;
      });
}

SupercellSymInfoFormatter<jsonParser> integral_site_coordinates() {
  return SupercellSymInfoFormatter<jsonParser>(
      "integral_site_coordinates",
      "Integer coordinates `(b, i, j, k)` of sites in the supercell, where `b` "
      "is the sublattice index of the site, and `(i,j,k)` are the integral "
      "coordinates of the unit cell containing the site.",
      [](SupercellSymInfo const &supercell_sym_info) -> jsonParser {
        jsonParser json = jsonParser::array();
        auto f = supercell_sym_info.unitcellcoord_index_converter();
        for (Index l = 0; l < f.total_sites(); l++) {
          jsonParser l_json;
          to_json(f(l), l_json);
          json.push_back(l_json);
        }
        return json;
      });
}

SupercellSymInfoFormatter<jsonParser> translation_permutations() {
  return SupercellSymInfoFormatter<jsonParser>(
      "translation_permutations",
      "Describes how sites permute due to translation. Outer array (with size "
      "equal to the supercell volume) corresponds to possible translations "
      "within the supercell, and the inner array (with size is equal to the "
      "total number of sites, prim basis size * supercell volume) describes "
      "how sites permute for the `i`th translation, according to `after[l] = "
      "before[translation_permutation[i][l]]`.",
      [](SupercellSymInfo const &supercell_sym_info) -> jsonParser {
        jsonParser json = jsonParser::array();
        to_json(supercell_sym_info.translation_permutations(), json);
        return json;
      });
}

SupercellSymInfoFormatter<jsonParser> factor_group() {
  return SupercellSymInfoFormatter<jsonParser>(
      "factor_group",
      "The supercell factor group is the subgroup of the prim factor group "
      "that keeps the supercell lattice invariant.",
      [](SupercellSymInfo const &supercell_sym_info) -> jsonParser {
        jsonParser json;
        write_symgroup(supercell_sym_info.factor_group(), json);
        return json;
      });
}

SupercellSymInfoFormatter<Index> factor_group_size() {
  return SupercellSymInfoFormatter<Index>(
      "factor_group_size",
      "The number of operations in the supercell factor group.",
      [](SupercellSymInfo const &supercell_sym_info) -> Index {
        return supercell_sym_info.factor_group().size();
      });
}

SupercellSymInfoFormatter<jsonParser> factor_group_permutations() {
  return SupercellSymInfoFormatter<jsonParser>(
      "factor_group_permutations",
      "Describes how sites permute due to supercell factor group operations. "
      "Outer array corresponds to supercell factor group operations (prim "
      "factor group operations that keep the supercell lattice invariant), and "
      "the inner array (with size is equal to the total number of sites, prim "
      "basis size * supercell volume) describes how sites permute for the "
      "`i`th supercell factor group operation, according to `after[l] = "
      "before[factor_group_permutations[i][l]]`.",
      [](SupercellSymInfo const &supercell_sym_info) -> jsonParser {
        jsonParser json = jsonParser::array();
        auto const &factor_group = supercell_sym_info.factor_group();
        for (int i = 0; i < factor_group.size(); ++i) {
          jsonParser permute_json;
          to_json(supercell_sym_info.factor_group_permute(i), permute_json);
          json.push_back(permute_json);
        }
        return json;
      });
}

SupercellSymInfoFormatter<std::string> point_group_name() {
  return SupercellSymInfoFormatter<std::string>(
      "point_group_name", "Supercell point group name.",
      [](SupercellSymInfo const &supercell_sym_info) -> std::string {
        return supercell_sym_info.factor_group().get_name();
      });
}

SupercellSymInfoFormatter<jsonParser> global_dof_reps() {
  return SupercellSymInfoFormatter<jsonParser>(
      "global_dof_reps",
      "Symmetry representations (matrices) that describe how values of global "
      "degrees of freedom transform under application of symmetry operations. "
      "The `i`th supercell factor group operation transforms vectorized DoF "
      "values according to `vector_after = global_dof_reps[dof_key][i] * "
      "vector_before`, where `dof_key` is the name of the DoF.",
      [](SupercellSymInfo const &supercell_sym_info) -> jsonParser {
        jsonParser json;
        for (auto const &pair : supercell_sym_info.global_dof_symreps()) {
          write_matrix_rep(pair.second, json[pair.first]);
        }
        return json;
      });
}

SupercellSymInfoFormatter<jsonParser> local_dof_reps() {
  return SupercellSymInfoFormatter<jsonParser>(
      "local_dof_reps",
      "Symmetry representations (matrices) that describe how values of local "
      "degrees of freedom transform under application of symmetry operations. "
      "The `i`th supercell factor group operation transforms vectorized local "
      "DoF values on sublattice `b` according to `vector_after = "
      "local_dof_reps[dof_key][b][i] * vector_before`, where `dof_key` is the "
      "name of the DoF. After transforming values at a site, values are "
      "permuted among sites.",
      [](SupercellSymInfo const &supercell_sym_info) -> jsonParser {
        jsonParser json;
        for (auto const &pair : supercell_sym_info.local_dof_symreps()) {
          json[pair.first] =
              jsonParser::array(pair.second.size(), jsonParser::object());
          Index b = 0;
          for (auto const &group_handle : pair.second) {
            write_matrix_rep(group_handle, json[pair.first][b]);
            ++b;
          }
        }
        return json;
      });
}

}  // namespace SupercellSymInfo_dataformatter_impl

template <>
DataFormatterDictionary<SupercellSymInfo>
make_attribute_dictionary<SupercellSymInfo>() {
  using namespace SupercellSymInfo_dataformatter_impl;
  DataFormatterDictionary<SupercellSymInfo> dict;
  dict.insert(transformation_matrix_to_super(), supercell_lattice(),
              supercell_lattice_column_matrix(), supercell_lattice_vectors(),
              supercell_lattice_params(), supercell_size(), supercell_volume(),
              unitcells(), integral_site_coordinates(),
              translation_permutations(), factor_group(), factor_group_size(),
              factor_group_permutations(), global_dof_reps(), local_dof_reps());

  return dict;
}

}  // namespace CASM
