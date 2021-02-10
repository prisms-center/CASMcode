#include "casm/crystallography/io/data/SharedPrim_data_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/io/json/SymGroup_json_io.hh"

namespace CASM {

namespace SharedPrim_dataformatter_impl {

typedef std::shared_ptr<const Structure> SharedPrim;

// use for ValueType=bool, int, double, std::string, jsonParser
template <typename ValueType>
using SharedPrimFormatter = GenericDatumFormatter<ValueType, SharedPrim>;

typedef Generic1DDatumFormatter<Eigen::VectorXd, SharedPrim>
    VectorXdSharedPrimFormatter;

typedef Generic1DDatumFormatter<Eigen::VectorXi, SharedPrim>
    VectorXiSharedPrimFormatter;

typedef Generic2DDatumFormatter<Eigen::MatrixXd, SharedPrim>
    MatrixXdSharedPrimFormatter;

typedef Generic2DDatumFormatter<Eigen::MatrixXi, SharedPrim>
    MatrixXiSharedPrimFormatter;

// -- lattice point group info

SharedPrimFormatter<jsonParser> lattice_point_group() {
  return SharedPrimFormatter<jsonParser>(
      "lattice_point_group", "Lattice point group, in JSON format.",
      [](SharedPrim const &shared_prim) -> jsonParser {
        SymGroup lattice_pg{
            SymGroup::lattice_point_group(shared_prim->lattice())};
        jsonParser json;
        write_symgroup(lattice_pg, json);
        return json;
      });
}

SharedPrimFormatter<std::string> lattice_point_group_name() {
  return SharedPrimFormatter<std::string>(
      "lattice_point_group_name", "Lattice point group name.",
      [](SharedPrim const &shared_prim) -> std::string {
        SymGroup lattice_pg{
            SymGroup::lattice_point_group(shared_prim->lattice())};
        return lattice_pg.get_name();
      });
}

SharedPrimFormatter<Index> lattice_point_group_size() {
  return SharedPrimFormatter<Index>(
      "lattice_point_group_size", "Lattice point group size.",
      [](SharedPrim const &shared_prim) -> Index {
        SymGroup lattice_pg{
            SymGroup::lattice_point_group(shared_prim->lattice())};
        return lattice_pg.size();
      });
}

// -- factor group info

SharedPrimFormatter<jsonParser> factor_group() {
  return SharedPrimFormatter<jsonParser>(
      "factor_group", "Factor group, in JSON format.",
      [](SharedPrim const &shared_prim) -> jsonParser {
        jsonParser json;
        write_symgroup(shared_prim->factor_group(), json);
        return json;
      });
}

SharedPrimFormatter<std::string> factor_group_name() {
  return SharedPrimFormatter<std::string>(
      "factor_group_name", "Factor group name.",
      [](SharedPrim const &shared_prim) -> std::string {
        return shared_prim->factor_group().get_name();
      });
}

SharedPrimFormatter<Index> factor_group_size() {
  return SharedPrimFormatter<Index>("factor_group_size", "Factor group size.",
                                    [](SharedPrim const &shared_prim) -> Index {
                                      return shared_prim->factor_group().size();
                                    });
}

// -- crystal point group info

SharedPrimFormatter<jsonParser> crystal_point_group() {
  return SharedPrimFormatter<jsonParser>(
      "crystal_point_group", "Crystal point group, in JSON format.",
      [](SharedPrim const &shared_prim) -> jsonParser {
        jsonParser json;
        write_symgroup(shared_prim->point_group(), json);
        return json;
      });
}

SharedPrimFormatter<std::string> crystal_point_group_name() {
  return SharedPrimFormatter<std::string>(
      "crystal_point_group_name", "Crystal point group name.",
      [](SharedPrim const &shared_prim) -> std::string {
        return shared_prim->point_group().get_name();
      });
}

SharedPrimFormatter<Index> crystal_point_group_size() {
  return SharedPrimFormatter<Index>("crystal_point_group_size",
                                    "Crystal point group size",
                                    [](SharedPrim const &shared_prim) -> Index {
                                      return shared_prim->point_group().size();
                                    });
}

// -- misc.

SharedPrimFormatter<double> volume() {
  return SharedPrimFormatter<double>(
      "volume", "Prim cell volume (length^3)",
      [](SharedPrim const &shared_prim) -> double {
        return shared_prim->lattice().volume();
      });
}

VectorXdSharedPrimFormatter lattice() {
  return VectorXdSharedPrimFormatter(
      "lattice", "Lattice vectors, unrolled: (a0, a1, a2, b0, ...)",
      [](SharedPrim const &shared_prim) -> Eigen::VectorXd {
        Eigen::Matrix<double, 3, 3, Eigen::ColMajor> L =
            shared_prim->lattice().lat_column_mat();
        return Eigen::Map<Eigen::VectorXd>(L.data(), L.size());
      });
}

MatrixXdSharedPrimFormatter lattice_column_matrix() {
  return MatrixXdSharedPrimFormatter(
      "lattice_column_matrix", "Lattice vectors, as column vector matrix",
      [](SharedPrim const &shared_prim) -> Eigen::MatrixXd {
        return shared_prim->lattice().lat_column_mat();
      });
}

VectorXdSharedPrimFormatter lattice_params() {
  return VectorXdSharedPrimFormatter(
      "lattice_params", "Lattice parameters, as: (a, b, c, alpha, beta, gamma)",
      [](SharedPrim const &shared_prim) -> Eigen::VectorXd {
        Eigen::VectorXd res(6);
        Lattice const &lattice = shared_prim->lattice();
        res << lattice.length(0), lattice.length(1), lattice.length(2),
            lattice.angle(0), lattice.angle(1), lattice.angle(2);
        return res;
      });
}

SharedPrimFormatter<jsonParser> asymmetric_unit() {
  return SharedPrimFormatter<jsonParser>(
      "asymmetric_unit", "Sets of indices of equivalent basis sites",
      [](SharedPrim const &shared_prim) -> jsonParser {
        adapter::Adapter<xtal::SymOpVector, SymGroup> adapter;
        auto symop_vector = adapter(shared_prim->factor_group().begin(),
                                    shared_prim->factor_group().end());
        auto asym_unit =
            xtal::make_asymmetric_unit(shared_prim->structure(), symop_vector);
        jsonParser json;
        to_json(asym_unit, json);
        return json;
      });
}

}  // namespace SharedPrim_dataformatter_impl

template <>
DataFormatterDictionary<std::shared_ptr<const Structure>>
make_attribute_dictionary<std::shared_ptr<const Structure>>() {
  using namespace SharedPrim_dataformatter_impl;
  DataFormatterDictionary<SharedPrim> dict;
  dict.insert(lattice_point_group(), lattice_point_group_name(),
              lattice_point_group_size(), factor_group(), factor_group_name(),
              factor_group_size(), crystal_point_group(),
              crystal_point_group_name(), crystal_point_group_size(), volume(),
              lattice(), lattice_column_matrix(), lattice_params(),
              asymmetric_unit());

  return dict;
}

}  // namespace CASM
