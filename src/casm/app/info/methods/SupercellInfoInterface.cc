#include "casm/app/info/methods/SupercellInfoInterface.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DatumFormatterAdapter.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/clex/io/stream/Configuration_stream_io.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/LatticeEnumEquivalents.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/io/data/SupercellSymInfo_data_io.hh"
#include "casm/symmetry/io/json/SymGroup_json_io.hh"

namespace CASM {

// declarations:
namespace {

std::vector<Lattice> parse_params_supercells(
    jsonParser const &params_json,
    std::shared_ptr<Structure const> shared_prim);

jsonParser &superlattice_info_to_json(
    const xtal::Superlattice &lat, jsonParser &json,
    std::shared_ptr<Structure const> shared_prim);

}  // namespace

namespace {

std::shared_ptr<Structure const> open_shared_prim(fs::path root) {
  ProjectSettings settings = open_project_settings(root);
  return std::make_shared<Structure const>(
      read_prim(settings.dir().prim(), settings.crystallography_tol()));
}

/// Data structure for formatting properties of prim and supercell_sym_info
struct SupercellInfoData {
  SupercellInfoData(std::shared_ptr<Structure const> _shared_prim,
                    SupercellSymInfo const &_supercell_sym_info,
                    jsonParser const &_params_json)
      : shared_prim(_shared_prim),
        supercell_sym_info(_supercell_sym_info),
        params_json(_params_json) {}

  std::shared_ptr<Structure const> shared_prim;

  SupercellSymInfo const &supercell_sym_info;

  jsonParser const &params_json;
};

// use for ValueType=bool, int, double, std::string, jsonParser
template <typename ValueType>
using SupercellInfoFormatter =
    GenericDatumFormatter<ValueType, SupercellInfoData>;

typedef Generic1DDatumFormatter<Eigen::VectorXd, SupercellInfoData>
    VectorXdSupercellInfoFormatter;

typedef Generic1DDatumFormatter<Eigen::VectorXi, SupercellInfoData>
    VectorXiSupercellInfoFormatter;

typedef Generic2DDatumFormatter<Eigen::MatrixXd, SupercellInfoData>
    MatrixXdSupercellInfoFormatter;

typedef Generic2DDatumFormatter<Eigen::MatrixXi, SupercellInfoData>
    MatrixXiSupercellInfoFormatter;

SupercellInfoFormatter<std::string> supercell_name() {
  return SupercellInfoFormatter<std::string>(
      "supercell_name",
      "A name given to all equivalent super lattices of the prim lattice. See "
      "the input option `supercell_name` for a more detailed description.",
      [](SupercellInfoData const &data) -> std::string {
        return make_supercell_name(data.shared_prim->point_group(),
                                   data.supercell_sym_info.prim_lattice(),
                                   data.supercell_sym_info.supercell_lattice());
      });
}

SupercellInfoFormatter<std::string> canonical_supercell_name() {
  return SupercellInfoFormatter<std::string>(
      "canonical_supercell_name",
      "Name of the canonical equivalent supercell. See the input option "
      "`supercell_name` for a more detailed description.",
      [](SupercellInfoData const &data) -> std::string {
        return make_canonical_supercell_name(
            data.shared_prim->point_group(),
            data.supercell_sym_info.prim_lattice(),
            data.supercell_sym_info.supercell_lattice());
      });
}

SupercellInfoFormatter<bool> is_canonical() {
  return SupercellInfoFormatter<bool>(
      "is_canonical", "Returns true if supercell lattice is in canonical form.",
      [](SupercellInfoData const &data) -> bool {
        return xtal::canonical::check(
            data.supercell_sym_info.supercell_lattice(),
            data.shared_prim->point_group());
      });
}

SupercellInfoFormatter<Index> multiplicity() {
  return SupercellInfoFormatter<Index>(
      "multiplicity", "Number of equivalent supercells",
      [](SupercellInfoData const &data) -> Index {
        return data.shared_prim->factor_group().size() /
               data.supercell_sym_info.factor_group().size();
      });
}

SupercellInfoFormatter<jsonParser> orbit() {
  return SupercellInfoFormatter<jsonParser>(
      "orbit",
      "Output an array of objects containings equivalent supercell lattices, "
      "specified as lattice row vector matrices.",
      [](SupercellInfoData const &data) -> jsonParser {
        jsonParser orbit = jsonParser::array();

        xtal::Lattice canonical_lattice = xtal::canonical::equivalent(
            data.supercell_sym_info.supercell_lattice(),
            data.shared_prim->point_group());

        LatticeEnumEquivalents enumerator(canonical_lattice,
                                          data.shared_prim->factor_group());

        for (xtal::Lattice const &supercell_lattice : enumerator) {
          xtal::Superlattice superlattice{data.shared_prim->lattice(),
                                          supercell_lattice};
          jsonParser tjson;
          to_json(supercell_lattice.lat_column_mat(),
                  tjson["supercell_lattice_column_matrix"]);
          to_json(superlattice.transformation_matrix_to_super(),
                  tjson["transformation_matrix_to_super"]);
          orbit.push_back(tjson);
        }
        return orbit;
      });
}

SupercellInfoFormatter<jsonParser> frac_coordinate() {
  return SupercellInfoFormatter<jsonParser>(
      "frac_coordinate",
      "Fractional coordinate (with respect to the supercell lattice vectors) "
      "of every site in the supercell. For coordinate conversions, the order "
      "in which sites are listed is the same as for `cart_coordinate` and "
      "`integral_site_coordinates`. The site order is determined by "
      "`linear_site_index = linear_unitcell_index + supercell_volume * "
      "sublattice_index`, where `linear_unitcell_index` is an index into the "
      "`unitcells` list.",
      [](SupercellInfoData const &data) -> jsonParser {
        auto const &prim = *data.shared_prim;
        jsonParser json = jsonParser::array();
        auto f = data.supercell_sym_info.unitcellcoord_index_converter();
        for (Index l = 0; l < f.total_sites(); l++) {
          jsonParser l_json;
          xtal::Coordinate coord = f(l).coordinate(prim);
          coord.set_lattice(data.supercell_sym_info.supercell_lattice(), CART);
          to_json_array(coord.const_frac(), l_json);
          json.push_back(l_json);
        }
        return json;
      });
}

SupercellInfoFormatter<jsonParser> cart_coordinate() {
  return SupercellInfoFormatter<jsonParser>(
      "cart_coordinate",
      "Cartesian coordinate of every site in the supercell. For coordinate "
      "conversions, the order in which sites are listed is the same as for "
      "`frac_coordinate` and `integral_site_coordinates`. The site order is "
      "determined by `linear_site_index = linear_unitcell_index + "
      "supercell_volume * sublattice_index`, where `linear_unitcell_index` is "
      "an index into the `unitcells` list.",
      [](SupercellInfoData const &data) -> jsonParser {
        auto const &prim = *data.shared_prim;
        jsonParser json = jsonParser::array();
        auto f = data.supercell_sym_info.unitcellcoord_index_converter();
        for (Index l = 0; l < f.total_sites(); l++) {
          jsonParser l_json;
          to_json_array(f(l).coordinate(prim).const_cart(), l_json);
          json.push_back(l_json);
        }
        return json;
      });
}

SupercellInfoFormatter<jsonParser> is_supercell_of() {
  return SupercellInfoFormatter<jsonParser>(
      "is_supercell_of",
      "For each supercell lattice, L, in the `params/supercells` array, "
      "determine if the input supercell lattice, S, is a supercell of L. "
      "Returns an array of JSON objects with the following results for each "
      "input supercell lattice: `value` (bool), and the values of `S`, `L` and "
      "`T` in `S = L * T`, where S and L are column vector matrices and T is "
      "an integer transformation matrix if S is a supercell of L and T is a "
      "floating point matrix if S is not a supercell of L.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::vector<xtal::Lattice> supercells =
            parse_params_supercells(data.params_json, data.shared_prim);
        double tol = TOL;
        xtal::Lattice const &S = data.supercell_sym_info.supercell_lattice();
        jsonParser result = jsonParser::array();
        for (xtal::Lattice const &L : supercells) {
          auto pair = xtal::is_superlattice(S, L, tol);
          jsonParser json;
          json["value"] = pair.first;
          json["S"] = S.lat_column_mat();
          json["L"] = L.lat_column_mat();
          if (pair.first) {
            json["T"] = lround(pair.second);
          } else {
            json["T"] = pair.second;
          }
          result.push_back(json);
        }
        return result;
      });
}

SupercellInfoFormatter<jsonParser> is_unitcell_of() {
  return SupercellInfoFormatter<jsonParser>(
      "is_unitcell_of",
      "For each supercell lattice, L, in the `params/supercells` array, "
      "determine if the input supercell lattice, S, is a unit cell of L. "
      "Returns an array of JSON objects with the following results for each "
      "input supercell lattice: `value` (bool), and the values of `S`, `L` "
      "and `T` in `L = S * T`, where S and L are column vector matrices and T "
      "is an integer transformation matrix if S is a unit cell of L and T is a "
      "floating point matrix if S is not a unit cell of L.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::vector<xtal::Lattice> supercells =
            parse_params_supercells(data.params_json, data.shared_prim);
        double tol = TOL;
        xtal::Lattice const &S = data.supercell_sym_info.supercell_lattice();
        jsonParser result = jsonParser::array();
        for (xtal::Lattice const &L : supercells) {
          auto pair = xtal::is_superlattice(L, S, tol);
          jsonParser json;
          json["value"] = pair.first;
          json["S"] = S.lat_column_mat();
          json["L"] = L.lat_column_mat();
          if (pair.first) {
            json["T"] = lround(pair.second);
          } else {
            json["T"] = pair.second;
          }
          result.push_back(json);
        }
        return result;
      });
}

SupercellInfoFormatter<jsonParser> is_equivalent_to() {
  return SupercellInfoFormatter<jsonParser>(
      "is_equivalent_to",
      "For each supercell lattice, L, in the `params/supercells` array, "
      "determine if the input supercell lattice, S, is equivalent to L. "
      "Returns an array of JSON objects with the following results for each "
      "input supercell lattice: `value` (bool), and the values of `S`, `L` "
      "and `U` in `S = L * U`, where S and L are column vector matrices, and U "
      "is a unimodular matrix (integer transformation matrix with "
      "determinant==1) if S is equivalent to L and U is a floating point "
      "matrix if S is not equivalent to L.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::vector<Lattice> supercells =
            parse_params_supercells(data.params_json, data.shared_prim);
        xtal::Lattice const &S = data.supercell_sym_info.supercell_lattice();
        xtal::LatticeIsEquivalent f{S};
        jsonParser result = jsonParser::array();
        for (xtal::Lattice const &L : supercells) {
          bool value = f(L);
          jsonParser json;
          json["value"] = value;
          json["S"] = S.lat_column_mat();
          json["L"] = L.lat_column_mat();
          if (value) {
            json["U"] = lround(f.U());
          } else {
            json["U"] = f.U();
          }
          result.push_back(json);
        }
        return result;
      });
}

SupercellInfoFormatter<jsonParser> is_equivalent_to_a_supercell_of() {
  return SupercellInfoFormatter<jsonParser>(
      "is_equivalent_to_a_supercell_of",
      "For each supercell lattice, L, in the `params/supercells` array, "
      "determine if the input supercell lattice, S, is equivalent to a "
      "supercell of L. Returns an array of JSON objects with the following "
      "results for each input supercell lattice: `value` (bool), and the "
      "values of `S`, `R`, `L` and `T` in `S = R * L * T`, where S and L are "
      "column vector matrices, R is a Cartesian symmetry matrix of a prim "
      "point group operation, and T is an integer transformation matrix if S "
      "is equivalent to a supercell of L and T is a floating point matrix if S "
      "is not equivalent to a supercell of L.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::vector<xtal::Lattice> supercells =
            parse_params_supercells(data.params_json, data.shared_prim);
        double tol = TOL;
        xtal::Lattice const &S = data.supercell_sym_info.supercell_lattice();
        auto begin = data.shared_prim->point_group().begin();
        auto end = data.shared_prim->point_group().end();

        jsonParser result = jsonParser::array();
        for (xtal::Lattice const &L : supercells) {
          auto pair = xtal::is_equivalent_superlattice(S, L, begin, end, tol);
          bool value = (pair.first != end);
          jsonParser json;
          json["value"] = value;
          json["S"] = S.lat_column_mat();
          json["L"] = L.lat_column_mat();
          if (value) {
            json["R"] = pair.first->matrix();
            json["T"] = lround(pair.second);
          } else {
            json["R"] = Eigen::Matrix3d::Identity();
            json["T"] = is_superlattice(S, L, tol).second;
          }
          result.push_back(json);
        }
        return result;
      });
}

SupercellInfoFormatter<jsonParser> commensurate_superlattice_of() {
  return SupercellInfoFormatter<jsonParser>(
      "commensurate_superlattice_of",
      "Calculate a superlattice, S, which is a supercell of all supercell "
      "lattices, L(i), in the `params/supercells` array. Returns a JSON object "
      "describing S such that `S = L_i*T_i` for all i, where S and L_i are "
      "column vector matrices, and T_i are integer transformation matrices.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::vector<Lattice> supercells =
            parse_params_supercells(data.params_json, data.shared_prim);
        Lattice S = make_commensurate_superduperlattice(supercells.begin(),
                                                        supercells.end());
        // uses point group of S, so remains equivalent to original S,
        // and it may not be canonical w.r.t prim point group
        S = xtal::canonical::equivalent(S);
        xtal::Superlattice superlattice{data.shared_prim->lattice(), S};
        jsonParser json;
        superlattice_info_to_json(superlattice, json, data.shared_prim);
        return json;
      });
}

SupercellInfoFormatter<jsonParser> minimal_commensurate_superlattice_of() {
  return SupercellInfoFormatter<jsonParser>(
      "minimal_commensurate_superlattice_of",
      "Calculate a superlattice, S, which the minimum volume superlattice "
      "that is equivalent to a supercell of all supercell lattices, L(i), in "
      "the `params/supercells` array. Returns a JSON object describing the "
      "minimum volume superlattice, S, such that `S = R_i*L_i*T_i` for all i, "
      "where S and L_i are column vector matrices, R_i is a Cartesian symmetry "
      "matrix of a prim point group operation, and T_i are integer "
      "transformation matrices.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::vector<Lattice> supercells =
            parse_params_supercells(data.params_json, data.shared_prim);
        auto begin = data.shared_prim->point_group().begin();
        auto end = data.shared_prim->point_group().end();
        Lattice S = make_minimal_commensurate_superduperlattice(
            supercells.begin(), supercells.end(), begin, end);

        // make canonical w.r.t prim point group
        S = xtal::canonical::equivalent(S, data.shared_prim->point_group());
        xtal::Superlattice superlattice{data.shared_prim->lattice(), S};
        jsonParser json;
        superlattice_info_to_json(superlattice, json, data.shared_prim);
        return json;
      });
}

SupercellInfoFormatter<jsonParser> fully_commensurate_superlattice_of() {
  return SupercellInfoFormatter<jsonParser>(
      "fully_commensurate_superlattice_of",
      "Calculate a superlattice, S, which is a supercell of all equivalents "
      "of all input supercell lattices, L(i), in the `params/supercells` "
      "array. Returns a JSON object describing the superlattice S such that "
      "`S = R_j*L_i*T_{i,j}` for all (i, j), where S and L_i are column vector "
      "matrices, R_j is a Cartesian symmetry matrix of a prim point group "
      "operation, and T_{i,j} are integer transformation matrices.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::vector<Lattice> supercells =
            parse_params_supercells(data.params_json, data.shared_prim);
        auto begin = data.shared_prim->point_group().begin();
        auto end = data.shared_prim->point_group().end();
        Lattice S = make_fully_commensurate_superduperlattice(
            supercells.begin(), supercells.end(), begin, end);

        // make canonical w.r.t prim point group
        S = xtal::canonical::equivalent(S, data.shared_prim->point_group());
        xtal::Superlattice superlattice{data.shared_prim->lattice(), S};
        jsonParser json;
        superlattice_info_to_json(superlattice, json, data.shared_prim);
        return json;
      });
}

SupercellInfoFormatter<jsonParser> enforce_minimum_size() {
  return SupercellInfoFormatter<jsonParser>(
      "enforce_minimum_size",
      "Calculate the minimum size superlattice, S, of the prim lattice, P, "
      "which is a supercell of the input supercell lattice, S_init, and "
      "greater or equal to some specified minimum size (calculated as "
      "integer multiples of the prim cell). Returns a JSON object describing "
      "the superlattice S, such that `S = P * T = S_init * M`, and the "
      "supercell_size of S is greater than or equal to `params/minimum_size` "
      "(integer, default == 1), where S, S_init, and P are column vector "
      "matrices, and T and M are integer transformation matrices. If "
      "`params/fixed_shape` (bool) is `true` (default == false), then M is "
      "fixed to be a multiple of the identity matrix, `M = m * I`, where m is "
      "an integer, and I is the identity matrix.",
      [](SupercellInfoData const &data) -> jsonParser {
        Log &log = CASM::log();
        ParentInputParser parser{data.params_json};
        std::runtime_error error_if_invalid{
            "Error reading `enforce_minimum_size` input"};

        Index minimum_size = 1;
        parser.optional(minimum_size, "minimum_size");

        bool fixed_shape = false;
        parser.optional(fixed_shape, "fixed_shape");

        report_and_throw_if_invalid(parser, log, error_if_invalid);

        auto begin = data.shared_prim->point_group().begin();
        auto end = data.shared_prim->point_group().end();
        Eigen::Matrix3l T_init =
            data.supercell_sym_info.transformation_matrix_to_super();
        Eigen::Matrix3l M =
            enforce_min_volume(begin, end, data.shared_prim->lattice(),
                               T_init.cast<int>(), minimum_size, fixed_shape)
                .cast<long>();

        Eigen::Matrix3l T = T_init * M;
        xtal::Lattice S = make_superlattice(data.shared_prim->lattice(), T);

        if (!fixed_shape) {
          // uses point group of S, so remains equivalent to original S,
          // and it may not be canonical w.r.t prim point group
          S = xtal::canonical::equivalent(S);
        }

        xtal::Superlattice superlattice{data.shared_prim->lattice(), S};
        jsonParser json;
        superlattice_info_to_json(superlattice, json, data.shared_prim);
        return json;
      });
}

SupercellInfoFormatter<jsonParser> default_configuration() {
  return SupercellInfoFormatter<jsonParser>(
      "default_configuration",
      "All degrees of freedom (DoF) of the default configuration in the "
      "supercell, formatted as JSON.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::shared_ptr<Supercell const> shared_supercell =
            std::make_shared<Supercell const>(
                data.shared_prim, data.supercell_sym_info.supercell_lattice());
        jsonParser json = jsonParser::object();
        Configuration configuration{shared_supercell};
        to_json(configuration, json);
        return json;
      });
}

SupercellInfoFormatter<jsonParser> default_structure() {
  return SupercellInfoFormatter<jsonParser>(
      "default_structure",
      "Structure resulting from the default configuration in the supercell, "
      "formatted as JSON.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::shared_ptr<Supercell const> shared_supercell =
            std::make_shared<Supercell const>(
                data.shared_prim, data.supercell_sym_info.supercell_lattice());
        jsonParser json = jsonParser::object();
        Configuration configuration{shared_supercell};
        to_json(make_simple_structure(configuration), json);
        return json;
      });
}

SupercellInfoFormatter<jsonParser> default_structure_with_vacancies() {
  return SupercellInfoFormatter<jsonParser>(
      "default_structure_with_vacancies",
      "Structure resulting from the default configuration in the supercell, "
      "including vacancies, formatted as JSON.",
      [](SupercellInfoData const &data) -> jsonParser {
        std::shared_ptr<Supercell const> shared_supercell =
            std::make_shared<Supercell const>(
                data.shared_prim, data.supercell_sym_info.supercell_lattice());
        jsonParser json = jsonParser::object();
        Configuration configuration{shared_supercell};
        std::set<std::string> const &excluded_species = {};
        to_json(make_simple_structure(configuration), json, excluded_species);
        return json;
      });
}

SupercellInfoFormatter<std::string> default_poscar() {
  return SupercellInfoFormatter<std::string>(
      "default_poscar",
      "Structure resulting from the default configuration in the supercell, "
      "formatted as a VASP POSCAR.",
      [](SupercellInfoData const &data) -> std::string {
        std::shared_ptr<Supercell const> shared_supercell =
            std::make_shared<Supercell const>(
                data.shared_prim, data.supercell_sym_info.supercell_lattice());
        std::string supercell_name =
            make_supercell_name(data.shared_prim->point_group(),
                                data.supercell_sym_info.prim_lattice(),
                                data.supercell_sym_info.supercell_lattice());
        Configuration configuration{shared_supercell};
        std::stringstream ss;
        VaspIO::PrintPOSCAR p{make_simple_structure(configuration),
                              supercell_name};
        // Va are ignored by default
        p.sort();
        p.print(ss);
        return ss.str();
      });
}

SupercellInfoFormatter<std::string> default_poscar_with_vacancies() {
  return SupercellInfoFormatter<std::string>(
      "default_poscar_with_vacancies",
      "Structure resulting from the default configuration in the supercell, "
      "formatted as a VASP POSCAR.",
      [](SupercellInfoData const &data) -> std::string {
        std::shared_ptr<Supercell const> shared_supercell =
            std::make_shared<Supercell const>(
                data.shared_prim, data.supercell_sym_info.supercell_lattice());
        std::string supercell_name =
            make_supercell_name(data.shared_prim->point_group(),
                                data.supercell_sym_info.prim_lattice(),
                                data.supercell_sym_info.supercell_lattice());
        Configuration configuration{shared_supercell};
        std::stringstream ss;
        VaspIO::PrintPOSCAR p{make_simple_structure(configuration),
                              supercell_name};
        // Va are ignored by default
        p.ignore() = {};  // do not ignore vacancies
        p.sort();
        p.print(ss);
        return ss.str();
      });
}

SupercellInfoFormatter<jsonParser> fill_supercell() {
  return SupercellInfoFormatter<jsonParser>(
      "fill_supercell",
      "Return a configuration constructed by tiling a sub-configuration into "
      "the specified supercell. Requires the following params: "
      "\"sub_configuration\", (object) the JSON representation of a "
      "sub-configuration. The output is the JSON representation of a "
      "configuration.",
      [](SupercellInfoData const &data) -> jsonParser {
        ParentInputParser parser{data.params_json};
        std::unique_ptr<Configuration> sub_configuration_ptr =
            parser.require<Configuration>("sub_configuration",
                                          data.shared_prim);
        std::runtime_error error_if_invalid{
            "Error reading \"fill_supercell\" input"};
        report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

        std::shared_ptr<Supercell const> shared_supercell =
            std::make_shared<Supercell const>(
                data.shared_prim, data.supercell_sym_info.supercell_lattice());

        jsonParser json = jsonParser::object();
        Configuration configuration =
            fill_supercell(*sub_configuration_ptr, shared_supercell);
        to_json(configuration, json);
        return json;
      });
}

SupercellInfoFormatter<jsonParser> configuration_from_normal_coordinate() {
  return SupercellInfoFormatter<jsonParser>(
      "configuration_from_normal_coordinate",
      "Return a configuration with DoF values corresponding to a normal "
      "coordinate in a DoF space basis, filled into the specified supercell. "
      "This method may be used for global DoF or continuous local DoF, not occ "
      "DoF. A symmetry operation may be applied if necessary to fill the "
      "requested supercell. Requires the following params: \"dof_space\", (see "
      "\"normal_coordinate\" property in `casm info --desc "
      "ConfigurationInfo`); \"normal_coordinate\" (vector) a normal coordinate "
      "of size equal to the DoF space basis dimension. May include params: "
      "\"default_configuration\" (object) JSON representation of a "
      "configuration in the DoF space  to use as the default configuration "
      "before applying the normal coordinates, this can be used to set other "
      "types of DoF. The output is: \"configuration\" the JSON "
      "representation of a configuration; \"symop\", the symmetry operation "
      "used to fill the supercell. If the requested supercell can not be tiled "
      "by the DoF space, even upon symmetry application, then a \"Could not "
      "tile supercell\" error message is returned.",
      [](SupercellInfoData const &data) -> jsonParser {
        ParentInputParser parser{data.params_json};
        std::unique_ptr<DoFSpace> dof_space_ptr =
            parser.require<DoFSpace>("dof_space", data.shared_prim);
        Eigen::VectorXd normal_coordinate;
        parser.require(normal_coordinate, "normal_coordinate");
        std::unique_ptr<Configuration> default_configuration =
            parser.optional<Configuration>("default_configuration",
                                           data.shared_prim);
        std::runtime_error error_if_invalid{
            "Error reading \"configuration_from_normal_coordinate\" input"};
        report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

        DoFSpace const &dof_space = *dof_space_ptr;
        Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
        if (dof_space.transformation_matrix_to_super().has_value()) {
          T = *dof_space.transformation_matrix_to_super();
        }

        auto shared_supercell =
            std::make_shared<Supercell const>(data.shared_prim, T);
        Configuration init_config(shared_supercell);

        if (default_configuration != nullptr) {
          if (default_configuration->supercell()
                  .sym_info()
                  .transformation_matrix_to_super() != T) {
            jsonParser json;
            json =
                "`default_configuration` supercell does not match `dof_space` "
                "supercell";
            return json;
          }
          init_config = *default_configuration;
        }

        set_dof_value(init_config, dof_space, normal_coordinate);

        std::shared_ptr<Supercell const> final_shared_supercell =
            std::make_shared<Supercell const>(
                data.shared_prim, data.supercell_sym_info.supercell_lattice());

        FillSupercell f(final_shared_supercell);
        SymOp const *symop_ptr =
            f.find_symop(init_config.supercell().lattice());
        if (symop_ptr == nullptr) {
          jsonParser json;
          json = "Could not tile supercell";
          return json;
        }
        Configuration final_config = f(init_config);

        jsonParser json;
        to_json(final_config, json["configuration"]);
        write_symop(data.shared_prim->factor_group(),
                    symop_ptr->master_group_index(), json["symop"]);
        return json;
      });
}

}  // namespace

namespace adapter {

template <typename ToType, typename FromType>
struct Adapter;

template <>
struct Adapter<SupercellSymInfo, SupercellInfoData> {
  SupercellSymInfo const &operator()(SupercellInfoData const &adaptable) const {
    return adaptable.supercell_sym_info;
  }
};

}  // namespace adapter

namespace {

DataFormatterDictionary<SupercellInfoData> make_supercell_info_dict() {
  DataFormatterDictionary<SupercellInfoData> supercell_info_dict;

  // properties that require prim and supercell_sym_info
  supercell_info_dict.insert(
      supercell_name(), canonical_supercell_name(), is_canonical(),
      multiplicity(), orbit(), frac_coordinate(), cart_coordinate(),
      is_supercell_of(), is_unitcell_of(), is_equivalent_to(),
      is_equivalent_to_a_supercell_of(), commensurate_superlattice_of(),
      minimal_commensurate_superlattice_of(),
      fully_commensurate_superlattice_of(), enforce_minimum_size(),
      default_configuration(), default_structure(),
      default_structure_with_vacancies(), default_poscar(),
      default_poscar_with_vacancies(), fill_supercell(),
      configuration_from_normal_coordinate());

  // properties that only require supercell_sym_info
  auto sym_info_dict = make_dictionary<SupercellSymInfo>();
  for (auto it = sym_info_dict.begin(); it != sym_info_dict.end(); ++it) {
    if (it->type() == DatumFormatterClass::Property) {
      supercell_info_dict.insert(
          make_datum_formatter_adapter<SupercellInfoData, SupercellSymInfo>(
              *it));
    }
  }

  return supercell_info_dict;
}

/// Local struct allows for InputParser of Lattice using a local specification
struct SupercellInfoLatticeIO {
  SupercellInfoLatticeIO(Lattice const &_lattice) : lattice(_lattice) {}
  Lattice lattice;
};

/// Parse SupercellInfoLatticeIO (essentially Lattice) from a JSON object). See
/// the description of allowed options in `SupercellInfoInterface::desc`.
void parse(InputParser<SupercellInfoLatticeIO> &parser,
           std::shared_ptr<Structure const> shared_prim) {
  Log &log = CASM::log();

  // read "transformation_matrix_to_super"
  Eigen::Matrix3l T;
  if (parser.self.contains("transformation_matrix_to_super")) {
    parser.optional(T, "transformation_matrix_to_super");

    // or read "supercell_lattice_row_vectors"
  } else if (parser.self.contains("supercell_lattice_row_vectors")) {
    Eigen::Matrix3d L_transpose;
    parser.optional(L_transpose, "supercell_lattice_row_vectors");
    Lattice super_lattice{L_transpose.transpose()};
    try {
      T = make_transformation_matrix_to_super(shared_prim->lattice(),
                                              super_lattice, TOL);
    } catch (std::exception &e) {
      auto pair = is_superlattice(super_lattice, shared_prim->lattice(), TOL);
      log << "The transformation_matrix_to_super determined from "
             "\"supercell_lattice_row_vectors\" is not approximately integer. "
             "Found:\n"
          << pair.second << std::endl;
      parser.insert_error("supercell_lattice_row_vectors", e.what());
    }

    // or read "supercell_lattice_column_matrix"
  } else if (parser.self.contains("supercell_lattice_column_matrix")) {
    Eigen::Matrix3d L;
    parser.optional(L, "supercell_lattice_column_matrix");
    Lattice super_lattice{L};
    try {
      T = make_transformation_matrix_to_super(shared_prim->lattice(),
                                              super_lattice, TOL);
    } catch (std::exception &e) {
      auto pair = is_superlattice(super_lattice, shared_prim->lattice(), TOL);
      log << "The transformation_matrix_to_super determined from "
             "\"supercell_lattice_column_matrix\" is not approximately "
             "integer. Found:\n"
          << pair.second << std::endl;
      parser.insert_error("supercell_lattice_column_matrix", e.what());
    }
    // or read "supercell_name"
  } else if (parser.self.contains("supercell_name")) {
    std::string supercell_name;
    parser.optional(supercell_name, "supercell_name");
    try {
      xtal::Superlattice superlattice = make_superlattice_from_supercell_name(
          shared_prim->factor_group(), shared_prim->lattice(), supercell_name);
      T = superlattice.transformation_matrix_to_super();
    } catch (std::exception &e) {
      parser.insert_error("supercell_name", e.what());
    }
    // else use Identity (prim cell)
  } else {
    T = Eigen::Matrix3l::Identity();
  }

  // read "make_canonical"
  bool make_canonical = false;
  parser.optional(make_canonical, "make_canonical");

  if (!parser.valid()) {
    return;
  }

  // construct SupercellSymInfo
  Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
  if (make_canonical) {
    super_lattice =
        xtal::canonical::equivalent(super_lattice, shared_prim->point_group());
  }

  parser.value = notstd::make_unique<SupercellInfoLatticeIO>(super_lattice);
}

/// Parse "params"/"supercells" (JSON array) as vector<Lattice>, using the
/// SupercellInfoLatticeIO parser for each element of the array. See
/// `SupercellInfoInterface::desc`.
std::vector<Lattice> parse_params_supercells(
    jsonParser const &params_json,
    std::shared_ptr<Structure const> shared_prim) {
  Log &log = CASM::log();
  ParentInputParser parser{params_json};
  std::runtime_error error_if_invalid{"Error reading params/supercells input"};

  if (!parser.self.contains("supercells")) {
    parser.error.insert("Required params/supercells does not exist");
    report_and_throw_if_invalid(parser, log, error_if_invalid);
  }

  if (!parser.self["supercells"].is_array()) {
    parser.error.insert("Input params/supercells must be an array.");
    report_and_throw_if_invalid(parser, log, error_if_invalid);
  }

  std::vector<Lattice> supercells;
  for (int i = 0; i < parser.self["supercells"].size(); ++i) {
    fs::path option{"supercells"};
    option /= std::to_string(i);
    auto shared_lattice_parser =
        parser.subparse<SupercellInfoLatticeIO>(option, shared_prim);
    if (!parser.valid()) {
      report_and_throw_if_invalid(parser, log, error_if_invalid);
    }
    supercells.push_back(shared_lattice_parser->value->lattice);
  }
  return supercells;
}

jsonParser &superlattice_info_to_json(
    const xtal::Superlattice &superlattice, jsonParser &json,
    std::shared_ptr<Structure const> shared_prim) {
  Lattice const &S = superlattice.superlattice();
  Eigen::Matrix3l const &T = superlattice.transformation_matrix_to_super();

  std::string supercell_name = make_supercell_name(shared_prim->point_group(),
                                                   superlattice.prim_lattice(),
                                                   superlattice.superlattice());

  std::string canonical_supercell_name = make_canonical_supercell_name(
      shared_prim->point_group(), superlattice.prim_lattice(),
      superlattice.superlattice());

  Eigen::VectorXd supercell_lattice_params(6);
  supercell_lattice_params << S.length(0), S.length(1), S.length(2), S.angle(0),
      S.angle(1), S.angle(2);

  json["transformation_matrix_to_super"] = T;
  json["supercell_lattice_row_vectors"] = S.lat_column_mat().transpose();
  json["supercell_lattice_column_matrix"] = S.lat_column_mat();
  json["supercell_name"] = supercell_name;
  json["canonical_supercell_name"] = canonical_supercell_name;
  json["supercell_lattice_params"] = supercell_lattice_params;
  json["supercell_size"] = superlattice.size();
  return json;
}

}  // namespace

std::string SupercellInfoInterface::desc() const {
  std::string description =
      "Get information about a supercell. The supercell is specified by   \n"
      "the prim and one of the following (else the primitive cell is      \n"
      "used):                                                             \n"
      "- transformation_matrix_to_super                                   \n"
      "- supercell_lattice_row_vectors                                    \n"
      "- supercell_lattice_column_matrix                                  \n"
      "- supercell_name                                                   \n\n"

      "Additionally, if the `make_canonical` option can be used to specify\n"
      "that the canonical equivalent supercell should be used.            \n\n"

      "From all super lattices which are equivalent via crystal point     \n"
      "group operations of prim, the canonical super lattice is the       \n"
      "equivalent lattice in niggli form which has the most favorable     \n"
      "orientation according to the CASM orientation comparision criteria.\n"
      "These criteria compare the lattices as column vector matrices and  \n"
      "first favor symmetric matrices and then favor non-negative lattice \n"
      "vectors which are aligned nearest to the Cartesian axes.           \n\n";

  std::string custom_options =
      "  prim: JSON object (optional, default=prim of current project)    \n"
      "    See `casm format --prim` for details on the prim format.       \n\n"

      "  transformation_matrix_to_super: 3x3 array of integer (optional)  \n"
      "    Transformation matrix T, defining the supercell lattice vectors\n"
      "    S, in terms of the prim lattice vectors, P: `S = P * T`, where \n"
      "    S and P are column vector matrices.                            \n\n"

      "  supercell_lattice_row_vectors: 3x3 array of integer (optional)   \n"
      "    Supercell lattice vectors, as a row vector matrix.             \n\n"

      "  supercell_lattice_column_matrix: 3x3 array of integer (optional) \n"
      "    Supercell lattice vectors, as a column vector matrix.          \n\n"

      "  supercell_name: string (optional)                                \n"
      "    A name given to all equivalent super lattices of the prim      \n"
      "    lattice. For the canonical super lattice, the name is          \n"
      "    constructed from the hermite normal form of                    \n"
      "    `transformation_matrix_to_super`. For a non-canonical super    \n"
      "    lattice, the name is the constructed from the name of the      \n"
      "    canonical super lattice and the index of the prim factor group \n"
      "    operation that (excluding the shift) transforms the canonical  \n"
      "    super lattice into this super lattice.                         \n"
      "                                                                   \n"
      "    Example 1: Canonical supercell name                            \n"
      "                                                                   \n"
      "        \"SCEL8_4_2_1_1_3_2\"                                      \n"
      "                                                                   \n"
      "    Example 2: Non-canonical supercell name representing a         \n"
      "    re-orientation by application of prim factor group operation   \n"
      "    with index 2 (indexing starting at 0).                         \n"
      "                                                                   \n"
      "        \"SCEL8_4_2_1_1_3_2.2\"                                    \n\n"

      "  make_canonical: bool (optional, default=false)                   \n"
      "    If \"true\", the canonical equivalent supercell is used.       \n\n"

      "  params: JSON object (optional, default={})                       \n"
      "    JSON object containing additional parameters for methods that  \n"
      "    require them. See method descriptions for details of how they  \n"
      "    are used or other specific parameters not documented here.     \n"
      "    Options include:                                               \n\n"

      "      supercells: array of JSON object                             \n"
      "        Each JSON object specifies one lattice, using the same     \n"
      "        options allowed for specifying the supercell:              \n"
      "        - transformation_matrix_to_super                           \n"
      "        - supercell_lattice_row_vectors                            \n"
      "        - supercell_lattice_column_matrix                          \n"
      "        - supercell_name                                           \n"
      "        - make_canonical                                           \n\n"

      "  properties: array of string                                      \n"
      "    An array of strings specifying which supercell properties to   \n"
      "    output. The allowed options are:                               \n\n";

  std::stringstream ss;
  auto dict = make_supercell_info_dict();
  dict.print_help(ss, DatumFormatterClass::Property);

  return name() + ": \n\n" + description + custom_options + ss.str();
}

std::string SupercellInfoInterface::name() const { return "SupercellInfo"; }

/// Run `prim` info method
void SupercellInfoInterface::run(jsonParser const &json_options,
                                 PrimClex const *primclex,
                                 fs::path root) const {
  Log &log = CASM::log();

  ParentInputParser parser{json_options};
  std::runtime_error error_if_invalid{"Error reading SupercellInfo input"};

  // read "prim"
  std::shared_ptr<Structure const> shared_prim;
  if (parser.self.contains("prim")) {
    // prim provided in input
    xtal::BasicStructure basic_structure;
    parser.optional<xtal::BasicStructure>(basic_structure, "prim", TOL);
    if (parser.valid()) {
      shared_prim = std::make_shared<Structure const>(basic_structure);
    }
  } else if (primclex != nullptr) {
    // if project provided via api
    shared_prim = primclex->shared_prim();
  } else {
    // if project contains current working directory
    if (!root.empty()) {
      try {
        shared_prim = open_shared_prim(root);
      } catch (std::exception &e) {
        parser.insert_error("prim", e.what());
      }
    } else {
      std::stringstream msg;
      msg << "Error in SupercellInfo: No \"prim\" in input and no project "
             "provided "
             "or found.";
      parser.insert_error("prim", msg.str());
    }
  }
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  auto shared_lattice_parser =
      parser.parse_as<SupercellInfoLatticeIO>(shared_prim);

  // read "properties"
  std::vector<std::string> properties;
  parser.require(properties, "properties");

  // read "params"
  jsonParser params_json = jsonParser::object();
  parser.optional(params_json, "params");
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  SupercellSymInfo supercell_sym_info = make_supercell_sym_info(
      *shared_prim, shared_lattice_parser->value->lattice);

  DataFormatterDictionary<SupercellInfoData> supercell_info_dict =
      make_supercell_info_dict();

  // format
  auto formatter = supercell_info_dict.parse(properties);
  jsonParser json;
  SupercellInfoData data{shared_prim, supercell_sym_info, params_json};
  formatter.to_json(data, json);
  log << json << std::endl;
}

}  // namespace CASM
