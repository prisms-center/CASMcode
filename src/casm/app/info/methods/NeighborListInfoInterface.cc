#include "casm/app/info/methods/NeighborListInfoInterface.hh"

#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/io/json/ClusterSpecs_json_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/symmetry/SupercellSymInfo.hh"

namespace CASM {

namespace {

/// Data structure for formatting properties of prim and supercell_sym_info
struct NeighborListInfoData {
  NeighborListInfoData(std::shared_ptr<Structure const> _shared_prim,
                       SupercellSymInfo const &_supercell_sym_info,
                       PrimNeighborList const &_prim_neighbor_list,
                       SuperNeighborList const &_supercell_neighbor_list)
      : shared_prim(_shared_prim),
        supercell_sym_info(_supercell_sym_info),
        prim_neighbor_list(_prim_neighbor_list),
        supercell_neighbor_list(_supercell_neighbor_list) {}

  std::shared_ptr<Structure const> shared_prim;
  SupercellSymInfo const &supercell_sym_info;
  PrimNeighborList const &prim_neighbor_list;
  SuperNeighborList const &supercell_neighbor_list;
};

// use for ValueType=bool, int, double, std::string, jsonParser
template <typename ValueType>
using NeighborListInfoFormatter =
    GenericDatumFormatter<ValueType, NeighborListInfoData>;

NeighborListInfoFormatter<jsonParser> prim_neighbor_list() {
  return NeighborListInfoFormatter<jsonParser>(
      "prim_neighbor_list",
      "Contains an array of unitcells, the ordered set of unitcells the make "
      "up the neighborhood of the origin unit cell, along with the "
      "nlist_weight_matrix and nlist_sublat_indices.",
      [](NeighborListInfoData const &data) -> jsonParser {
        jsonParser json;
        std::vector<UnitCell> unitcells{data.prim_neighbor_list.begin(),
                                        data.prim_neighbor_list.end()};
        to_json(unitcells, json["unitcells"], jsonParser::as_flattest());
        json["weight_matrix"] = data.prim_neighbor_list.weight_matrix();
        json["sublat_indices"] = data.prim_neighbor_list.sublat_indices();
        return json;
      });
}

NeighborListInfoFormatter<jsonParser> supercell_neighbor_list() {
  return NeighborListInfoFormatter<jsonParser>(
      "supercell_neighbor_list",
      "Contains, for each unitcell in the supercell, "
      "`linear_unitcell_indices`, the neighborhood of unitcells, and "
      "`linear_site_indices` the neighborhood of sites. The linear indices can "
      "be used to lookup coordinates with the `unitcells` and "
      "`integral_site_coordinates` properties of the supercell.",
      [](NeighborListInfoData const &data) -> jsonParser {
        auto const &sym_info = data.supercell_sym_info;
        auto const &T = sym_info.transformation_matrix_to_super();
        Index volume = T.determinant();
        auto const &supercell_lattice = sym_info.supercell_lattice();
        auto const &L = supercell_lattice.lat_column_mat();
        auto const &scel_nlist = data.supercell_neighbor_list;

        jsonParser json;
        json["transformation_matrix_to_super"] = T;
        json["supercell_lattice_column_matrix"] = L;
        json["supercell_lattice_vectors"] = L.transpose();
        json["supercell_name"] = make_supercell_name(
            data.shared_prim->point_group(), sym_info.prim_lattice(),
            sym_info.supercell_lattice());
        json["supercell_volume"] = volume;

        json["n_neighbor_sites"] = scel_nlist.sites(0).size();
        json["n_neighbor_unitcells"] = scel_nlist.unitcells(0).size();
        json["periodic_images_of_neighborhood_overlap"] = scel_nlist.overlaps();

        json["linear_site_indices"] = jsonParser::array();
        json["linear_unitcell_indices"] = jsonParser::array();
        for (Index l = 0; l < volume; ++l) {
          json["linear_site_indices"].push_back(scel_nlist.sites(l));
          json["linear_unitcell_indices"].push_back(scel_nlist.unitcells(l));
        }
        return json;
      });
}

DataFormatterDictionary<NeighborListInfoData> make_neighbor_list_info_dict() {
  DataFormatterDictionary<NeighborListInfoData> neighbor_list_info_dict;

  // properties that require prim and supercell_sym_info
  neighbor_list_info_dict.insert(prim_neighbor_list(),
                                 supercell_neighbor_list());
  return neighbor_list_info_dict;
}

// used with `for_all_orbits` to expand the neighbor list with the orbits
// contructed by a ClusterSpecs object
struct ExpandPrimNeighborList {
  ExpandPrimNeighborList(PrimNeighborList &_prim_neighbor_list)
      : prim_neighbor_list(_prim_neighbor_list) {}

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) const {
    for (const auto &orbit : orbits) {
      for (const auto &equiv : orbit) {
        for (const auto &site : equiv) {
          prim_neighbor_list.expand(site.unitcell());
        }
      }
    }
  }

  PrimNeighborList &prim_neighbor_list;
};

}  // namespace

std::string NeighborListInfoInterface::desc() const {
  std::string description =
      "Get prim and supercell neighbor list information. The supercell is \n"
      "specified by the prim and one of the following (else the primitive \n"
      "cell is used):                                                     \n"
      "- transformation_matrix_to_super                                   \n"
      "- supercell_lattice_vectors                                        \n"
      "- supercell_lattice_column_matrix                                  \n"
      "- supercell_name                                                   \n\n";

  std::string custom_options =
      "  prim: JSON object (optional, default=prim of current project)    \n"
      "    See `casm format --prim` for details on the prim format.       \n\n"

      "  nlist_weight_matrix: 3x3 array of integer (optional)             \n"
      "    The neighbor list weight matrix, W, defines the canonical order\n"
      "    of neighboring UnitCell through lexicographically sorting      \n"
      "    [r, i, j, k], where r = (i,j,k).transpose() * W * (i,j,k). If  \n"
      "    not provided, it is obtained from 1) the settings of the       \n"
      "    current project, or 2) an appropriate default for the prim     \n"
      "    lattice so that unit cells are added to the neighborhood in an \n"
      "    approximate sphere around the origin.                          \n\n"

      "  nlist_sublat_indices: 3x3 array of integer (optional)            \n"
      "    The indices of sublattices that should be included in the      \n"
      "    supercell neighbor list. If not provided, it is obtained from  \n"
      "    1) the settings of the current project, or 2) indices of the   \n"
      "    sites that have >= 2 occupant DoF, or continuous DoF.          \n\n"

      "  cluster_specs: array of JSON objects (optional)                  \n"
      "    The `cluster_specs` array holds one or more JSON descriptions  \n"
      "    of cluster orbits, as used in `bspecs.json`. The prim neighbor \n"
      "    is expanded to include the sites in all cluster orbits.        \n\n"

      "  basis_set_names: array of string (optional)                      \n"
      "    Names of basis sets in the current project whose orbits are    \n"
      "    all added to the prim neighbor list.                           \n\n"

      "  unitcells: array of [i,j,k] (optional)                           \n"
      "    Array of unit cells (specified by [i, j, k] multiples of the   \n"
      "    prim lattice vectors) that should be added to the prim neighbor\n"
      "    list.                                                          \n\n"

      "  transformation_matrix_to_super: 3x3 array of integer (optional)  \n"
      "    Transformation matrix T, defining the supercell lattice vectors\n"
      "    S, in terms of the prim lattice vectors, P: `S = P * T`, where \n"
      "    S and P are column vector matrices.                            \n\n"

      "  supercell_lattice_vectors: 3x3 array of integer (optional)       \n"
      "    Supercell lattice vectors, as a row vector matrix.             \n\n"

      "  supercell_lattice_column_matrix: 3x3 array of integer (optional) \n"
      "    Supercell lattice vectors, as a column vector matrix.          \n\n"

      "  supercell_name: string (optional)                                \n"
      "    Unique name given to a supercell, based on the hermite normal  \n"
      "    form, of the transformation_matrix_to_super and, if not        \n"
      "    canonical, the index of the prim factor group operation that   \n"
      "    transforms the canonical supercell into this supercell.        \n\n"

      "  properties: array of string (optional, default=[])               \n"
      "    An array of strings specifying which neighbor list properties  \n"
      "    to output. The default value of an empty array will print all  \n"
      "    properties. The allowed options are:                           \n\n";

  std::stringstream ss;
  auto dict = make_neighbor_list_info_dict();
  dict.print_help(ss, DatumFormatterClass::Property);

  return name() + ": \n\n" + description + custom_options + ss.str();
}

std::string NeighborListInfoInterface::name() const {
  return "NeighborListInfo";
}

/// Run `NeighborListInfo` info method
void NeighborListInfoInterface::run(jsonParser const &json_options,
                                    PrimClex const *primclex) const {
  Log &log = CASM::log();

  ParentInputParser parser{json_options};
  std::runtime_error error_if_invalid{"Error reading NeighborListInfo input"};

  std::unique_ptr<PrimClex> primclex_ptr;
  if (primclex == nullptr) {
    fs::path root = find_casmroot(fs::current_path());
    if (!root.empty()) {
      primclex_ptr = notstd::make_unique<PrimClex>(root);
      primclex = primclex_ptr.get();
    }
  }

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
    std::stringstream msg;
    msg << "Error in NeighborListInfo: No \"prim\" in input and no project "
           "provided or found.";
    parser.insert_error("prim", msg.str());
  }
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // read "nlist_weight_matrix" and "nlist_sublat_indices"
  Eigen::Matrix3l nlist_weight_matrix;
  std::set<int> nlist_sublat_indices;
  if (parser.self.contains("nlist_weight_matrix") ||
      parser.self.contains("nlist_sublat_indices")) {
    parser.optional(nlist_weight_matrix, "nlist_weight_matrix");
    parser.optional(nlist_sublat_indices, "nlist_sublat_indices");

  } else if (primclex != nullptr) {
    // if project provided via api
    try {
      nlist_weight_matrix = primclex->settings().nlist_weight_matrix();
      nlist_sublat_indices = primclex->settings().nlist_sublat_indices();
    } catch (std::exception &e) {
      parser.error.insert(e.what());
    }
  } else {
    nlist_weight_matrix =
        default_nlist_weight_matrix(*shared_prim, shared_prim->lattice().tol());
    nlist_sublat_indices = default_nlist_sublat_indices(*shared_prim);
  }
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  PrimNeighborList prim_neighbor_list{nlist_weight_matrix,
                                      nlist_sublat_indices.begin(),
                                      nlist_sublat_indices.end()};
  ExpandPrimNeighborList neighbor_list_expander{prim_neighbor_list};

  // read "cluster_specs" (optional)
  if (parser.self.contains("cluster_specs")) {
    std::vector<jsonParser> cluster_specs_json_vec;
    parser.optional(cluster_specs_json_vec, "cluster_specs");

    for (jsonParser const &cluster_specs_json : cluster_specs_json_vec) {
      InputParser<ClusterSpecs> parser{cluster_specs_json, shared_prim};
      if (!parser.valid()) {
        continue;
      }
      for_all_orbits(*parser.value, err_log(), neighbor_list_expander);
    }
  }

  // read "basis_set_names" (optional)
  if (primclex != nullptr && parser.self.contains("basis_set_names")) {
    std::vector<std::string> basis_set_names;
    parser.optional(basis_set_names, "basis_set_names");
    for (auto basis_set_name : basis_set_names) {
      if (!primclex->has_basis_set_specs(basis_set_name)) {
        std::stringstream msg;
        msg << "No basis set named: " << basis_set_name;
        parser.insert_error("basis_set_names", msg.str());
      }

      ClexBasisSpecs const &basis_set_specs =
          primclex->basis_set_specs(basis_set_name);
      for_all_orbits(*basis_set_specs.cluster_specs, err_log(),
                     neighbor_list_expander);
    }
  }

  // read "unitcells" (optional)
  if (parser.self.contains("unitcells")) {
    std::vector<UnitCell> unitcells;
    parser.optional(unitcells, "unitcells");
    prim_neighbor_list.expand(unitcells.begin(), unitcells.end());
  }

  // read "transformation_matrix_to_super"
  Eigen::Matrix3l T;
  if (parser.self.contains("transformation_matrix_to_super")) {
    parser.optional(T, "transformation_matrix_to_super");

    // or read "supercell_lattice_vectors"
  } else if (parser.self.contains("supercell_lattice_vectors")) {
    Eigen::Matrix3d L_transpose;
    parser.optional(L_transpose, "supercell_lattice_vectors");
    Lattice super_lattice{L_transpose.transpose()};
    T = make_transformation_matrix_to_super(shared_prim->lattice(),
                                            super_lattice, TOL);

    // or read "supercell_lattice_column_matrix"
  } else if (parser.self.contains("supercell_lattice_column_matrix")) {
    Eigen::Matrix3d L;
    parser.optional(L, "supercell_lattice_column_matrix");
    Lattice super_lattice{L};
    T = make_transformation_matrix_to_super(shared_prim->lattice(),
                                            super_lattice, TOL);

    // or read "supercell_name"
  } else if (parser.self.contains("supercell_name")) {
    std::string supercell_name;
    parser.optional(supercell_name, "supercell_name");
    xtal::Superlattice superlattice = make_superlattice_from_supercell_name(
        shared_prim->factor_group(), shared_prim->lattice(), supercell_name);
    T = superlattice.transformation_matrix_to_super();

    // else use Identity (prim cell)
  } else {
    T = Eigen::Matrix3l::Identity();
  }

  // read "properties"
  std::vector<std::string> properties;
  parser.optional(properties, "properties");
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  Lattice supercell_lattice = make_superlattice(shared_prim->lattice(), T);
  SupercellSymInfo supercell_sym_info =
      make_supercell_sym_info(*shared_prim, supercell_lattice);
  SuperNeighborList supercell_neighbor_list{T, prim_neighbor_list};

  auto dict = make_neighbor_list_info_dict();
  // output all properties if empty
  if (properties.empty()) {
    for (auto it = dict.begin(); it != dict.end(); ++it) {
      if (it->type() == DatumFormatterClass::Property) {
        properties.push_back(it->name());
      }
    }
  }

  // format
  auto formatter = dict.parse(properties);
  jsonParser json;
  NeighborListInfoData data{shared_prim, supercell_sym_info, prim_neighbor_list,
                            supercell_neighbor_list};
  formatter.to_json(data, json);
  log << json << std::endl;
}

}  // namespace CASM
