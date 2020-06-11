#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/global/enum/json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

  /// \brief Write IntegralCluster to JSON object
  ///
  /// Format:
  /// \code
  /// {
  ///   "min_length" : number,
  ///   "max_length" : number,
  ///   "sites" : [
  ///     [b, i, j, k],
  ///     ...
  ///   ]
  /// }
  /// \endcode
  jsonParser &to_json(const IntegralCluster &clust, jsonParser &json) {
    json.put_obj();
    ClusterInvariants invariants {clust};
    json["min_length"] = invariants.displacement().front();
    json["max_length"] = invariants.displacement().back();
    json["sites"].put_array(clust.begin(), clust.end());
    return json;
  }

  /// \brief Read IntegralCluster from JSON
  ///
  /// Format:
  /// \code
  /// {
  ///   "min_length" : number,
  ///   "max_length" : number,
  ///   "coordinate_mode" : ("FRAC", "CART", "INT" (default)) (optional)
  ///   "sites" : [
  ///     [b, i, j, k],
  ///     ...
  ///   ]
  /// }
  /// \endcode
  ///
  /// - Also accepts "prototype" in place of "sites"
  void from_json(IntegralCluster &clust, const jsonParser &json, double xtal_tol) {

    std::string name;
    if(json.contains("sites")) {
      name = "sites";
    }
    else if(json.contains("prototype")) {
      name = "prototype";
    }
    else {
      Log &err_log = default_err_log();
      err_log.error("Reading IntegralCluster from JSON");
      err_log << "Expected 'sites' or 'prototype' containing a list of cluster site coordinates\n";
      err_log << "Input:\n" << json << "\n" << std::endl;

      throw std::runtime_error("Error reading IntegralCluster from JSON");
    }


    CASM::COORD_TYPE coord_type = INTEGRAL;
    json.get_if(coord_type, "coordinate_mode");

    if(coord_type == INTEGRAL) {
      from_json(clust.elements(), json[name]);
    }
    else {
      clust.elements().clear();
      for(auto it = json[name].begin(); it != json[name].end(); ++it) {
        Eigen::Vector3d vcoord;
        from_json(vcoord, *it);

        xtal::Coordinate tcoord(vcoord, clust.prim().lattice(), coord_type);

        clust.elements().emplace_back(xtal::UnitCellCoord::from_coordinate(clust.prim(), tcoord, xtal_tol));
      }
    }
    return;
  }

  IntegralCluster jsonConstructor<IntegralCluster>::from_json(
    const jsonParser &json,
    Structure const &prim,
    double xtal_tol) {
    IntegralCluster clust(prim);
    CASM::from_json(clust, json, xtal_tol);
    return clust;
  }

  /// \brief Parse IntegralCluster from JSON
  ///
  /// Format:
  /// \code
  /// {
  ///   "coordinate_mode" : ("FRAC", "CART", "INT" (default)) (optional)
  ///   "sites" : [
  ///     [b, i, j, k],
  ///     ...
  ///   ]
  /// }
  /// \endcode
  ///
  /// - Also accepts "prototype" in place of "sites"
  void parse(
    InputParser<IntegralCluster> &parser,
    const std::shared_ptr<Structure const> &shared_prim) {

    std::string name;
    const jsonParser &json = parser.self;
    if(json.contains("sites")) {
      name = "sites";
    }
    else if(json.contains("prototype")) {
      name = "prototype";
    }
    else {
      parser.error.insert("Error reading IntegralCluster from JSON: Expected 'sites' or 'prototype' containing a list of cluster site coordinates.");
      return;
    }

    double xtal_tol = shared_prim->lattice().tol();
    CASM::COORD_TYPE coord_type;
    parser.optional_else(coord_type, "coordinate_mode", INTEGRAL);

    parser.value = notstd::make_unique<IntegralCluster>(*shared_prim);
    auto &clust = *parser.value;

    if(coord_type == INTEGRAL) {
      parser.require(clust.elements(), name);
    }
    else {
      std::vector<Eigen::Vector3d> coord_vec;
      parser.require(coord_vec, name);

      try {
        for(const auto &coord : coord_vec) {
          xtal::Coordinate tcoord {coord, shared_prim->lattice(), coord_type};
          clust.elements().emplace_back(
            xtal::UnitCellCoord::from_coordinate(*shared_prim, tcoord, xtal_tol));
        }
      }
      catch(std::exception &e) {
        parser.error.insert("Error: could not read coordinates from '" + name + "'");
        return;
      }
    }
    return;
  }
}
