#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  // const std::string traits<IntegralCluster>::name = "IntegralCluster";


  IntegralCluster::IntegralCluster(PrimType const &prim):
    m_prim_ptr(&prim) {}

  typename IntegralCluster::PrimType const &IntegralCluster::prim() const {
    return *m_prim_ptr;
  }

  /// \brief Access vector of elements
  std::vector<xtal::UnitCellCoord> &IntegralCluster::elements() {
    return m_element;
  }

  /// \brief const Access vector of elements
  const std::vector<xtal::UnitCellCoord> &IntegralCluster::elements() const {
    return m_element;
  }

  /// \brief Return the coordinate corresponding to element(i)
  xtal::Coordinate IntegralCluster::coordinate(size_type i) const {
    return this->element(i).coordinate(prim());
  }

  /// \brief Translate the cluster by a UnitCell translation
  IntegralCluster &IntegralCluster::operator+=(xtal::UnitCell trans) {
    for(auto it = this->begin(); it != this->end(); ++it) {
      *it += trans;
    }
    return *this;
  }
}

#include "casm/symmetry/SymTools.hh"
namespace CASM {

  namespace sym {
    /// Apply SymOp to IntegralCluster
    ///
    /// Specialization of template from casm/symmetry/Symtools.hh:
    ///   template <typename Transform, typename Object, typename... Args>
    ///   Object &apply(const Transform &transformation, Object &obj, const Args &... args);
    template<>
    IntegralCluster &apply(SymOp const &op, IntegralCluster &clust, Structure const &prim) {
      for(auto &e : clust) {
        sym::apply(op, e, prim);
      }
      return clust;
    }
  }
}


#include "casm/app/AppIO.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/global/enum/json_io.hh"

namespace CASM {

  // TODO: remove
  // /// \brief Print IntegralCluster to stream, using default Printer<IntegralCluster>
  // std::ostream &operator<<(std::ostream &sout, IntegralCluster const &clust) {
  //   OrbitPrinterOptions opt;
  //   opt.coord_type = INTEGRAL;
  //   SitesPrinter printer {opt};
  //   Log log(sout);
  //   printer.print(clust, log);
  //   return sout;
  // }
  //
  // /// \brief Write IntegralCluster to JSON object
  // ///
  // /// Format:
  // /// \code
  // /// {
  // ///   "min_length" : number,
  // ///   "max_length" : number,
  // ///   "sites" : [
  // ///     [b, i, j, k],
  // ///     ...
  // ///   ]
  // /// }
  // /// \endcode
  // jsonParser &to_json(IntegralCluster const &clust, jsonParser &json) {
  //   json.put_obj();
  //   ClusterInvariants invariants {clust};
  //   json["min_length"] = invariants.displacement().front();
  //   json["max_length"] = invariants.displacement().back();
  //   json["sites"].put_array(clust.begin(), clust.end());
  //   return json;
  // }
  //
  // /// \brief Read IntegralCluster from JSON
  // ///
  // /// Format:
  // /// \code
  // /// {
  // ///   "min_length" : number,
  // ///   "max_length" : number,
  // ///   "coordinate_mode" : ("FRAC", "CART", "INT" (default)) (optional)
  // ///   "sites" : [
  // ///     [b, i, j, k],
  // ///     ...
  // ///   ]
  // /// }
  // /// \endcode
  // ///
  // /// - Also accepts "prototype" in place of "sites"
  // void from_json(IntegralCluster &clust, jsonParser const &json, double xtal_tol) {
  //
  //   std::string name;
  //   if(json.contains("sites")) {
  //     name = "sites";
  //   }
  //   else if(json.contains("prototype")) {
  //     name = "prototype";
  //   }
  //   else {
  //     Log &err_log = default_err_log();
  //     err_log.error("Reading IntegralCluster from JSON");
  //     err_log << "Expected 'sites' or 'prototype' containing a list of cluster site coordinates\n";
  //     err_log << "Input:\n" << json << "\n" << std::endl;
  //
  //     throw std::runtime_error("Error reading IntegralCluster from JSON");
  //   }
  //
  //
  //   CASM::COORD_TYPE coord_type = INTEGRAL;
  //   json.get_if(coord_type, "coordinate_mode");
  //
  //   if(coord_type == INTEGRAL) {
  //     from_json(clust.elements(), json[name]);
  //   }
  //   else {
  //     clust.elements().clear();
  //     for(auto it = json[name].begin(); it != json[name].end(); ++it) {
  //       Eigen::Vector3d vcoord;
  //       from_json(vcoord, *it);
  //
  //       Coordinate tcoord(vcoord, clust.prim().lattice(), coord_type);
  //
  //       clust.elements().emplace_back(UnitCellCoord::from_coordinate(clust.prim(), tcoord, xtal_tol));
  //     }
  //   }
  //   return;
  // }
  //
  // IntegralCluster jsonConstructor<IntegralCluster>::from_json(
  //   jsonParser const &json,
  //   Structure const &prim,
  //   double xtal_tol) {
  //   IntegralCluster clust(prim);
  //   CASM::from_json(clust, json, xtal_tol);
  //   return clust;
  // }

}
