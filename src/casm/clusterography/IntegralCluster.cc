#include "casm/clusterography/IntegralCluster.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/global/enum/json_io.hh"

namespace CASM {

  template class ClusterSymCompare<SymCompare<CRTPBase<AperiodicSymCompare<IntegralCluster> > > >;
  template class AperiodicSymCompare<IntegralCluster>;
  /*
  template bool ClusterSymCompare<AperiodicSymCompare<IntegralCluster> >::compare_impl(
    IntegralCluster const &,
    IntegralCluster const &) const;
  template bool ClusterSymCompare<AperiodicSymCompare<IntegralCluster> >::invariants_compare_impl(
    ClusterInvariants<IntegralCluster> const &,
    ClusterInvariants<IntegralCluster> const &) const;
  */

  template class ClusterSymCompare<SymCompare<CRTPBase<PrimPeriodicSymCompare<IntegralCluster> > > >;
  template class PrimPeriodicSymCompare<IntegralCluster>;
  /*
  template bool ClusterSymCompare<PrimPeriodicSymCompare<IntegralCluster> >::compare_impl(
    IntegralCluster const &,
    IntegralCluster const &) const;
  template bool ClusterSymCompare<PrimPeriodicSymCompare<IntegralCluster> >::invariants_compare_impl(
    ClusterInvariants<IntegralCluster> const &,
    ClusterInvariants<IntegralCluster> const &) const;
  */

  template class ClusterSymCompare<SymCompare<CRTPBase<ScelPeriodicSymCompare<IntegralCluster> > > >;
  template class ScelPeriodicSymCompare<IntegralCluster>;
  /*
  template bool ClusterSymCompare<ScelPeriodicSymCompare<IntegralCluster> >::compare_impl(
    IntegralCluster const &,
    IntegralCluster const &) const;
  template bool ClusterSymCompare<ScelPeriodicSymCompare<IntegralCluster> >::invariants_compare_impl(
    ClusterInvariants<IntegralCluster> const &,
    ClusterInvariants<IntegralCluster> const &) const;
  */

  template class ClusterSymCompare<SymCompare<CRTPBase<WithinScelSymCompare<IntegralCluster> > > >;
  template class WithinScelSymCompare<IntegralCluster>;
  /*
  template bool ClusterSymCompare<WithinScelSymCompare<IntegralCluster> >::compare_impl(
    IntegralCluster const &,
    IntegralCluster const &) const;
  template bool ClusterSymCompare<WithinScelSymCompare<IntegralCluster> >::invariants_compare_impl(
    ClusterInvariants<IntegralCluster> const &,
    ClusterInvariants<IntegralCluster> const &) const;
  */

  /// \brief Print IntegralCluster to stream, using default Printer<IntegralCluster>
  std::ostream &operator<<(std::ostream &sout, const IntegralCluster &clust) {
    OrbitPrinterOptions opt;
    opt.coord_type = INTEGRAL;
    SitesPrinter printer {opt};
    Log log(sout);
    printer.print(clust, log);
    return sout;
  }

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
    json["min_length"] = clust.min_length();
    json["max_length"] = clust.max_length();
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
      //UnitCellCoord coord(clust.prim());
      //  !!TODO!!  from_json(clust.elements(), json[name], clust.prim());
      from_json(clust.elements(), json[name]);
    }
    else {
      clust.elements().clear();
      for(auto it = json[name].begin(); it != json[name].end(); ++it) {
        Eigen::Vector3d vcoord;
        from_json(vcoord, *it);

        Coordinate tcoord(vcoord, clust.prim().lattice(), coord_type);

        clust.elements().emplace_back(tcoord, xtal_tol);
      }
    }
    return;
  }

  IntegralCluster jsonConstructor<IntegralCluster>::from_json(
    const jsonParser &json,
    const Structure &prim,
    double xtal_tol) {
    IntegralCluster clust(prim);
    CASM::from_json(clust, json, xtal_tol);
    return clust;
  }

  IntegralCluster jsonConstructor<IntegralCluster>::from_json(
    const jsonParser &json,
    const PrimClex &primclex) {
    IntegralCluster clust(primclex.prim());
    CASM::from_json(clust, json, primclex.crystallography_tol());
    return clust;
  }

}
