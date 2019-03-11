#include "casm/kinetics/DiffusionTransformationTraits.hh"

#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"

#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"

namespace CASM {

  template class DatabaseTypeOrbit<PrimPeriodicDiffTransSymCompare>;

  UnitCellCoord traits<Kinetics::DiffusionTransformation>::position(
    const Kinetics::DiffusionTransformation &diff_trans) {
    return diff_trans.occ_transform()[0].uccoord;
  }

  void PrimPeriodicDiffTransOrbitTraits::write_pos(const OrbitType &orbit, std::ostream &sout) {
    std::vector<PrimPeriodicDiffTransOrbit> container;
    container.push_back(orbit);
    PrototypePrinter<Kinetics::DiffusionTransformation> printer;
    Log out(sout);
    print_clust(container.begin(), container.end(), out, printer);
  };

  std::string PrimPeriodicDiffTransOrbitTraits::generate_name_impl(const OrbitType &orbit) {
    if(orbit.id() == "none") {
      const auto &db = orbit.primclex().db<PrimPeriodicDiffTransOrbit>();
      auto find_it = db.search(orbit);
      if(find_it != db.end()) {
        return find_it->name();
      }
    }
    return traits<PrimPeriodicDiffTransOrbit>::orbit_type_name + "/" + orbit.id();
  };


  const std::string traits<PrimPeriodicDiffTransOrbit>::name = "PrimPeriodicDiffTransOrbit";

  const std::string traits<PrimPeriodicDiffTransOrbit>::short_name = "diff_trans";

  const std::string traits<PrimPeriodicDiffTransOrbit>::orbit_type_name = "diff_trans";

  /// does lexicographical comparison
  bool traits<PrimPeriodicDiffTransOrbit>::name_compare(std::string A, std::string B) {
    return A < B;
  };
}

