#include "casm/kinetics/PrimPeriodicDiffTransOrbitTraits.hh"

#include <string>

namespace CASM {

  const std::string traits<PrimPeriodicDiffTransOrbit>::name = "PrimPeriodicDiffTransOrbit";

  const std::string traits<PrimPeriodicDiffTransOrbit>::short_name = "diff_trans";

  const std::string traits<PrimPeriodicDiffTransOrbit>::orbit_type_name = "diff_trans";
  /// does lexicographical comparison
  bool traits<PrimPeriodicDiffTransOrbit>::name_compare(std::string A, std::string B) {
    return A < B;
  };
}
