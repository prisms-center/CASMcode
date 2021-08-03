#include "casm/app/monte2/standard_monte2_method_interfaces.hh"

#include <vector>

#include "casm/app/monte2/Monte2Interface.hh"
#include "casm/app/monte2/methods/CanonicalMonteCarloInterface.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

/// A map containing interfaces that allow `casm enum` to run the enumeration
/// methods in libcasm
Monte2InterfaceVector make_standard_monte2_method_interfaces() {
  Monte2InterfaceVector vec;
  vec.emplace_back(notstd::make_cloneable<CanonicalMonteCarloInterface>());

  return vec;
}

}  // namespace CASM
