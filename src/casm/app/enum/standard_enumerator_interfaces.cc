#include "casm/app/enum/standard_enumerator_interfaces.hh"

#include <vector>

#include "casm/app/enum/EnumInterface.hh"
#include "casm/app/enum/methods/ConfigEnumAllOccupationsInterface.hh"
#include "casm/app/enum/methods/ConfigEnumRandomLocalInterface.hh"
#include "casm/app/enum/methods/ConfigEnumRandomOccupationsInterface.hh"
#include "casm/app/enum/methods/ConfigEnumSiteDoFsInterface.hh"
#include "casm/app/enum/methods/ConfigEnumStrainInterface.hh"
#include "casm/app/enum/methods/ScelEnumInterface.hh"
#include "casm/app/enum/methods/SuperConfigEnumInterface.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

/// A map containing interfaces that allow `casm enum` to run the enumeration
/// methods in libcasm
EnumInterfaceVector make_standard_enumerator_interfaces() {
  EnumInterfaceVector vec;
  vec.emplace_back(notstd::make_cloneable<ConfigEnumAllOccupationsInterface>());
  vec.emplace_back(notstd::make_cloneable<ConfigEnumRandomLocalInterface>());
  vec.emplace_back(
      notstd::make_cloneable<ConfigEnumRandomOccupationsInterface>());
  vec.emplace_back(notstd::make_cloneable<ConfigEnumSiteDoFsInterface>());
  vec.emplace_back(notstd::make_cloneable<ConfigEnumStrainInterface>());
  vec.emplace_back(notstd::make_cloneable<ScelEnumInterface>());
  vec.emplace_back(notstd::make_cloneable<SuperConfigEnumInterface>());

  return vec;
}

}  // namespace CASM
