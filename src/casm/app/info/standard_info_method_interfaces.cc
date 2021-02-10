#include "casm/app/info/standard_info_method_interfaces.hh"

#include <vector>

#include "casm/app/info/InfoInterface.hh"
#include "casm/app/info/methods/PrimInfoInterface.hh"
#include "casm/app/info/methods/SupercellInfoInterface.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

/// A vector containing `casm info` method interfaces
InfoInterfaceVector make_standard_info_method_interfaces() {
  InfoInterfaceVector vec;
  vec.emplace_back(notstd::make_cloneable<PrimInfoInterface>());
  vec.emplace_back(notstd::make_cloneable<SupercellInfoInterface>());
  return vec;
}

}  // namespace CASM
