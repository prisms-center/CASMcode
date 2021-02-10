#ifndef CASM_info_StandardInfoMethodInterfaces
#define CASM_info_StandardInfoMethodInterfaces

#include <vector>

#include "casm/app/info/InfoInterface.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class InfoInterfaceBase;
typedef std::vector<notstd::cloneable_ptr<InfoInterfaceBase>>
    InfoInterfaceVector;

/// A vector containing `casm info` method interfaces
InfoInterfaceVector make_standard_info_method_interfaces();

}  // namespace CASM

#endif
