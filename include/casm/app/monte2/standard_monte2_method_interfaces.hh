#ifndef CASM_monte2_StandardMethodInterfaces
#define CASM_monte2_StandardMethodInterfaces

#include <vector>

#include "casm/app/monte2/Monte2Interface.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class Monte2InterfaceBase;
typedef std::vector<notstd::cloneable_ptr<Monte2InterfaceBase>>
    Monte2InterfaceVector;

/// A map containing interfaces that allow `casm monte2` to run the enumeration
/// methods in libcasm
Monte2InterfaceVector make_standard_monte2_method_interfaces();

}  // namespace CASM

#endif
