#ifndef CASM_enum_StandardEnumerationInterfaces
#define CASM_enum_StandardEnumerationInterfaces

#include <vector>
#include "casm/app/enum/EnumInterface.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  class EnumInterfaceBase;
  typedef std::vector<notstd::cloneable_ptr<EnumInterfaceBase>> EnumInterfaceVector;

  /// A map containing interfaces that allow `casm enum` to run the enumeration methods in libcasm
  EnumInterfaceVector make_standard_enumerator_interfaces();

}

#endif
