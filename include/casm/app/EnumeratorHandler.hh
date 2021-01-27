#ifndef CASM_EnumeratorHandler
#define CASM_EnumeratorHandler

#include <map>
#include <memory>

#include "casm/app/casm_functions.hh"
#include "casm/app/enum/standard_enumerator_interfaces.hh"

namespace CASM {

namespace Completer {
class EnumOption;
}

class RuntimeLibrary;
class ProjectSettings;

class EnumeratorHandler : public notstd::Cloneable {
  CLONEABLE_NEEDS_DESTRUCTOR_DEF(EnumeratorHandler)
 public:
  EnumeratorHandler(ProjectSettings const &set);

  EnumInterfaceVector &get() { return m_enumerator; }

  EnumInterfaceVector const &get() const { return m_enumerator; }

 private:
  ProjectSettings const *m_set;

  EnumInterfaceVector m_enumerator;

  std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_lib;
};

/// \brief Load enumerator plugins from a CASM project
template <typename EnumInterfaceVectorInserter, typename RuntimeLibInserter>
std::pair<EnumInterfaceVectorInserter, RuntimeLibInserter>
load_enumerator_plugins(ProjectSettings const &set,
                        EnumInterfaceVectorInserter enum_it,
                        RuntimeLibInserter lib_it);

}  // namespace CASM

#endif
