#ifndef CASM_EnumeratorHandler
#define CASM_EnumeratorHandler

#include <map>
#include <memory>

#include "casm/app/casm_functions.hh"

namespace CASM {

  namespace Completer {
    class EnumOption;
  }

  typedef InterfaceMap<Completer::EnumOption> EnumeratorMap;
  class RuntimeLibrary;
  class ProjectSettings;

  class EnumeratorHandler {

  public:

    EnumeratorHandler(const ProjectSettings &set);

    ~EnumeratorHandler() {
      // order of deletion matters
      m_enumerator->clear();
      m_lib.clear();
    }

    EnumeratorMap &map() {
      return *m_enumerator;
    }

    const EnumeratorMap &map() const {
      return *m_enumerator;
    }

  private:

    const ProjectSettings *m_set;

    notstd::cloneable_ptr<EnumeratorMap> m_enumerator;

    std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_lib;

  };

  /// \brief Load enumerator plugins from a CASM project
  template<typename EnumeratorMapInserter, typename RuntimeLibInserter>
  std::pair<EnumeratorMapInserter, RuntimeLibInserter>
  load_enumerator_plugins(
    const ProjectSettings &set,
    EnumeratorMapInserter enum_it,
    RuntimeLibInserter lib_it);

}

#endif
