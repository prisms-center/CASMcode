#ifndef CASM_EnumeratorHandler
#define CASM_EnumeratorHandler

#include <map>
#include <memory>

#include "casm/container/Enumerator.hh"
#include "casm/system/RuntimeLibrary.hh"

namespace CASM {

  class ProjectSettings;

  class EnumeratorHandler {

  public:

    EnumeratorHandler(const ProjectSettings &set);

    ~EnumeratorHandler() {
      // order of deletion matters
      m_enumerator.clear();
      m_lib.clear();
    }

    EnumeratorMap &map() {
      return m_enumerator;
    }

    const EnumeratorMap &map() const {
      return m_enumerator;
    }

  private:

    const ProjectSettings *m_set;

    EnumeratorMap m_enumerator;

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
