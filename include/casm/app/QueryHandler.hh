#ifndef CASM_ConfigQueryHandler
#define CASM_ConfigQueryHandler

#include <map>
#include <memory>

#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/DataFormatter.hh"

namespace CASM {

  class PrimClex;

  template<typename DataObject>
  class QueryHandler {

  public:

    QueryHandler(const PrimClex &primclex) {};

    ~QueryHandler() {
      // order of deletion matters
      m_dict.clear();
      m_lib.clear();
    }

    DataFormatterDictionary<DataObject> &dict() {
      return m_dict;
    }

    const DataFormatterDictionary<DataObject> &dict() const {
      return m_dict;
    }

  private:

    PrimClex const *m_primclex;

    DataFormatterDictionary<DataObject> m_dict;

    std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_lib;

  };

  /// \brief Load enumerator plugins from a CASM project
  template<typename DataFormatterDictInserter, typename RuntimeLibInserter>
  std::pair<DataFormatterDictInserter, RuntimeLibInserter>
  load_query_plugins(
    const PrimClex &primclex,
    DataFormatterDictInserter dict_it,
    RuntimeLibInserter lib_it);

}

#endif
