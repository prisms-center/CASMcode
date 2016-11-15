#ifndef CASM_QueryHandler
#define CASM_QueryHandler

#include <map>
#include <memory>

#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/DataFormatter.hh"

namespace CASM {

  class ProjectSettings;

  template<typename DataObject>
  struct QueryTraits {};

  template<typename _DataObject>
  class QueryHandler : public notstd::Cloneable {

  public:

    typedef _DataObject DataObject;

    QueryHandler(const ProjectSettings &set);

    ~QueryHandler();

    DataFormatterDictionary<DataObject> &dict();

    const DataFormatterDictionary<DataObject> &dict() const;

    /// \brief Set the selection to be used for the 'selected' column
    ///
    /// - ToDo: generalize ConfigIO::Selected
    void set_selected(const typename QueryTraits<DataObject>::Selected &selection);

    /// \brief Set the selection to be used for the 'selected' column
    ///
    /// - ToDo: generalize ConstConfigSelection
    void set_selected(const typename QueryTraits<DataObject>::Selection &selection);

    /// \brief Add user-defined query alias
    ///
    /// - Aliases are added to memory, but not saved to file until ProjectSettings
    ///   is saved
    void add_alias(const std::string &alias_name, const std::string &alias_command);

    /// \brief const Access aliases map
    ///
    /// - key: alias name
    /// - mapped value: alias command
    const std::map<std::string, std::string> &aliases() const {
      return m_aliases;
    }

    std::unique_ptr<QueryHandler<DataObject> > clone() const {
      return std::unique_ptr<QueryHandler<DataObject> >(this->_clone());
    }

  private:

    /// \brief Access aliases map
    ///
    /// - key: alias name
    /// - value: alias command
    std::map<std::string, std::string> &_aliases() {
      return m_aliases;
    }

    QueryHandler<DataObject> *_clone() const override {
      return new QueryHandler<DataObject>(*this);
    }

    const ProjectSettings *m_set;

    std::map<std::string, std::string> m_aliases;

    DataFormatterDictionary<DataObject> m_dict;

    std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_lib;

  };

  /// \brief Load enumerator plugins from a CASM project
  template<typename DataFormatterDictInserter, typename RuntimeLibInserter>
  std::pair<DataFormatterDictInserter, RuntimeLibInserter>
  load_query_plugins(
    const ProjectSettings &set,
    DataFormatterDictInserter dict_it,
    RuntimeLibInserter lib_it);

}

#endif
