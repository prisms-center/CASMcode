#ifndef CASM_QueryHandler
#define CASM_QueryHandler

#include <map>
#include <memory>

#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class RuntimeLibrary;
class ProjectSettings;

namespace DB {
template <typename DataObject>
class Selected;
template <typename DataObject>
class Selection;
}  // namespace DB

template <typename _DataObject>
class QueryHandler : public notstd::Cloneable {
  CLONEABLE_NEEDS_DESTRUCTOR_DEF(QueryHandler)
 public:
  typedef _DataObject DataObject;

  QueryHandler(const ProjectSettings &set);

  DataFormatterDictionary<DataObject> &dict();

  const DataFormatterDictionary<DataObject> &dict() const;

  /// \brief Set the selection to be used for the 'selected' column
  void set_selected(const DB::Selected<DataObject> &selection);

  /// \brief Set the selection to be used for the 'selected' column
  void set_selected(const DB::Selection<DataObject> &selection);

  /// \brief Add user-defined query alias
  ///
  /// - Aliases are added to memory, but not saved to file until ProjectSettings
  ///   is saved
  void add_alias(const std::string &alias_name,
                 const std::string &alias_command);

  /// \brief const Access aliases map
  ///
  /// - key: alias name
  /// - mapped value: alias command
  const std::map<std::string, std::string> &aliases() const {
    return m_aliases;
  }

 private:
  /// \brief Access aliases map
  ///
  /// - key: alias name
  /// - value: alias command
  std::map<std::string, std::string> &_aliases() { return m_aliases; }

  const ProjectSettings *m_set;

  std::map<std::string, std::string> m_aliases;

  notstd::cloneable_ptr<DataFormatterDictionary<DataObject> > m_dict;

  std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_lib;
};

/// \brief Load enumerator plugins from a CASM project
template <typename DataFormatterDictInserter, typename RuntimeLibInserter>
std::pair<DataFormatterDictInserter, RuntimeLibInserter> load_query_plugins(
    const ProjectSettings &set, DataFormatterDictInserter dict_it,
    RuntimeLibInserter lib_it);

}  // namespace CASM

#endif
