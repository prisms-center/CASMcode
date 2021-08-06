#ifndef CASM_DB_Remove
#define CASM_DB_Remove

#include "casm/database/ConfigData.hh"

namespace CASM {
namespace Completer {
class RmOption;
}
}  // namespace CASM

// To be specialized for calculable 'ConfigType' classes:
//   StructureMap<ConfigType>::map(std::string, DatabaseIterator<ConfigType>
//   hint, inserter result) Import<ConfigType>::desc Import<ConfigType>::run
//   Import<ConfigType>::_import_formatter
//   Update<ConfigType>::desc
//   Update<ConfigType>::run
//   Update<ConfigType>::_update_formatter
//   Remove<ConfigType>::desc
//   Remove<ConfigType>::run

namespace CASM {
namespace DB {

/// Generic ConfigType-dependent part of Remove
template <typename _ConfigType>
class RemoveT : protected ConfigData {
 public:
  typedef _ConfigType ConfigType;

  RemoveT(const PrimClex &primclex, std::string report_dir);

  /// \brief Erase properties via structure origin file path
  void erase_properties_via_origin(std::vector<std::string> const &origins,
                                   bool dry_run);

  /// \brief Erase Configurations that have no data
  void erase(const DB::Selection<ConfigType> &selection, bool dry_run);

  /// \brief Erase data and files (permanently), but not Configuration
  void erase_data(const DB::Selection<ConfigType> &selection, bool dry_run);

  /// \brief Removes Configurations and data and files (permanently)
  ///
  /// - Data are always associated with one 'from' configuration, so the
  ///   selection here indicates 'from' configurations
  /// - The 'to' configurations are updated with the new best mapping properties
  void erase_all(const DB::Selection<ConfigType> &selection, bool dry_run);

 private:
  void _erase_report(const std::vector<std::string> &fail);

  std::string m_report_dir;
};

// To be specialized for ConfigType (default implementation exists for
// ConfigTypes)
template <typename ConfigType>
class Remove : public RemoveT<ConfigType> {
 public:
  Remove(const PrimClex &primclex, std::string report_dir);
  static std::string desc();
  static int run(const PrimClex &, const Completer::RmOption &opt);
};
}  // namespace DB
}  // namespace CASM

#endif
