#ifndef CASM_DB_Update
#define CASM_DB_Update

#include "casm/database/ConfigData.hh"

namespace CASM {
namespace Completer {
class UpdateOption;
}
}  // namespace CASM

// To be specialized for calculable 'ConfigType' classes:
//   StructureMap<ConfigType>::map(fs::path, DatabaseIterator<ConfigType> hint,
//   inserter result) Import<ConfigType>::desc Import<ConfigType>::run
//   Import<ConfigType>::_import_formatter
//   Update<ConfigType>::desc
//   Update<ConfigType>::run
//   Update<ConfigType>::_update_formatter
//   Remove<ConfigType>::desc
//   Remove<ConfigType>::run

namespace CASM {
namespace DB {

template <typename T>
class Selection;

struct UpdateSettings {
  UpdateSettings(bool _output_as_json = true)
      : output_as_json(_output_as_json) {}

  void set_default() { *this = UpdateSettings(); }

  // Output reports as JSON instead of columns
  bool output_as_json;
};

/// Generic ConfigType-dependent part of Import
template <typename _ConfigType>
class UpdateT : protected ConfigData {
 public:
  typedef _ConfigType ConfigType;

  /// \brief Constructor
  UpdateT(const PrimClex &primclex, const StructureMap<ConfigType> &mapper,
          UpdateSettings const &_set, std::string report_dir);

  /// \brief Re-parse calculations 'from' all selected configurations
  void update(const DB::Selection<ConfigType> &selection, bool force);

  UpdateSettings const &settings() const { return m_set; }

 protected:
  // Allow ConfigType to specialize the report formatting for 'update'
  virtual DataFormatter<ConfigIO::Result> _update_formatter() const = 0;

  void _update_report(std::vector<ConfigIO::Result> &results,
                      const DB::Selection<ConfigType> &selection) const;

 private:
  const StructureMap<ConfigType> &m_structure_mapper;

  UpdateSettings m_set;

  std::string m_report_dir;
};

// To be specialized for ConfigType (no default implemenation exists)
template <typename ConfigType>
class Update;
/*: public UpdateT<ConfigType> {
public:
  static const std::string desc;
  int run(const PrimClex &, const jsonParser &input, const
Completer::UpdateOption &opt);

private:
  // Allow ConfigType to specialize the report formatting for 'update'
  DataFormatter<ConfigIO::Result> _update_formatter() const override;
  };*/

}  // namespace DB
}  // namespace CASM

#endif
