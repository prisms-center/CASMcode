#ifndef CASM_DB_ConfigData
#define CASM_DB_ConfigData

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/dataformatter/DataFormatter.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools.hh"
#include "casm/clex/MappedProperties.hh"
#include "casm/global/definitions.hh"

// To be specialized for calculable 'ConfigType' classes:
//   Import<ConfigType>::desc
//   Import<ConfigType>::run
//   Import<ConfigType>::_import_formatter
//   Update<ConfigType>::desc
//   Update<ConfigType>::run
//   Update<ConfigType>::_update_formatter
//   Remove<ConfigType>::desc
//   Remove<ConfigType>::run

namespace CASM {

class PrimClex;
class Supercell;

namespace DB {

template <typename T>
class Selection;
template <typename ValueType>
class Database;
template <typename ValueType>
class DatabaseIterator;
class PropertiesDatabase;

/// Create a new report directory to avoid overwriting existing results
std::string create_report_dir(std::string report_dir);

namespace ConfigIO {

/// Data structure for data import results
struct ImportData {
  ImportData()
      : preexisting(false),
        copy_data(false),
        copy_more(false),
        is_best(false) {}

  // base responsibility:
  bool preexisting;
  bool copy_data;
  bool copy_more;
  bool is_best;
};

/// Data structure for mapping / import results
struct Result {
  Result() : has_data(false), has_complete_data(false), is_new_config(false) {}

  // Set 'to'/'from' as empty strings if no mapping possible
  MappedProperties properties;

  // Path to properties.calc.json or POS file that was imported
  std::string pos_path;

  // If a properties.calc.json file is found in standard locations
  bool has_data;

  // If all PrimClex required properties exist
  bool has_complete_data;

  // If 'to' Configuration did not exist in the database prior to mapping
  bool is_new_config;

  // If mapping failed, stores any error message that may be generated
  std::string fail_msg;

  ImportData import_data;
};

GenericDatumFormatter<std::string, Result> initial_path();

GenericDatumFormatter<std::string, Result> final_path();

GenericDatumFormatter<std::string, Result> fail_msg();

GenericDatumFormatter<std::string, Result> data_origin();

GenericDatumFormatter<std::string, Result> to_configname();

GenericDatumFormatter<bool, Result> has_data();

GenericDatumFormatter<bool, Result> has_complete_data();

GenericDatumFormatter<bool, Result> preexisting_data();

GenericDatumFormatter<bool, Result> import_data();

GenericDatumFormatter<bool, Result> import_additional_files();

GenericDatumFormatter<double, Result> lattice_deformation_cost();

GenericDatumFormatter<double, Result> atomic_deformation_cost();

GenericDatumFormatter<double, Result> energy();

GenericDatumFormatter<double, Result> score();

GenericDatumFormatter<double, Result> best_score();

GenericDatumFormatter<bool, Result> is_best();

GenericDatumFormatter<bool, Result> is_new_config();

GenericDatumFormatter<bool, Result> selected();

/// Insert default formatters to dictionary, for 'casm import'
void default_import_formatters(DataFormatterDictionary<Result> &dict,
                               PropertiesDatabase &db_props);

/// Insert default formatters to dictionary, for 'casm update'
void default_update_formatters(DataFormatterDictionary<Result> &dict,
                               PropertiesDatabase &db_props);

}  // namespace ConfigIO

template <typename T>
struct TypeTag {
  using type = T;
};

class ConfigData {
 public:
  /// \brief Checks if pos_path can be used to resolve a properties.calc.json,
  /// and return its path Return path to properties.calc.json or POSCAR-type
  /// file that will be imported checking a couple possible locations relative
  /// to pos_path
  ///
  /// checks:
  /// 1) is a JSON file? is pos_path ends in ".json" or ".JSON", return pos_path
  /// 2) assume pos_path is /path/to/POS, checks for
  /// /path/to/calctype.current/properties.calc.json 3) assume pos_path is
  /// /path/to/POS, checks for /path/to/properties.calc.json else returns
  /// pos_path

  // Otherwise return pos_path
  static std::string resolve_struc_path(std::string pos_path,
                                        PrimClex const &_pclex);

  template <typename ConfigType>
  ConfigData(const PrimClex &_primclex, TypeTag<ConfigType>);

  const PrimClex &primclex() const { return m_primclex; }

  Database<Supercell> &db_supercell() const;

  template <typename ConfigType>
  Database<ConfigType> &db_config() const;

  /// Uses primclex().settings().default_clex().calctype
  PropertiesDatabase &db_props() const;

  /// \brief Path to default calctype training_data directory for config
  std::string calc_dir(const std::string configname) const;

  /// \brief Return true if there are existing files in the traning_data
  /// directory
  ///        for a particular configuration
  bool has_existing_files(const std::string &to_configname) const;

  /// \brief Return true if there are existing files in the traning_data
  /// directory
  ///        for a particular configuration
  bool has_existing_data(const std::string &to_configname) const;

  bool has_existing_data_or_files(const std::string &to_configname) const;

  /// Check if 'properties.calc.json' file has not changed since last read
  bool no_change(const std::string &configname) const;

  /// \brief Remove existing files in the traning_data directory for a
  /// particular
  ///        configuration
  void rm_files(const std::string &configname, bool dry_run) const;

  /// \brief Copy files in the same directory as properties.calc.json into the
  ///        traning_data directory for a particular configuration
  void cp_files(ConfigIO::Result &res, bool dry_run,
                bool copy_additional_files) const;

 private:
  const PrimClex &m_primclex;
  std::function<PropertiesDatabase &()> m_db_props_func;
};

// To be specialized for ConfigType (no default implemenation exists)
template <typename ConfigType>
class StructureMap;

/*
  template<typename ConfigType>
  class StructureMap {

public:

  typedef std::back_insert_iterator<std::vector<ConfigIO::Result> >
map_result_inserter;

  /// \brief Specialized import method for ConfigType
  ///
  /// \param p Path to structure or properties.calc.json file. Not guaranteed to
exist or be valid.
  /// \param hint Iterator to 'from' config for 'casm update', or 'end' if
unknown as with 'casm import'.
  /// \param result Insert iterator of Result objects to output mapping results
  ///
  /// - Should output one or more mapping results from the structure located at
specied path
  /// - >1 result handles case of non-primitive configurations
  /// - responsible for filling in Result data structure
  /// - If 'hint' is not nullptr, use hint as 'from' config, else 'from' == 'to'
  map_result_inserter map(
    std::string p,
    DatabaseIterator<ConfigType> hint,
    map_result_inserter result) const;
};
*/
}  // namespace DB
}  // namespace CASM

#endif
