#ifndef CASM_DB_ConfigData
#define CASM_DB_ConfigData

#include <boost/filesystem.hpp>
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/database/MappedProperties.hh"


// To be specialized for calculable 'ConfigType' classes:
//   StructureMap<ConfigType>::map(fs::path, DatabaseIterator<ConfigType> hint, inserter result)
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

    template<typename T> class Selection;
    template<typename ValueType> class Database;
    template<typename ValueType> class DatabaseIterator;
    class PropertiesDatabase;

    /// Create a new report directory to avoid overwriting existing results
    fs::path create_report_dir(fs::path report_dir);

    namespace ConfigIO {

      /// Data structure for mapping / import results
      struct Result {

        Result() :
          pos(""),
          has_data(false),
          has_complete_data(false),
          is_new_config(false),
          fail_msg("") {}


        // structure or properties.calc.json location as input
        fs::path pos;

        // Set 'to'/'from' as empty strings if no mapping possible
        MappedProperties mapped_props;

        // If a properties.calc.json file is found in standard locations
        bool has_data;

        // If all PrimClex required properties exist
        bool has_complete_data;

        // If 'to' Configuration did not exist in the database prior to mapping
        bool is_new_config;

        // If mapping failed, stores any error message that may be generated
        std::string fail_msg;
      };

      /// Data structure for data import results
      struct ImportData {

        ImportData() :
          preexisting(false),
          copy_data(false),
          copy_more(false),
          last_insert("") {}

        // base responsibility:
        bool preexisting;
        bool copy_data;
        bool copy_more;
        fs::path last_insert;
      };

      GenericDatumFormatter<std::string, Result> pos();

      GenericDatumFormatter<std::string, Result> fail_msg();

      /// Use 'from_configname' as 'configname'
      GenericDatumFormatter<std::string, Result> configname();

      GenericDatumFormatter<std::string, Result> from_configname();

      GenericDatumFormatter<std::string, Result> to_configname();

      GenericDatumFormatter<bool, Result> has_data();

      GenericDatumFormatter<bool, Result> has_complete_data();

      GenericDatumFormatter<bool, Result> preexisting_data(const std::map<std::string, ImportData> &data_results);

      GenericDatumFormatter<bool, Result> import_data(const std::map<std::string, ImportData> &data_results);

      GenericDatumFormatter<bool, Result> import_additional_files(const std::map<std::string, ImportData> &data_results);

      GenericDatumFormatter<double, Result> lattice_deformation_cost();

      GenericDatumFormatter<double, Result> basis_deformation_cost();

      GenericDatumFormatter<double, Result> relaxed_energy();

      GenericDatumFormatter<double, Result> score();

      GenericDatumFormatter<double, Result> best_score();

      GenericDatumFormatter<bool, Result> is_best();

      GenericDatumFormatter<bool, Result> selected();

      /// Insert default formatters to dictionary, for 'casm import'
      void default_import_formatters(
        DataFormatterDictionary<Result> &dict,
        PropertiesDatabase &db_props,
        const std::map<std::string, ImportData> &data_results);

      /// Insert default formatters to dictionary, for 'casm update'
      void default_update_formatters(
        DataFormatterDictionary<Result> &dict,
        PropertiesDatabase &db_props);

    }


    class ConfigDataGeneric : public Logging {

    public:

      ConfigDataGeneric(const PrimClex &_primclex, Log &_file_log);

      const PrimClex &primclex() const {
        return m_primclex;
      }

      Log &file_log() const {
        return m_file_log;
      }

      Database<Supercell> &db_supercell() const;

      /// Uses primclex().settings().default_clex().calctype
      virtual PropertiesDatabase &db_props() const = 0;


      /// \brief Path to default calctype training_data directory for config
      fs::path calc_dir(const std::string configname) const;

      /// \brief Return path to properties.calc.json that will be imported
      ///        checking a couple possible locations relative to pos_path
      ///
      /// checks:
      /// 1) is a JSON file? is pos_path ends in ".json" or ".JSON", return pos_path
      /// 2) assume pos_path is /path/to/POS, checks for /path/to/calctype.current/properties.calc.json
      /// 3) assume pos_path is /path/to/POS, checks for /path/to/properties.calc.json
      /// else returns empty path
      ///
      fs::path calc_properties_path(fs::path pos_path) const;

      /// \brief Return true if there are existing files in the traning_data directory
      ///        for a particular configuration
      bool has_existing_files(const std::string &from_configname) const;

      /// \brief Return true if there are existing files in the traning_data directory
      ///        for a particular configuration
      bool has_existing_data(const std::string &from_configname) const;

      bool has_existing_data_or_files(const std::string &from_configname) const;

      /// Check if 'properties.calc.json' file has not changed since last read
      bool no_change(const std::string &configname) const;

      /// \brief Remove existing files in the traning_data directory for a particular
      ///        configuration
      void rm_files(const std::string &configname, bool dry_run) const;

      /// \brief Copy files in the same directory as properties.calc.json into the
      ///        traning_data directory for a particular configuration
      std::pair<bool, bool> cp_files(
        const fs::path &pos_path,
        const std::string &configname,
        bool dry_run,
        bool copy_additional_files) const;


    private:
      const PrimClex &m_primclex;
      Log &m_file_log;
    };


    template<typename _ConfigType>
    class ConfigData : public ConfigDataGeneric {
    public:
      typedef _ConfigType ConfigType;

      ConfigData(const PrimClex &_primclex, Log &_file_log) :
        ConfigDataGeneric(_primclex, _file_log) {}

      Database<ConfigType> &db_config() const;

      PropertiesDatabase &db_props() const override;

    };


    // To be specialized for ConfigType (no default implemenation exists)
    template<typename ConfigType>
    class StructureMap {

    public:

      typedef std::back_insert_iterator<std::vector<ConfigIO::Result> > map_result_inserter;

      /// \brief Specialized import method for ConfigType
      ///
      /// \param p Path to structure or properties.calc.json file. Not guaranteed to exist or be valid.
      /// \param hint Iterator to 'from' config for 'casm update', or 'end' if unknown as with 'casm import'.
      /// \param result Insert iterator of Result objects to output mapping results
      ///
      /// - Should output one or more mapping results from the structure located at specied path
      /// - >1 result handles case of non-primitive configurations
      /// - responsible for filling in Result data structure
      /// - If 'hint' is not nullptr, use hint as 'from' config, else 'from' == 'to'
      map_result_inserter map(
        fs::path p,
        DatabaseIterator<ConfigType> hint,
        map_result_inserter result) const;
    };

  }
}

#endif
