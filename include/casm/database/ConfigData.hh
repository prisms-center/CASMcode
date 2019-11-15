#ifndef CASM_DB_ConfigData
#define CASM_DB_ConfigData

#include <boost/filesystem.hpp>
#include "casm/global/definitions.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/dataformatter/DataFormatter.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools.hh"



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

    /// \brief Struct with optional parameters for Config Mapping
    /// Specifies default parameters for all values, in order to simplify
    /// parsing from JSON
    struct MappingSettings {
      MappingSettings(double _lattice_weight = 0.5,
                      bool _ideal = false,
                      bool _strict = false,
                      bool _primitive_only = false,
                      std::vector<std::string> _forced_lattices = {},
                      std::string _filter = "",
                      double _cost_tol = CASM::TOL,
                      double _min_va_frac = 0.,
                      double _max_va_frac = 0.5,
                      double _max_vol_change = 0.3) :
        lattice_weight(_lattice_weight),
        ideal(_ideal),
        strict(_strict),
        primitive_only(_primitive_only),
        forced_lattices(_forced_lattices),
        filter(_filter),
        cost_tol(_cost_tol),
        min_va_frac(_min_va_frac),
        max_va_frac(_max_va_frac),
        max_vol_change(_max_vol_change) {}

      void set_default() {
        *this = MappingSettings();
      }

      double lattice_weight;
      bool ideal;
      bool strict;
      bool primitive_only;
      std::vector<std::string> forced_lattices;
      std::string filter;
      double cost_tol;
      double min_va_frac;
      double max_va_frac;
      double max_vol_change;

    };

    jsonParser &to_json(MappingSettings const &_set, jsonParser &_json);

    jsonParser const &from_json(MappingSettings &_set, jsonParser const &_json);

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
          map_result(""),
          has_data(false),
          has_complete_data(false),
          is_new_config(false),
          fail_msg("") {}

        // structure or properties.calc.json location as input
        fs::path pos;

        // Set 'to'/'from' as empty strings if no mapping possible
        ConfigMapperResult::MapData map_result;

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

    template<typename T>
    struct TypeTag {
      using type = T;
    };


    class ConfigData {
    public:
      /// \brief Return path to properties.calc.json that will be imported
      ///        checking a couple possible locations relative to pos_path
      ///
      /// checks:
      /// 1) is a JSON file? is pos_path ends in ".json" or ".JSON", return pos_path
      /// 2) assume pos_path is /path/to/POS, checks for /path/to/calctype.current/properties.calc.json
      /// 3) assume pos_path is /path/to/POS, checks for /path/to/properties.calc.json
      /// else returns empty path
      ///
      static fs::path calc_properties_path(fs::path pos_path, PrimClex const &_pclex);

      // If pos_path can be used to resolve a properties.calc.json, return its path.
      // Otherwise return pos_path
      static fs::path resolve_struc_path(fs::path pos_path, PrimClex const &_pclex);

      template<typename ConfigType>
      ConfigData(const PrimClex &_primclex, Log &_file_log, TypeTag<ConfigType>);

      const PrimClex &primclex() const {
        return m_primclex;
      }

      Log &file_log() const {
        return m_file_log;
      }

      Database<Supercell> &db_supercell() const;

      template<typename ConfigType>
      Database<ConfigType> &db_config() const;

      /// Uses primclex().settings().default_clex().calctype
      PropertiesDatabase &db_props() const;

      /// \brief Path to default calctype training_data directory for config
      fs::path calc_dir(const std::string configname) const;

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
      std::function<PropertiesDatabase&()> m_db_props_func;
    };


    // To be specialized for ConfigType (no default implemenation exists)
    template<typename ConfigType>
    class StructureMap;
    /*
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
    */
  }
}

#endif
