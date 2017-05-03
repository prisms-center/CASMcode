#ifndef CASM_DB_Import
#define CASM_DB_Import

#include "casm/database/PropertiesDatabase.hh"
#include "casm/app/casm_functions.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"

namespace CASM {
  namespace Completer {
    class ImportOption;
    class UpdateOption;
    class RmOption;
  }
}

namespace CASM {
  class Supercell;
}

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
  namespace DB {

    template<typename T> class Selection;

    /// Create a new report directory to avoid overwriting existing results
    fs::path create_report_dir(fs::path report_dir);

    /// Construct pos_paths from input args --pos && --batch
    template<typename OutputIterator>
    std::pair<OutputIterator, int> construct_pos_paths(
      const PrimClex &primclex,
      const Completer::ImportOption &import_opt,
      OutputIterator result);

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

      /// Insert default formatters to dictionary
      void default_formatters(
        DataFormatterDictionary<Result> &dict,
        PropertiesDatabase &db_props,
        const std::map<std::string, ImportData> &data_results);

    }


    class ConfigDataGeneric {

    public:

      ConfigDataGeneric(const PrimClex &_primclex, Log &_file_log);

      const PrimClex &primclex() const {
        return m_primclex;
      }

      Log &file_log() const {
        return m_file_log;
      }

    protected:

      Database<Supercell> &db_supercell() const;

      virtual PropertiesDatabase &db_props() const = 0;

      /// \brief Path to default calctype training_data directory for config
      virtual fs::path calc_dir(const std::string configname) const = 0;

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
      ///
      /// - First: calc_props_path = _calc_properties_path(pos_path) to get properties.calc.json location
      /// - If calc_props_path.empty(), return
      /// - else if !copy_additional_files copy properties.calc.json file only and return
      /// - else, recursively copy all files from calc_props_path.remove_filename()
      ///   to the training data directory for the current calctype
      void cp_files(
        const fs::path &pos_path,
        const std::string &configname,
        bool dry_run,
        bool copy_additional_files) const;


    protected:

      void _import_report(
        std::vector<ConfigIO::Result> &results,
        const std::map<std::string, ConfigIO::ImportData> &data_results);

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

      /// \brief Path to default calctype training_data directory for config
      fs::path calc_dir(const std::string configname) const override;
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


    /// Generic ConfigType-dependent part of Import
    template<typename _ConfigType>
    class ImportT : protected ConfigData<_ConfigType> {

    public:

      typedef _ConfigType ConfigType;

      /// \brief Constructor
      ImportT(
        const PrimClex &primclex,
        const StructureMap<ConfigType> &mapper,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir,
        Log &_file_log) :
        ConfigData<_ConfigType>(primclex, _file_log),
        m_structure_mapper(mapper),
        m_import_data(import_data),
        m_copy_additional_files(copy_additional_files),
        m_overwrite(overwrite),
        m_report_dir(report_dir) {}

      template<typename PathIterator>
      void import(PathIterator begin, PathIterator end);

    private:

      /// Allow ConfigType to specialize the report formatting for 'import'
      virtual DataFormatter<ConfigIO::Result> _import_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const = 0;


      const StructureMap<ConfigType> &m_structure_mapper;

      // attempt to import calculation results into database, else just insert
      // configurations w/out data
      bool m_import_data;

      // attempt to copy extra files from the directory where the structure is
      // being imported from to the training_data directory
      bool m_copy_additional_files;

      // Allow overwriting of existing data by 'casm import'
      bool m_overwrite;

      fs::path m_report_dir;
    };


    /// Generic ConfigType-dependent part of Import
    template<typename _ConfigType>
    class UpdateT : protected ConfigData<_ConfigType> {

    public:

      typedef _ConfigType ConfigType;

      /// \brief Constructor
      UpdateT(
        const PrimClex &primclex,
        const StructureMap<ConfigType> &mapper,
        fs::path report_dir) :
        ConfigData<_ConfigType>(primclex, null_log()),
        m_structure_mapper(mapper),
        m_report_dir(report_dir) {}

      /// \brief Re-parse calculations 'from' all selected configurations
      void update(const DB::Selection<ConfigType> &selection, bool force);

    protected:

      // Allow ConfigType to specialize the report formatting for 'update'
      virtual DataFormatter<ConfigIO::Result> _update_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const = 0;

      void _update_report(std::vector<ConfigIO::Result> &results, const DB::Selection<ConfigType> &selection) const;

    private:

      const StructureMap<ConfigType> &m_structure_mapper;

      fs::path m_report_dir;
    };


    /// Generic ConfigType-dependent part of Remove
    template<typename _ConfigType>
    class RemoveT : protected ConfigData<_ConfigType> {

    public:

      typedef _ConfigType ConfigType;

      RemoveT(const PrimClex &primclex, fs::path report_dir, Log &_file_log) :
        ConfigData<_ConfigType>(primclex, _file_log),
        m_report_dir(report_dir) {}

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

      void _erase_fail_report(std::vector<std::string> fail);

      fs::path m_report_dir;
    };

    // To be specialized for ConfigType (no default implemenation exists)
    template<typename ConfigType>
    class Import : public ImportT<ConfigType> {
    public:
      static const std::string desc;
      int run(const PrimClex &, const jsonParser &input, const Completer::ImportOption &opt);

    private:
      // Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<ConfigIO::Result> _import_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const override;
    };

    // To be specialized for ConfigType (no default implemenation exists)
    template<typename ConfigType>
    class Update : public UpdateT<ConfigType> {
    public:
      static const std::string desc;
      int run(const PrimClex &, const jsonParser &input, const Completer::UpdateOption &opt);

    private:
      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<ConfigIO::Result> _update_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const override;
    };

    // To be specialized for ConfigType (no default implemenation exists)
    template<typename ConfigType>
    class Remove : public RemoveT<ConfigType> {
    public:
      static const std::string desc;
      int run(const PrimClex &, const Completer::RmOption &opt);
    };
  }
}

#endif
