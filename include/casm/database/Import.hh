#ifndef CASM_Import
#define CASM_Import

#include "casm/database/PropertiesDatabase.hh"

namespace CASM {

  namespace DB {

    /// Import/update/erase configurations and data, non-templated components
    class ImportBase : public Logging {

    public:

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

      typedef std::back_insert_iterator<std::vector<Result> > import_inserter;

      /// \brief Constructor
      ImportBase(
        const PrimClex &primclex,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir);

      /// Construct pos_paths from input args
      template<typename OutputIterator>
      static std::pair<OutputIterator, int> construct_pos_paths(
        const PrimClex &primclex,
        const Completer::ImportOption &import_opt,
        OutputIterator result);

      /// Create a new report directory to avoid overwriting existing results
      static fs::path create_report_dir(fs::path report_dir);


    private:

      DataBase<Supercell> &db_supercell() {
        return m_primclex.db<Supercell>();
      }

      /// \brief Return path to properties.calc.json that will be imported
      ///        checking a couple possible locations relative to pos_path
      ///
      /// checks:
      /// 1) is a JSON file? is pos_path ends in ".json" or ".JSON", return pos_path
      /// 2) assume pos_path is /path/to/POS, checks for /path/to/calctype.current/properties.calc.json
      /// 3) assume pos_path is /path/to/POS, checks for /path/to/properties.calc.json
      /// else returns empty path
      ///
      fs::path _calc_properties_path(fs::path pos_path) const;

      /// \brief Return true if there are existing files in the traning_data directory
      ///        for a particular configuration
      bool _has_existing_files(const std::string &from_configname) const;

      /// \brief Return true if there are existing files in the traning_data directory
      ///        for a particular configuration
      bool _has_existing_data(const std::string &from_configname) const;

      bool _has_existing_data_or_files(const std::string &from_configname) const;

      /// Check if 'properties.calc.json' file has not changed since last read
      bool _no_change(const std::string &configname) const;

      /// \brief Remove existing files in the traning_data directory for a particular
      ///        configuration
      void _rm_files(const std::string &configname) const;

      void _recurs_rm_files(const fs::path &p) const;

      /// \brief Copy files in the same directory as properties.calc.json into the
      ///        traning_data directory for a particular configuration
      ///
      /// - First: calc_props_path = _calc_properties_path(pos_path) to get properties.calc.json location
      /// - If calc_props_path.empty(), return
      /// - else if !m_copy_additional_files copy properties.calc.json file only and return
      /// - else, recursively copy all files from calc_props_path.remove_filename()
      ///   to the training data directory for the current calctype
      void _cp_files(const fs::path &pos_path, const std::string &configname) const;

      void _recurs_cp_files(const fs::path &from_dir, const fs::path &to_dir) const;

      void _import_report(std::vector<Result> &results);


      GenericDatumFormatter<std::string, Result> _pos() const;

      GenericDatumFormatter<std::string, Result> _fail_msg() const;

      /// Use 'from_configname' as 'configname'
      GenericDatumFormatter<std::string, Result> _configname() const;

      GenericDatumFormatter<std::string, Result> _from_configname() const;

      GenericDatumFormatter<std::string, Result> _to_configname() const;

      GenericDatumFormatter<bool, Result> _has_data() const;

      GenericDatumFormatter<bool, Result> _has_complete_data() const;

      GenericDatumFormatter<bool, Result> _preexisting_data(const std::map<std::string, ImportData> &data_results) const;

      GenericDatumFormatter<bool, Result> _import_data(const std::map<std::string, ImportData> &data_results) const;

      GenericDatumFormatter<bool, Result> _import_additional_files(const std::map<std::string, ImportData> &data_results) const;

      GenericDatumFormatter<double, Result> _lattice_deformation_cost() const;

      GenericDatumFormatter<double, Result> _basis_deformation_cost() const;

      GenericDatumFormatter<double, Result> _relaxed_energy() const;

      GenericDatumFormatter<double, Result> _score() const;

      GenericDatumFormatter<double, Result> _best_score() const;

      GenericDatumFormatter<std::string, Result> _selected() const;

      /// Insert default formatters to dictionary
      void _default_formatters(
        DataFormatterDictionary<Result> &dict,
        std::map<std::string, ImportResult> &data_results) const;


      // --- member variables ----------------

      const PrimClex &m_primclex;

      // path to directory where files are written
      fs::path m_report_dir;

      // log file to show file rm and cp
      fs::ofstream m_file_str;

      // formatting m_file_str;
      Log m_file_log;

      // attempt to import calculation results into database, else just insert
      // configurations w/out data
      bool m_import_data;

      // attempt to copy extra files from the directory where the structure is
      // being imported from to the training_data directory
      bool m_copy_additional_files;

      // Allow overwriting of existing data by 'casm import'
      bool m_overwrite;
    };


    /// Generic ConfigType-dependent part of Import
    template<typename _ConfigType>
    class ImportT: public ImportBase {

    public:

      typedef _ConfigType ConfigType;

      /// \brief Constructor
      ImportT(
        const PrimClex &primclex,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir) :
        Import(primclex, import_data, copy_additional_files, overwrite, report_dir) {}

      template<typename PathIterator>
      void import(PathIterator begin, PathIterator end);

      /// \brief Re-parse calculations 'from' all selected configurations
      void update(const DB::Selection<ConfigType> &selection, bool force);

      /// \brief Erase Configurations that have no data
      void erase(const DB::Selection<ConfigType> &selection);

      /// \brief Erase data and files (permanently), but not Configuration
      void erase_data(const DB::Selection<ConfigType> &selection);

      /// \brief Removes Configurations and data and files (permanently)
      ///
      /// - Data are always associated with one 'from' configuration, so the
      ///   selection here indicates 'from' configurations
      /// - The 'to' configurations are updated with the new best mapping properties
      void erase_all(const DB::Selection<ConfigType> &selection);

    private:

      DataBase<ConfigType> &db_config() {
        return m_primclex.db<ConfigType>();
      }

      PropertiesDatabase &db_props() {
        return m_primclex.db_props<ConfigType>();
      }

      /// \brief Path to default calctype training_data directory for config
      fs::path calc_dir(const std::string configname) const;

      void _update_report(std::vector<Result> &results, const DB::Selection<ConfigType> &selection) const;

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
      virtual import_inserter _import(
        fs::path p,
        DataBaseIterator<ConfigType> hint,
        import_inserter result) = 0;

      /// Allow ConfigType to specialize the report formatting for 'import'
      virtual DataFormatter<Result> _import_formatter(
        std::map<std::string, ImportResult> &data_results) const = 0;

      // Allow ConfigType to specialize the report formatting for 'update'
      virtual DataFormatter<Result> _update_formatter(
        std::map<std::string, ImportResult> &data_results) const = 0;
    };


    template<typename ConfigType>
    class Import {}

    /// \brief Access Import<ConfigType>::import via the API
    template<typename ConfigType>
    class ImportInterface : public InterfaceBase<Completer::ImportOption> {

    public:

      std::string help() const override {
        return Import<ConfigType>::import_help;
      }

      std::string name() const override {
        return QueryTraits<ConfigType>::short_name;
      }

      int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::ImportOption &import_opt) const override {
        return Import<ConfigType>::import(primclex, kwargs, import_opt);
      }

      std::unique_ptr<InterfaceBase<Completer::ImportOption> > clone() const {
        return std::unique_ptr<InterfaceBase<Completer::ImportOption> >(this->_clone());
      }

    private:

      ImportInterface<ConfigType> *_clone() const override {
        return new ImportInterface<ConfigType>(*this);
      }
    };

    /// \brief Access Import<ConfigType>::update via the API
    template<typename ConfigType>
    class UpdateInterface : public InterfaceBase<Completer::UpdateOption> {

    public:

      std::string help() const override {
        return Import<ConfigType>::update_help;
      }

      std::string name() const override {
        return QueryTraits<ConfigType>::short_name;
      }

      int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &update_opt) const override {
        return Import<ConfigType>::update(primclex, kwargs, update_opt);
      }

      std::unique_ptr<InterfaceBase<Completer::UpdateOption> > clone() const {
        return std::unique_ptr<InterfaceBase<Completer::UpdateOption> >(this->_clone());
      }

    private:

      UpdateInterface<ConfigType> *_clone() const override {
        return new UpdateInterface<ConfigType>(*this);
      }
    };

    /*
    /// \brief Access Import<ConfigType>::remove via the API
    template<typename ConfigType>
    class RemoveInterface : public InterfaceBase<Completer::RemoveOption> {

    public:

      std::string help() const override {
        return Import<ConfigType>::remove_help;
      }

      std::string name() const override {
        return QueryTraits<ConfigType>::short_name;
      }

      int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::RemoveOption &remove_opt) const override {
        return Import<ConfigType>::remove(primclex, kwargs, update_opt);
      }

      std::unique_ptr<InterfaceBase<Completer::RemoveOption> > clone() const {
        return std::unique_ptr<InterfaceBase<Completer::RemoveOption> >(this->_clone());
      }

    private:

      RemoveInterface<ConfigType> *_clone() const override {
        return new RemoveInterface<ConfigType>(*this);
      }
    };
    */


    /// Configuration-specialized Import
    template<>
    class Import<Configuration> : public ImportT<Configuration> {
    public:

      typedef ConfigMapper MapperType;

      /// Construct ConfigMapper from input args
      std::unique_ptr<ConfigMapper> make(
        const PrimClex &primclex,
        const Completer::ImportOption &import_opt);

      /// \brief Constructor
      Import(
        const PrimClex &primclex,
        const ConfigMapper &configmapper,
        std::vector<std::string> dof,
        bool primitive_only,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir);

      static const std::string import_help;
      static int import(PrimClex &primclex, const jsonParser &kwargs, const Completer::ImportOption &import_opt);

      static const std::string update_help;
      static int update(PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt);

      static const std::string remove_help;
      static int remove(PrimClex &primclex, const jsonParser &kwargs, const Completer::RemoveOption &import_opt);

    private:

      /// Construct ConfigMapper from input args
      ConfigMapper _make_configmapper(const PrimClex &primclex, const jsonParser &kwargs);

      /// \brief Specialized import method for ConfigType
      import_inserter _import(
        fs::path p,
        DataBaseIterator<Configuration> hint,
        import_inserter result) override;

      /// Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<Result> _import_formatter(
        std::map<std::string, ImportResult> &data_results) const override;

      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<Result> _update_formatter(
        std::map<std::string, ImportResult> &data_results) const override;

      /// \brief Read BasicStructure<Site> to be imported
      BasicStructure<Site> _make_structure(const fs::path &p);

      /// \brief Import Configuration with only occupation DoF
      bool _occupation_only() const;


      const MapperType &m_configmapper;
      std::vector<std::string> m_dof;
      bool m_primitive_only;
    };

    /*
    /// DiffTransConfiguration-specialized Import
    template<>
    class Import<DiffTransConfiguration> : public ImportT<DiffTransConfiguration> {
    public:

      /// \brief Constructor
      Import(
        const PrimClex& primclex,
        const DiffTransConfigMapper& configmapper,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir = primclex.root_dir() / "import_report");

      static const std::string import_help;
      static int import(PrimClex &primclex, const jsonParser &kwargs, const Completer::ImportOption &import_opt);

      static const std::string update_help;
      static int update(PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt);

      static const std::string remove_help;
      static int remove(PrimClex &primclex, const jsonParser &kwargs, const Completer::RemoveOption &import_opt);

    private:

      /// \brief Specialized import method for ConfigType
      import_inserter _import(
        fs::path p,
        DataBaseIterator<DiffTransConfiguration> hint,
        import_inserter result) override;

      /// Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<Result> _import_formatter(
        std::map<std::string, ImportResult>& data_results) const override;

      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<Result> _update_formatter(
        std::map<std::string, ImportResult>& data_results) const override;

    private:
      const DiffTransConfigMapper& m_configmapper;

    };
    */

  }

  namespace Completer {
    class ImportOption;
  }

  typedef InterfaceBase<Completer::ImportOption> ImportInterfaceBase;
  typedef InterfaceMap<Completer::ImportOption> ImporterMap;


}
