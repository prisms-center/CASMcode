#include "casm/database/DiffTransConfigImport.hh"

namespace CASM {
  namespace DB {

    // --- DiffTransConfiguration specializations --------------------------------


    /// \brief Constructor
    Import<Kinetics::DiffTransConfiguration>::Import(
      const PrimClex &primclex,
      const StructureMap<Kinetics::DiffTransConfiguration> &mapper,
      bool import_data,
      bool copy_additional_files,
      bool overwrite,
      fs::path report_dir,
      Log &file_log) :

      ImportT(primclex, mapper, import_data, copy_additional_files, overwrite, report_dir, file_log) {}

    const std::string Import<Kinetics::DiffTransConfiguration>::desc = "ToDo";

    int Import<Kinetics::DiffTransConfiguration>::run(
      const PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::ImportOption &import_opt) {
      return 0;
    }

    DataFormatter<ConfigIO::Result> Import<Kinetics::DiffTransConfiguration>::_import_formatter(
      const std::map<std::string, ConfigIO::ImportData> &data_results) const {

      // todo

      DataFormatterDictionary<ConfigIO::Result> dict;
      ConfigIO::default_import_formatters(dict, db_props(), data_results);



      std::vector<std::string> col = {
        "configname", "selected", "pos", "has_data", "has_complete_data",
        "import_data", "import_additional_files", "score", "best_score"
      };

      return dict.parse(col);
    }

    /// \brief Constructor
    Update<Kinetics::DiffTransConfiguration>::Update(
      const PrimClex &primclex,
      const StructureMap<Kinetics::DiffTransConfiguration> &mapper,
      fs::path report_dir) :
      UpdateT(primclex, mapper, report_dir) {}

    const std::string Update<Kinetics::DiffTransConfiguration>::desc = "ToDo";

    int Update<Kinetics::DiffTransConfiguration>::run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt) {
      return 0;
    }

    DataFormatter<ConfigIO::Result> Update<Kinetics::DiffTransConfiguration>::_update_formatter() const {

      // todo

      DataFormatterDictionary<ConfigIO::Result> dict;
      ConfigIO::default_update_formatters(dict, db_props());

      std::vector<std::string> col = {
        "configname", "selected", "to_configname", "has_data", "has_complete_data",
        "score", "best_score", "is_best",
        "lattice_deformation_cost", "basis_deformation_cost", "deformation_cost",
        "relaxed_energy"
      };

      return dict.parse(col);
    }

    /*
    int Import<DiffTransConfiguration>::import(
      PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::ImportOption &import_opt) {

    }

    const std::string Import<DiffTransConfiguration>::update_desc = "ToDo";

    int Import<DiffTransConfiguration>::update(
      PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::UpdateOption &update_opt) {

    }

    const std::string Import<DiffTransConfiguration>::remove_desc = "ToDo";

    int Import<DiffTransConfiguration>::remove(
      PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::RemoveOption &remove_opt) {

    }

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
    import_inserter Import<DiffTransConfiguration>::_import(
      fs::path p,
      DataBaseIterator<DiffTransConfiguration> hint,
      import_inserter result) override {

      // todo

      return import_inserter;
    }

    /// Allow ConfigType to specialize the report formatting for 'import'
    DataFormatter<Result> _import_formatter(
      const std::map<std::string, ImportData>& data_results) const {

      // todo

      DataFormatterDictionary<Result> dict;
      _default_formatters(dict, data_results);

      std::vector<std::string> col = {
        "configname", "selected", "pos", "has_data", "has_complete_data",
        "import_data", "import_additional_files", "score", "best_score"};

      return m_dict.parse(col);
    }

    // Allow ConfigType to specialize the report formatting for 'update'
    DataFormatter<Result> Import<DiffTransConfiguration>::_update_formatter(
      const std::map<std::string, ImportData>& data_results) const {

      // todo

      DataFormatterDictionary<Result> dict;
      _default_formatters(dict, data_results);

      std::vector<std::string> col = {
        "configname", "selected", "to_configname", "has_data", "has_complete_data",
        "score", "best_score"};

      return m_dict.parse(col);
    }
    */

  }
}
