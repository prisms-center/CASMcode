#include "casm/database/DiffTransConfigImport.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/kinetics/DiffTransConfigMapping.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace CASM {
  namespace DB {

    // --- DiffTransConfiguration specializations --------------------------------

    StructureMap<Kinetics::DiffTransConfiguration>::StructureMap(
      const PrimClex &primclex,
      std::unique_ptr<Kinetics::DiffTransConfigMapper> mapper) :
      ConfigData<Kinetics::DiffTransConfiguration>(primclex, null_log()),
      m_difftransconfigmapper(std::move(mapper)) {}

    StructureMap<Kinetics::DiffTransConfiguration>::StructureMap(
      const PrimClex &primclex,
      const jsonParser &kwargs) :
      ConfigData<Kinetics::DiffTransConfiguration>(primclex, null_log()) {

      // -- read settings --
      bool rotate = true;
      bool strict = false;
      bool ideal;
      kwargs.get_else(ideal, "ideal", false);

      double lattice_weight;
      kwargs.get_else(lattice_weight, "lattice_weight", 0.5);

      double max_vol_change;
      kwargs.get_else(max_vol_change, "max_vol_change", 0.3);

      double min_va_frac;
      kwargs.get_else(min_va_frac, "min_va_frac", 0.0);

      double max_va_frac;
      kwargs.get_else(max_va_frac, "max_va_frac", 0.5);

      // -- collect settings used --
      m_used.put_obj();
      m_used["ideal"] = ideal;
      m_used["lattice_weight"] = lattice_weight;
      m_used["max_vol_change"] = max_vol_change;
      m_used["min_va_frac"] = min_va_frac;
      m_used["max_va_frac"] = max_va_frac;

      // -- construct ConfigMapper --
      int map_opt = Kinetics::DiffTransConfigMapper::Options::none;
      if(rotate) map_opt |= Kinetics::DiffTransConfigMapper::Options::rotate;
      if(strict) map_opt |= Kinetics::DiffTransConfigMapper::Options::strict;
      if(!ideal) map_opt |= Kinetics::DiffTransConfigMapper::Options::robust;

      m_difftransconfigmapper.reset(new Kinetics::DiffTransConfigMapper(
                                      primclex,
                                      lattice_weight,
                                      max_vol_change,
                                      map_opt,
                                      primclex.crystallography_tol()));
      m_difftransconfigmapper->set_min_va_frac(min_va_frac);
      m_difftransconfigmapper->set_max_va_frac(max_va_frac);

    }

    const jsonParser &StructureMap<Kinetics::DiffTransConfiguration>::used() const {
      return m_used;
    }

    StructureMap<Kinetics::DiffTransConfiguration>::map_result_inserter StructureMap<Kinetics::DiffTransConfiguration>::map(
      fs::path p,
      DatabaseIterator<Kinetics::DiffTransConfiguration> hint,
      map_result_inserter result) const {
      //todo
      fs::path prop_path = this->calc_properties(p);

      ConfigIO::Result res;
      res.pos = (prop_path.empty() ? p : prop_path);

      std::unique_ptr<Kinetics::DiffTransConfiguration> hint_config;
      if(hint != db_difftransconfig().end()) {
        hint_config = notstd::make_unique<DiffTransConfiguration>(*hint);
        res.mapped_props.from = hint_config->name();
      }
      return result;
    }


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

    const std::string Import<Kinetics::DiffTransConfiguration>::desc =

      "Import DiffTransConfiguration: \n\n"

      "  'casm import' of DiffTransConfiguration proceeds in two steps: \n\n"

      "  1) For each set of files: \n"
      "     - Read structures from VASP POSCAR type files or CASM properties.calc.json \n"
      "       file. \n"
      "     - Map the structure onto a DiffTransConfiguration of the primitive crystal \n"
      "       structure. \n"
      "     - Record relaxation data (lattice & basis deformation cost). \n\n"

      "  2) If data import is requested, iterate over each import record and do \n"
      "     the following:\n"
      "     - If data is imported, the corresponding properties.calc.json file is\n"
      "       copied into the directory of the mapped configuration. Optionally, \n"
      "       additional files in the directory of the imported structure file may\n"
      "       also be copied. \n"
      "     - Reports are generated detailing the results of the import: \n"
      "       - import_map_fail: Structures that could not be mapped onto the     \n"
      "         primitive crystal structure. \n"
      "       - import_map_success: Configurations that were successfully mapped  \n"
      "         and imported into the Configuration database (or already existed).\n"
      "       - import_data_fail: Structures with data that would be imported     \n"
      "         except preexisting data prevents it. \n"
      "       - import_conflict: Configurations that were mapped to by multiple   \n"
      "         structures. \n\n"


      ;

    int Import<Kinetics::DiffTransConfiguration>::run(
      const PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::ImportOption &import_opt) {
      //todo
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

    const std::string Update<Kinetics::DiffTransConfiguration>::desc =


      "Update DiffTransConfiguration calculation results: \n\n"

      "  'casm update' of DiffTransConfiguration calculation results proceeds as follows: \n\n"
      "  For each DiffTransConfiguration in the input selection: \n"
      "   - Read properties.calc.json file from training_data directory.        \n"
      "   - Map the relaxed structure onto a DiffTransConfiguration of the primitive crystal\n"
      "     structure. \n"
      "   - Record relaxation data: \n"
      "     - Lattice & basis deformation cost \n"
      "     - Initial difftransconfiguration and relaxed difftransconfiguration \n\n"
      "   - If multiple difftransconfigurations relax onto a difftransconfiguration for which there \n"
      "     is no calculation data, the calculation data from the with the lowest \n"
      "     conflict resolution score is used for the relaxed difftransconfiguration.\n"
      "   - Both default and diiftransconfiguration-specific conflict resolution scoring\n"
      "     method can be set via: \n"
      "       'casm update --set-default-conflict-score -i <JSON>'\n"
      "       'casm update --set-default-conflict-score -s <JSON filename>'\n"
      "       'casm update --set-conflict-score configname -i <JSON>'\n"
      "       'casm update --set-conflict-score configname -s <JSON filename>'\n\n";

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
