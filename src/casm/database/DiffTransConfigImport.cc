#include "casm/database/DiffTransConfigImport.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/kinetics/DiffTransConfigMapping.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/import.hh"
#include "casm/app/update.hh"
#include "casm/app/rm.hh"
#include "casm/database/DiffTransConfigDatabase.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/Import_impl.hh"
#include "casm/database/Update_impl.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/casm_io/DataFormatter_impl.hh"
#include "casm/basis_set/DoF.hh"


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

      //These are names of the lattices you want to restrict searches to
      std::vector<std::string> forced_lattice_names;
      kwargs.get_else(forced_lattice_names, "forced_lattices", std::vector<std::string>());

      bool lattices_forced_in_settings = (forced_lattice_names.size() > 0);
      if(lattices_forced_in_settings) {
        m_used["forced_lattices"] = forced_lattice_names;
      }
      bool restricted;
      kwargs.get_else(restricted, "restricted", false);
      if(restricted) {
        m_used["restricted"] = restricted;
      }
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

      //If the settings specified at least one lattice, then force that on the configmapper
      if(lattices_forced_in_settings) {
        m_difftransconfigmapper->set_forced_lattices(forced_lattice_names);
      }

      //If the settings specified use boxiness, then force that on the configmapper
      if(restricted) {
        m_difftransconfigmapper->restricted();
      }
    }

    const jsonParser &StructureMap<Kinetics::DiffTransConfiguration>::used() const {
      return m_used;
    }

    StructureMap<Kinetics::DiffTransConfiguration>::map_result_inserter StructureMap<Kinetics::DiffTransConfiguration>::map(
      fs::path p,
      DatabaseIterator<Kinetics::DiffTransConfiguration> hint,
      map_result_inserter result) const {
      //todo
      fs::path prop_path = this->calc_properties_path(p);

      ConfigIO::Result res;
      res.pos = (prop_path.empty() ? p : prop_path);

      std::unique_ptr<Kinetics::DiffTransConfiguration> hint_config;
      if(hint != db_config().end()) {
        hint_config = notstd::make_unique<Kinetics::DiffTransConfiguration>(*hint);
        res.mapped_props.from = hint_config->name();
      }

      // do mapping
      DiffTransConfigMapperResult map_result;
      map_result = m_difftransconfigmapper->import_structure_occupation(res.pos, hint_config.get());
      if(!map_result.success) {
        res.fail_msg = map_result.fail_msg;
        *result++ = res;
        return result;
      }
      // if the result was a success, need to populate proper fields in
      // map_result.relaxation_properties["best_mapping"]
      if(res.pos.extension() == ".json" || res.pos.extension() == ".JSON") {
        jsonParser json(res.pos);
        if(json.contains("kra")) {
          map_result.kra = json["kra"].get<double>();
        }
        //Maybe move kra calculation from DiffTransConfigMapping to here
      }
      // insert in database (note that this also/only inserts primitive)
      if(!map_result.config->has_valid_from_occ()) {
        throw std::runtime_error("You forgot to set m_from_config_is_A");

      }
      Kinetics::DiffTransConfigInsertResult insert_result = map_result.config->insert();
      res.is_new_config = insert_result.insert_canonical;
      res.mapped_props.to = insert_result.canonical_it->name();
      // read in raw unmapped data and put into res
      if(!prop_path.empty()) {
        std::tie(res.mapped_props.unmapped, res.has_data, res.has_complete_data) = read_calc_properties<Kinetics::DiffTransConfiguration>(primclex(), prop_path);
      }



      // copy relaxation properties from best config mapping into 'mapped' props
      res.mapped_props.mapped[insert_result.canonical_it.name()] = map_result.relaxation_properties;
      //These two aren't really being used yet
      //res.mapped_props.mapped["cart_op"] = map_result.cart_op;
      //This is a hack right now because default conflict score looks for minimum
      // relaxed_energy which doesn't make sense for diff_trans_config
      res.mapped_props.mapped["relaxed_energy"] = map_result.kra;
      res.mapped_props.mapped["kra"] = map_result.kra;



      *result++ = res;

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

      "  'casm import' of DiffTransDiffTransConfiguration proceeds in two steps: \n\n"

      "  1) For each set of files: \n"
      "     - Read structures from folders VASP POSCAR type files or CASM properties.calc.json \n"
      "       file. \n"
      "	    - Note: for the paths for the --pos or --batch options please give the path to the\n"
      "	      parent folder containing image numbers. i.e. POSCARS for each endpoint are stored\n"
      "	      poscars/ <--give path to here"
      "		00/"
      " 	  POSCAR"
      "		01/"
      " 	  POSCAR"
      "		02/"
      " 	  POSCAR"
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
      // -- collect input settings --

      const po::variables_map &vm = import_opt.vm();
      jsonParser used;
      jsonParser _default;


      // get input report_dir, check if exists, and create new report_dir.i if necessary
      fs::path report_dir = primclex.dir().reports_dir() / "import_report";
      report_dir = create_report_dir(report_dir);

      // 'mapping' subsettings are used to construct ConfigMapper, and also returns
      // the 'used' settings
      jsonParser map_json;
      kwargs.get_else(map_json, "mapping", jsonParser());
      StructureMap<Kinetics::DiffTransConfiguration> mapper(primclex, map_json);
      used["mapping"] = mapper.used();

      // 'data' subsettings
      jsonParser data;
      kwargs.get_if(data, "data");

      bool import_data;
      if(vm.count("data")) {
        import_data = true;
      }
      else {
        data.get_else(import_data, "import", false);
      }
      used["data"]["import"] = import_data;

      bool copy_additional_files;
      data.get_else(copy_additional_files, "copy_additional_files", false);
      used["data"]["copy_additional_files"] = copy_additional_files;

      bool overwrite;
      data.get_else(overwrite, "overwrite", false);
      used["data"]["overwrite"] = overwrite;

      // -- print used settings --
      Log &log = primclex.log();
      log.read("Settings");
      log << used << std::endl << std::endl;

      // -- construct Import --
      Import<Kinetics::DiffTransConfiguration> f(
        primclex,
        mapper,
        import_data,
        copy_additional_files,
        overwrite,
        report_dir,
        primclex.log());

      // -- read structure file paths --
      std::set<fs::path> pos;
      auto res = construct_pos_paths(primclex, import_opt, std::inserter(pos, pos.end()));
      if(res.second) {
        return res.second;
      }

      // -- read structure file paths --
      f.import(pos.begin(), pos.end());

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

    int Update<Kinetics::DiffTransConfiguration>::run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &update_opt) {
      // -- collect input settings --

      const po::variables_map &vm = update_opt.vm();
      jsonParser used;

      bool force;
      if(vm.count("force")) {
        force = true;
      }
      else {
        kwargs.get_else(force, "force", false);
      }
      used["force"] = force;

      // get input report_dir, check if exists, and create new report_dir.i if necessary
      fs::path report_dir = primclex.dir().reports_dir() / "update_report";
      report_dir = create_report_dir(report_dir);

      // 'mapping' subsettings are used to construct ConfigMapper and return 'used' settings values
      // still need to figure out how to specify this in general
      jsonParser map_json;
      kwargs.get_else(map_json, "mapping", jsonParser());
      StructureMap<Kinetics::DiffTransConfiguration> mapper(primclex, map_json);
      used["mapping"] = mapper.used();

      // 'data' subsettings
      // bool import_data = true;
      // bool import_additional_files = false;
      // bool overwrite = false;

      // -- print used settings --
      Log &log = primclex.log();
      log.read("Settings");
      log << used << std::endl << std::endl;

      // -- construct Update --
      Update<Kinetics::DiffTransConfiguration> f(
        primclex,
        mapper,
        report_dir);

      // -- read selection --
      DB::Selection<Kinetics::DiffTransConfiguration> sel(primclex, update_opt.selection_path());

      // -- update --
      f.update(sel, force);
      return 0;
    }

    DataFormatter<ConfigIO::Result> Update<Kinetics::DiffTransConfiguration>::_update_formatter() const {

      // todo

      DataFormatterDictionary<ConfigIO::Result> dict;
      ConfigIO::default_update_formatters(dict, db_props());

      std::vector<std::string> col = {
        "configname", "selected", "to_configname", "has_data", "has_complete_data",
        "score", "best_score", "is_best",
        "lattice_deformation_cost", "basis_deformation_cost",
        "relaxed_energy"
      };

      return dict.parse(col);
    }


  }
}
