#include "casm/database/ConfigImport.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/jsonStruc.hh"
#include "casm/clex/Calculable.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/import.hh"
#include "casm/app/update.hh"
#include "casm/app/rm.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/Selection.hh"
#include "casm/database/Import_impl.hh"
#include "casm/database/Update_impl.hh"
#include "casm/casm_io/DataFormatter_impl.hh"
#include "casm/basis_set/DoF.hh"

namespace CASM {

  template class DataFormatter<DB::ConfigIO::Result>;

  namespace DB {

    // --- Configuration specializations ---------------------------------------

    // --- Import<Configuration> ---

    /// Construct with PrimClex and by moving a ConfigMapper
    StructureMap<Configuration>::StructureMap(
      const PrimClex &primclex,
      std::unique_ptr<ConfigMapper> mapper,
      bool primitive_only,
      std::vector<std::string> dof) :
      ConfigData<Configuration>(primclex, null_log()),
      m_configmapper(std::move(mapper)),
      m_primitive_only(primitive_only),
      m_dof(dof) {}

    /// Construct with PrimClex and JSON settings (see Import / Update desc)
    StructureMap<Configuration>::StructureMap(
      const PrimClex &primclex,
      const jsonParser &kwargs,
      bool primitive_only,
      std::vector<std::string> dof) :
      ConfigData<Configuration>(primclex, null_log()),
      m_primitive_only(primitive_only),
      m_dof(dof) {

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
      int map_opt = ConfigMapper::none;
      if(rotate) map_opt |= ConfigMapper::rotate;
      if(strict) map_opt |= ConfigMapper::strict;
      if(!ideal) map_opt |= ConfigMapper::robust;

      m_configmapper.reset(new ConfigMapper(
                             primclex,
                             lattice_weight,
                             max_vol_change,
                             map_opt,
                             primclex.crystallography_tol()));
      m_configmapper->set_min_va_frac(min_va_frac);
      m_configmapper->set_max_va_frac(max_va_frac);

    }

    const jsonParser &StructureMap<Configuration>::used() const {
      return m_used;
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
    StructureMap<Configuration>::map_result_inserter StructureMap<Configuration>::map(
      fs::path p,
      DatabaseIterator<Configuration> hint,
      map_result_inserter result) const {
      // need to set Result data (w/ defaults):
      // - fs::path pos = "";
      // - MappedProperties mapped_props {from:"", to:"", unmapped:{}, mapped:{}};
      // - bool has_data = false;
      // - bool has_complete_data = false;
      // - bool is_new_config = false;
      // - std::string fail_msg = "";

      // get path to properties.calc.json that will be imported
      //   (checking a couple possible locations relative to pos_path),
      //   or empty if none could be found
      fs::path prop_path = this->calc_properties_path(p);

      ConfigIO::Result res;
      res.pos = (prop_path.empty() ? p : prop_path);

      std::unique_ptr<Configuration> hint_config;
      if(hint != db_config().end()) {
        hint_config = notstd::make_unique<Configuration>(*hint);
        res.mapped_props.from = hint_config->name();
      }

      // read from structure file or properties.calc.json file (if exists)
      BasicStructure<Site> struc = this->_make_structure(res.pos);

      // do mapping
      ConfigMapperResult map_result;
      if(_occupation_only()) {
        map_result = m_configmapper->import_structure_occupation(struc, hint_config.get());
      }
      else {
        map_result =  m_configmapper->import_structure(struc);
      }

      if(!map_result.success) {
        res.fail_msg = map_result.fail_msg;
        *result++ = res;
        return result;
      }
      // if the result was a success, need to populate relaxed energy in
      // map_result.relaxation_properties["best_mapping"]["relaxed_energy"]
      if(res.pos.extension() == ".json" || res.pos.extension() == ".JSON") {
        jsonParser json(res.pos);
        if(json.contains("relaxed_energy")) {
          map_result.relaxation_properties["best_mapping"]["relaxed_energy"] = json["relaxed_energy"];
        }
      }
      // insert in database (note that this also/only inserts primitive)
      ConfigInsertResult insert_result = map_result.config->insert(m_primitive_only);

      res.is_new_config = insert_result.insert_canonical;

      res.mapped_props.to = insert_result.canonical_it.name();

      // check for and read raw 'unmapped' data, adds 'data_timestamp'
      if(!prop_path.empty()) {
        std::tie(res.mapped_props.unmapped, res.has_data, res.has_complete_data) =
          read_calc_properties<Configuration>(primclex(), prop_path);
      }

      // copy relaxation properties from best config mapping into 'mapped' props
      auto it = map_result.relaxation_properties["best_mapping"].begin();
      auto end = map_result.relaxation_properties["best_mapping"].end();
      for(; it != end; ++it) {
        res.mapped_props.mapped[it.name()] = *it;
      }
      res.mapped_props.mapped["best_assignment"] = map_result.best_assignment;
      res.mapped_props.mapped["cart_op"] = map_result.cart_op;

      // at this point, the mapped structure result is complete
      *result++ = res;

      // it may be the structure was not primitive:
      // - in which case we need to create a result indicating that the primitive
      //   was also inserted in the database,
      // - but don't try to scale the data for the primitive
      if(insert_result.canonical_it != insert_result.primitive_it) {
        ConfigIO::Result prim_res;
        prim_res.pos = res.pos;
        prim_res.mapped_props.from = insert_result.primitive_it.name();
        prim_res.mapped_props.to = insert_result.primitive_it.name();
        prim_res.is_new_config = insert_result.insert_primitive;
        *result++ = prim_res;
      }

      return result;
    }

    /// \brief Read BasicStructure<Site> to be imported
    ///
    /// If 'p.extension()' == ".json" or ".JSON", read as properties.calc.json
    /// Else, read as VASP POSCAR
    BasicStructure<Site> StructureMap<Configuration>::_make_structure(const fs::path &p) const {

      BasicStructure<Site> struc;
      if(p.extension() == ".json" || p.extension() == ".JSON") {
        jsonParser json(p);
        from_json(simple_json(struc, "relaxed_"), json);
      }
      else {
        fs::ifstream struc_stream(p);
        struc.read(struc_stream);
      }
      return struc;
    }

    /// \brief Import Configuration with only occupation DoF
    bool StructureMap<Configuration>::_occupation_only() const {
      return m_dof.size() == 1 && m_dof[0] == "occupation";
    }


    // --- Import<Configuration> ---

    /// \brief Constructor
    Import<Configuration>::Import(
      const PrimClex &primclex,
      const StructureMap<Configuration> &mapper,
      bool import_data,
      bool copy_additional_files,
      bool overwrite,
      fs::path report_dir,
      Log &file_log) :

      ImportT(primclex, mapper, import_data, copy_additional_files, overwrite,
              report_dir, file_log) {}


    const std::string Import<Configuration>::desc =

      "Import Configuration: \n\n"

      "  'casm import' of Configuration proceeds in two steps: \n\n"

      "  1) For each file: \n"
      "     - Read structure from VASP POSCAR type file or CASM properties.calc.json \n"
      "       file. \n"
      "     - Map the structure onto a Configuration of the primitive crystal \n"
      "       structure. \n"
      "     - Record relaxation data (lattice & basis deformation cost). \n\n"

      "  2) If data import is requested, iterate over each import record and do \n"
      "     the following:\n"
      "     - If multiple imported structures map onto a configuration for which \n"
      "       there is no calculation data, import calculation data from the     \n"
      "       structure with the lowest \"conflict_score\".                      \n"
      "     - If one or more imported structures map onto a configuration for    \n"
      "       which calculation data already exist, do not import any new data   \n"
      "       unless the \"overwrite\" option is given, in which case data and   \n"
      "       files are imported if the \"conflict_score\" will be improved. \n"
      "     - If data is imported, the corresponding properties.calc.json file is\n"
      "       copied into the directory of the mapped configuration. Optionally, \n"
      "       additional files in the directory of the imported structure file may\n"
      "       also by copided. \n"
      "     - Reports are generated detailing the results of the import: \n"
      "       - import_map_fail: Structures that could not be mapped onto the     \n"
      "         primitive crystal structure. \n"
      "       - import_map_success: Configurations that were successfully mapped  \n"
      "         and imported into the Configuration database (or already existed).\n"
      "       - import_data_fail: Structures with data that would be imported     \n"
      "         except preexisting data prevents it. \n"
      "       - import_conflict: Configurations that were mapped to by multiple   \n"
      "         structures. \n\n"

      "Settings: \n\n"

      "  primitive_only: bool (optional, default=false)\n"
      "      By convention, primitive configurations are always imported along with \n"
      "      non-primitive configurations. If false, only the primitive configuration\n"
      "      will be imported. Note: data from non-primitive configurations is never\n"
      "      used for primitive configurations.\n\n"

      "  mapping: JSON object (optional)\n"
      "      A JSON object containing the following options controlling the structure-\n"
      "      mapping algorithm:\n"

      "    lattice_weight: number in range [0.0, 1.0] (optional, default=0.5) \n"
      "        Candidate configurations are compared using \"deformation_cost\" to \n"
      "        determine the best mapping of the import structure to a             \n"
      "        configuration. The \"lattice_weight\" determines the relative weight\n"
      "        of the \"lattice_deformation_cost\" and \"basis_deformation_cost\"  \n"
      "        when calculating the total \"deformation_cost\". See _____ for      \n"
      "        details. \n\n"

      "    max_vol_change: number (optional, default=0.3)\n"
      "        Adjusts range of SCEL volumes searched while mapping imported \n"
      "        structure onto ideal crystal (only necessary if the presence of \n"
      "        vacancies makes the volume ambiguous). Default is +/- 30% of the rounded\n"
      "        integer volume (relaxed volume / primitive unit cell volume) of the \n"
      "        structure . Smaller values yield faster execution, larger values may \n"
      "        yield more accurate mapping.\n\n"

      "    max_va_frac: number (optional, default=0.5)\n"
      "        Places upper bound on the fraction of sites that are allowed to be \n"
      "        vacant after relaxed structure is mapped onto the ideal crystal. \n"
      "        Smaller values yield faster execution, larger values may yield more \n"
      "        accurate mapping. Has no effect if supercell volume can be inferred \n"
      "        from the number of atoms in the structure. Default value allows up to \n"
      "        50% of sites to be vacant.\n\n"

      "    min_va_frac: number (optional, default=0.0)\n"
      "        Places lower bound on the fraction of sites that are allowed to be \n"
      "        vacant after relaxed structure is mapped onto the ideal crystal. \n"
      "        Nonzero values may yield faster execution if updating configurations \n"
      "        that are known to have a large number of vacancies, at potential \n"
      "        sacrifice of mapping accuracy. Has no effect if supercell volume can \n"
      "        be inferred from the number of atoms in the structure. Default value \n"
      "        allows as few as 0% of sites to be vacant.\n\n"

      "    ideal: bool (optional, default=false)\n"
      "        Assume imported structures are unstrained (ideal) for faster importing.\n"
      "        Can be slower if used on deformed structures, in which case more \n"
      "        robust methods will be used\n\n"

      "  data: JSON object (optional)\n"
      "      A JSON object containing the following options controlling when calculated \n"
      "      properties are updated.\n\n"

      "    import: bool (optional, default=false)\n"
      "        If true (default), attempt to import calculation results. Results are \n"
      "        added from a \"properties.calc.json\" file which is checked for in the\n"
      "        following locations, relative to the input file path, 'pos': \n"

      "        1) Is 'pos' a JSON file? If 'pos' ends in \".json\" or \".JSON\", then\n"
      "           it is assumed to be a 'properties.calc.json' file.\n"
      "        2) If / path / to / pos, checks for / path / to / calctype.current / properties.calc.json\n"
      "        3) If / path / to / pos, checks for / path / to / properties.calc.json \n\n"

      "        If false, only configuration structure is imported.\n\n"

      "    copy_additional_files: bool (optional, default = false)\n"
      "        If true, attempt to copy all files & directories in the same directory\n"
      "        as the structure file or , if it exists, the properties.calc.json file.\n"
      "        Files & directories will only by copied if there are no existing files\n"
      "        or directories in the 'training_data' directory for the configuration \n"
      "        the structure as been mapped to or \"overwrite\"=true.\n\n"

      "    overwrite: bool (optional, default=false)\n"
      "        If true, data and files will be imported that overwrite existing\n"
      "        data and files, if the score calculated by the \"conflict_score\"\n"
      "        for the configuration being mapped to will be improved.\n\n";


    int Import<Configuration>::run(
      const PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::ImportOption &import_opt) {

      // -- collect input settings --

      const po::variables_map &vm = import_opt.vm();
      jsonParser used;
      jsonParser _default;

      // general settings
      bool primitive_only;
      kwargs.get_else(primitive_only, "primitive_only", false);
      used["primitive_only"] = primitive_only;

      // get input report_dir, check if exists, and create new report_dir.i if necessary
      fs::path report_dir = primclex.dir().root_dir() / "import_report";
      report_dir = create_report_dir(report_dir);

      // 'mapping' subsettings are used to construct ConfigMapper, and also returns
      // the 'used' settings
      std::vector<std::string> dof {"occupation"};
      jsonParser map_json;
      kwargs.get_else(map_json, "mapping", jsonParser());
      StructureMap<Configuration> mapper(primclex, map_json, primitive_only, dof);
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
      Import<Configuration> f(
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


    const std::string Update<Configuration>::desc =

      "Update Configuration calculation results: \n\n"

      "  'casm update' of Configuration calculation results proceeds as follows: \n\n"

      "  For each Configuration in the input selection: \n"
      "   - Read properties.calc.json file from training_data directory.        \n"
      "   - Map the relaxed structure onto a Configuration of the primitive crystal\n"
      "     structure. \n"
      "   - Record relaxation data: \n"
      "     - Lattice & basis deformation cost \n"
      "     - Initial configuration and relaxed configuration \n\n"
      "   - If multiple configurations relax onto a configuration for which there \n"
      "     is no calculation data, the calculation data from the with the lowest \n"
      "     conflict resolution score is used for the relaxed configuration.\n"
      "   - Both default and configuration-specific conflict resolution scoring\n"
      "     method can be set via: \n"
      "       'casm update --set-default-conflict-score -i <JSON>'\n"
      "       'casm update --set-default-conflict-score -s <JSON filename>'\n"
      "       'casm update --set-conflict-score configname -i <JSON>'\n"
      "       'casm update --set-conflict-score configname -s <JSON filename>'\n"

      "Settings: \n\n"

      "  force: bool (optional, default=false) \n"
      "    Force update all specified Configuration, else use timestamps to       \n"
      "    determine which to update. \n"

      "  mapping: JSON object (optional)\n"
      "      A JSON object containing the following options controlling the structure-\n"
      "      mapping algorithm:\n"

      "    lattice_weight: number in range [0.0, 1.0] (optional, default=0.5) \n"
      "        Candidate configurations are compared using \"deformation_cost\" to \n"
      "        determine the best mapping of the import structure to a             \n"
      "        configuration. The \"lattice_weight\" determines the relative weight\n"
      "        of the \"lattice_deformation_cost\" and \"basis_deformation_cost\"  \n"
      "        when calculating the total \"deformation_cost\". See _____ for      \n"
      "        details. \n\n"

      "    max_vol_change: number (optional, default=0.3)\n"
      "        Adjusts range of SCEL volumes searched while mapping imported \n"
      "        structure onto ideal crystal (only necessary if the presence of \n"
      "        vacancies makes the volume ambiguous). Default is +/- 30% of the rounded\n"
      "        integer volume (relaxed volume / primitive unit cell volume) of the \n"
      "        structure . Smaller values yield faster execution, larger values may \n"
      "        yield more accurate mapping.\n\n"

      "    max_va_frac: number (optional, default=0.5)\n"
      "        Places upper bound on the fraction of sites that are allowed to be \n"
      "        vacant after relaxed structure is mapped onto the ideal crystal. \n"
      "        Smaller values yield faster execution, larger values may yield more \n"
      "        accurate mapping. Has no effect if supercell volume can be inferred \n"
      "        from the number of atoms in the structure. Default value allows up to \n"
      "        50% of sites to be vacant.\n\n"

      "    min_va_frac: number (optional, default=0.0)\n"
      "        Places lower bound on the fraction of sites that are allowed to be \n"
      "        vacant after relaxed structure is mapped onto the ideal crystal. \n"
      "        Nonzero values may yield faster execution if updating configurations \n"
      "        that are known to have a large number of vacancies, at potential \n"
      "        sacrifice of mapping accuracy. Has no effect if supercell volume can \n"
      "        be inferred from the number of atoms in the structure. Default value \n"
      "        allows as few as 0% of sites to be vacant.\n\n"

      "    ideal: bool (optional, default=false)\n"
      "        Assume imported structures are unstrained (ideal) for faster importing.\n"
      "        Can be slower if used on deformed structures, in which case more \n"
      "        robust methods will be used\n\n"

      "Conflict Resolution: \n\n"

      "  Which metric should be used to determine which calculation results \n"
      "  should be used for a particular configuration if multiple results\n"
      "  map to the same configuration and the self-mapping result is not \n"
      "  available. Should consist of a \"method\" and method-dependent \n"
      "  parameters. The \"method\" options and associated parameters are: \n"

      "    \"deformation_cost\":\n"
      "       \"lattice_weight\": number, in range [0, 1.0]\n"

      "       Uses a weighted sum of cost functions for lattice and basis \n"
      "       deformation. See below for complete definition. Ex: \n"
      "         {\"method\":\"deformation_cost\":, \"lattice_weight\":0.5} \n\n"

      "    \"minimum\":\n"
      "       \"property\": property name (ex: \"relaxed_energy\")\n"

      "       Reads the specified property from the mapped properties and\n"
      "       selects the minimum to be the best mapping. Ex: \n"
      "         {\"method\":\"minimum\":, \"property\": \"relaxed_energy\"} \n\n"

      "    \"maximum\":\n"
      "       \"property\": property name (ex: \"relaxed_energy\")\n"

      "       Reads the specified property from the mapped properties and\n"
      "       selects the maximum to be the best mapping. Ex: \n"
      "         {\"method\":\"maximum\":, \"property\": \"relaxed_energy\"} \n\n"

      "    \"direct_selection\":\n"
      "       \"name\": configname to force as 'best' (ex: \"SCEL3_1_1_3_0_0_0/4\")\n"

      "       Directly specify which configuration's properties should be used. Ex: \n"
      "         {\"method\":\"direct_selection\":, \"name\": \"SCEL3_1_1_3_0_0_0/4\"} \n\n"

      "  The default value used by CASM is:\n"
      "    {\"method\": \"minimum\", \"property\": \"relaxed_energy\"}\n\n"

      "Deformation cost:\n"

      "  The \"deformation_cost\" is:\n\n"

      "      deformation_cost = w*lattice_deformation_cost + \n"
      "                         (1.0-w)*basis_deformation_cost,\n\n"

      "  where \"w\" is the \"lattice_weight\" factor (default=0.5),\n\n"

      "  the \"lattice_deformation_cost\" is the mean-square displacement of a \n"
      "  points on the surface of a sphere of volume equal to the atomic volume \n"
      "  when  it is deformed by the volume preserving deviatoric deformation \n"
      "  matrix, F_deviatoric:\n\n"

      "      V = relaxed_atomic_volume;\n"
      "      F = deformation matrix (lattice_relaxed = F*lattice_ideal);\n"
      "      F_deviatoric = F/pow(F.determinant(), 1./3.);\n"
      "      I = 3x3 identity matrix;\n\n"

      "      lattice_deformation_cost = pow( 3.*V / (4.*pi), 2.0/3.0) / 3.0 * \n"
      "          (0.5 * (F.t * F / pow(std::abs(F.determinant()), 2.0/3.0) - I)).squaredNorm()\n\n"

      "  and the \"basis_deformation_cost\" is a cost function for the amount of\n"
      "  basis site relaxation:\n\n"

      "      D = 3xN matrix of basis site displacements (displacements are applied \n"
      "          before strain)\n"
      "      Natoms = number of atoms in configuration\n"
      "      basis_deformation_cost = (F*D * D.transpose() * F.transpose()).trace() \n"
      "          / (max(Natoms, 1.))\n\n";



    /// Allow ConfigType to specialize the report formatting for 'import'
    DataFormatter<ConfigIO::Result> Import<Configuration>::_import_formatter(
      const std::map<std::string, ConfigIO::ImportData> &data_results) const {

      DataFormatterDictionary<ConfigIO::Result> dict;
      ConfigIO::default_import_formatters(dict, db_props(), data_results);

      std::vector<std::string> col = {
        "configname", "selected", "pos", "has_data", "has_complete_data",
        "preexisting_data", "import_data", "import_additional_files",
        "score", "best_score", "is_best",
        "lattice_deformation_cost", "basis_deformation_cost",
        "relaxed_energy"
      };

      return dict.parse(col);
    }


    // --- Update<Configuration> ---

    /// \brief Constructor
    Update<Configuration>::Update(
      const PrimClex &primclex,
      const StructureMap<Configuration> &mapper,
      fs::path report_dir) :
      UpdateT(primclex, mapper, report_dir) {}

    int Update<Configuration>::run(
      const PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::UpdateOption &update_opt) {

      // -- collect input settings --

      const po::variables_map &vm = update_opt.vm();
      jsonParser used;

      // general settings
      bool primitive_only = false;

      bool force;
      if(vm.count("force")) {
        force = true;
      }
      else {
        kwargs.get_else(force, "force", false);
      }
      used["force"] = force;

      // get input report_dir, check if exists, and create new report_dir.i if necessary
      fs::path report_dir = primclex.dir().root_dir() / "update_report";
      report_dir = create_report_dir(report_dir);

      // 'mapping' subsettings are used to construct ConfigMapper and return 'used' settings values
      // still need to figure out how to specify this in general
      std::vector<std::string> dof {"occupation"};
      jsonParser map_json;
      kwargs.get_else(map_json, "mapping", jsonParser());
      StructureMap<Configuration> mapper(primclex, map_json, primitive_only, dof);
      used["mapping"] = mapper.used();

      // 'data' subsettings
      bool import_data = true;
      bool import_additional_files = false;
      bool overwrite = false;

      // -- print used settings --
      Log &log = primclex.log();
      log.read("Settings");
      log << used << std::endl << std::endl;

      // -- construct Update --
      Update<Configuration> f(
        primclex,
        mapper,
        report_dir);

      // -- read selection --
      DB::Selection<Configuration> sel(primclex, update_opt.selection_path());

      // -- update --
      f.update(sel, force);

      return 0;
    }

    // Allow ConfigType to specialize the report formatting for 'update'
    DataFormatter<ConfigIO::Result> Update<Configuration>::_update_formatter() const {

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
