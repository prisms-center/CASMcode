#include "casm/database/ConfigImport.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/clex/io/json/ConfigMapping_json_io.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/Import_impl.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/Update_impl.hh"

namespace CASM {

template class DataFormatter<DB::ConfigIO::Result>;

namespace Local {

/// Construct MappedProperties from a SimpleStructure solution of the
/// configuration mapping algorithm
///
/// \param simple_structure A SimpleStructure solution of the mapping algorithm.
/// This represents
///        what the configuration mapping algorithm considered to be the best
///        way to map the calculated SimpleStructure to the particular
///        Configuration under consideration.
/// \param dof_managed_properties A list of the SimpleStructure properties, both
/// site and global,
///        that correspond to BasicStructure degrees of freedom and thus should
///        not be included in MappedProperties
/// \param lattice_deformation_cost The lattice mapping score. Expected from
/// xtal::LatticeNode::cost. \param atomic_deformation_cost The basis mapping
/// score. Expected from xtal::AssignmentNode::cost. \param total_cost Total
/// mapping score. Expected from xtal::MappingNode::cost. Depending on
///        mapping method options may not be a linear combination of
///        lattice_deformation_cost and atomic_deformation_cost.
///
/// Note:
/// - Property names in "simple_structure" and "dof_managed_properties" must
/// follow CASM property
///   naming conventions as documented for AnisoValTraits.
MappedProperties make_mapped_properties(
    SimpleStructure const &simple_structure,
    std::set<std::string> const &dof_managed_properties,
    double lattice_deformation_cost, double atomic_deformation_cost,
    double total_cost) {
  MappedProperties result;

  for (auto const &prop : simple_structure.properties) {
    if (!dof_managed_properties.count(prop.first)) {
      result.global[prop.first] = prop.second;
      // If "*strain" is a property, rather than a DoF, we will also store the
      // lattice
      if (prop.first.find("strain") != std::string::npos) {
        result.global["latvec"] = simple_structure.lat_column_mat;
      }
    }
  }

  for (auto const &prop : simple_structure.mol_info.properties) {
    if (!dof_managed_properties.count(prop.first)) {
      result.site[prop.first] = prop.second;
      // If "disp" is a property, rather than a DoF, we will also store the
      // coordinates
      if (prop.first == "disp") {
        result.site["coordinate"] = simple_structure.mol_info.coords;
      }
    }
  }

  result.scalar("lattice_deformation_cost") = lattice_deformation_cost;
  result.scalar("atomic_deformation_cost") = atomic_deformation_cost;
  result.scalar("total_cost") = total_cost;
  return result;
}

static MappedProperties _make_mapped_properties(
    MappingNode const &_node, ConfigMapperResult::Individual const &_map) {
  return make_mapped_properties(
      _map.resolved_struc, _map.dof_managed_properties, _node.lattice_node.cost,
      _node.atomic_node.cost, _node.cost);
}
}  // namespace Local

namespace DB {

// --- Configuration specializations ---------------------------------------

// --- Import<Configuration> ---

ConfigMapping::Settings const &StructureMap<Configuration>::settings() const {
  return m_configmapper->settings();
}

/// Construct with PrimClex and ConfigMapping::Settings (see Import / Update
/// desc)
StructureMap<Configuration>::StructureMap(ConfigMapping::Settings const &_set,
                                          const PrimClex &primclex)
    : m_primclex_ptr(&primclex) {
  // -- construct ConfigMapper --
  m_configmapper.reset(
      new ConfigMapper(primclex, _set, primclex.crystallography_tol()));
}

/// \brief Specialized import method for ConfigType
///
/// \param p Path to structure or properties.calc.json file. Not guaranteed to
/// exist or be valid. \param hint std::unique_ptr<Configuration> to 'from'
/// config for 'casm update', or null if unknown as with 'casm import'. \param
/// result Insert iterator of Result objects to output mapping results
///
/// - Should output one or more mapping results from the structure located at
/// specied path
/// - >1 result handles case of non-primitive configurations
/// - responsible for filling in Result data structure
StructureMap<Configuration>::map_result_inserter
StructureMap<Configuration>::map(
    fs::path p, std::vector<std::string> const &req_properties,
    std::unique_ptr<Configuration> const &hint_config,
    map_result_inserter result) const {
  // need to set Result data (w/ defaults):
  // - std::string pos = "";
  // - MappedProperties mapped_props {origin:"", to:"", unmapped:{}, mapped:{}};
  // - bool has_data = false;
  // - bool has_complete_data = false;
  // - bool is_new_config = false;
  // - std::string fail_msg = "";
  ConfigIO::Result res;
  res.pos_path = p.string();

  if (!fs::exists(res.pos_path)) {
    res.has_files = false;
    res.fail_msg = "Specified file does not exist!";
  } else {
    res.has_files = true;
  }
  // read from structure file or properties.calc.json file (if exists)
  SimpleStructure sstruc = this->_make_structure(res.pos_path);

  // do mapping
  ConfigMapperResult map_result =
      m_configmapper->import_structure(sstruc, hint_config.get());

  if (!map_result.success()) {
    res.fail_msg = map_result.fail_msg;
    *result++ = std::move(res);
    return result;
  }

  if (map_result.n_optimal() > 1) {
    res.fail_msg = "There were " + std::to_string(map_result.n_optimal()) +
                   " optimal mappings, when only one was expected.";
    std::cout << "#maps: " << map_result.maps.size() << std::endl;
    for (auto const &map : map_result.maps) {
      std::cout << "---" << std::endl;
      jsonParser json;
      to_json(map.second.resolved_struc, json["resolved_struc"]);
      to_json(Local::_make_mapped_properties(map.first, map.second),
              json["mapped_properties"]);
      std::cout << "---" << std::endl;
    }
    *result++ = std::move(res);
    return result;
  }

  for (auto const &map : map_result.maps) {
    // insert in database (note that this also/only inserts primitive)
    ConfigInsertResult insert_result =
        map.second.config.insert(settings().primitive_only);

    if (insert_result.canonical_it !=
        m_primclex_ptr->db<Configuration>().end()) {
      res.is_new_config = insert_result.insert_canonical;

      res.properties = Local::_make_mapped_properties(map.first, map.second);
      res.properties.file_data = p.string();
      res.properties.to = insert_result.canonical_it.name();
      res.properties.origin = p.string();

      if (hint_config) {
        res.properties.init_config = hint_config->name();
      }

      res.has_data = false;
      res.has_complete_data = true;
      for (std::string const &propname : req_properties) {
        if (res.properties.global.count(propname) ||
            res.properties.site.count(propname)) {
          res.has_data = true;
        } else {
          res.has_complete_data = false;
        }
      }

      // at this point, the mapped structure result is complete
      *result++ = res;
    }

    // it may be the structure was not primitive:
    // - in which case we need to create a result indicating that the primitive
    //   was also inserted in the database,
    // - but don't try to scale the data for the primitive
    if ((insert_result.primitive_it !=
         m_primclex_ptr->db<Configuration>().end()) &&
        (insert_result.primitive_it != insert_result.canonical_it)) {
      ConfigIO::Result prim_res;
      prim_res.pos_path = res.pos_path;
      prim_res.has_files = false;
      prim_res.properties.file_data = p.string();
      prim_res.properties.origin =
          "prim:" +
          res.properties.origin;  // insert_result.primitive_it.name();
      prim_res.properties.to = insert_result.primitive_it.name();
      prim_res.is_new_config = insert_result.insert_primitive;
      // at this point, the mapped structure result is complete
      *result++ = prim_res;
    }
  }
  return result;
}

/// \brief Read BasicStructure to be imported
///
/// If 'p.extension()' == ".json" or ".JSON", read as properties.calc.json
/// Else, read as VASP POSCAR
SimpleStructure StructureMap<Configuration>::_make_structure(
    const fs::path &p) const {
  SimpleStructure sstruc;
  if (p.extension() == ".json" || p.extension() == ".JSON") {
    jsonParser json(p);
    from_json(sstruc, json);
  } else {
    fs::ifstream struc_stream(p);
    BasicStructure struc = BasicStructure::from_poscar_stream(struc_stream);
    sstruc = make_simple_structure(struc);
  }
  return sstruc;
}

// --- Import<Configuration> ---

/// \brief Constructor
Import<Configuration>::Import(const PrimClex &primclex,
                              const StructureMap<Configuration> &mapper,
                              ImportSettings const &_set,
                              std::string const &report_dir)
    :

      ImportT(primclex, mapper, _set, report_dir) {}

const std::string Import<Configuration>::desc =

    "Import Configuration: \n\n"

    "  'casm import' of Configuration proceeds in two steps: \n\n"

    "  1) For each file: \n"
    "     - Read structure from a CASM structure file (typically named \n"
    "       structure.json or properties.calc.json) or VASP POSCAR type file.\n"
    "     - Map the structure onto a Configuration of the primitive crystal \n"
    "       structure. \n"
    "     - Record relaxation data (lattice & basis deformation cost). \n\n"

    "  2) If properties or files import is requested (--properties, \n"
    "     --copy-structure-files, or --copy-additional-files), iterate over \n"
    "     each import record and do the following:\n"
    "     - If property import is requested, properties are imported from \n"
    "       all successfully mapped structures. \n"
    "     - If multiple imported structures map onto the same configuration, \n"
    "       the properties from the structure with the best value for \n"
    "       \"conflict_score\" are used when querying properties of the \n"
    "       configuration. By default, property \"energy\" is used for the \n"
    "       conflict score and the minimum energy is best.\n"
    "     - If it is requested that the structure file or additional files \n"
    "       are copied into the training_data directories, copying \n"
    "       occurs if either the configuration being mapped to has no \n"
    "       properties or files, or if the structure being imported has the \n"
    "       best conflict score and the \"overwrite\" option is given. \n"
    "     - Reports are generated detailing the results of the import: \n"
    "       - map_fail: Structures that could not be mapped onto the \n"
    "         primitive crystal structure. \n"
    "       - map_success: Configurations that were successfully mapped and\n"
    "         imported into the Configuration database (or already existed).\n"
    "       - import_data_fail: Structures with data that would be imported\n"
    "         except preexisting data prevents it. \n"
    "       - import_conflict: Configurations that were mapped to by multiple\n"
    "         structures. \n\n"

    "Settings: \n\n"

    "  mapping: JSON object (optional)\n"
    "      A JSON object containing the following options controlling the "
    "structure-\n"
    "      mapping algorithm:\n"

    "    primitive_only: bool (optional, default=false)\n"
    "        By convention, primitive configurations are always imported along "
    "with \n"
    "        non-primitive configurations. If false, only the primitive "
    "configuration\n"
    "        will be imported. Note: data from non-primitive configurations is "
    "never\n"
    "        used for primitive configurations.\n\n"

    "    lattice_weight: number in range [0.0, 1.0] (optional, default=0.5) \n"
    "        Candidate configurations are compared using \"deformation_cost\" "
    "to \n"
    "        determine the best mapping of the import structure to a           "
    "  \n"
    "        configuration. The \"lattice_weight\" determines the relative "
    "weight\n"
    "        of the \"lattice_deformation_cost\" and "
    "\"atomic_deformation_cost\"  \n"
    "        when calculating the total \"deformation_cost\". See description \n"
    "        of Deformation Cost in \"casm update --desc\" for details. \n\n"

    "    max_vol_change: number (optional, default=0.3)\n"
    "        Adjusts range of SCEL volumes searched while mapping imported \n"
    "        structure onto ideal crystal (only necessary if the presence of \n"
    "        vacancies makes the volume ambiguous). Default is +/- 30% of the "
    "rounded\n"
    "        integer volume (relaxed volume / primitive unit cell volume) of "
    "the \n"
    "        structure . Smaller values yield faster execution, larger values "
    "may \n"
    "        yield more accurate mapping.\n\n"

    "    max_va_frac: number (optional, default=0.5)\n"
    "        Places upper bound on the fraction of sites that are allowed to "
    "be \n"
    "        vacant after relaxed structure is mapped onto the ideal crystal. "
    "\n"
    "        Smaller values yield faster execution, larger values may yield "
    "more \n"
    "        accurate mapping. Has no effect if supercell volume can be "
    "inferred \n"
    "        from the number of atoms in the structure. Default value allows "
    "up to \n"
    "        50% of sites to be vacant.\n\n"

    "    min_va_frac: number (optional, default=0.0)\n"
    "        Places lower bound on the fraction of sites that are allowed to "
    "be \n"
    "        vacant after relaxed structure is mapped onto the ideal crystal. "
    "\n"
    "        Nonzero values may yield faster execution if updating "
    "configurations \n"
    "        that are known to have a large number of vacancies, at potential "
    "\n"
    "        sacrifice of mapping accuracy. Has no effect if supercell volume "
    "can \n"
    "        be inferred from the number of atoms in the structure. Default "
    "value \n"
    "        allows as few as 0% of sites to be vacant.\n\n"

    "    ideal: bool (optional, default=false)\n"
    "        Assume imported structures are in the setting of the ideal "
    "crytal.\n"
    "        This results in faster mapping, but may not identify the ideal "
    "mapping.\n"
    "        If large mapping costs are encountered, try re-running with ideal "
    ": false\n\n"

    "    robust: bool (optional, default=false)\n"
    "        Perform additional checks to determine if mapping is degenerate "
    "in cost\n"
    "        to other mappings, which can occur if the imported structure has "
    "symmetry\n"
    "        that is incompatible with prim.json. Results in slower "
    "execution.\n\n"

    "    filter: string (optional) \n"
    "        Restricts the import to only consider supercells that match a "
    "provided\n"
    "        casm query expression.\n\n"

    "    forced_lattices: array of strings (optional) \n"
    "        Restricts the import to only consider supercells provided via a "
    "list of\n"
    "        a list of their conventional names (i.e., "
    "\"SCEL2_2_1_1_1_1_0\").\n\n"

    "  data: JSON object (optional)\n"
    "      A JSON object containing the following options controlling when "
    "calculated \n"
    "      properties are updated.\n\n"

    "    import_properties: bool (optional, default=false)\n"
    "        If true, attempt to import structure properties. \n"
    "        Properties are read from the structure files indicated by --pos \n"
    "        or --batch. If the file ends in \".json\" or \".JSON\" it is \n"
    "        read as a CASM structure file, otherwise it is read as a VASP \n"
    "        POSCAR.\n\n"

    "    copy_structure_files: bool (optional, default=false)\n"
    "        If true, attempt to copy structure files into the project as \n"
    "        \"properties.calc.json\". Files will only by copied if there \n"
    "        are no existing files or data for the configuration the      \n"
    "        structure has been mapped to, or it is the best scoring      \n"
    "        mapping and \"overwrite\"=true. Copying \"properties.calc.json\" \n"
    "        implies \"import_properties\"=true.\n\n"

    "    copy_additional_files: bool (optional, default=false)\n"
    "        If true, attempt to copy all files & directories in the same \n"
    "        directory as the structure file or, if it exists, the \n"
    "        properties.calc.json file. Files & directories will only be \n"
    "        copied if there are no existing files or directories for the \n"
    "        configuration the structure has been mapped to, or it is the \n"
    "        best scoring mapping and \"overwrite\"=true. Note that this \n"
    "        flag also requires setting \"copy_structure_files\"=true.\n\n"

    "    overwrite: bool (optional, default=false)\n"
    "        If true, data and files will be imported that overwrite existing\n"
    "        data and files, if the score calculated by the \n"
    "        \"conflict_score\" for the configuration being mapped to will be\n"
    "        improved.\n\n";

int Import<Configuration>::run(const PrimClex &primclex,
                               const jsonParser &kwargs,
                               const Completer::ImportOption &import_opt) {
  // -- collect input settings --

  ConfigMapping::Settings map_settings;
  if (kwargs.contains("mapping")) from_json(map_settings, kwargs["mapping"]);

  ImportSettings import_settings;
  if (kwargs.contains("data")) from_json(import_settings, kwargs["data"]);

  // get input report_dir, check if exists, and create new report_dir.i if
  // necessary
  std::string report_dir =
      (fs::path(primclex.dir().reports_dir()) / "import_report").string();
  report_dir = create_report_dir(report_dir);

  // 'mapping' subsettings are used to construct ConfigMapper, and also returns
  // the 'used' settings
  StructureMap<Configuration> mapper(map_settings, primclex);

  jsonParser used_settings;
  used_settings["mapping"] = mapper.settings();
  used_settings["data"] = import_settings;

  // -- print used settings --
  Log &log = CASM::log();
  log.read("Settings");
  log << used_settings << std::endl << std::endl;

  // -- construct Import --
  Import<Configuration> f(primclex, mapper, import_settings, report_dir);

  // -- read structure file paths --
  std::set<fs::path> pos;
  auto res =
      construct_pos_paths(primclex, import_opt, std::inserter(pos, pos.end()));
  if (res.second) {
    return res.second;
  }

  // -- read structure file paths --
  f.import(pos.begin(), pos.end());

  return 0;
}

const std::string Update<Configuration>::desc =

    "Update Configuration calculation results: \n\n"

    "  'casm update' of Configuration calculation results proceeds as follows: "
    "\n\n"

    "  For each Configuration in the input selection: \n"
    "   - Read properties.calc.json file from training_data directory.        "
    "\n"
    "   - Map the relaxed structure onto a Configuration of the primitive "
    "crystal\n"
    "     structure. \n"
    "   - Record relaxation data: \n"
    "     - Lattice & basis deformation cost \n"
    "     - Initial configuration and relaxed configuration \n\n"
    "   - If multiple configurations relax onto a configuration for which "
    "there \n"
    "     is no calculation data, the calculation data from the with the "
    "lowest \n"
    "     conflict resolution score is used for the relaxed configuration.\n\n"

    "Settings: \n\n"

    "  force: bool (optional, default=false) \n"
    "    Force update all specified Configuration, else use timestamps to      "
    " \n"
    "    determine which to update. \n"

    "Settings: \n\n"

    "  mapping: JSON object (optional)\n"
    "      A JSON object containing the following options controlling the "
    "structure-\n"
    "      mapping algorithm:\n"

    "    primitive_only: bool (optional, default=false)\n"
    "        By convention, primitive configurations are always imported along "
    "with \n"
    "        non-primitive configurations. If false, only the primitive "
    "configuration\n"
    "        will be imported. Note: data from non-primitive configurations is "
    "never\n"
    "        used for primitive configurations.\n\n"

    "    lattice_weight: number in range [0.0, 1.0] (optional, default=0.5) \n"
    "        Candidate configurations are compared using \"deformation_cost\" "
    "to \n"
    "        determine the best mapping of the import structure to a           "
    "  \n"
    "        configuration. The \"lattice_weight\" determines the relative "
    "weight\n"
    "        of the \"lattice_deformation_cost\" and "
    "\"atomic_deformation_cost\"  \n"
    "        when calculating the total \"deformation_cost\". See description below. \n\n"

    "    max_vol_change: number (optional, default=0.3)\n"
    "        Adjusts range of SCEL volumes searched while mapping imported \n"
    "        structure onto ideal crystal (only necessary if the presence of \n"
    "        vacancies makes the volume ambiguous). Default is +/- 30% of the "
    "rounded\n"
    "        integer volume (relaxed volume / primitive unit cell volume) of "
    "the \n"
    "        structure . Smaller values yield faster execution, larger values "
    "may \n"
    "        yield more accurate mapping.\n\n"

    "    max_va_frac: number (optional, default=0.5)\n"
    "        Places upper bound on the fraction of sites that are allowed to "
    "be \n"
    "        vacant after relaxed structure is mapped onto the ideal crystal. "
    "\n"
    "        Smaller values yield faster execution, larger values may yield "
    "more \n"
    "        accurate mapping. Has no effect if supercell volume can be "
    "inferred \n"
    "        from the number of atoms in the structure. Default value allows "
    "up to \n"
    "        50% of sites to be vacant.\n\n"

    "    min_va_frac: number (optional, default=0.0)\n"
    "        Places lower bound on the fraction of sites that are allowed to "
    "be \n"
    "        vacant after relaxed structure is mapped onto the ideal crystal. "
    "\n"
    "        Nonzero values may yield faster execution if updating "
    "configurations \n"
    "        that are known to have a large number of vacancies, at potential "
    "\n"
    "        sacrifice of mapping accuracy. Has no effect if supercell volume "
    "can \n"
    "        be inferred from the number of atoms in the structure. Default "
    "value \n"
    "        allows as few as 0% of sites to be vacant.\n\n"

    "    ideal: bool (optional, default=false)\n"
    "        Assume imported structures are in the setting of the ideal "
    "crytal.\n"
    "        This results in faster mapping, but may not identify the ideal "
    "mapping.\n"
    "        If large mapping costs are encountered, try re-running with ideal "
    ": false\n\n"

    "    fix_volume: bool (optional, default=false)\n"
    "        Assume imported structures have the same integer volume as the "
    "starting\n"
    "        configuration. All supercells of this volume are considered for "
    "mapping.\n"
    "        This assumption may fail for systems that allow vacancies, when "
    "vacancy\n"
    "        concentration is high. Increases execution speed in these "
    "cases.\n\n"

    "    fix_lattice: bool (optional, default=false)\n"
    "        Assume imported structures have the same lattice as the starting "
    "configuration\n"
    "        Only this supercell will be considered for mapping, but all "
    "orientational\n"
    "        relationships will still be considered. Increases execution "
    "speed, especially\n"
    "        at large supercell volume, but cannot detect relaxation to a "
    "different supercell.\n\n"

    "    robust: bool (optional, default=false)\n"
    "        Perform additional checks to determine if mapping is degenerate "
    "in cost\n"
    "        to other mappings, which can occur if the imported structure has "
    "symmetry\n"
    "        that is incompatible with prim.json. Results in slower "
    "execution.\n\n"

    "  data: JSON object (optional)\n"
    "      A JSON object containing the following options controlling how "
    "calculated \n"
    "      properties are updated. After structural relaxation, a structural "
    "that \n"
    "      began as structure 'A' may map more closely onto a different "
    "structure, 'B'\n"
    "      In some cases, multiple structures may map onto the same 'B', and a "
    "scoring\n"
    "      metric is used to specify which set of calculation data is "
    "associated with \n"
    "      configuration 'B'. The following values determine the scoring "
    "metric:\n\n"

    "    \"deformation_cost\":\n"
    "       \"lattice_weight\": number, in range [0, 1.0]\n"

    "       Uses a weighted sum of cost functions for lattice and basis \n"
    "       deformation. See below for complete definition. Ex: \n"
    "         {\"method\":\"deformation_cost\":, \"lattice_weight\":0.5} \n\n"

    "    \"minimum\":\n"
    "       \"property\": property name (ex: \"energy\")\n"

    "       Reads the specified property from the mapped properties and\n"
    "       selects the minimum to be the best mapping. Ex: \n"
    "         {\"method\":\"minimum\":, \"property\": \"energy\"} \n\n"

    "    \"maximum\":\n"
    "       \"property\": property name (ex: \"energy\")\n"

    "       Reads the specified property from the mapped properties and\n"
    "       selects the maximum to be the best mapping. Ex: \n"
    "         {\"method\":\"maximum\":, \"property\": \"energy\"} \n\n"

    "    \"direct_selection\":\n"
    "       \"name\": configname to force as 'best' (ex: "
    "\"SCEL3_1_1_3_0_0_0/4\")\n"

    "       Directly specify which configuration's properties should be used. "
    "Ex: \n"
    "         {\"method\":\"direct_selection\":, \"name\": "
    "\"SCEL3_1_1_3_0_0_0/4\"} \n\n"

    "  The default value used by CASM is:\n"
    "    {\"method\": \"minimum\", \"property\": \"energy\"}\n\n"

    "Deformation cost:\n"

    "  The \"deformation_cost\" is:\n\n"

    "      deformation_cost = w*lattice_deformation_cost + \n"
    "                         (1.0-w)*atomic_deformation_cost,\n\n"

    "  where \"w\" is the \"lattice_weight\" factor (default=0.5),\n\n"

    "  the \"lattice_deformation_cost\" is the mean-square displacement of a \n"
    "  points on the surface of a sphere of volume equal to the atomic volume "
    "\n"
    "  when  it is deformed by the volume preserving deviatoric deformation \n"
    "  matrix, F_deviatoric:\n\n"

    "      V = relaxed_atomic_volume;\n"
    "      F = deformation matrix (lattice_relaxed = F*lattice_ideal);\n"
    "      F_deviatoric = F/pow(F.determinant(), 1./3.);\n"
    "      I = 3x3 identity matrix;\n\n"

    "      lattice_deformation_cost = pow( 3.*V / (4.*pi), 2.0/3.0) / 3.0 * \n"
    "          (0.5 * (F.t * F / pow(std::abs(F.determinant()), 2.0/3.0) - "
    "I)).squaredNorm()\n\n"

    "  and the \"atomic_deformation_cost\" is a cost function for the amount "
    "of\n"
    "  basis site relaxation:\n\n"

    "      D = 3xN matrix of basis site displacements (displacements are "
    "applied \n"
    "          before strain)\n"
    "      Natoms = number of atoms in configuration\n"
    "      atomic_deformation_cost = (F*D * D.transpose() * "
    "F.transpose()).trace() \n"
    "          / (max(Natoms, 1.))\n\n";

/// Allow ConfigType to specialize the report formatting for 'import'
DataFormatter<ConfigIO::Result> Import<Configuration>::_import_formatter()
    const {
  DataFormatterDictionary<ConfigIO::Result> dict;
  ConfigIO::default_import_formatters(dict, db_props());

  std::vector<std::string> col = {"initial_path",
                                  "selected",
                                  "fail_msg",
                                  "to_configname",
                                  "final_path",
                                  "is_new_config",
                                  "has_any_required_properties",
                                  "has_all_required_properties",
                                  "preexisting_data",
                                  "preexisting_files",
                                  "import_properties",
                                  "properties_origin",
                                  "import_structure_file",
                                  "import_additional_files",
                                  "score",
                                  "best_score",
                                  "is_new_best",
                                  "lattice_deformation_cost",
                                  "atomic_deformation_cost",
                                  "energy"};

  return dict.parse(col);
}

// --- Update<Configuration> ---

/// \brief Constructor
Update<Configuration>::Update(const PrimClex &primclex,
                              const StructureMap<Configuration> &mapper,
                              UpdateSettings const &set,
                              std::string const &report_dir)
    : UpdateT(primclex, mapper, set, report_dir) {}

int Update<Configuration>::run(const PrimClex &primclex,
                               const jsonParser &kwargs,
                               const Completer::UpdateOption &update_opt) {
  // -- collect input settings --

  jsonParser used;

  // general settings
  bool force;
  kwargs.get_else(force, "force", false);
  used["force"] = force;

  // TODO: this could take more settings, for now output_as_json fixed true
  UpdateSettings update_settings;

  // 'mapping' subsettings are used to construct ConfigMapper and return 'used'
  // settings values still need to figure out how to specify this in general
  jsonParser map_json;
  kwargs.get_else(map_json, "mapping", jsonParser());
  ConfigMapping::Settings map_settings =
      map_json.get<ConfigMapping::Settings>();

  StructureMap<Configuration> mapper(map_settings, primclex);
  used["mapping"] = mapper.settings();

  // 'data' subsettings
  // bool import_data = true;
  // bool import_additional_files = false;
  // bool overwrite = false;

  // -- print used settings --
  Log &log = CASM::log();
  log.read("Settings");
  log << used << std::endl << std::endl;

  // get input report_dir, check if exists, and create new report_dir.i if
  // necessary
  std::string report_dir =
      (fs::path(primclex.dir().reports_dir()) / "update_report").string();
  report_dir = create_report_dir(report_dir);

  // -- construct Update --
  Update<Configuration> f(primclex, mapper, update_settings, report_dir);

  // -- read selection --
  DB::Selection<Configuration> sel(primclex, update_opt.selection_path());

  // -- update --
  f.update(sel, force);

  return 0;
}

// Allow ConfigType to specialize the report formatting for 'update'
DataFormatter<ConfigIO::Result> Update<Configuration>::_update_formatter()
    const {
  DataFormatterDictionary<ConfigIO::Result> dict;
  ConfigIO::default_update_formatters(dict, db_props());

  std::vector<std::string> col = {"initial_path",
                                  "properties_origin",
                                  "selected",
                                  "fail_msg",
                                  "to_configname",
                                  "has_any_required_properties",
                                  "has_all_required_properties",
                                  "score",
                                  "best_score",
                                  "is_new_best",
                                  "lattice_deformation_cost",
                                  "atomic_deformation_cost",
                                  "energy"};

  return dict.parse(col);
}

}  // namespace DB
}  // namespace CASM
