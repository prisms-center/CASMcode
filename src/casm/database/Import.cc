
namespace CASM {
  namespace DB {

    /// \brief Constructor
    ImportBase::ImportBase(
      const PrimClex &primclex,
      bool import_data,
      bool copy_additional_files,
      bool overwrite,
      fs::path report_dir) :

      Logging(primclex),
      m_report_dir(report_dir),
      m_file_str(fs::open(m_report_dir / "files")),
      m_file_log(m_file_str),
      m_primclex(primclex),
      m_import_data(import_data),
      m_copy_additional_files(copy_additional_files),
      m_overwrite(overwrite) {}

    /// Create a new report directory to avoid overwriting existing results
    fs::path ImportBase::create_report_dir(fs::path report_dir) {
      Index i = 0;
      while(fs::exists(report_dir + "." + std::to_string(i))) {
        ++i;
      }
      report_dir = report_dir + "." + std::to_string(i);
      fs::create_directory(report_dir);
      return report_dir;
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
    fs::path ImportBase::_calc_properties_path(fs::path pos_path) const {

      // check 1: is a JSON file
      if(pos_path.extension() == ".json" || pos_path.extension() == ".JSON") {
        return pos_path;
      }

      // check 2: /path/to/POS -> /path/to/calctype.current/properties.calc.json
      {
        fs::path dft_path = pos_path;
        dft_path.remove_filename();
        (dft_path /= ("calctype." + primclex().settings().default_clex().calctype)) /= "properties.calc.json";
        if(fs::exists(dft_path)) {
          return dft_path;
        }
      }

      // check 3: /path/to/POS -> /path/to/properties.calc.json
      {
        fs::path dft_path = pos_path;
        dft_path.remove_filename();
        dft_path /= "properties.calc.json";
        if(fs::exists(dft_path)) {
          return dft_path;
        }
      }

      // not found, return empty path
      return fs::path();
    }

    /// \brief Return true if there are existing files in the traning_data directory
    ///        for a particular configuration
    bool ImportBase::_has_existing_files(const std::string &from_configname) const {
      fs::path p = _calc_dir(from_configname);
      if(!fs::exists(p)) {
        return false;
      }
      return std::distance(fs::directory_iterator(p), directory_iterator());
    }

    /// \brief Return true if there are existing files in the traning_data directory
    ///        for a particular configuration
    bool ImportBase::_has_existing_data(const std::string &from_configname) const {
      return db_props().find_via_from(from_configname) != db_props().end();
    }

    bool ImportBase::_has_existing_data_or_files(const std::string &from_configname) const {
      return _has_existing_data(from_configname) || _has_existing_files(from_configname);
    }

    /// Check if 'properties.calc.json' file has not changed since last read
    ///
    /// - Compares 'data_timestamp' && fs::last_write_time
    bool ImportBase::_no_change(const std::string &configname) const {
      fs::path prop_path = calc_properties_path(configname);
      if(!prop_path.empty()) {
        auto it = db_props().find_via_from(configname);
        if(it != db_props().end()) {
          auto json_it = it->unmapped.find("data_timestamp");
          if(json_it != it->unmapped.end() &&
             json_it->get<time_t>() == fs::last_write_time(prop_path)) {

            return true;
          }
        }
      }
      return false;
    }

    /// \brief Remove existing files in the traning_data directory for a particular
    ///        configuration
    void ImportBase::_rm_files(const std::string &configname) const {
      fs::path p = _calc_dir(configname);
      if(!fs::exists(p)) {
        return;
      }
      m_file_log.custom(std::string("Remove calculation files: ") + configname);
      _recurs_rm_files(p);
      m_file_log << std::endl;
    }

    void ImportBase::_recurs_rm_files(const fs::path &p) const {
      auto it = fs::directory_iterator(p);
      auto end = directory_iterator()
      for(; it != end; ++it) {
        if(fs::is_regular_file(*it)) {
          m_file_log << "rm " << *it << std::endl;
          fs::remove(*it);
        }
        else {
          _recurs_rm_file(*it);
          m_file_log << "rm " << *it << std::endl;
        }
      }
      m_file_log << "rm " << p << std::endl;
    }

    /// \brief Copy files in the same directory as properties.calc.json into the
    ///        traning_data directory for a particular configuration
    ///
    /// - First: calc_props_path = _calc_properties_path(pos_path) to get properties.calc.json location
    /// - If calc_props_path.empty(), return
    /// - else if !m_copy_additional_files copy properties.calc.json file only and return
    /// - else, recursively copy all files from calc_props_path.remove_filename()
    ///   to the training data directory for the current calctype
    void ImportBase::_cp_files(const fs::path &pos_path, const std::string &configname) const {
      fs::path p = _calc_dir(configname);
      if(!fs::exists(p)) {
        fs::create_directory(p);
      }

      fs::path calc_props_path = _calc_properties_path(pos_path);
      if(calc_props_path.empty()) {
        return;
      }

      m_file_log.custom(std::string("Copy calculation files: ") + configname);
      if(!m_copy_additional_files) {
        m_file_log << "cp " << calc_props_path << " " << p << std::endl;
        fs::copy_file(calc_props_path, p);
      }
      else {
        _recurs_cp_files(calc_props_path.remove_filename(), p);
      }
      m_file_log << std::endl;
    }

    void ImportBase::_recurs_cp_files(const fs::path &from_dir, const fs::path &to_dir) const {
      auto it = fs::directory_iterator(from_dir);
      auto end = directory_iterator()
      for(; it != end; ++it) {
        if(fs::is_regular_file(*it)) {
          m_file_log << "cp " << *it << " " << to_dir << std::endl;
          fs::copy_file(*it, to_dir);
        }
        else {
          fs::path new_to_dir = to_dir / it->filename();
          fs::create_directory(new_to_dir);
          _recurs_cp_file(*it, new_to_dir);
        }
      }
    }

    void ImportBase::_import_report(std::vector<Result> &results) {

      // map_fail: could not map
      // map_success: could map
      // import_data_fail: would import but couldn't (score < best_score && !data_results.count(from))
      // import_data_conflicts: conflicts with other in import batch && preexisting
      // - pos, config, score_method, chosen?, overwrite?, import data?, import additional files?, score, best_score, is_preexisting?


      // list of structures that could not be mapped
      std::vector<Result> map_fail;

      // list of structures that could be mapped
      std::vector<Result> map_success;

      // list of structures that would be imported except preexisting data prevents it
      std::vector<Result> import_data_fail;

      for(long i = 0; i < results.size(); ++i) {
        const auto &res = results[i];
        if(res.mapped_props.to.empty()) {
          map_fail.push_back(res);
        }
        else {
          map_success.push_back(res);
          if(res.has_data && db_props().score(res.mapped_props) < db_props().best_score(res.mapped_props.to)) {
            import_data_fail.push_back(res);
          }
        }
      }

      // list of conflicts (multiple config with same 'from')
      std::map<std::string, std::vector<long> > conflict_count;
      std::vector<Result> conflict;

      for(long i = 0; i < results.size(); ++i) {
        auto it = conflict_count.find(res.mapped_props.from);
        if(it == conflict_count.end()) {
          conflict_count[res.mapped_props.from] = std::vector<long>(1, i);
        }
        else {
          it->second.push_back(i);
        }
      }
      for(const auto &val : conflict_count) {
        if(val.second.size() > 1) {
          for(const auto &i : val.second) {
            conflict.push_back(results[i]);
          }
        }
      }


      // output a 'batch' file with paths to structures that could not be imported
      if(map_fail.size()) {

        fs::path p = m_report_dir / "import_map_fail";
        fs::ostream sout(p);

        log() << "WARNING: Could not import " << map_fail.size() << " structures." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        DataFormatterDictionary<Result> dict(_pos(), _fail_msg());
        auto formatter = m_dict.parse({"pos", "fail_msg"})
                         sout << formatter(map_fail.begin(), map_fail.end());
      }

      // - pos, config, score_method, import data?, import additional files?, score, best_score, is_preexisting?
      auto formatter = this->_import_formatter(dict, data_results);

      if(map_success.size()) {

        fs::path p = m_report_dir / "import_map_success";
        fs::ostream sout(p);

        log() << "Successfully imported " << map_success.size() << " structures." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(map_success.begin(), map_successs.end());
      }

      if(import_data_fail.size()) {

        fs::path p = m_report_dir / "import_data_fail";
        fs::ostream sout(p);

        log() << "WARNING: Did not import data from "
              << import_data_fail.size() << " structures which have are a mapping score"
              " better than the existing data." << std::endl;
        log() << "  You may wish to inspect these structures and allow overwriting "
              "or remove existing data manually." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(import_data_fail.begin(), import_data_fail.end());
      }

      if(conflict.size()) {
        fs::path p = m_report_dir / "import_conflict";
        fs::ostream sout(p);

        log() << "WARNING: Imported data from structures that mapped to the same configuration." << std::endl
              << "  Data can only be imported from one of the conflicting structures." << std::endl
              << "  Based on the current conflict resolution method the 'best' result was automatically chosen, " << std::endl
              << "  but you may wish to inspect these results and manually select which structures to import." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(conflict.begin(), conflict.end());
      }
    }


    GenericDatumFormatter<std::string, Result> ImportBase::_pos() const {
      return GenericDatumFormatter<bool, Result>("pos", "", [&](const Result & res) {
        return res.pos.string();
      });
    }

    GenericDatumFormatter<std::string, Result> ImportBase::_fail_msg() const {
      return GenericDatumFormatter<bool, Result>("fail_msg", "", [&](const Result & res) {
        return res.fail_msg;
      });
    }

    /// Use 'from_configname' as 'configname'
    GenericDatumFormatter<std::string, Result> ImportBase::_configname() const {
      return GenericDatumFormatter<bool, Result>("configname", "", [&](const Result & res) {
        return res.mapped_props.from;
      });
    }

    GenericDatumFormatter<std::string, Result> ImportBase::_from_configname() const {
      return GenericDatumFormatter<bool, Result>("from_configname", "", [&](const Result & res) {
        return res.mapped_props.from;
      });
    }

    GenericDatumFormatter<std::string, Result> ImportBase::_to_configname() const {
      return GenericDatumFormatter<bool, Result>("to_configname", "", [&](const Result & res) {
        return res.mapped_props.to;
      });
    }

    GenericDatumFormatter<bool, Result> ImportBase::_has_data() const {
      return GenericDatumFormatter<bool, Result>("has_data", "", [&](const Result & res) {
        return res.has_data;
      });
    }

    GenericDatumFormatter<bool, Result> ImportBase::_has_complete_data() const {
      return GenericDatumFormatter<bool, Result>("has_complete_data", "", [&](const Result & res) {
        return res.has_data;
      });
    }

    GenericDatumFormatter<bool, Result> ImportBase::_preexisting_data(const std::map<std::string, ImportData> &data_results) const {
      return GenericDatumFormatter<bool, Result>(
               "preexisting_data",
               "",
      [&](const Result & res) {
        return data_results.find(res.mapped_props.from)->preexisting;
      },
      [&](const Result & res) {
        return data_results.find(res.mapped_props.from) != data_results.end();
      });
    }

    GenericDatumFormatter<bool, Result> ImportBase::_import_data(const std::map<std::string, ImportData> &data_results) const {
      return GenericDatumFormatter<bool, Result>(
               "import_data",
               "",
      [&](const Result & res) {
        return data_results.find(res.mapped_props.from)->last_insert == res.pos;
      },
      [&](const Result & res) {
        return data_results.count(res.mapped_props.from) != 0;
      });
    }

    GenericDatumFormatter<bool, Result> ImportBase::_import_additional_files(const std::map<std::string, ImportData> &data_results) const {
      return GenericDatumFormatter<bool, Result>(
               "import_additional_files",
               "",
      [&](const Result & res) {
        auto it = data_results.find(res.mapped_props.from);
        if(it != data_results.end()) {
          return it->copy_more;
        }
        return false;
      });
    }

    GenericDatumFormatter<double, Result> ImportBase::_lattice_deformation_cost() const {
      return GenericDatumFormatter<bool, Result>(
               "lattice_deformation_cost", "",
      [&](const Result & res) {
        return res.mapped_props.mapped["lattice_deformation_cost"];
      },
      [&](const Result & res) {
        return res.mapped_props.mapped.contains("lattice_deformation_cost");
      });
    }

    GenericDatumFormatter<double, Result> ImportBase::_basis_deformation_cost() const {
      return GenericDatumFormatter<bool, Result>(
               "basis_deformation_cost", "",
      [&](const Result & res) {
        return res.mapped_props.mapped["basis_deformation_cost"];
      },
      [&](const Result & res) {
        return res.mapped_props.mapped.contains("basis_deformation_cost");
      });
    }

    GenericDatumFormatter<double, Result> ImportBase::_relaxed_energy() const {
      return GenericDatumFormatter<bool, Result>(
               "relaxed_energy", "",
      [&](const Result & res) {
        return res.mapped_props.mapped["relaxed_energy"];
      },
      [&](const Result & res) {
        return res.mapped_props.mapped.contains("relaxed_energy");
      });
    }

    GenericDatumFormatter<double, Result> ImportBase::_score() const {
      return GenericDatumFormatter<bool, Result>(
               "score", "",
      [&](const Result & res) {
        return this->db_props().score(res.mapped_props);
      },
      [&](const Result & res) {
        return res.has_data;
      });
    }

    GenericDatumFormatter<double, Result> ImportBase::_best_score() const {
      return GenericDatumFormatter<bool, Result>(
               "score", "",
      [&](const Result & res) {
        return this->db_props().best_score(res.mapped_props.to);
      },
      [&](const Result & res) {
        return this->db_props().find_via_to(res.mapped_props.to) != this->db_props().end();
      });
    }

    GenericDatumFormatter<double, Result> ImportBase::_is_best() const {
      return GenericDatumFormatter<bool, Result>(
               "is_best", "",
      [&](const Result & res) {
        return res.mapped_props.from == this->db_props().relaxed_from(res.mapped_props.to);
      },
      [&](const Result & res) {
        return this->db_props().find_via_to(res.mapped_props.to) != this->db_props().end();
      });
    }

    GenericDatumFormatter<std::string, Result> ImportBase::_selected() const {
      return GenericDatumFormatter<std::string, Result>("selected", "", [&](const Result & res) {
        return false;
      });
    }

    /// Insert default formatters to dictionary
    void ImportBase::_default_formatters(
      DataFormatterDictionary<Result> &dict,
      std::map<std::string, ImportResult> &data_results) {

      dict.insert(
        _pos(),
        _fail_msg(),
        _configname(),
        _from_configname(),
        _to_configname(),
        _has_data(),
        _has_complete_data(),
        _prexisting_data(data_results),
        _import_data(data_results),
        _import_additional_files(data_results),
        _lattice_deformation_cost(),
        _basis_deformation_cost(),
        _relaxed_energy(),
        _score(),
        _best_score(),
        _is_best(),
        _selected());
    }


    // --- Configuration specializations ---------------------------------------

    /*

    still need to devide interface for setting conflict score methods

    "    conflict_score: JSON object (optional, default= see below)\n"
    "        Which metric should be used to determine which calculation results \n"
    "        should be used for a particular configuration if multiple results\n"
    "        map to the same configuration and the self-mapping result is not \n"
    "        available. Should consist of a \"method\" and method-dependent \n"
    "        parameters. The \"method\" options and associated parameters are: \n"

    "          \"deformation_cost\":\n"
    "             \"lattice_weight\": number, in range [0, 1.0]\n"

    "             Uses a weighted sum of cost functions for lattice and basis \n"
    "             deformation. See below for complete definition.\n"

    "          \"minimum\":\n"
    "             \"property\": property name (ex: \"relaxed_energy\")\n"

    "             Reads the specified property from the mapped properties and\n"
    "             selects the minimum to be the best mapping.\n"

    "          \"maximum\":\n"
    "             \"property\": property name (ex: \"some property\")\n"

    "             Reads the specified property from the mapped properties and\n"
    "             selects the maximum to be the best mapping.\n"

    "          \"direct_selection\":\n"
    "             \"name\": configname to force as 'best' (ex: \"SCEL3_1_1_3_0_0_0/4\")\n"

    "             Directly specify which configuration's properties should\n"

    "        The default value is:\n"
    "          \"conflict_score\" : {\n"
    "             \"method\": \"minimum\",\n"
    "             \"property\": \"relaxed_energy\"\n"
    "          }\n\n"

    "  Deformation cost:\n"

    "        The \"deformation_cost\" is:\n\n"

    "            deformation_cost = w*lattice_deformation_cost + \n"
    "                               (1.0-w)*basis_deformation_cost,\n\n"

    "        where \"w\" is the \"lattice_weight\" factor (default=0.5),\n\n"

    "        the \"lattice_deformation_cost\" is the mean-square displacement of a \n"
    "        points on the surface of a sphere of volume equal to the atomic volume \n"
    "        when  it is deformed by the volume preserving deviatoric deformation \n"
    "        matrix, F_deviatoric:\n\n"

    "            V = relaxed_atomic_volume;\n"
    "            F = deformation matrix (lattice_relaxed = F*lattice_ideal);\n"
    "            F_deviatoric = F/pow(F.determinant(), 1./3.);\n"
    "            I = 3x3 identity matrix;\n\n"

    "            lattice_deformation_cost = pow( 3.*V / (4.*pi), 2.0/3.0) / 3.0 * \n"
    "                (0.5 * (F.t * F / pow(std::abs(F.determinant()), 2.0/3.0) - I)).squaredNorm()\n\n"

    "        and the \"basis_deformation_cost\" is a cost function for the amount of\n"
    "        basis site relaxation:\n\n"

    "            D = 3xN matrix of basis site displacements (displacements are applied \n"
    "                before strain)\n"
    "            Natoms = number of atoms in configuration\n"
    "            basis_deformation_cost = (F*D * D.transpose() * F.transpose()).trace() \n"
    "                / (max(Natoms, 1.))\n\n"

    */

    /// \brief Constructor
    Import<Configuration>::Import(
      const PrimClex &primclex,
      const ConfigMapper &configmapper,
      std::vector<std::string> dof,
      bool primitive_only,
      bool import_data,
      bool copy_additional_files,
      bool overwrite,
      fs::path report_dir) :

      ImportT(primclex, import_data, copy_additional_files, overwrite, report_dir),
      m_configmapper(configmapper),
      m_dof(dof),
      m_primitive_only(primitive_only) {}

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
    import_inserter Import<Configuration>::_import(
      fs::path p,
      DataBaseIterator<Configuration> hint,
      import_inserter result,
      bool force) override {

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
      fs::path prop_path = _calc_properties_path(pos);

      Result res;
      res.pos = pos;

      if(hint != db_config().end()) {
        res.mapped_props.from = hint.name();
      }

      // read from structure file or properties.calc.json file (if exists)
      BasicStructure<Site> struc = prop_path.empty() ? _make_structure(pos) : _make_structure(prop_path);

      // do mapping
      ConfigMapperResult map_result;
      if(_occupation_only()) {
        map_result = m_configmapper.import_structure_occupation(struc), hint);
      }
      else {
        map_result =  m_configmapper.import_structure(struc), hint);
      }

      if(!map_result.success) {
        res.fail_msg = map_result.fail_msg;
        *result++ = res;
        return result;
      }

      // insert in database (note that this also/only inserts primitive)
      ConfigInsertResult insert_result = map_result.config->insert(m_primitive_only);
      res.is_new_config = insert_result.insert_canonical;
      res.mapped_props.to = insert_result.canonical_it.name();

      // check for and read raw 'unmapped' data, adds 'data_timestamp'
      if(!prop_path.empty()) {
        std::tie(res.mapped_props.unmapped, res.has_data, res.has_complete_data) =
          Configuration::read_calc_properties(primclex(), prop_path);
      }

      // copy relaxation properties from best config mapping into 'mapped' props
      auto it = map_result.relaxation_properties["best_mapping"].begin();
      auto end = map_result.relaxation_properties["best_mapping"].end();
      for(; it != end; ++it) {
        res.mapped_props.mapped[it.name()] = *it;
      }
      res.mapped_props["best_assignment"] = map_result.best_assignment;
      res.mapped_props["cart_op"] = map_result.cart_op;

      // at this point, the mapped structure result is complete
      *result++ = res;

      // it may be the structure was not primitive:
      // - in which case we need to create a result indicating that the primitive
      //   was also inserted in the database,
      // - but don't try to scale the data for the primitive
      if(insert_result.canonical_it != insert_result.primitive_it) {
        Result prim_res;
        prim_res.pos = pos;
        prim_res.mapped_props.from = insert_result.primitive_it.name();
        prim_res.mapped_props.to = insert_result.primitive_it.name();
        prim_res.is_new_config = insert_result.insert_primitive;
        *result++ = prim_res;
      }

      return result;
    }


    const std::string Import<Configuration>::import_help =

      "Import Configuration: \n\n"

      "  'casm import' of Configuration proceeds in two steps: \n\n"

      "  1) For each file: \n"
      "     - Read structure from VASP POSCAR type file or CASM properties.calc.json \n"
      "       file. \n"
      "     - Map the structure onto a Configuration of the primitive crystal \n"
      "       structure. \n"
      "     - Record relaxation data (lattice & basis deformation cost). \n\n"
      "
      "  2) If data import is requested, iterate over each import record and do \n"
"     the following:
      \n"
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
      "           it is assumed to be a 'properties.calc.json' file.
      "        2) If / path / to / pos, checks for / path / to / calctype.current / properties.calc.json\n"
        "        3) If / path / to / pos, checks for / path / to / properties.calc.json \n\n"

          "        If false, only configuration structure is imported.\n\n"

"    import_additional_files :
          bool (optional, default = false)\n"
          "        If true, attempt to copy all files & directories in the same directory\n"
          "        as the structure file or , if it exists, the properties.calc.json file.\n"
            "        Files & directories will only by copied if there are no existing files\n"
              "        or directories in the 'training_data' directory for the configuration \n"
                "        the structure as been mapped to or \"overwrite\"=true.\n\n"

                "    overwrite: bool (optional, default=false)\n"
                "        If true, data and files will be imported that overwrite existing\n"
                "        data and files, if the score calculated by the \"conflict_score\"\n"
                "        for the configuration being mapped to will be improved.\n\n";


                int Import<Configuration>::import(
                  PrimClex &primclex,
                  const jsonParser &kwargs,
                  const Completer::ImportOption &import_opt) {

                // -- collect input settings --

                po::variables_map &vm = import_opt.vm();
                  ConfigMapper configmapper;
                  jsonParser used;
                  jsonParser _default;

                  // general settings
                  bool primitive_only;
                  kwargs.get_else(primitive_only, "primitive_only", false);
                  used["primitive_only"] = primitive_only;

                  // get input report_dir, check if exists, and create new report_dir.i if necessary
                  fs::path report_dir = primclex.dir().root() / "import_report");
                  report_dir = ImportBase::create_report_dir(report_dir);

                  // 'mapping' subsettings are used to
                  std::tie(configmapper, used["mapping"]) = Import<Configuration>::_make_configmapper(primclex, kwargs);

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
                  used["data"]["import"] = import;

                  bool import_additional_files;
                  data.get_else(copy_additional_files, "import_additional_files", false);
                  used["data"]["import_additional_files"] = import_additional_files;

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
                    configmapper,
                  {"occupation"}, // still need to figure out how to specify this in general
                  primitive_only,
                  import_data,
                  import_additional_files,
                  overwrite,
                  report_dir);

                  // -- read structure file paths --
                  std::vector<fs::path> pos;
                  auto res = ImportBase::construct_pos_paths(primclex, import_opt, std::back_inserter(pos));
                  if(res.second) {
                  return res.second;
                }

                // -- read structure file paths --
                f.import(pos.begin(), pos.end());

                  return 0;
                }

    const std::string Import<Configuration>::update_help =

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
      "     \"conflict_score\" is used for the relaxed configuration.\n\n"

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
      "        robust methods will be used\n\n";


    int Import<Configuration>::update(
      PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::UpdateOption &update_opt) {

      // -- collect input settings --

      po::variables_map &vm = update_opt.vm();
      ConfigMapper configmapper;
      jsonParser used;
      jsonParser _default;

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
      fs::path report_dir = primclex.dir().root() / "update_report");
      report_dir = ImportBase::create_report_dir(report_dir);

      // 'mapping' subsettings are used to
      std::tie(configmapper, used["mapping"]) = Import<Configuration>::_make_configmapper(primclex, kwargs);

      // 'data' subsettings
      bool import_data = true;
      bool import_additional_files = false;
      bool overwrite = false;

      // -- print used settings --
      Log &log = primclex.log();
      log.read("Settings");
      log << used << std::endl << std::endl;

      // -- construct Import --
      Import<Configuration> f(
        primclex,
        configmapper,
      {"occupation"}, // still need to figure out how to specify this in general
      primitive_only,
      import_data,
      import_additional_files,
      overwrite,
      report_dir);

      // -- read selection --
      DB::Selection<Configuration> sel(primclex, update_opt.selection_path());

      // -- update --
      f.update(sel, force);

      return 0;
    }

    const std::string Import<Configuration>::remove_help = "ToDo";

    int Import<Configuration>::remove(
      PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::RemoveOption &remove_opt) {
      return 0;
    }

    /// Construct ConfigMapper from input args
    std::pair<ConfigMapper, jsonParser> Import<Configuration>::_make_configmapper(
      const PrimClex &primclex,
      const jsonParser &kwargs) {

      // -- read settings --
      bool rotate = true;
      bool strict = false;
      bool ideal;
      kwargs.get_else(ideal, kwargs["ideal"], false);

      double lattice_weight;
      kwargs.get_else(lattice_weight, kwargs["lattice_weight"], 0.5);

      double max_vol_change;
      kwargs.get_else(max_vol_change, kwargs["max_vol_change"], 0.3);

      double min_va_frac;
      kwargs.get_else(min_va_frac, kwargs["min_va_frac"], 0.0);

      double max_va_frac;
      kwargs.get_else(max_va_frac, kwargs["max_va_frac"], 0.5);

      // -- collect settings used --
      jsonParser used;
      used["ideal"] = ideal;
      used["lattice_weight"] = lattice_weight;
      used["max_vol_change"] = max_vol_change;
      used["min_va_frac"] = min_va_frac;
      used["max_va_frac"] = max_va_frac;

      // -- construct ConfigMapper --
      int map_opt = ConfigMapper::none;
      if(rotate) map_opt |= ConfigMapper::rotate;
      if(strict) map_opt |= ConfigMapper::strict;
      if(!ideal) map_opt |= ConfigMapper::robust;

      ConfigMapper configmapper(
        primclex,
        lattice_weight,
        max_vol_change,
        map_opt,
        tol);
      configmapper.set_min_va_frac(min_va_frac);
      configmapper.set_max_va_frac(max_va_frac);

      return std::make_pair(configmapper, used);
    }

    /// Allow ConfigType to specialize the report formatting for 'import'
    DataFormatter<Result> Import<Configuration>::_import_formatter(
      std::map<std::string, ImportResult> &data_results) const {

      DataFormatterDictionary<Result> dict;
      _default_formatters(dict, data_results);

      std::vector<std::string> col = {
        "configname", "selected", "pos", "has_data", "has_complete_data",
        "import_data", "import_additional_files", "score", "best_score", "is_best",
        "lattice_deformation_cost", "basis_deformation_cost", "deformation_cost",
        "relaxed_energy"
      };

      return m_dict.parse(col);
    }

    // Allow ConfigType to specialize the report formatting for 'update'
    DataFormatter<Result> Import<Configuration>::_update_formatter(
      std::map<std::string, ImportResult> &data_results) const {

      DataFormatterDictionary<Result> dict;
      _default_formatters(dict, data_results);

      std::vector<std::string> col = {
        "configname", "selected", "to_configname", "has_data", "has_complete_data",
        "score", "best_score", "is_best",
        "lattice_deformation_cost", "basis_deformation_cost", "deformation_cost",
        "relaxed_energy"
      };

      return m_dict.parse(col);
    }

    /// \brief Read BasicStructure<Site> to be imported
    ///
    /// If 'p.extension()' == ".json" or ".JSON", read as properties.calc.json
    /// Else, read as VASP POSCAR
    BasicStructure<Site> Import<Configuration>::_make_structure(const fs::path &p) {

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
    bool Import<Configuration>::_occupation_only() const {
      return m_dof.size() == 1 && m_dof[0] == "occupation";
    }

    // --- DiffTransConfiguration specializations --------------------------------

    /*
    /// \brief Constructor
    Import<DiffTransConfiguration>::Import(
      const PrimClex& primclex,
      const DiffTransConfigMapper& configmapper,
      bool import_data,
      bool copy_additional_files,
      bool overwrite,
      fs::path report_dir) :

      ImportT(primclex, import_data, copy_additional_files, overwrite, report_dir) :
      m_configmapper(configmapper) {}

    const std::string Import<DiffTransConfiguration>::import_help = "ToDo";

    int Import<DiffTransConfiguration>::import(
      PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::ImportOption &import_opt) {

    }

    const std::string Import<DiffTransConfiguration>::update_help = "ToDo";

    int Import<DiffTransConfiguration>::update(
      PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::UpdateOption &update_opt) {

    }

    const std::string Import<DiffTransConfiguration>::remove_help = "ToDo";

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
      std::map<std::string, ImportResult>& data_results) const {

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
      std::map<std::string, ImportResult>& data_results) const {

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
