#include <cstring>
#include <tuple>

#include "casm/app/casm_functions.hh"
#include "casm/crystallography/jsonStruc.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/casm_io/VaspIO.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {
  namespace Import_impl {
    // record some pieces of data from each import.  These are used to resolve conflicts at the end
    // the 'relaxjson' datarecord stores relaxation properties that will be merged into Configuration::calc_properties() during the final step
    enum DataType {path = 0, relaxjson = 1, contcar = 2, energy = 3};
    using Data = std::tuple<fs::path, jsonParser, std::string, std::pair<bool, double> >;
  }

  namespace Completer {
    ImportOption::ImportOption(): OptionHandlerBase("import") {}

    double ImportOption::lattice_weight() const {
      return m_lattice_weight;
    }

    double ImportOption::vol_tolerance() const {
      return m_vol_tolerance;
    }

    double ImportOption::min_va_frac() const {
      return m_min_va_frac;
    }

    double ImportOption::max_va_frac() const {
      return m_max_va_frac;
    }
    const std::vector<fs::path> &ImportOption::pos_vec() const {
      return m_pos_vec;
    }

    const fs::path &ImportOption::batch_path() const {
      return m_batch_path;
    }

    void ImportOption::initialize() {
      add_help_suboption();

      m_desc.add_options()
      ("pos,p", po::value<std::vector<fs::path> >(&m_pos_vec)->multitoken()->value_name(ArgHandler::path()), "Path(s) to structure(s) being imported (multiple allowed, but no wild-card matching)")
      ("cost-weight,w", po::value<double>(&m_lattice_weight)->default_value(0.5),
       "Adjusts cost function for mapping optimization (cost=w*lattice_deformation+(1-w)*basis_deformation).")
      ("min-energy", "Resolve mapping conflicts based on energy rather than deformation.")
      ("max-vol-change", po::value<double>(&m_vol_tolerance)->default_value(0.25),
       "Adjusts range of SCEL volumes searched while mapping imported structure onto ideal crystal (only necessary if the presence of vacancies makes the volume ambiguous). Default is +/- 25% of relaxed_vol/prim_vol. Smaller values yield faster import, larger values may yield more accurate mapping.")
      ("max-va-frac", po::value<double>(&m_max_va_frac)->default_value(0.5),
       "Places upper bound on the fraction of sites that are allowed to be vacant after imported structure is mapped onto the ideal crystal. Smaller values yield faster execution, larger values may yield more accurate mapping. Has no effect if supercell volume can be inferred from the number of atoms in the structure. Default value allows up to 50% of sites to be vacant.")
      ("min-va-frac", po::value<double>(&m_min_va_frac)->default_value(0.),
       "Places lower bound on the fraction of sites that are allowed to be vacant after imported structure is mapped onto the ideal crystal. Nonzero values may yield faster execution if updating configurations that are known to have a large number of vacancies, at potential sacrifice of mapping accuracy.  Has no effect if supercell volume can be inferred from the number of atoms in the structure. Default value allows as few as 0% of sites to be vacant.")
      ("batch,b", po::value<fs::path>(&m_batch_path)->value_name(ArgHandler::path()), "Path to batch file, which should list one structure file path per line (can be used in combination with --pos)")
      ("rotate,r", "Rotate structure to be consistent with setting of PRIM")
      ("ideal,i", "Assume imported structures are unstrained (ideal) for faster importing. Can be slower if used on deformed structures, in which case more robust methods will be used")
      //("strict,s", "Request that symmetrically equivalent configurations be treated as distinct.")
      ("data,d", "Attempt to extract calculation data (properties.calc.json file) from the enclosing directory of the structure files, if it is available")
      ("copy-additional-files", "Recursively copy other files from the same directory as the properties.calc.json file.");

      return;
    }

  }

  bool _has_existing_files(fs::path p);
  fs::path _calc_properties_path(const PrimClex &primclex, fs::path pos_path);
  void _cp_files(const fs::path &pos_path, const Configuration &config, bool copy_additional_files, Log &log);
  void _recurs_cp_files(const fs::path &from_dir, const fs::path &to_dir, Log &log);

  // ///////////////////////////////////////
  // 'import' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  /// Import proceeds in two steps.
  ///   1) read each file, map it onto a Configuration of the PrimClex
  ///       - record relaxation data for each one
  ///
  ///   2) If data import was requested, iterate over each import record and do the following:
  ///       - if multiple imported structures map onto a configuration for which there is no calculation data,
  ///         import calculation data from the structure with the lowest mapping cost
  ///       - if one or more imported structuress map onto a configuration for which calculation
  ///         data already exist, do not import any new data
  ///       - if data is imported, the corresponding properties.calc.json file is copied into the directory of the
  ///         mapped configuration. A structure file, relaxed_structure.vasp is also written to the directory.
  ///       - relaxed_structure.vasp gives the relaxed structure in a setting and orientation that matches the
  ///         generated POS file
  ///
  int import_command(const CommandArgs &args) {

    double tol = TOL;
    COORD_TYPE coordtype = FRAC;

    double vol_tol;
    double lattice_weight;
    std::vector<fs::path> pos_paths;
    //fs::path dft_path, batch_path;
    fs::path batch_path;
    po::variables_map vm;
    /// Set command line options using boost program_options
    Completer::ImportOption import_opt;

    try {

      po::store(po::parse_command_line(args.argc, args.argv, import_opt.desc()), vm); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        args.log << std::endl;
        args.log << import_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << import_opt.desc() << std::endl;

        args.log << "DESCRIPTION" << std::endl;
        args.log << "    Import structure specified by --pos. If it doesn't exist make a directory for it and copy data over" << std::endl;
        args.log << "    If a *.json file is specified, it will be interpreted as a 'calc.properties.json' file." << std::endl;
        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      vol_tol = import_opt.vol_tolerance();
      lattice_weight = import_opt.lattice_weight();
      pos_paths = import_opt.pos_vec();
      batch_path = import_opt.batch_path();

    }
    catch(po::error &e) {
      args.err_log << import_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return 3;
    }
    catch(std::exception &e) {
      args.err_log << import_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return 4;

    }

    if(!vm.count("pos") && !vm.count("batch")) {
      args.err_log << import_opt.desc() << std::endl;
      args.err_log << "No structures specified for import (specify structures using --pos or --batch)." << std::endl;
      return 5;
    }

    //read all the import paths
    if(vm.count("batch")) {
      if(!fs::exists(batch_path)) {
        args.err_log << "ERROR: Batch import file does not exist at " << batch_path << "\n";
        return 6;
      }
      fs::ifstream batchfile(batch_path);
      fs::path tpath;
      while(batchfile >> tpath) {
        pos_paths.push_back(tpath);
        batchfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      if(tpath != pos_paths.back() && tpath.string().size() != 0 && fs::exists(tpath))
        pos_paths.push_back(tpath);
    }
    {
      if(pos_paths.size() == 0) {
        args.err_log <<   "ERROR: No files specified for import.\n";
        if(vm.count("batch"))
          args.err_log << "       Check batch file for errors.\n";
        return 7;
      }
      bool missing_files(false);
      for(auto it = pos_paths.cbegin(); it != pos_paths.cend(); ++it) {
        if(!fs::exists(*it)) {
          if(!missing_files)
            args.err_log << "*** ERROR: Missing file(s):\n";
          missing_files = true;
          args.err_log   << "           " << fs::absolute(*it) << "\n";
        }
      }
      if(missing_files)
        return 8;
    }

    COORD_MODE C(coordtype);

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);
    fs::path pwd = fs::current_path();

    int map_opt = ConfigMapper::none;
    if(vm.count("rotate")) map_opt |= ConfigMapper::rotate;
    if(vm.count("strict")) map_opt |= ConfigMapper::strict;
    if(!vm.count("ideal")) map_opt |= ConfigMapper::robust;
    ConfigMapper configmapper(primclex, lattice_weight, vol_tol, map_opt, tol);
    configmapper.set_min_va_frac(import_opt.min_va_frac());
    configmapper.set_max_va_frac(import_opt.max_va_frac());


    // import_map keeps track of mapping collisions -- only used if vm.count("data")
    // import_map[config_name] gives a list all the configuration paths that mapped onto configuration 'config_name' :  import_map[config_name][i].first
    //                         along with a list of the mapping properties {lattice_deformation, basis_deformation}  :  import_map[config_name][i].second
    std::map<Configuration *, std::vector<Import_impl::Data> > import_map;
    std::vector<std::string > error_log;
    Index n_unique(0);
    // iterate over structure files
    args.log << "  Beginning import of " << pos_paths.size() << " configuration" << (pos_paths.size() > 1 ? "s" : "") << "...\n" << std::endl;
    for(auto it = pos_paths.begin(); it != pos_paths.end(); ++it) {
      if(it != pos_paths.begin())
        args.log << "\n***************************\n" << std::endl;

      fs::path pos_path, import_path;

      pos_path = fs::absolute(*it, pwd);
      std::string imported_name;
      BasicStructure<Site> import_struc;

      // If user requested data import, try to get structural data from properties.calc.json, instead of POS, etc.
      // Since properties.calc.json would be used during 'casm update' to validate relaxation
      if(vm.count("data")) {
        fs::path dft_path = _calc_properties_path(primclex, pos_path);
        if(!dft_path.empty()) {
          pos_path = dft_path;
        }
      }

      //Import structure and make note of path
      jsonParser relax_data;
      std::pair<bool, double> checkenergy(false, 0.0);

      try {
        if(pos_path.extension() == ".json" || pos_path.extension() == ".JSON") {
          jsonParser datajson(pos_path);
          if(datajson.contains("relaxed_energy")) {
            checkenergy = std::pair<bool, double>(true, datajson["relaxed_energy"].get<double>());
          }
          from_json(simple_json(import_struc, "relaxed_"), datajson);
        }
        else {
          fs::ifstream struc_stream(pos_path);
          import_struc.read(struc_stream);
        }

        Eigen::Matrix3d cart_op;
        std::vector<Index> best_assignment;
        jsonParser fullrelax_data;
        if(configmapper.import_structure_occupation(import_struc, imported_name, fullrelax_data, best_assignment, cart_op, true)) {
          args.log << "  " << pos_path << "\n  was imported successfully as " << imported_name << std::endl << std::endl;
          n_unique++;
        }
        else {
          args.log << "  " << pos_path << "\n  mapped onto pre-existing equivalent structure " << imported_name << std::endl << std::endl;
        }
        relax_data = fullrelax_data["best_mapping"];
        args.log << "  Relaxation stats -> lattice_deformation = " << relax_data["lattice_deformation"].get<double>()
                 << "      basis_deformation = " << relax_data["basis_deformation"].get<double>() << std::endl << std::endl;;
      }
      catch(std::exception &e) {
        args.err_log << "  ERROR: Unable to import " << pos_path << " because \n"
                     << "    -> " << e.what() << "\n\n";
        error_log.push_back(it->string() + "\n     -> " + e.what());
        if(it != pos_paths.cend()) {
          args.log << "  Continuing...\n";
        }
        continue;
      }

      Configuration &imported_config(primclex.configuration(imported_name));

      imported_config.push_back_source(jsonParser(std::make_pair("import_mapped", pos_path.string())));



      //Exit if user did not request not to copy data
      if(!vm.count("data"))
        continue;

      std::stringstream contcar_ss;
      VaspIO::PrintPOSCAR(import_struc).print(contcar_ss);
      import_map[&imported_config].push_back(Import_impl::Data(pos_path, relax_data, contcar_ss.str(), checkenergy));

    }

    // All the mapping is finished; now we migrate data, if requested
    std::stringstream conflict_log;
    if(vm.count("data")) {
      args.log << "  Attempting to import data..." << std::endl;
      auto it(import_map.begin()), end_it(import_map.end());
      for(; it != end_it; ++it) {
        Configuration &imported_config = *(it->first);
        std::vector<Import_impl::Data> &data_vec(it->second);

        fs::path import_path = fs::absolute(imported_config.get_path(), pwd);
        bool preexisting(false);
        if(_has_existing_files(imported_config.calc_dir())) {
          preexisting = true;
        }
        Index mult = data_vec.size() + Index(preexisting);
        double best_weight(1e19);
        double best_energy(1e19);
        Index best_conflict(0), best_ind(0);
        if(mult > 1) {
          conflict_log <<  "  CONFLICT -> " <<  mult << " matching structures for config " << imported_config.name() << ": " <<  std::endl;
          double w = lattice_weight;
          Index conflict_ind(1);
          if(preexisting) {
            conflict_ind++;
            conflict_log << "           1) Pre-existing data for " << imported_config.name() << " before import." << std::endl
                         << "              Relaxation stats:" << std::endl;
            if(imported_config.calc_properties().contains("basis_deformation") && imported_config.calc_properties().contains("lattice_deformation")) {
              double ld = imported_config.calc_properties()["lattice_deformation"].get<double>();
              double bd = imported_config.calc_properties()["basis_deformation"].get<double>();
              conflict_log << "                -- lattice_deformation = " << ld << ";  basis_deformation = " << bd << ";  weighted avg = " << w *ld + (1.0 - w)*bd << std::endl;
              if(!vm.count("min-energy")) {
                best_weight = w * ld + (1.0 - w) * bd;
                best_conflict = 0;
                best_ind = -1;
              }
            }
            else {
              conflict_log << "                -- lattice_deformation = unknown;  basis_deformation = unknown;  weighted avg = unknown" << std::endl;
            }
            if(imported_config.calc_properties().contains("relaxed_energy")) {
              conflict_log << "                -- relaxed_energy = " << imported_config.calc_properties()["relaxed_energy"].get<double>() << std::endl;
              if(vm.count("min-energy")) {
                best_energy = imported_config.calc_properties()["relaxed_energy"].get<double>();
                best_conflict = 0;
                best_ind = -1;
              }
            }
            else
              conflict_log << "                -- relaxed_energy = unknown" << std::endl;
            conflict_log << std::endl;
          }
          for(Index i = 0; i < data_vec.size(); i++) {
            fs::path pos_path = std::get<Import_impl::path>(data_vec[i]);
            conflict_log << "           " << conflict_ind++ << ") Structure imported from " << pos_path << "." << std::endl
                         << "              Relaxation stats:" << std::endl;
            jsonParser &relaxjson = std::get<Import_impl::relaxjson>(data_vec[i]);
            double ld = relaxjson["lattice_deformation"].get<double>();
            double bd = relaxjson["basis_deformation"].get<double>();
            conflict_log << "                -- lattice_deformation = " << ld << ";  basis_deformation = " << bd << ";  weighted avg = " << w *ld + (1.0 - w)*bd << std::endl;
            if(std::get<Import_impl::energy>(data_vec[i]).first)
              conflict_log << "                -- relaxed_energy = " << std::get<Import_impl::energy>(data_vec[i]).second << std::endl;
            else
              conflict_log << "                -- relaxed_energy = unknown" << std::endl;
            conflict_log << std::endl;
            if(vm.count("min-energy")) {
              if(std::get<Import_impl::energy>(data_vec[i]).first) {
                if(std::get<Import_impl::energy>(data_vec[i]).second < best_energy) {
                  best_energy = std::get<Import_impl::energy>(data_vec[i]).second;
                  best_conflict = conflict_ind - 1;
                  best_ind = i;
                }
              }
            }
            else {
              if(w * ld + (1.0 - w)*bd < best_weight) {
                best_weight = w * ld + (1.0 - w) * bd;
                best_conflict = conflict_ind - 1;
                best_ind = i;
              }
            }
          }
          if(preexisting) {
            conflict_log << "          ==> Resolution: No data will be imported since data already exists" << std::endl;
            if(valid_index(best_ind)) {
              if(!vm.count("min-energy"))
                conflict_log << "          *** WARNING: Conflicting config #" << best_conflict << " maps more closely onto ideal crystal! ***" << std::endl;
              else
                conflict_log << "          *** WARNING: Conflicting config #" << best_conflict << " has a lower energy! ***" << std::endl;
            }
          }
          else {
            if(!vm.count("min-energy"))
              conflict_log << "          ==> Resolution: Import data from closest match, structure #" << best_conflict  << std::endl;
            else
              conflict_log << "          ==> Resolution: Import data from lowest energy config, structure #" << best_conflict  << std::endl;
          }
          conflict_log << "\n          ----------------------------------------------\n" << std::endl;
        }
        if(preexisting) {
          continue;
        }

        fs::path pos_path = std::get<Import_impl::path>(data_vec[best_ind]);
        if(pos_path.extension() != ".json" && pos_path.extension() != ".JSON") {
          args.log << "  No calculation data was found in the enclosing directory of \n"
                   << "    " << pos_path << std::endl
                   << "  Continuing..." << std::endl;
          continue;
        }

        fs::path import_target = import_path / ("calctype." + primclex.settings().default_clex().calctype);

        _cp_files(pos_path, imported_config, vm.count("copy-additional-files"), args.log);

        if(!fs::exists(imported_config.get_pos_path()))
          imported_config.write_pos();

        {
          fs::ofstream contcar_out(import_target / "relaxed_structure.vasp");
          contcar_out << std::get<Import_impl::contcar>(data_vec[best_ind]);
        }

        jsonParser calc_data;
        if(!imported_config.read_calc_properties(calc_data)) {
          args.log << "  WARNING: Some properties from " << pos_path << " were not valid. Viable values will still be recorded.\n";
        }

        jsonParser &relaxjson = std::get<Import_impl::relaxjson>(data_vec[best_ind]);
        //append relaxjson onto calc_data
        auto jit = relaxjson.cbegin(), jit_end = relaxjson.cend();
        for(; jit != jit_end; ++jit) {
          calc_data[jit.name()] = *jit;
        }
        imported_config.set_calc_properties(calc_data);

        imported_config.push_back_source(jsonParser(std::make_pair("data_inferred_from_mapping", pos_path.string())));

      }
    }
    args.log << "\n***************************\n" << std::endl;

    args.log << "  Finished importing " << pos_paths.size() <<  " structures";
    if(n_unique == 0)
      args.log << " (none of these are new or unique)";
    else if(n_unique < pos_paths.size())
      args.log << " (only " << n_unique << " of these " << (n_unique == 1 ? "is" : "are") << " new and unique)";
    args.log << "." <<  std::endl;

    //Update directories
    args.log << "  Writing SCEL..." << std::endl;
    primclex.print_supercells();
    args.log << "  Writing config_list..." << std::endl << std::endl;
    primclex.write_config_list();
    args.log << "  DONE" << std::endl << std::endl;

    if(error_log.size() > 0) {
      args.log << "  WARNING: --The following paths could not be imported due to errors:\n";
      for(auto it = error_log.cbegin(); it != error_log.cend(); ++it) {
        args.log << *it
                 << "\n        ----------------------------------------------\n" << std::endl;
      }
      args.log << "\n" << std::endl;
    }
    if(conflict_log.str().size()) {
      args.log << "  WARNING: -- The following conflicts were found\n" << std::endl
               << conflict_log.str() << std::endl;

      args.log << "  Please review these conflicts.  A different resolution can be obtained by removing datafiles from\n"
               << "  the training_data directory and performing an import using a manually reduced set of files.\n";
    }
    args.log << "  DONE" << std::endl << std::endl;

    args.log << std::endl;

    return 0;
  };

  bool _has_existing_files(fs::path p) {
    if(!fs::exists(p)) {
      return false;
    }
    return std::distance(fs::directory_iterator(p), fs::directory_iterator());
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
  fs::path _calc_properties_path(const PrimClex &primclex, fs::path pos_path) {

    // check 1: is a JSON file
    if(pos_path.extension() == ".json" || pos_path.extension() == ".JSON") {
      return pos_path;
    }

    // check 2: /path/to/POS -> /path/to/calctype.current/properties.calc.json
    {
      fs::path dft_path = pos_path;
      dft_path.remove_filename();
      (dft_path /= ("calctype." + primclex.settings().default_clex().calctype)) /= "properties.calc.json";
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

  /// \brief Copy files in the same directory as properties.calc.json into the
  ///        traning_data directory for a particular configuration
  ///
  /// - First: calc_props_path = _calc_properties_path(pos_path) to get properties.calc.json location
  /// - If calc_props_path.empty(), return
  /// - else if !copy_additional_files copy properties.calc.json file only and return
  /// - else, recursively copy all files from calc_props_path.remove_filename()
  ///   to the training data directory for the current calctype
  void _cp_files(const fs::path &pos_path, const Configuration &config, bool copy_additional_files, Log &log) {
    fs::path p = config.calc_dir();
    if(!fs::exists(p)) {
      fs::create_directories(p);
    }

    fs::path calc_props_path = _calc_properties_path(config.get_primclex(), pos_path);
    if(calc_props_path.empty()) {
      return;
    }

    log.custom(std::string("Copy calculation files: ") + config.name());
    if(!copy_additional_files) {
      log << "cp " << calc_props_path << " " << p << std::endl;
      fs::copy_file(calc_props_path, p / calc_props_path.filename());
    }
    else {
      _recurs_cp_files(calc_props_path.remove_filename(), p, log);
    }
    log << std::endl;
  }

  void _recurs_cp_files(const fs::path &from_dir, const fs::path &to_dir, Log &log) {
    auto it = fs::directory_iterator(from_dir);
    auto end = fs::directory_iterator();
    for(; it != end; ++it) {
      if(fs::is_regular_file(*it)) {
        log << "cp " << *it << " " << to_dir << std::endl;
        fs::copy_file(*it, to_dir / it->path().filename());
      }
      else {
        fs::path new_to_dir = to_dir / it->path().filename();
        fs::create_directories(new_to_dir);
        _recurs_cp_files(*it, new_to_dir, log);
      }
    }
  }
}
