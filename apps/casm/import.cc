#include "import.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/casm_io/FileSystemInterface.hh"

namespace CASM {

  //I don't like declaring these functions here, but getting it to compile is something I don't want
  //to waste time on right now and it didn't work the normal ways
  //bool copyDirectory(fs::path const &source, fs::path const &destination);
  //bool makeDirectory(boost::filesystem::path dirPath, std::ostream &stream = std::cerr, bool errorWarn = true);

  // ///////////////////////////////////////
  // 'import' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int import_command(int argc, char *argv[]) {

    double tol(TOL);
    COORD_TYPE coordtype = FRAC;
    double vol_tol(0.25);
    double lattice_weight(0.5);
    std::vector<fs::path> pos_paths;
    fs::path dft_path, batch_path;
    bool same_dir(false), no_import(true);
    po::variables_map vm;
    /// Set command line options using boost program_options
    po::options_description desc("'casm import' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    ("pos,p", po::value<std::vector<fs::path> >(&pos_paths)->multitoken(), "Path(s) to structure(s) being imported (multiple allowed, but no wild-card matching)")
    ("cost-weight,w", po::value<double>(&lattice_weight)->default_value(0.5),
     "Adjusts cost function for mapping optimization (cost=w*lattice_deformation+(1-w)*basis_deformation).")
    ("max-vol-change", po::value<double>(&vol_tol)->default_value(0.25),
     "Adjusts range of SCEL volumes searched while mapping imported structure onto ideal crystal (only necessary if the presence of vacancies makes the volume ambiguous). Default is +/- 25% of relaxed_vol/prim_vol. Smaller values yield faster import, larger values may yield more accurate mapping.")
    ("batch,b", po::value<fs::path>(&batch_path), "Path to batch file, which should list one structure file path per line (can be used in combination with --pos)")
    ("rotate,r", "Rotate structure to be consistent with setting of PRIM")
    ("ideal,i", "Assume imported structures are unstrained (ideal) for faster importing. Can be slower if used on deformed structures, in which case more robust methods will be used")
    //("strict,s", "Request that symmetrically equivalent configurations be treated as distinct.")
    ("data,d", "Attempt to extract calculation data from the enclosing directory of the structure files.");

    try {

      po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        std::cout << std::endl;
        std::cout << desc << std::endl;

        std::cout << "DESCRIPTION" << std::endl;
        std::cout << "    Import structure specified by --pos. If it doesn't exist make a directory for it and copy data over" << std::endl;
        std::cout << "    If a *.json file is specified, it will be interpreted as a 'calc.properties.json' file." << std::endl;
        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

    }
    catch(po::error &e) {
      std::cerr << desc << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return 3;
    }
    catch(std::exception &e) {
      std::cerr << desc << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return 4;

    }

    if(!vm.count("pos") && !vm.count("batch")) {
      std::cerr << desc << std::endl;
      std::cerr << "No structures specified for import (specify structures using --pos or --batch)." << std::endl;
      return 5;
    }

    //read all the import paths
    if(vm.count("batch")) {
      if(!fs::exists(batch_path)) {
        std::cerr << "ERROR: Batch import file does not exist at " << batch_path << "\n";
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
        std::cerr <<   "ERROR: No files specified for import.\n";
        if(vm.count("batch"))
          std::cerr << "       Check batch file for errors.\n";
        return 7;
      }
      bool missing_files(false);
      for(auto it = pos_paths.cbegin(); it != pos_paths.cend(); ++it) {
        if(!fs::exists(*it)) {
          if(!missing_files)
            std::cerr << "*** ERROR: Missing file(s):\n";
          missing_files = true;
          std::cerr   << "           " << *it << "\n";
        }
      }
      if(missing_files)
        return 8;
    }

    COORD_MODE C(coordtype);

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cerr << "ERROR in 'casm import': No casm project found." << std::endl;
      return 7;
    }
    fs::current_path(root);


    std::cout << "\n***************************" << std::endl << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;


    int map_opt = ConfigMapper::none;
    if(vm.count("rotate")) map_opt |= ConfigMapper::rotate;
    if(vm.count("strict")) map_opt |= ConfigMapper::strict;
    if(!vm.count("ideal")) map_opt |= ConfigMapper::robust;
    ConfigMapper configmapper(primclex, lattice_weight, vol_tol, map_opt, tol);

    // import_map keeps track of mapping collisions -- only used if vm.count("data")
    // import_map[config_name] gives a list all the configuration paths that mapped onto configuration 'config_name' :  import_map[config_name][i].first
    //                         along with a list of the mapping properties {lattice_deformation, basis_deformation}  :  import_map[config_name][i].second
    std::map<std::string, std::vector<std::pair<std::string, std::vector<double> > > > import_map;
    std::vector<std::string > error_log;
    Index n_unique(0);
    // iterate over structure files
    std::cout << "  Beginning import of " << pos_paths.size() << " configuration" << (pos_paths.size() > 1 ? "s" : "") << "...\n" << std::endl;
    for(auto it = pos_paths.begin(); it != pos_paths.end(); ++it) {
      if(it != pos_paths.begin())
        std::cout << "\n***************************\n" << std::endl;

      fs::path pos_path, import_path;

      pos_path = fs::absolute(*it);
      std::string imported_name;
      BasicStructure<Site> import_struc;

      // If user requested data import, try to get structural data from properties.calc.json, instead of POS, etc.
      // Since properties.calc.json would be used during 'casm update' to validate relaxation
      if(vm.count("data") && pos_path.extension() != ".json" && pos_path.extension() != ".JSON") {
        fs::path dft_path = pos_path;
        dft_path.remove_filename();
        (dft_path /= ("calctype." + primclex.get_curr_calctype())) /= "properties.calc.json";
        if(!fs::exists(dft_path)) {
          dft_path = pos_path;
          dft_path.remove_filename();
          dft_path /= "properties.calc.json";
          if(!fs::exists(dft_path)) {
            dft_path = fs::path();
          }
        }
        if(!dft_path.empty())
          pos_path = dft_path;
      }


      //Import structure and make note of path
      bool new_import = false;
      jsonParser relax_data;
      try {

        if(pos_path.extension() == ".json" || pos_path.extension() == ".JSON") {
          from_json(simple_json(import_struc, "relaxed_"), jsonParser(pos_path));
        }
        else {
          fs::ifstream struc_stream(pos_path);
          import_struc.read(struc_stream);
        }

        Eigen::Matrix3d cart_op;
        std::vector<Index> best_assignment;
        if(configmapper.import_structure_occupation(import_struc, imported_name, relax_data, best_assignment, cart_op)) {
          std::cout << "  " << pos_path << "\nwas imported successfully as " << imported_name << std::endl << std::endl;
          n_unique++;
          new_import = true;
        }
        else {
          std::cout << "  " << pos_path << "\n  mapped onto pre-existing equivalent structure " << imported_name << std::endl << std::endl;
        }
        std::cout << "  Relaxation stats -> lattice_deformation = " << relax_data["lattice_deformation"].get<double>()
                  << "      basis_deformation = " << relax_data["basis_deformation"].get<double>() << std::endl << std::endl;;
      }
      catch(std::exception &e) {
        std::cerr << "  ERROR: Unable to import " << pos_path << " because \n"
                  << "    -> " << e.what() << "\n\n";
        error_log.push_back(it->string() + "\n     -> " + e.what());
        if(it != pos_paths.cend()) {
          std::cout << "  Continuing...\n";
        }
        continue;
      }

      Configuration &imported_config(primclex.configuration(imported_name));

      imported_config.push_back_source(json_pair("import_mapped", pos_path.string()));

      import_path = imported_config.get_path();
      import_path = fs::absolute(import_path);


      //Exit if user did not request not to copy data
      if(!vm.count("data"))
        continue;


      bool prevent_import(false);
      std::vector<double> conflict_stats;
      if(fs::exists(imported_config.calc_properties_path())) {
        prevent_import = true;
        if(import_map[imported_name].size() == 0) {

          if(imported_config.calc_properties().contains("basis_deformation") && imported_config.calc_properties().contains("lattice_deformation")) {
            conflict_stats.push_back(imported_config.calc_properties()["lattice_deformation"].get<double>());
            conflict_stats.push_back(imported_config.calc_properties()["basis_deformation"].get<double>());
            //if(imported_config.calc_properties().contains("relaxed_energy"))
            //conflict_stats.push_back(imported_config.calc_properties()["relaxed_energy"].get<double>());
          }
          import_map[imported_name].push_back(std::make_pair("A pre-existing copy", conflict_stats));
        }
        //std::cout << "Calculation data already exists for configuration " << imported_name << ". I will not overwrite it.\n Continuing..." << std::endl;
        //continue;
      }


      std::cout << "  Attempting to import data..." << std::endl;
      if(pos_path.extension() != ".json" && pos_path.extension() != ".JSON") {
        std::cout << "  No calculation data was found in the enclosing directory. Continuing..." << std::endl;
        continue;
      }

      if(relax_data.contains("basis_deformation") && relax_data.contains("lattice_deformation")) {
        conflict_stats.clear();
        conflict_stats.push_back(relax_data["lattice_deformation"].get<double>());
        conflict_stats.push_back(relax_data["basis_deformation"].get<double>());
        import_map[imported_name].push_back(std::make_pair(pos_path.string(), conflict_stats));
      }

      if(prevent_import)
        continue;

      fs::path import_target = import_path / ("calctype." + primclex.get_curr_calctype());
      if(!fs::exists(import_target))
        fs::create_directories(import_target);

      fs::copy_file(pos_path, imported_config.calc_properties_path());

      if(!fs::exists(imported_config.get_pos_path()))
        primclex.configuration(imported_name).write_pos();

      jsonParser calc_data;
      if(!imported_config.read_calc_properties(calc_data)) {
        std::cout << "  WARNING: Some properties from " << *it << " were not valid. Viable values will still be recorded.\n";
      }
      //append relax_data onto calc_data
      auto jit = relax_data.cbegin(), jit_end = relax_data.cend();
      for(; jit != jit_end; ++jit) {
        calc_data[jit.name()] = *jit;
      }
      imported_config.set_calc_properties(calc_data);

      imported_config.push_back_source(json_pair("data_inferred_from_mapping", pos_path.string()));

    }
    std::cout << "\n***************************\n" << std::endl;

    std::cout << "  Finished importing " << pos_paths.size() <<  " structures";
    if(n_unique == 0)
      std::cout << " (none of these are new or unique)";
    else if(n_unique < pos_paths.size())
      std::cout << " (only " << n_unique << " of these " << (n_unique == 1 ? "is" : "are") << " new and unique)";
    std::cout << "." <<  std::endl;

    //Update directories
    std::cout << "  Writing SCEL..." << std::endl;
    primclex.print_supercells();
    std::cout << "  Writing config_list..." << std::endl << std::endl;
    primclex.write_config_list();
    std::cout << "  DONE" << std::endl << std::endl;

    if(error_log.size() > 0) {
      std::cout << "  WARNING: --The following paths could not be imported due to errors:\n";
      for(auto it = error_log.cbegin(); it != error_log.cend(); ++it) {
        std::cout << *it
                  << "\n        ----------------------------------------------\n" << std::endl;
      }
      std::cout << "\n" << std::endl;
    }
    auto it = import_map.cbegin(), it_end = import_map.cend();
    bool no_conflict = true;
    for(; it != it_end; ++it) {
      if(it->second.size() > 1) {
        if(no_conflict) {
          no_conflict = false;
          std::cout << "  WARNING: --The following conflicts were found\n"
                    << "           --for each conflict, calculation data was committed for the first duplicate listed\n";
        }
        std::cout   << "        " << it->second.size() << " copies of configuration " << it->first << " were found:\n";
        for(auto it1 = it->second.cbegin(); it1 != it->second.cend(); ++it1) {
          std::cout << "            " << it1->first << std::endl
                    << "              -- lattice_deformation = ";
          if(it1->second.size())
            std::cout << (it1->second)[0];
          else
            std::cout << "unknown";
          std::cout << " and basis_deformation = ";
          if(it1->second.size())
            std::cout << (it1->second)[1] << std::endl;
          else
            std::cout << "unknown" << std::endl;
        }
        std::cout << "\n        ----------------------------------------------\n" << std::endl;
      }

    }
    if(!no_conflict)
      std::cout << "  Please review these conflicts.  A different resolution can be obtained by removing datafiles from\n"
                << "  the training_data directory and performing an import using a manually reduced set of files.\n";
    std::cout << "  DONE" << std::endl << std::endl;

    std::cout << std::endl;

    return 0;
  };

}

