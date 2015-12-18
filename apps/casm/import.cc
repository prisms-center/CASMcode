#include "import.hh"

#include <cstring>
#include <tuple>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"
//#include "casm/clex/ConfigMapping.hh"
//#include "casm/casm_io/FileSystemInterface.hh"

namespace CASM {
  namespace Import_impl {
    // record some pieces of data from each import.  These are used to resolve conflicts at the end
    // the 'relaxjson' datarecord stores relaxation properties that will be merged into Configuration::calc_properties() during the final step
    enum DataType {path = 0, relaxjson = 1, contcar = 2, energy = 3};
    using Data = std::tuple<fs::path, jsonParser, std::string, std::pair<bool, double> >;
  }

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
    ("min-energy", "Resolve mapping conflicts based on energy rather than deformation.")
    ("max-vol-change", po::value<double>(&vol_tol)->default_value(0.25),
     "Adjusts range of SCEL volumes searched while mapping imported structure onto ideal crystal (only necessary if the presence of vacancies makes the volume ambiguous). Default is +/- 25% of relaxed_vol/prim_vol. Smaller values yield faster import, larger values may yield more accurate mapping.")
    ("batch,b", po::value<fs::path>(&batch_path), "Path to batch file, which should list one structure file path per line (can be used in combination with --pos)")
    ("rotate,r", "Rotate structure to be consistent with setting of PRIM")
    ("ideal,i", "Assume imported structures are unstrained (ideal) for faster importing. Can be slower if used on deformed structures, in which case more robust methods will be used")
    //("strict,s", "Request that symmetrically equivalent configurations be treated as distinct.")
    ("data,d", "Attempt to extract calculation data from the enclosing directory of the structure files, if it is available");

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
      return 9;
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
    std::map<Configuration *, std::vector<Import_impl::Data> > import_map;
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
      std::pair<bool, double> checkenergy(false, 0.0);
      double energy;
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
          std::cout << "  " << pos_path << "\n  was imported successfully as " << imported_name << std::endl << std::endl;
          n_unique++;
          new_import = true;
        }
        else {
          std::cout << "  " << pos_path << "\n  mapped onto pre-existing equivalent structure " << imported_name << std::endl << std::endl;
        }
        relax_data = fullrelax_data["best_mapping"];
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



      //Exit if user did not request not to copy data
      if(!vm.count("data"))
        continue;

      std::stringstream contcar_ss;
      import_struc.print5_occ(contcar_ss);
      import_map[&imported_config].push_back(Import_impl::Data(pos_path, relax_data, contcar_ss.str(), checkenergy));

    }

    // All the mapping is finished; now we migrate data, if requested
    std::stringstream conflict_log;
    if(vm.count("data")) {
      std::cout << "  Attempting to import data..." << std::endl;
      auto it(import_map.begin()), end_it(import_map.end());
      for(; it != end_it; ++it) {
        Configuration &imported_config = *(it->first);
        std::vector<Import_impl::Data> &data_vec(it->second);

        fs::path import_path = fs::absolute(imported_config.get_path());
        bool preexisting(false);
        if(fs::exists(imported_config.calc_properties_path()))
          preexisting = true;
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
          std::cout << "  No calculation data was found in the enclosing directory of \n"
                    << "    " << pos_path << std::endl
                    << "  Continuing..." << std::endl;
          continue;
        }

        fs::path import_target = import_path / ("calctype." + primclex.get_curr_calctype());
        if(!fs::exists(import_target))
          fs::create_directories(import_target);

        fs::copy_file(pos_path, imported_config.calc_properties_path());

        if(!fs::exists(imported_config.get_pos_path()))
          imported_config.write_pos();

        {
          fs::ofstream contcar_out(import_target / "relaxed_structure.vasp");
          contcar_out << std::get<Import_impl::contcar>(data_vec[best_ind]);
        }

        jsonParser calc_data;
        if(!imported_config.read_calc_properties(calc_data)) {
          std::cout << "  WARNING: Some properties from " << pos_path << " were not valid. Viable values will still be recorded.\n";
        }

        jsonParser &relaxjson = std::get<Import_impl::relaxjson>(data_vec[best_ind]);
        //append relaxjson onto calc_data
        auto jit = relaxjson.cbegin(), jit_end = relaxjson.cend();
        for(; jit != jit_end; ++jit) {
          calc_data[jit.name()] = *jit;
        }
        imported_config.set_calc_properties(calc_data);

        imported_config.push_back_source(json_pair("data_inferred_from_mapping", pos_path.string()));

      }
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
    if(conflict_log.str().size()) {
      std::cout << "  WARNING: -- The following conflicts were found\n" << std::endl
                << conflict_log.str() << std::endl;

      std::cout << "  Please review these conflicts.  A different resolution can be obtained by removing datafiles from\n"
                << "  the training_data directory and performing an import using a manually reduced set of files.\n";
    }
    std::cout << "  DONE" << std::endl << std::endl;

    std::cout << std::endl;

    return 0;
  };

}

