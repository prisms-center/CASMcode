#include "update.hh"

#include <cstring>

#include "casm/crystallography/jsonStruc.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/CASM_classes.hh"
//#include "casm/misc/Time.hh"
#include "casm_functions.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'update' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int update_command(int argc, char *argv[]) {

    po::variables_map vm;
    int choice;
    std::string scellname;
    int refid, configid;
    double tol(TOL);
    double vol_tol(0.25);
    double lattice_weight(0.5);
    po::options_description desc("'casm update' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    ("cost-weight,w", po::value<double>(&lattice_weight)->default_value(0.5),
     "Adjusts cost function for mapping optimization (cost=w*lattice_deformation+(1-w)*basis_deformation)")
    ("max-vol-change", po::value<double>(&vol_tol)->default_value(0.25),
     "Adjusts range of SCEL volumes searched while mapping imported structure onto ideal crystal (only necessary if the presence of vacancies makes the volume ambiguous). Default is +/- 25% of relaxed_vol/prim_vol. Smaller values yield faster import, larger values may yield more accurate mapping.")
    ("force,f", "Force all configurations to update (otherwise, use timestamps to determine which configurations to update)");

    try {
      po::store(po::parse_command_line(argc, argv, desc), vm);

      /** --help option
       */
      if(vm.count("help")) {
        std::cout << std::endl;
        std::cout << desc << std::endl;
        std::cout << "DESCRIPTION" << std::endl;
        std::cout << "    Updates all values and files after manual changes or configuration \n";
        std::cout << "    calculations.\n";
        std::cout << "\n";

        return 0;
      }

      po::notify(vm);

    }
    catch(po::error &e) {
      std::cerr << desc << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return 1;
    }

    catch(std::exception &e) {
      std::cerr << desc << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return 1;
    }

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cerr << "Error in 'casm update': No casm project found." << std::endl;
      return 1;
    }
    fs::current_path(root);

    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;


    ConfigMapper configmapper(primclex, lattice_weight, vol_tol, ConfigMapper::rotate | ConfigMapper::robust, tol);
    std::cout << "Reading calculation data... " << std::endl << std::endl;
    std::vector<std::string> bad_config_report;
    std::vector<std::string> prop_names = primclex.get_curr_property();
    PrimClex::config_iterator it = primclex.config_begin();
    Index num_updated(0);
    for(; it != primclex.config_end(); ++it) {
      /// Read properties.calc.json file containing externally calculated properties
      ///   location: casmroot/supercells/SCEL_NAME/CONFIG_ID/CURR_CALCTYPE/properties.calc.json
      ///
      ///   Will read as many curr_property as found in properties.calc.json


      //std::cout << "begin Configuration::read_calculated()" << std::endl;

      /// properties.calc.json: contains calculated properties
      ///   Currently only loading those properties that have references
      fs::path filepath = it->calc_properties_path();
      // determine if there is fresh data to read and put it in 'calc_properties'
      jsonParser parsed_props;
      if(fs::exists(filepath)) {
        time_t datatime, filetime;
        // Compare 'datatime', from config_list database to 'filetime', from filesystem timestamp
        it->calc_properties().get_if(datatime, "data_timestamp");
        filetime = fs::last_write_time(filepath);

        if(!vm.count("force") && filetime == datatime) {
          continue;
        }
        num_updated++;
        std::cout << std::endl << "***************************" << std::endl << std::endl;
        std::cout << "Working on " << filepath.string() << "\n";

        //json relax_data;
        it->read_calc_properties(parsed_props);
        bool new_config_flag;
        std::string imported_name;

        {
          //Convert relaxed structure into a configuration, merge calculation data
          BasicStructure<Site> relaxed_struc;
          jsonParser json(filepath);
          from_json(simple_json(relaxed_struc, "relaxed_"), json);
          std::vector<Index> best_assignment;
          Eigen::Matrix3d cart_op;
          json.put_obj();
          try {
            new_config_flag = configmapper.import_structure_occupation(relaxed_struc,
                                                                       &(*it),
                                                                       imported_name,
                                                                       json,
                                                                       best_assignment,
                                                                       cart_op);
          }
          catch(std::exception &e) {
            std::cerr << "\nError: Unable to map relaxed structure data contained in " << filepath << " onto PRIM.\n"
                      << "       " << e.what() << std::endl;
            //throw std::runtime_error(std::string("Unable to map relaxed structure data contained in ") + filepath.string() + " onto PRIM.\n");
            return 1;
          }

          //copy data over
          for(auto jit = json.cbegin(); jit != json.cend(); ++jit) {
            parsed_props[jit.name()] = *jit;
          }
        }
        if(imported_name == it->name()) {
          it->set_calc_properties(parsed_props);
          continue;
        }

        Configuration &imported_config = primclex.configuration(imported_name);
        // Structure is mechanically unstable!
        std::cout << "Configuration " << it->name() << " appears to be mechanically unstable!\n"
                  << "After relaxation, it most closely maps onto " << " configuration " << imported_name << ", which"
                  << (new_config_flag ?
                      " has been automatically added to"
                      : " already exists in")
                  << " your project.\n";
        // Note the instability:
        it->push_back_source(json_unit("mechanically_unstable"));
        it->push_back_source(json_pair("relaxed_to", imported_name));
        imported_config.push_back_source(json_pair("relaxation_of", it->name()));
        bad_config_report.push_back(std::string("  - ") + it->name() + " relaxed to " + imported_name);

        // if imported_config has no properties, copy them over
        if(!fs::exists(imported_config.calc_properties_path())
           && (imported_config.calc_properties().is_null() || imported_config.calc_properties().size() == 0)) {
          std::cout << "Because no calculation data exists for configuration " << imported_name << ",\n"
                    << "it will assume the data from " << filepath << "\n";

          if(!fs::exists(imported_config.get_pos_path()))
            primclex.configuration(imported_name).write_pos();

          fs::path import_target = imported_config.calc_properties_path();
          import_target.remove_filename();
          if(!fs::exists(import_target))
            fs::create_directories(import_target);

          std::cout << "New file path is " << import_target << "\n";
          fs::copy_file(filepath, imported_config.calc_properties_path());


          parsed_props["data_timestamp"] = fs::last_write_time(imported_config.calc_properties_path());

          imported_config.set_calc_properties(parsed_props);
          imported_config.push_back_source(json_pair("data_inferred_from_mapping", it->name()));

          continue;
        }
        // prior data exist -- we won't copy data over, but do some validation to see if data are compatible
        else {
          const jsonParser &extant_props = imported_config.calc_properties();
          bool data_mismatch = false;

          auto prop_it = extant_props.cbegin(), prop_end = extant_props.cend();
          for(; prop_it != prop_end; ++prop_it) {
            if(!parsed_props.contains(prop_it.name())) {
              data_mismatch = true;
              continue;
            }

            // Try to ensure that the property is some sort of energy and convertible to scalar
            if(!prop_it->is_number() || (prop_it.name()).find("energy") == std::string::npos)
              continue;

            if(!almost_equal(prop_it->get<double>(), parsed_props[prop_it.name()].get<double>(), 1e-4)) {
              if(parsed_props[prop_it.name()].get<double>() < prop_it->get<double>()) {
                std::cout << "\nWARNING: Mapped configuration " << imported_name << " has \n"
                          << "                    " << prop_it.name() << "=" << prop_it->get<double>() << "\n"
                          << "         Which is higher than the relaxed value for mechanically unstable configuration " << it->name() << " which is\n"
                          << "                    " << prop_it.name() << "=" << parsed_props[prop_it.name()].get<double>() << "\n"
                          << "         This suggests that " << imported_name << " may be a metastable minimum.  Please investigate further.\n";
                bad_config_report.back() += ", **which may be metastable**";
                continue;
              }
              data_mismatch = true;
            }
          }
          if(data_mismatch)
            std::cout << "WARNING: The data parsed from \n"
                      << "             " << filepath << "\n"
                      << "         is incompatible with existing data for configuration " << imported_name << "\n"
                      << "         even though " << it->name() << " was found to relax to " << imported_name << "\n";
        }

        std::cout << std::endl;
      }
    }
    std::cout << std::endl << "***************************" << std::endl << std::endl;
    std::cout << "  DONE: ";
    if(num_updated == 0)
      std::cout <<  "No new data were detected." << std::endl << std::endl;
    else {
      std::cout <<  "Analyzed new data for " << num_updated << " configurations." << std::endl << std::endl;
      std::cout << "Generating references... " << std::endl << std::endl;
      /// This also re-writes the config_list.json
      primclex.generate_references();
      std::cout << "  DONE" << std::endl << std::endl;
      if(bad_config_report.size() > 0) {
        std::cout << "Some abonormal relaxations were detected:" << std::endl << std::endl
                  << "           *** Final Relaxation Report ***" << std::endl;
        for(auto report_it = bad_config_report.cbegin(); report_it != bad_config_report.cend(); ++report_it)
          std::cout << *report_it << std::endl;
        std::cout << "\nIt is recommended that you review these configurations more carefully.\n" << std::endl;
      }
    }
    return 0;
  }

}

