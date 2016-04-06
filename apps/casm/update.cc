#include "update.hh"

#include <cstring>

#include "casm/crystallography/jsonStruc.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/CASM_classes.hh"
//#include "casm/misc/Time.hh"
#include "casm_functions.hh"

namespace CASM {
  namespace Update_impl {
    // record some pieces of data from each import.  These are used to resolve conflicts at the end
    // the 'relaxjson' datarecord stores relaxation properties that will be merged into Configuration::calc_properties() during the final step
    enum DataType {relaxjson = 0, self_map = 1, new_config = 2};
    using Data = std::tuple<jsonParser, bool, bool>;
  }


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

    std::map<Configuration *, std::map<Configuration *, Update_impl::Data> > update_map;

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
      if(!fs::exists(filepath)) {
        if(it->calc_properties().size() > 0 && !(it->calc_properties().is_null())) {
          //clear the calculated properties if the data file is missing -- this probably means that the user has deleted it
          it->set_calc_properties(jsonParser());
        }
      }
      else {
        time_t datatime(0), filetime(fs::last_write_time(filepath));
        // Compare 'datatime', from config_list database to 'filetime', from filesystem timestamp
        it->calc_properties().get_if(datatime, "data_timestamp");

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
          //NOTE: reusing jsonParser to collect relaxation data
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
          if(imported_name != it->name() && json.contains("suggested_mapping")) {
            auto j_it(json["suggested_mapping"].cbegin()), j_end(json["suggested_mapping"].cend());
            for(; j_it != j_end; ++j_it) {
              parsed_props[j_it.name()] = *j_it;
            }
            update_map[&(*it)][&(*it)] = Update_impl::Data(parsed_props, false, new_config_flag);
          }

          Configuration &imported_config = primclex.configuration(imported_name);
          auto j_it(json["best_mapping"].cbegin()), j_end(json["best_mapping"].cend());
          for(; j_it != j_end; ++j_it) {
            parsed_props[j_it.name()] = *j_it;
          }

          update_map[&imported_config][&(*it)] = Update_impl::Data(parsed_props, it->name() == imported_name, new_config_flag);
        }


      }

    }


    // Now loop over update_map to construct a relaxation report and record calculation data appropriately
    double w = lattice_weight;
    std::stringstream relax_log;
    for(auto it = update_map.begin(); it != update_map.end(); ++it) {
      Configuration &imported_config = *(it->first);
      std::map<Configuration *, Update_impl::Data> &datamap(it->second);
      bool self_mapped(false);

      double lowest_energy(1e20), self_energy(1e20);
      std::map<Configuration *, Update_impl::Data>::iterator best_it;
      double best_cost(1e20);

      //enum DataType {name = 0, relaxjson = 1, self_map = 2};
      std::vector<std::string> report;
      for(auto it2 = datamap.begin(); it2 != datamap.end(); ++it2) {
        std::stringstream t_ss;
        //if(imported_config.name() == it->name()) {
        //it->set_calc_properties(parsed_props);
        //continue;
        //}
        Configuration &source_config(*(it2->first));
        const jsonParser &source_data(std::get<Update_impl::relaxjson>(it2->second));
        if(&source_config == &imported_config) {
          self_mapped = std::get<Update_impl::self_map>(it2->second);
        }
        else {
          // Note the instability:
          source_config.push_back_source(jsonParser("mechanically_unstable"));
          source_config.push_back_source(jsonParser(std::make_pair("relaxed_to", imported_config.name())));
          imported_config.push_back_source(jsonParser(std::make_pair("relaxation_of", source_config.name())));
        }
        double bd = source_data["basis_deformation"].get<double>();
        double ld = source_data["lattice_deformation"].get<double>();
        double tcost = w * ld + (1.0 - w) * bd;

        if(tcost < best_cost) {
          best_cost = tcost;
          best_it = it2;
        }

        t_ss << "Starting configuration " << source_config.name() << ":" << std::endl;
        t_ss << "                  Result of mapping relaxation onto " << imported_config.name() << std::endl;
        t_ss << "                    -- lattice_deformation = " << ld << ";  basis_deformation = " << bd << ";  weighted_avg = " << tcost << std::endl;

        if(&source_config != &imported_config) {
          jsonParser const &alt_data(std::get<Update_impl::relaxjson>(update_map[&source_config][&source_config]));
          bd = alt_data["basis_deformation"].get<double>();
          ld = alt_data["lattice_deformation"].get<double>();
          tcost = w * ld + (1.0 - w) * bd;
          t_ss << "                  Result of mapping relaxation onto " <<  source_config.name() << std::endl;
          t_ss << "                    -- lattice_deformation = " << ld << ";  basis_deformation = " << bd << ";  weighted_avg = " << tcost << std::endl;
        }
        if(source_data.contains("relaxed_energy")) {
          double energy = source_data["relaxed_energy"].get<double>();
          if(it2->first == it->first)
            self_energy = energy;
          if(energy < lowest_energy)
            lowest_energy = energy;
          t_ss << "                    -- relaxed_energy = " << source_data["relaxed_energy"].get<double>() << std::endl;
        }
        else
          t_ss << "                    -- relaxed_energy = unknown" << std::endl;
        report.push_back(t_ss.str());
      }
      Index best_ind(std::distance(datamap.begin(), best_it));
      if(datamap.size() > 1 || (!self_mapped && datamap.find(&imported_config) == datamap.end())) {
        relax_log << " " << datamap.size() << " structure" << (datamap.size() > 1 ? "s have" : " has") <<  " mapped onto Configuration " << imported_config.name()
                  << ", which " << (std::get<Update_impl::new_config>(datamap.cbegin()->second)
                                    ? "has been automatically added to your project"
                                    : "already existed in your project")
                  << std::endl << std::endl;

        for(Index i = 0; i < report.size(); i++) {
          if(i == best_ind) {
            relax_log << "    best ==> " << i + 1 << ") " << report[i] << std::endl;
          }
          else {
            relax_log << "              " << i + 1 << ") " << report[i] << std::endl;
          }
        }

        // if imported_config has no properties, copy properties from best config
        fs::path filepath((best_it->first)->calc_properties_path());
        jsonParser &parsed_props(std::get<Update_impl::relaxjson>(best_it->second));

        if(!fs::exists(imported_config.calc_properties_path())
           && (imported_config.calc_properties().is_null() || imported_config.calc_properties().size() == 0)) {
          relax_log << "  Because no calculation data exists for configuration " << imported_config.name() << ",\n"
                    << "  it will acquire the data from " << filepath << "\n";

          if(!fs::exists(imported_config.get_pos_path()))
            primclex.configuration(imported_config.name()).write_pos();

          fs::path import_target = imported_config.calc_properties_path();
          import_target.remove_filename();
          if(!fs::exists(import_target))
            fs::create_directories(import_target);

          //relax_log << "New file path is " << import_target << "\n";
          //fs::copy_file(filepath, imported_config.calc_properties_path());

          //parsed_props["data_timestamp"] = fs::last_write_time(imported_config.calc_properties_path());

          imported_config.set_calc_properties(parsed_props);
          imported_config.push_back_source(jsonParser(std::make_pair("data_inferred_from_mapping", (best_it->first)->name())));
        }

        // prior data exist -- we won't copy data over, but do some validation to see if data are compatible
        else {
          relax_log << "  There is already data associated with " << imported_config.name() << " -- it will not be overwritten.\n\n";
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
                relax_log << "\n     WARNING: Mapped configuration " << imported_config.name() << " has \n"
                          << "                    " << prop_it.name() << "=" << prop_it->get<double>() << "\n"
                          << "         Which is higher than the relaxed value for mechanically unstable configuration " << (best_it->first)->name() << " which is\n"
                          << "                    " << prop_it.name() << "=" << parsed_props[prop_it.name()].get<double>() << "\n"
                          << "         This suggests that " << imported_config.name() << " may be a metastable minimum.  Please investigate further.\n";
                continue;
              }
              data_mismatch = true;
            }
          }
          if(data_mismatch)
            relax_log << "       WARNING: There are calculated properties from \n"
                      << "             " << filepath << "\n"
                      << "         that do not match existing calculations for configuration " << imported_config.name() << " (they differ by more than 1E-04)\n"
                      << "         even though " << (best_it->first)->name() << " was found to relax to " << imported_config.name() << "\n";
        }
        relax_log << "\n          ----------------------------------------------\n" << std::endl;
        if(self_mapped && datamap.find(&imported_config) != datamap.end())
          imported_config.set_calc_properties(std::get<Update_impl::relaxjson>(datamap[&imported_config]));
      }
      else if(self_mapped) { //don't report anything, just record the data
        imported_config.set_calc_properties(std::get<Update_impl::relaxjson>(best_it->second));
      }
      else {
        //relax_log << "Starting configuration " << source_config.name() << " is mechanically unstable.\n";
        //relax_log << "\n          ----------------------------------------------\n" << std::endl;
      }
    }

    std::cout << std::endl << "***************************" << std::endl << std::endl;
    std::cout << "  DONE: ";
    if(num_updated == 0)
      std::cout <<  "No new data were detected." << std::endl << std::endl;
    else {
      std::cout << "Analyzed new data for " << num_updated << " configurations." << std::endl << std::endl;

      if(relax_log.str().size() > 0) {
        std::cout << "WARNING: Abnormal relaxations were detected:\n" << std::endl
                  << "           *** Final Relaxation Report ***" << std::endl
                  << relax_log.str()
                  << std::endl << std::endl;
        std::cout << "\nIt is recommended that you review these configurations more carefully.\n" << std::endl;
      }

    }

    if(num_updated > 0) {
      std::cout << "Writing to configuration database..." << std::endl << std::endl;
      primclex.write_config_list();
    }

    return 0;
  }
}

