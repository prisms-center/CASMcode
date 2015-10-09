#include "select.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"

namespace CASM {

  bool write_selection(const ConfigSelection<true> &config_select, bool force, const fs::path &out_path, bool write_json, bool only_selected) {
    if(fs::exists(out_path) && !force) {
      std::cerr << "File " << out_path << " already exists. Use --force to force overwrite." << std::endl;
      return 1;
    }

    std::cout  << "Writing selection: " << out_path << std::endl;
    if(write_json || out_path.extension() == ".json" || out_path.extension() == ".JSON") {
      jsonParser json;
      config_select.to_json(json);
      SafeOfstream sout;
      sout.open(out_path);
      json.print(sout.ofstream(), only_selected);
      sout.close();
    }
    else {
      SafeOfstream sout;
      sout.open(out_path);
      config_select.print(sout.ofstream(), only_selected);
      sout.close();
    }
    std::cout << "  DONE." << std::endl;
    return 0;
  }

  // ///////////////////////////////////////
  // 'select' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int select_command(int argc, char *argv[]) {

    //casm enum [—supercell min max] [—config supercell ] [—hopconfigs hop.background]
    //- enumerate supercells and configs and hop local configurations

    std::vector<std::string> setcom;
    std::vector<std::string> selection;
    fs::path out_path;
    COORD_TYPE coordtype;
    po::variables_map vm;
    bool force;

    // casm select --set <RPN> -c ‘config_list_A’ (-i force in-place? -back for backup?)
    // - Change ‘is_selected’ status in place
    // - Use ‘master’, or just don’t provide -c to indicate master list
    // casm select --set <RPN> -c ‘config_list_A’ -o ‘config_list_C’
    // - Apply ‘--set’, then write all configurations (selected and unselected)
    // - Use ‘master’, or just don’t provide -c to indicate master list
    // casm select --subset -c ‘config_list_A’ -o ‘config_list_C’
    // - Write selected configurations in A into C
    // - Use ‘master’, or just don’t provide -c to indicate master list
    // casm select --union -c ‘config_list_A’ ‘config_list_B’ ... -o ‘config_list_C’
    // - Write configurations selected in any of A, B, etc. into C
    // - Use ‘master’ to indicate master list
    // casm select --intersection -c ‘config_list_A’, ‘config_list_B’ -o ‘config_list_C’
    // - Write configurations selected in all of A, B, ... into C
    // - Use ‘master’ to indicate master list


    try {

      /// Set command line options using boost program_options
      po::options_description desc("'casm select' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("config,c", po::value<std::vector<std::string> >(&selection)->multitoken(),
       "One or more configuration files to operate on. If not given, or if given the keyword \"MASTER\" the master list is used.")
      ("output,o", po::value<fs::path>(&out_path), "Name for output file")
      ("json", "Write JSON output (otherwise CSV, unless output extension in .json or .JSON)")
      ("subset", "Write selected configurations in input list")
      ("union", "Write configurations selected in any of the input lists")
      ("intersection", "Write configurations selected in all of the input lists")
      ("set", po::value<std::vector<std::string> >(&setcom)->multitoken(), "Criteria for selecting configurations.  Call using --set [OPT ...]")
      ("force,f", po::value(&force)->zero_tokens(), "Overwrite output file");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw



        bool call_help = false;

        if(!vm.count("help")) {
          if(vm.count("set") + vm.count("subset") + vm.count("union") + vm.count("intersection") != 1) {
            std::cout << desc << std::endl;
            std::cout << "Error in 'casm select'. Must use one of --set, --subset, --union, or --intersection." << std::endl;
            //call_help = true;
            return 1;
          }

          if(vm.count("subset") && !vm.count("output")) {
            std::cout << "Error in 'casm select --subset'. Please specify an --output file." << std::endl;
            //call_help = true;
            return 1;
          }

          if(vm.count("union") && !vm.count("output")) {
            std::cout << "Error in 'casm select --union'. Please specify an --output file." << std::endl;
            //call_help = true;
            return 1;
          }

          if(vm.count("intersection") && !vm.count("output")) {
            std::cout << "Error in 'casm select --intersection'. Please specify an --output file." << std::endl;
            //call_help = true;
            return 1;
          }
        }

        /** --help option
        */
        if(vm.count("help") || call_help) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION" << std::endl;
          std::cout << "\n";
          std::cout << "    Use '--set [on | off] [options]' for setting 'selected'.\n";
          std::cout << "    Options are listed using reverse polish notation.\n";
          std::cout << "    Operators include: re (regex match whole string),\n";
          std::cout << "    rs (regex match anywhere in string), eq (==), ne (!=),\n";
          std::cout << "    lt (<), le (<=), gt (>), ge (>=), add, sub, mult, div, pow, \n";
          std::cout << "    AND, OR, NOT, XOR. Variable keywords are:                   \n" <<
                    "      'scelname': ex. SCELV_A_B_C_D_E_F                         \n" <<
                    "      'configname': ex. SCELV_A_B_C_D_E_F/I                     \n" <<
                    "      'scel_size': supercell volume, as number of primitive cells  \n" <<
                    "      'is_calculated': \"1\" if properties have been calculated,\n" <<
                    "         \"0\" otherwise                                        \n" <<
                    "      'is_groundstate': \"1\" if ground state, \"0\" if not,    \n" <<
                    "         \"unknown\" otherwise                                  \n" <<
                    "      'dist_from_hull': distance of formation energy from convex\n" <<
                    "         hull, or \"unknown\"                                   \n" <<
                    "      'formation_energy': calculated formation energy, or \"unknown\"\n" <<
                    "      'clex(formation_energy)': cluster expansion predicted     \n" <<
                    "         formation energy                                       \n" <<
                    "      'comp(x)': composition, where x is one of 'a', 'b', ...   \n" <<
                    "      'atom_frac(X)': atomic fraction, where X is the name of an\n" <<
                    "         atom in the prim. Vacancies are not included.          \n" <<
                    "      'site_frac(X)': fraction of all sites occupied by X, where\n" <<
                    "         X is the name of an atom in the prim. Vacancies are    \n" <<
                    "         included.                                              \n";
          std::cout << "\n";
          std::cout << "    Examples:\n";
          std::cout << "      casm select --set on scelname SCEL3 rs\n";
          std::cout << "      - set 'selected'=true for all volume 3 configurations\n";
          std::cout << "\n";
          std::cout << "      casm select --set on scelname SCEL3 rs true_comp(A) 0.5 lt AND\n";
          std::cout << "      - set 'selected'=true for all volume 3 configurations\n";
          std::cout << "        with A true_composition less than 0.5.\n";
          std::cout << "\n";
          std::cout << "      casm select --set on comp(a) comp(b) add 0.75 lt\n";
          std::cout << "      - set 'selected'=true for all configurations with\n";
          std::cout << "        'a' + 'b' composition less than 0.75.\n";
          std::cout << "\n";
          std::cout << "    Use '--subset', '--union', and '--intersection' to operate on\n" <<
                    "    configuration lists.\n\n";

          std::cout << "    Note: It is currently necessary to use ' -X', with quotes and \n" <<
                    "    a leading space to parse a negative number.\n\n";

          if(call_help)
            return 1;

          return 0;
        }

        po::notify(vm); // throws on error, so do after help in case
        // there are any problems

        if(vm.count("set") && vm.count("config") && selection.size() != 1) {
          std::cout << "Error in 'casm select --set'. Zero (implies using Master list) or one lists should be input." << std::endl;
          return 1;
        }

        if(vm.count("subset") && vm.count("config") && selection.size() != 1) {
          std::cout << "Error in 'casm select --subset'. Zero (implies using Master list) or one lists should be input." << std::endl;
          return 1;
        }

      }
      catch(po::error &e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << desc << std::endl;
        return 1;
      }
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return 1;

    }

    for(int i = 0; i < selection.size(); i++) {
      if(selection[i] != "MASTER") {
        selection[i] = fs::absolute(fs::path(selection[i])).string();
      }
    }
    out_path = fs::absolute(out_path);


    // switch to root directory
    fs::path orig = fs::current_path();
    fs::path root = find_casmroot(orig);
    if(root.empty()) {
      std::cout << "Error: No casm project found." << std::endl;
      return 1;
    }
    fs::current_path(root);


    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    if(vm.count("set")) {

      bool only_selected = false;

      Array<std::string> criteria;
      for(int i = 0; i < setcom.size(); i++)
        criteria.push_back(setcom[i]);

      std::cout << "Set selection: " << criteria << std::endl << std::endl;

      /// Prepare for calculating correlations. Maybe this should get put into Clexulator.
      const DirectoryStructure &dir = primclex.dir();
      const ProjectSettings &set = primclex.settings();
      if(fs::exists(dir.clexulator_src(set.name(), set.bset()))) {
        primclex.read_global_orbitree(dir.clust(set.bset()));
        primclex.generate_full_nlist();
        primclex.generate_supercell_nlists();

      }

      if(!vm.count("config") || (selection.size() == 1 && selection[0] == "MASTER")) {

        if(!vm.count("output")) {
          for(auto it = primclex.config_begin(); it != primclex.config_end(); ++it) {
            it->set_selected(get_selection(criteria, *it, it->selected()));
          }

          std::cout << "  DONE." << std::endl << std::endl;

          std::cout << "Writing config_list..." << std::endl;
          primclex.write_config_list();
          std::cout << "  DONE." << std::endl;

        }
        else {
          ConfigSelection<true> config_select(primclex);
          for(auto it = config_select.config_begin(); it != config_select.config_end(); ++it) {
            it.set_selected(get_selection(criteria, *it, it.selected()));
          }

          std::cout << "  DONE." << std::endl << std::endl;

          return write_selection(config_select, vm.count("force"), out_path, vm.count("json"), only_selected);

        }
      }
      else {
        ConfigSelection<true> config_select(primclex, selection[0]);
        for(auto it = config_select.config_begin(); it != config_select.config_end(); ++it) {
          it.set_selected(get_selection(criteria, *it, it.selected()));
        }

        bool force = vm.count("force");
        if(!vm.count("output")) {
          out_path = selection[0];
          force = true;
        }

        return write_selection(config_select, force, out_path, vm.count("json"), only_selected);
      }
    }

    if(vm.count("subset")) {

      bool only_selected = true;

      if(!vm.count("config") || (selection.size() == 1 && selection[0] == "MASTER")) {
        ConfigSelection<true> config_select(primclex);
        return write_selection(config_select, force, out_path, vm.count("json"), only_selected);
      }
      else {
        ConfigSelection<true> config_select(primclex, selection[0]);
        return write_selection(config_select, force, out_path, vm.count("json"), only_selected);
      }
    }

    if(vm.count("union")) {

      if(selection.size() == 0) {
        selection.push_back("MASTER");
      }

      // load initial selection into config_select
      ConfigSelection<true> config_select;
      if(selection[0] == "MASTER") {
        config_select = ConfigSelection<true>(primclex);
      }
      else {
        config_select = ConfigSelection<true>(primclex, selection[0]);
      }

      // remove unselected configurations
      auto it = config_select.config_cbegin();
      while(it != config_select.config_cend()) {
        if(!it.selected()) {
          it = config_select.erase(it);
        }
        else {
          ++it;
        }
      }

      // loop through other lists, inserting all selected configurations
      for(int i = 1; i < selection.size(); i++) {
        ConfigSelection<true> tselect;
        if(selection[i] == "MASTER") {
          tselect = ConfigSelection<true>(primclex);
        }
        else {
          tselect = ConfigSelection<true>(primclex, selection[i]);
        }
        for(auto it = tselect.config_cbegin(); it != tselect.config_cend(); ++it) {
          if(it.selected()) {
            config_select.insert(std::make_pair(it.name(), it.selected()));
          }
        }
      }

      bool only_selected = true;
      return write_selection(config_select, force, out_path, vm.count("json"), only_selected);

    }

    if(vm.count("intersection")) {

      if(selection.size() == 0) {
        selection.push_back("MASTER");
      }

      // load initial selection into config_select
      ConfigSelection<true> config_select;
      if(selection[0] == "MASTER") {
        config_select = ConfigSelection<true>(primclex);
      }
      else {
        config_select = ConfigSelection<true>(primclex, selection[0]);
      }

      // remove unselected configurations
      auto it = config_select.config_cbegin();
      while(it != config_select.config_cend()) {
        if(!it.selected()) {
          it = config_select.erase(it);
        }
        else {
          ++it;
        }
      }

      // loop through other lists, keeping only configurations selected in the other lists
      for(int i = 1; i < selection.size(); i++) {
        ConfigSelection<true> tselect;
        if(selection[i] == "MASTER") {
          tselect = ConfigSelection<true>(primclex);
        }
        else {
          tselect = ConfigSelection<true>(primclex, selection[i]);
        }
        auto it = config_select.config_cbegin();
        while(it != config_select.config_cend()) {
          auto find_it = tselect.find(it.name());
          if(find_it == tselect.config_end()) {
            it = config_select.erase(it);
          }
          else if(!find_it.selected()) {
            it = config_select.erase(it);
          }
          else {
            ++it;
          }
        }
      }

      bool only_selected = true;
      return write_selection(config_select, force, out_path, vm.count("json"), only_selected);

    }

    std::cout << "\n***************************\n" << std::endl;

    std::cout << std::endl;

    return 0;
  };

}


