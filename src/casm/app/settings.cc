#include "casm/app/casm_functions.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  class ClexDescription;
  bool clex_exists(const DirectoryStructure &dir, const ClexDescription &desc);

  namespace {
    std::string _wdefaultval(std::string name, std::pair<std::string, std::string> val) {
      return name + ": '" + val.first + "' (" + val.second + ")\n";
    }

    std::string _wdefaultval(std::string name, std::pair<fs::path, std::string> val) {
      return _wdefaultval(name, std::make_pair(val.first.string(), val.second));
    }
  }

  // ///////////////////////////////////////
  // 'settings' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int settings_command(const CommandArgs &args) {

    std::string single_input;
    std::vector<std::string> multi_input;
    COORD_TYPE coordtype;
    po::variables_map vm;

    try {

      /// Set command line options using boost program_options
      po::options_description desc("'casm settings' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("list,l", "List project settings")
      ("new-property", po::value<std::string>(&single_input), "Create cluster expansions for a new property")
      ("new-bset", po::value<std::string>(&single_input), "Create a new basis set")
      ("new-calctype", po::value<std::vector<std::string> >(&multi_input)->multitoken(), "Create a new calculation type")
      ("new-ref", po::value<std::string>(&single_input), "Create a new calculation reference")
      ("new-eci", po::value<std::string>(&single_input), "Create a new set of effective cluster interactions (ECI)")
      ("set-formation-energy", po::value<std::string>(&single_input), "Specify the cluster expansion to use for formation energy")
      ("new-clex", po::value<std::string>(&single_input), "Create a new cluster expansion")
      ("set-default-clex", po::value<std::string>(&single_input), "Set the cluster expansion that CASM uses or acts on by default")
      ("erase-clex", po::value<std::string>(&single_input), "Erase the specified cluster expansion. Does not erase underlying bset, eci, etc.")
      ("clex", po::value<std::string>(), "The cluster expansion for which to set property, bset, calctype, ref, or eci")
      ("set-property", po::value<std::string>(&single_input), "Set the current basis set")
      ("set-bset", po::value<std::string>(&single_input), "Set the basis set")
      ("set-calctype", po::value<std::vector<std::string> >(&multi_input)->multitoken(), "Set the calculation type")
      ("set-ref", po::value<std::vector<std::string> >(&multi_input)->multitoken(), "Set the calculation reference")
      ("set-eci", po::value<std::string>(&single_input), "Set the effective cluster interactions (ECI)")
      ("set-view-command", po::value<std::string>(&single_input), "Set the command used by 'casm view'.")
      ("set-cxx", po::value<std::string>(&single_input), "Set the c++ compiler. Use '' to revert to default.")
      ("set-cxxflags", po::value<std::string>(&single_input), "Set the c++ compiler options. Use '' to revert to default.")
      ("set-soflags", po::value<std::string>(&single_input), "Set the shared library compilation options. Use '' to revert to default.")
      ("set-casm-prefix", po::value<std::string>(&single_input), "Set the casm prefix. Use '' to revert to default.")
      ("set-boost-prefix", po::value<std::string>(&single_input), "Set the boost prefix. Use '' to revert to default.");

      try {
        po::store(po::parse_command_line(args.argc, args.argv, desc), vm); // can throw

        bool call_help = false;

        std::vector<std::string> all_opt = {"list",
                                            "new-property", "new-bset", "new-calctype", "new-ref", "new-eci",
                                            "new-clex", "set-formation-energy", "erase-clex",
                                            "set-default-clex", "set-property", "set-bset", "set-calctype", "set-ref", "set-eci",
                                            "set-cxx", "set-cxxflags", "set-soflags", "set-casm-prefix", "set-boost-prefix",
                                            "set-view-command"
                                           };
        int option_count = 0;
        for(int i = 0; i < all_opt.size(); i++) {
          option_count += vm.count(all_opt[i]);
        }

        // must call one and only one option at a time:
        if(option_count == 0) {
          std::cout << "Error in 'casm settings'. No option selected." << std::endl;
          call_help = true;
        }
        else if(option_count > 1) {
          std::cout << "Error in 'casm settings'. Use one option (other than --clex) at a time." << std::endl;
          call_help = true;
        }

        // --help option
        if(vm.count("help") || call_help) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION" << std::endl;
          std::cout << "\n";
          std::cout << "    Often it is useful to try multiple different basis sets, \n" <<
                    "    calculation settings, references, or ECI fits in a single\n" <<
                    "    project. The 'casm settings' option helps to organize    \n" <<
                    "    these within a project and quickly switch between        \n" <<
                    "    different settings.                                      \n";
          std::cout << "\n";
          std::cout << "    Examples:\n";
          std::cout << "                                                             \n" <<
                    "      casm settings --list                                   \n" <<
                    "      - List all settings, with '*' for current settings     \n\n" <<

                    "      casm settings --new-clex 'my_newclex'                  \n" <<
                    "      casm settings --set-default-clex 'other_clex'          \n" <<
                    "      - Creates a new group of settings for a cluster        \n" <<
                    "        expansion                                            \n" <<
                    "      - Includes property, calctype, ref, bset, and eci      \n" <<
                    "      - Can be used in Monte Carlo input files, and as       \n" <<
                    "        arguments to 'casm select' and 'casm query'          \n" <<
                    "        properties such as 'clex' and 'corr' to specify which\n" <<
                    "        basis functions and eci to use.                      \n\n" <<

                    "      casm settings --set-formation-energy 'other_clex'      \n" <<
                    "      - The cluster expansion 'other_clex' is copied to one  \n" <<
                    "        named 'formation_energy' which is used as the default\n" <<
                    "        for grand canoncial Monte Carlo calculations and     \n" <<
                    "        cluster expansion convex hull calculations.          \n\n" <<

                    "      casm settings --new-property 'my_new_property'         \n" <<
                    "      casm settings --new-bset 'my_new_bset'                 \n" <<
                    "      casm settings --new-calctype 'my_new_calctype' ['my_new_ref']\n" <<
                    "      casm settings --new-ref 'my_new_ref'                   \n" <<
                    "      casm settings --new-eci 'my_new_eci'                   \n" <<
                    "      - Creates new settings directories with appropriate    \n" <<
                    "        names                                                \n" <<
                    "      - For --new-property, a new 'default' eci is created.  \n" <<
                    "      - For --new-calctype, a new reference may optionally be\n" <<
                    "        speficied. If it is not, a new 'default' reference is\n" <<
                    "        created.                                             \n" <<
                    "      - For --new-ref, a new reference is created for the    \n" <<
                    "        current calctype                                     \n" <<
                    "      - For --new-eci, a new eci directory is created for the\n" <<
                    "         current clex, calctype and ref.                     \n\n" <<


                    "      casm settings --set-property 'other_property'          \n" <<
                    "      casm settings --set-bset 'other_bset'                  \n" <<
                    "      casm settings --set-calctype 'other_calctype' ['other_ref'] \n" <<
                    "      casm settings --set-ref ['other_calctype'] 'other_ref' \n" <<
                    "      casm settings --set-eci 'other_eci'                    \n" <<
                    "      - Switch the current settings                          \n" <<
                    "      - For --set-property, standard options are:            \n" <<
                    "        - 'formation_energy' (the only standard option for now)\n" <<
                    "        - 'formation_energy' (the only standard option for now)\n" <<
                    "        given, 'ref_default' is used, or if that doesn't     \n" <<
                    "        exist the first one found is used.                   \n" <<
                    "      - For --set-calctype 'other_ref' is optional, if not   \n" <<
                    "        given, 'ref_default' is used, or if that doesn't     \n" <<
                    "        exist the first one found is used.                   \n" <<
                    "      - For --set-ref, 'other_calctype' is optional if among \n" <<
                    "        all calctype there is no other reference called      \n" <<
                    "        'other_ref'. Otherwise it is required.               \n" <<
                    "      - For --set-ref, 'other_calctype' is optional if among \n" <<
                    "        all calctype there is no other reference called      \n" <<
                    "        'other_ref'. Otherwise it is required.               \n\n" <<

                    "      casm settings --set-view-command 'casm.view \"open -a /Applications/VESTA/VESTA.app\"'\n" <<
                    "      - Sets the command used by 'casm view' to open         \n" <<
                    "        visualization software.                              \n" <<
                    "      - Will be executed with '/path/to/POSCAR' as an        \n" <<
                    "        argument, the location of a POSCAR for a configuration\n" <<
                    "        selected for visualization.                          \n" <<
                    "\n";

          if(call_help)
            return 1;

          return 0;
        }

        po::notify(vm); // throws on error, so do after help in case
        // there are any problems



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

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    DirectoryStructure dir(root);
    ProjectSettings set(root);
    ClexDescription clex_desc = set.default_clex();

    if(vm.count("clex")) {
      auto it = set.cluster_expansions().find(vm["clex"].as<std::string>());
      if(it == set.cluster_expansions().end()) {
        args.err_log.error("Invalid --clex value");
        args.err_log << vm["clex"].as<std::string>() << " not found. expected basis set specifications file at: " << dir.bspecs(clex_desc.bset) << "\n" << std::endl;
        return ERR_INVALID_ARG;
      }
      clex_desc = it->second;
    }

    // disply all settings
    if(vm.count("list")) {
      set.print_summary(args.log);
      return 0;
    }

    // create new cluster_expansions/clex.property directory and sub-directories
    else if(vm.count("new-property")) {
      if(set.new_clex_dir(single_input)) {
        std::cout << "Created new property '" << single_input << "'.\n\n";
        clex_desc.property = single_input;

        // and create matching eci
        if(set.new_eci_dir(clex_desc.property, clex_desc.calctype, clex_desc.ref, clex_desc.bset, "default")) {
          clex_desc.eci = "default";
          std::cout << "Created new eci 'default'.\n\n";
        }
        else {
          std::cout << "Could not create new eci 'default'.\n\n";
          return 1;
        }
      }
      else {
        std::cout << "Could not create new property '" << single_input << "'.\n\n";
        return 1;
      }

      clex_desc.name = single_input;
      if(set.new_clex(clex_desc)) {
        std::cout << "Created new cluster expansion named '" << single_input << "'\n\n";
      }
      else {
        std::cout << "Could not create new cluster expansion named '" << single_input << "'.\n\n";
        return 1;
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Set '" << clex_desc.name << "' as default cluster expansion.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not set '" << clex_desc.name << "' as default cluster expansion.\n\n";
        return 1;
      }
    }

    // create new bset, and default eci for current calctype and ref
    else if(vm.count("new-bset")) {
      if(set.new_bset_dir(single_input)) {
        clex_desc.bset = single_input;
        std::cout << "Created new bset '" << single_input << "'.\n\n";

        // and create matching eci
        if(set.new_eci_dir(clex_desc.property, clex_desc.calctype, clex_desc.ref, clex_desc.bset, "default")) {
          clex_desc.eci = "default";
          std::cout << "Created new eci 'default'.\n\n";
        }
        else {
          std::cout << "Could not create new eci 'default'.\n\n";
          return 1;
        }
      }
      else {
        std::cout << "Could not create new bset '" << single_input << "'.\n\n";
        return 1;
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, false, false, true);
        }
        std::cout << "Switched to new bset '" << clex_desc.bset << "' and 'default' eci.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch new bset '" << clex_desc.bset << "' and 'default' eci.\n\n";
        return 1;
      }
    }

    // create new calctype, and ref (either as specified, or default)
    else if(vm.count("new-calctype")) {

      // create new calctype
      if(set.new_calc_settings_dir(multi_input[0])) {
        clex_desc.calctype = multi_input[0];
        std::cout << "Created new calctype '" << multi_input[0] << "'.\n\n";

        // if ref given, create specified ref
        if(multi_input.size() > 1) {
          if(set.new_ref_dir(multi_input[0], multi_input[1])) {
            clex_desc.ref = multi_input[1];
            std::cout << "Created new ref '" << multi_input[1] << "'.\n\n";

            // and create matching eci
            if(set.new_eci_dir(clex_desc.property, clex_desc.calctype, clex_desc.ref, clex_desc.bset, "default")) {
              std::cout << "Created new eci 'default'.\n\n";
            }
            else {
              std::cout << "Could not create new eci 'default'.\n\n";
              return 1;
            }
          }
          else {
            std::cout << "Could not create new ref '" << multi_input[0] << "'.\n\n";
            return 1;
          }
        }
        // else create default ref
        else {
          if(set.new_ref_dir(multi_input[0], "default")) {
            clex_desc.ref = "default";
            std::cout << "Created new ref 'default'.\n\n";

            // and create matching eci
            if(set.new_eci_dir(clex_desc.property, clex_desc.calctype, clex_desc.ref, clex_desc.bset, "default")) {
              clex_desc.eci = "default";
              std::cout << "Created new eci 'default'.\n\n";
            }
            else {
              std::cout << "Could not create new eci 'default'.\n\n";
              return 1;
            }
          }
          else {
            std::cout << "Could not create new ref 'default'.\n\n";
            return 1;
          }
        }

        return 0;
      }
      else {
        std::cout << "Could not create new calctype '" << multi_input[0] << "'.\n\n";
        return 1;
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Switched to new calctype '" << clex_desc.calctype << "', ref '" << clex_desc.ref << "', and '" << clex_desc.eci << "' eci.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to new calctype '" << clex_desc.calctype << "', ref '" << clex_desc.ref << "', and '" << clex_desc.eci << "' eci.\n\n";
        return 1;
      }
    }

    // create new ref for current calctype
    else if(vm.count("new-ref")) {
      if(set.new_ref_dir(clex_desc.calctype, single_input)) {
        clex_desc.ref = single_input;
        std::cout << "Created new ref '" << single_input << "'.\n\n";

        // and create matching eci
        if(set.new_eci_dir(clex_desc.property, clex_desc.calctype, clex_desc.ref, clex_desc.bset, "default")) {
          clex_desc.eci = "default";
          std::cout << "Created new eci 'default'.\n\n";
        }
        else {
          std::cout << "Could not create new eci 'default'.\n\n";
          return 1;
        }
      }
      else {
        std::cout << "Could not create new ref '" << single_input << "'.\n\n";
        return 1;
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Switched to new ref '" << clex_desc.ref << "', and '" << clex_desc.eci << "' eci.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to new ref '" << clex_desc.ref << "', and '" << clex_desc.eci << "' eci.\n\n";
        return 1;
      }
    }

    // create new eci directory
    else if(vm.count("new-eci")) {
      if(set.new_eci_dir(clex_desc.property, clex_desc.calctype, clex_desc.ref, clex_desc.bset, single_input)) {
        clex_desc.eci = single_input;
        std::cout << "Created new eci '" << single_input << "'.\n\n";
      }
      else {
        std::cout << "Could not create new eci '" << single_input << "'.\n\n";
        return 1;
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        //try {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, false, false, true);
        }
        /*}
        catch(std::exception &e) {
          args.err_log.error("unknown");
          args.err_log << "something happened\n" << std::endl;
        }
        */
        std::cout << "Switched to new eci '" << clex_desc.eci << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to new eci '" << clex_desc.eci << "'.\n\n";
        return 1;
      }
    }

    // create new cluster expansion settings group
    else if(vm.count("new-clex")) {
      clex_desc.name = single_input;
      if(set.new_clex(clex_desc)) {
        std::cout << "Created new cluster expansion named '" << single_input << "'\n\n";
      }
      else {
        std::cout << "Could not create new cluster expansion named '" << single_input << "'.\n\n";
        return 1;
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Set '" << clex_desc.name << "' as default cluster expansion.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not set '" << clex_desc.name << "' as default cluster expansion.\n\n";
        return 1;
      }
    }

    // erase cluster expansion settings group
    else if(vm.count("erase-clex")) {
      if(set.default_clex().name == single_input) {
        std::cout << "Coult not erase the cluster expansion named '" << single_input << "' "
                  << "because it is the default cluster expansion.\n\n";
        return 1;
      }

      if(set.cluster_expansions().size() == 1) {
        std::cout << "Could not erase the cluster expansion named '" << single_input << "' "
                  << "because it is the only cluster expansion.\n\n";
        return 1;
      }

      if(set.erase_clex(set.clex(single_input))) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Erased cluster expansion named '" << single_input << "'\n\n";
        return 0;
      }
      else {
        std::cout << "Could not erase cluster expansion named '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set clex to use by default for formation energy
    else if(vm.count("set-formation-energy")) {

      if(clex_desc.property != "formation_energy") {
        std::cerr << "Attempting to use cluster expansion '" << clex_desc.name
                  << "' for formation_energy, but it has 'property' value '"
                  << clex_desc.property << "'.\n\n";
        return ERR_INVALID_ARG;
      }

      auto it = set.cluster_expansions().find("formation_energy");
      if(it != set.cluster_expansions().end()) {
        set.erase_clex(it->second);
      }
      clex_desc.name = "formation_energy";
      if(set.new_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Now using cluster expansion '" << single_input
                  << "' as the default for formation_energy.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not use cluster expansion '" << single_input
                  << "' as the default for formation_energy.\n\n";
        return 1;
      }
    }

    // set default clex
    else if(vm.count("set-default-clex")) {
      if(set.set_default_clex(single_input)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Switched to cluster expansion '" << single_input << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to cluster expansion '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set bset
    else if(vm.count("set-property")) {
      clex_desc.property = single_input;
      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Switched to property '" << single_input << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to property '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set bset
    else if(vm.count("set-bset")) {
      clex_desc.bset = single_input;
      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Switched to bset '" << single_input << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to bset '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set calctype. if ref given, set also; else select default, else select first found
    else if(vm.count("set-calctype")) {

      // if ref given
      if(multi_input.size() > 1) {
        clex_desc.calctype = multi_input[0];
        clex_desc.ref = multi_input[1];
      }
      // else figure out ref
      else {
        std::vector<std::string> all = dir.all_ref(clex_desc.calctype);

        // if no ref, create default
        if(!all.size()) {
          std::cout << "No ref found. Creating 'default'.\n\n";
          set.new_ref_dir(clex_desc.calctype, "default");
          clex_desc.ref = "default";
        }
        // if default exists, set to default
        else if(std::find(all.begin(), all.end(), "default") != all.end()) {
          clex_desc.ref = "default";
        }
        // else set first found
        else {
          clex_desc.ref = all[0];
        }
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Switched to calctype '" << clex_desc.calctype << "' and ref '" << clex_desc.ref << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to calctype '" << clex_desc.calctype << "' and ref '" << clex_desc.ref << "'.\n\n";
        return 1;
      }
    }

    // if given calctype and ref, set both; if only given one argument set ref in current calctype
    else if(vm.count("set-ref")) {
      if(multi_input.size() > 1) {
        clex_desc.calctype = multi_input[0];
        clex_desc.ref = multi_input[1];
      }
      else {
        clex_desc.ref = multi_input[0];
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        std::cout << "Switched to calctype '" << clex_desc.calctype << "' and ref '" << clex_desc.ref << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to calctype '" << clex_desc.calctype << "' and ref '" << clex_desc.ref << "'.\n\n";
        return 1;
      }
    }

    // set eci
    else if(vm.count("set-eci")) {
      clex_desc.eci = single_input;
      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, false, false, true);
        }
        std::cout << "Switched to eci '" << single_input << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to eci '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set cxx
    else if(vm.count("set-cxx")) {
      set.set_cxx(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      std::cout << "Set " << _wdefaultval("cxx", set.cxx());
      std::cout << "Compile command is now: '" << set.compile_options() << "'\n\n";
      std::cout << "Shard object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set cxxflags
    else if(vm.count("set-cxxflags")) {
      set.set_cxxflags(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      std::cout << "Set " << _wdefaultval("cxxflags", set.cxxflags());
      std::cout << "Compile command is now: '" << set.compile_options() << "'\n\n";

      return 0;
    }

    // set soflags
    else if(vm.count("set-soflags")) {
      set.set_soflags(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      std::cout << "Set " << _wdefaultval("soflags", set.soflags());
      std::cout << "Shard object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set casm prefix
    else if(vm.count("set-casm-prefix")) {
      set.set_casm_prefix(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      std::cout << "Set " << _wdefaultval("casm_prefix", set.casm_prefix());
      std::cout << "Compile command is now: '" << set.compile_options() << "'\n\n";
      std::cout << "Shard object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set boost prefix
    else if(vm.count("set-boost-prefix")) {
      set.set_casm_prefix(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      std::cout << "Set " << _wdefaultval("boost_prefix", set.boost_prefix());
      std::cout << "Compile command is now: '" << set.compile_options() << "'\n\n";
      std::cout << "Shard object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set 'casm view' command
    else if(vm.count("set-view-command")) {
      set.set_view_command(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, false);
      }

      std::cout << "Set view command to: '" << set.view_command() << "'\n\n";

      return 0;
    }

    std::cout << std::endl;

    return 0;
  };


}
