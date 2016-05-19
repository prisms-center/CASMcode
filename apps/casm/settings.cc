#include <cstring>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"

namespace CASM {

  void list_settings(const fs::path &root) {

    DirectoryStructure dir(root);
    ProjectSettings set(root);

    std::vector<std::string> all = dir.all_bset();
    std::cout << "Basis sets:\n";
    for(int i = 0; i < all.size(); i++) {
      std::cout << "  " << all[i];
      if(set.bset() == all[i]) {
        std::cout << "*";
      }
      std::cout << "\n";
    }
    std::cout << "\n";

    all = dir.all_calctype();
    std::cout << "Training Data & References:\n";
    for(int i = 0; i < all.size(); i++) {
      std::vector<std::string> all_ref = dir.all_ref(all[i]);
      for(int j = 0; j < all_ref.size(); j++) {
        std::cout << "  " << all[i] << " / " << all_ref[j];
        if(set.calctype() == all[i] && set.ref() == all_ref[j]) {
          std::cout << "*";
        }
        std::cout << "\n";
      }
      std::cout << "\n";
    }

    all = dir.all_eci(set.clex(), set.calctype(), set.ref(), set.bset());
    std::cout << "ECI for current settings:\n";
    for(int i = 0; i < all.size(); i++) {
      std::cout << "  " << all[i];
      if(set.eci() == all[i]) {
        std::cout << "*";
      }
      std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "Compiler settings: '" << set.compile_options() << "'\n\n";

    std::cout << "SO settings: '" << set.so_options() << "'\n\n";

    std::cout << "'casm view' command: '" << set.view_command() << "'\n\n";


  }

  // ///////////////////////////////////////
  // 'settings' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int settings_command(const CommandArgs& args) {

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
      ("new-bset", po::value<std::string>(&single_input), "Create a new basis set")
      ("new-calctype", po::value<std::vector<std::string> >(&multi_input)->multitoken(), "Create a new calculation type")
      ("new-ref", po::value<std::string>(&single_input), "Create a new calculation reference")
      ("new-eci", po::value<std::string>(&single_input), "Create a new set of effective clust interactions (ECI)")
      ("set-bset", po::value<std::string>(&single_input), "Set the current basis set")
      ("set-calctype", po::value<std::vector<std::string> >(&multi_input)->multitoken(), "Set the current calculation type")
      ("set-ref", po::value<std::vector<std::string> >(&multi_input)->multitoken(), "Set the current calculation reference")
      ("set-eci", po::value<std::string>(&single_input), "Set the current effective clust interactions (ECI)")
      ("set-view-command", po::value<std::string>(&single_input), "Set the command used by 'casm view'.")
      ("set-compile-options", po::value<std::string>(&single_input), "Set the compiler options.")
      ("unset-compile-options", "Use the default compiler options.")
      ("set-so-options", po::value<std::string>(&single_input), "Set the options for generating shared libraries.")
      ("unset-so-options", "Use the default options for generating shared libraries.");

      try {
        po::store(po::parse_command_line(args.argc, args.argv, desc), vm); // can throw

        bool call_help = false;

        std::vector<std::string> all_opt = {"list", "new-bset", "new-calctype", "new-ref", "new-eci",
                                            "set-bset", "set-calctype", "set-ref", "set-eci",
                                            "set-compile-options", "set-so-options", "set-view-command",
                                            "unset-compile-options", "unset-so-options",
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
          std::cout << "Error in 'casm settings'. Use one option at a time." << std::endl;
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

                    "      casm settings --new-bset 'my_new_bset'                 \n" <<
                    "      casm settings --new-calctype 'my_new_calctype' ['my_new_ref']\n" <<
                    "      casm settings --new-ref 'my_new_ref'                   \n" <<
                    "      casm settings --new-eci 'my_new_eci'                   \n" <<
                    "      - Creates new settings directories with appropriate    \n" <<
                    "        names                                                \n" <<
                    "      - For --new-calctype, a new reference may optionally be\n" <<
                    "        speficied. If it is not, a new 'default' reference is\n" <<
                    "        created.                                     \n" <<
                    "      - For --new-ref, a new reference is created for the    \n" <<
                    "        current calctype                                     \n" <<
                    "      - For --new-eci, a new eci directory is created for the\n" <<
                    "         current clex, calctype and ref.                     \n\n" <<

                    "      casm settings --set-bset 'other_bset'                  \n" <<
                    "      casm settings --set-calctype 'other_calctype' ['other_ref'] \n" <<
                    "      casm settings --set-ref ['other_calctype'] 'other_ref' \n" <<
                    "      casm settings --set-eci 'other_eci'                    \n" <<
                    "      - Switch the current settings                          \n" <<
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

    fs::path &root = args.root;
    DirectoryStructure dir(root);
    ProjectSettings set(root);


    // disply all settings
    if(vm.count("list")) {
      list_settings(root);
      return 0;
    }

    // create new bset, and default eci for current calctype and ref
    else if(vm.count("new-bset")) {
      if(set.new_bset_dir(single_input)) {
        set.set_bset(single_input);
        set.commit();
        std::cout << "Created and switched to new bset '" << single_input << "'.\n\n";

        // and create matching eci
        if(set.new_eci_dir(set.clex(), set.calctype(), set.ref(), set.bset(), "default")) {
          set.set_eci(set.clex(), set.calctype(), set.ref(), set.bset(), "default");
          set.commit();
          std::cout << "Created and switched to new eci '" << "default" << "'.\n\n";
          return 0;
        }
        else {
          std::cout << "Could not create new eci '" << "default" << "'.\n\n";
          return 1;
        }
      }
      else {
        std::cout << "Could not create new bset '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // create new calctype, and ref (either as specified, or default)
    else if(vm.count("new-calctype")) {

      // create new calctype
      if(set.new_calc_settings_dir(multi_input[0])) {
        set.set_calctype(multi_input[0]);
        set.commit();
        std::cout << "Created and switched to new calctype '" << multi_input[0] << "'.\n\n";

        // if ref given, create specified ref
        if(multi_input.size() > 1) {
          if(set.new_ref_dir(multi_input[0], multi_input[1])) {
            set.set_ref(multi_input[0], multi_input[1]);
            set.commit();
            std::cout << "Created and switched to new ref '" << multi_input[1] << "'.\n\n";

            // and create matching eci
            if(set.new_eci_dir(set.clex(), set.calctype(), set.ref(), set.bset(), "default")) {
              set.set_eci(set.clex(), set.calctype(), set.ref(), set.bset(), "default");
              set.commit();
              std::cout << "Created and switched to new eci '" << "default" << "'.\n\n";
              return 0;
            }
            else {
              std::cout << "Could not create new eci '" << "default" << "'.\n\n";
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
            set.set_ref(multi_input[0], "default");
            set.commit();
            std::cout << "Created and switched to new ref '" << "default" << "'.\n\n";

            // and create matching eci
            if(set.new_eci_dir(set.clex(), set.calctype(), set.ref(), set.bset(), "default")) {
              set.set_eci(set.clex(), set.calctype(), set.ref(), set.bset(), "default");
              set.commit();
              std::cout << "Created and switched to new eci '" << "default" << "'.\n\n";
              return 0;
            }
            else {
              std::cout << "Could not create new eci '" << "default" << "'.\n\n";
              return 1;
            }
          }
          else {
            std::cout << "Could not create new ref '" << "default" << "'.\n\n";
            return 1;
          }
        }

        return 0;
      }
      else {
        std::cout << "Could not create new calctype '" << multi_input[0] << "'.\n\n";
        return 1;
      }
    }

    // create new ref for current calctype
    else if(vm.count("new-ref")) {
      if(set.new_ref_dir(set.calctype(), single_input)) {
        set.set_ref(set.calctype(), single_input);
        set.commit();
        std::cout << "Created and switched to new ref '" << single_input << "'.\n\n";

        // and create matching eci
        if(set.new_eci_dir(set.clex(), set.calctype(), set.ref(), set.bset(), "default")) {
          set.set_eci(set.clex(), set.calctype(), set.ref(), set.bset(), "default");
          set.commit();
          std::cout << "Created and switched to new eci '" << "default" << "'.\n\n";
          return 0;
        }
        else {
          std::cout << "Could not create new eci '" << "default" << "'.\n\n";
          return 1;
        }
      }
      else {
        std::cout << "Could not create new ref '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // create new eci directory
    else if(vm.count("new-eci")) {
      if(set.new_eci_dir(set.clex(), set.calctype(), set.ref(), set.bset(), single_input)) {
        set.set_eci(set.clex(), set.calctype(), set.ref(), set.bset(), single_input);
        set.commit();
        std::cout << "Created and switched to new eci '" << single_input << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not create new eci '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set bset
    else if(vm.count("set-bset")) {
      if(set.set_bset(single_input)) {
        set.commit();
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
      // set calctype
      if(set.set_calctype(multi_input[0])) {
        set.commit();
        std::cout << "Switched to calctype '" << multi_input[0] << "'.\n\n";

        // if ref given
        if(multi_input.size() > 1) {

          // set ref
          if(set.set_ref(multi_input[0], multi_input[1])) {
            set.commit();
            std::cout << "Switched to ref '" << multi_input[1] << "'.\n\n";
            return 0;
          }
          else {
            std::cout << "Could not switch to ref '" << multi_input[1] << "'.\n\n";
            return 1;
          }

        }
        // else figure out
        else {
          std::vector<std::string> all = dir.all_ref(set.calctype());

          // if no ref, create default
          if(!all.size()) {
            if(set.new_ref_dir(set.calctype(), "default")) {
              set.set_ref(set.calctype(), "default");
              set.commit();
              std::cout << "Created and switched to new ref '" << "default" << "'.\n\n";
              return 0;
            }
            else {
              std::cout << "Found not ref, and could not create new ref '" << "default" << "'.\n\n";
              return 1;
            }
          }
          // if default exists, set to default
          else if(std::find(all.begin(), all.end(), "default") != all.end()) {
            if(set.set_ref(set.calctype(), "default")) {
              set.commit();
              std::cout << "Switched to ref '" << "default" << "'.\n\n";
              return 0;
            }
            else {
              std::cout << "Could not switch to bset '" << "default" << "'.\n\n";
              return 1;
            }
          }
          // else set first found
          else {
            if(set.set_ref(set.calctype(), all[0])) {
              set.commit();
              std::cout << "Switched to ref '" << all[0] << "'.\n\n";
              return 0;
            }
            else {
              std::cout << "Could not switch to bset '" << all[0] << "'.\n\n";
              return 1;
            }
          }

        }

        return 0;
      }
      else {
        std::cout << "Could not switch to calctype '" << multi_input[0] << "'.\n\n";
        return 1;
      }
    }

    // if given calctype and ref, set both; if only given one argument set ref in current calctype
    else if(vm.count("set-ref")) {
      if(multi_input.size() > 1) {
        if(set.set_calctype(multi_input[0])) {
          set.commit();
          std::cout << "Switched to calctype '" << multi_input[0] << "'.\n\n";

          if(set.set_ref(set.calctype(), multi_input[1])) {
            set.commit();
            std::cout << "Switched to ref '" << multi_input[1] << "'.\n\n";
            return 0;
          }
          else {
            std::cout << "Could not switch to ref '" << multi_input[1] << "'.\n\n";
            return 1;
          }
        }
        else {
          std::cout << "Could not switch to calctype '" << multi_input[0] << "'.\n\n";
          return 1;
        }
      }
      else {
        if(set.set_ref(set.calctype(), multi_input[0])) {
          set.commit();
          std::cout << "Switched to ref '" << multi_input[0] << "'.\n\n";
          return 0;
        }
        else {
          std::cout << "Could not switch to ref '" << multi_input[0] << "'.\n\n";
          return 1;
        }
      }
    }

    // set eci
    else if(vm.count("set-eci")) {
      if(set.set_eci(set.clex(), set.calctype(), set.ref(), set.bset(), single_input)) {
        set.commit();
        std::cout << "Switched to eci '" << single_input << "'.\n\n";
        return 0;
      }
      else {
        std::cout << "Could not switch to eci '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set compile options
    else if(vm.count("set-compile-options")) {
      set.set_compile_options(single_input);
      set.commit();

      std::cout << "Set compile options to: '" << set.compile_options() << "'\n\n";

      return 0;
    }

    // set compile options
    else if(vm.count("unset-compile-options")) {
      set.set_compile_options(RuntimeLibrary::default_compile_options());
      set.commit();

      std::cout << "Set compile options to: '" << set.compile_options() << "'\n\n";

      return 0;
    }

    // set compile options
    else if(vm.count("set-so-options")) {
      set.set_so_options(single_input);
      set.commit();

      std::cout << "Set so options to: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set compile options
    else if(vm.count("unset-so-options")) {
      set.set_so_options(RuntimeLibrary::default_so_options());
      set.commit();

      std::cout << "Set so options to: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set compile options
    else if(vm.count("set-view-command")) {
      set.set_view_command(single_input);
      set.commit();

      std::cout << "Set view command to: '" << set.view_command() << "'\n\n";

      return 0;
    }

    std::cout << std::endl;

    return 0;
  };


}