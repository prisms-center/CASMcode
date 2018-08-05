#include "casm/app/casm_functions.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/PrimClex.hh"

#include "casm/completer/Handlers.hh"

#include <functional>

namespace CASM {

  struct ClexDescription;
  bool clex_exists(const DirectoryStructure &dir, const ClexDescription &desc);

  namespace {
    std::string _wdefaultval(std::string name, std::pair<std::string, std::string> val) {
      return name + ": '" + val.first + "' (" + val.second + ")\n";
    }

    std::string _wdefaultval(std::string name, std::pair<fs::path, std::string> val) {
      return _wdefaultval(name, std::make_pair(val.first.string(), val.second));
    }
  }

  namespace Completer {
    SettingsOption::SettingsOption(): OptionHandlerBase("settings") {}

    const std::string &SettingsOption::input_str() const {
      return m_input_str;
    }

    const std::vector<std::string> &SettingsOption::input_vec() const {
      return m_input_vec;
    }

    void SettingsOption::initialize() {
      add_help_suboption();

      m_desc.add_options()
      ("list,l", "List project settings")
      ("new-property", po::value<std::string>(&m_input_str), "Create cluster expansions for a new property")
      ("new-bset", po::value<std::string>(&m_input_str), "Create a new basis set")
      ("new-calctype", po::value<std::string>(&m_input_str), "Create a new calculation type")
      ("new-ref", po::value<std::string>(&m_input_str), "Create a new calculation reference")
      ("new-eci", po::value<std::string>(&m_input_str), "Create a new set of effective cluster interactions (ECI)")
      ("set-formation-energy", po::value<std::string>(&m_input_str), "Specify the cluster expansion to use for formation energy")
      ("new-clex", po::value<std::string>(&m_input_str), "Create a new cluster expansion")
      ("set-default-clex", po::value<std::string>(&m_input_str), "Set the cluster expansion that CASM uses or acts on by default")
      ("erase-clex", po::value<std::string>(&m_input_str), "Erase the specified cluster expansion. Does not erase underlying bset, eci, etc.")
      ("clex", po::value<std::string>(), "The cluster expansion for which to set property, bset, calctype, ref, or eci")
      ("set-property", po::value<std::string>(&m_input_str), "Set the current basis set")
      ("set-bset", po::value<std::string>(&m_input_str), "Set the basis set")
      ("set-calctype", po::value<std::string>(&m_input_str), "Set the calculation type")
      ("set-ref", po::value<std::string>(&m_input_str), "Set the calculation reference")
      ("set-eci", po::value<std::string>(&m_input_str), "Set the effective cluster interactions (ECI)")
      ("set-all", po::value<std::vector<std::string> >(&m_input_vec)->multitoken(), "Set the current property, calctype, ref, bset, and eci all at once.")
      ("set-view-command", po::value<std::string>(&m_input_str), "Set the command used by 'casm view'.")
      ("set-cxx", po::value<std::string>(&m_input_str), "Set the c++ compiler. Use '' to revert to default.")
      ("set-cxxflags", po::value<std::string>(&m_input_str), "Set the c++ compiler options. Use '' to revert to default.")
      ("set-soflags", po::value<std::string>(&m_input_str), "Set the shared library compilation options. Use '' to revert to default.")
      ("set-casm-prefix", po::value<std::string>(&m_input_str), "Set the casm prefix. Use '' to revert to default.")
      ("set-casm-includedir", po::value<std::string>(&m_input_str), "Set the casm includedir. Use '' to revert to default.")
      ("set-casm-libdir", po::value<std::string>(&m_input_str), "Set the casm libdir. Use '' to revert to default.")
      ("set-boost-prefix", po::value<std::string>(&m_input_str), "Set the boost prefix. Use '' to revert to default.")
      ("set-boost-includedir", po::value<std::string>(&m_input_str), "Set the boost includedir. Use '' to revert to default.")
      ("set-boost-libdir", po::value<std::string>(&m_input_str), "Set the boost libdir. Use '' to revert to default.");
      return;
    }

  }

  namespace {

    /// "create": for 'new-X'
    /// "do_not_create": for 'set-X'
    enum create_mode {create, do_not_create};

    /// Just a local data structure to ease passing
    struct Data {

      typedef std::pair<std::string, create_mode> pair_type;

      Data(PrimClex *_primclex,
           DirectoryStructure &_dir,
           ProjectSettings &_set,
           ClexDescription &_desc,
           Log &_log,
           Log &_err_log) :
        primclex(_primclex),
        dir(_dir),
        set(_set),
        desc(_desc),
        log(_log),
        err_log(_err_log),
        property(pair_type(desc.property, do_not_create)),
        calctype(pair_type(desc.calctype, do_not_create)),
        ref(pair_type(desc.ref, do_not_create)),
        bset(pair_type(desc.bset, do_not_create)),
        eci(pair_type(desc.eci, do_not_create)) {}

      PrimClex *primclex;
      DirectoryStructure &dir;
      ProjectSettings &set;
      ClexDescription &desc;
      Log &log;
      Log &err_log;

      pair_type property;
      pair_type calctype;
      pair_type ref;
      pair_type bset;
      pair_type eci;

      int update() {

        // the resulting clex if successful
        ClexDescription tdesc(desc.name, property.first, calctype.first, ref.first, bset.first, eci.first);

        bool res = try_new("property", property,
        [&]() {
          return contains(dir.all_property(), tdesc.property);
        },
        [&]() {
          return set.new_clex_dir(tdesc.property);
        }) &&
        try_new("calctype", calctype,
        [&]() {
          return contains(dir.all_calctype(), tdesc.calctype);
        },
        [&]() {
          return set.new_calc_settings_dir(tdesc.calctype);
        }) &&
        try_new("ref", ref,
        [&]() {
          return contains(dir.all_ref(tdesc.calctype), tdesc.ref);
        },
        [&]() {
          return set.new_ref_dir(tdesc.calctype, tdesc.ref);
        }) &&
        try_new("bset", bset,
        [&]() {
          return contains(dir.all_bset(), tdesc.bset);
        },
        [&]() {
          return set.new_bset_dir(tdesc.bset);
        }) &&
        try_new("eci", eci,
        [&]() {
          return contains(dir.all_eci(tdesc.property, tdesc.calctype, tdesc.ref, tdesc.bset), tdesc.eci);
        },
        [&]() {
          return set.new_eci_dir(tdesc.property, tdesc.calctype, tdesc.ref, tdesc.bset, tdesc.eci);
        });

        // If could not create settings directories
        if(!res) {
          return 1;
        }

        if(clex_exists(dir, tdesc) && set.set_default_clex(tdesc)) {
          set.commit();

          bool read_settings = true;
          bool read_composition = false;

          bool read_chem_ref = false;
          bool read_configs = false;
          if((desc.property != tdesc.property) ||
             (desc.calctype != tdesc.calctype) ||
             (desc.ref != tdesc.ref)) {
            read_chem_ref = true;
            read_configs = true;
          }

          bool clear_clex = false;
          if((desc.property != tdesc.property) ||
             (desc.bset != tdesc.bset) ||
             (desc.eci != tdesc.eci)) {
            clear_clex = true;
          }

          desc = tdesc;
          log << "Updated default settings:\n";
          bool is_default = true;
          int indent = 0;
          desc.print(log, is_default, indent);

          if(primclex) {
            primclex->refresh(read_settings, read_composition, read_chem_ref, read_configs, clear_clex);
          }
          return 0;
        }
        else {
          err_log << "Unknown error: Could not switch settings." << std::endl;
          err_log << "Tried to switch to:\n";
          tdesc.print(err_log, false, 0);
          return 1;
        }
      }

      /// \brief Try to create new settings, if requested
      ///
      /// \param set_name: 'property', 'calctype', 'ref', etc.
      /// \param set: pair of, for example, <'default', create_mode>
      /// \param check_f: function returns true if setting exists, false if not
      /// \param new_f: function to create setting
      ///
      /// if create and check_f():
      ///   return success;
      /// if create and !check_f():
      ///   return new_f();
      /// if do_not_create and !check_f():
      ///   return fail;
      /// if do_not_create and check_f():
      ///   return success;
      bool try_new(std::string set_name, pair_type set, std::function<bool ()> check_f, std::function<bool ()> new_f) {

        create_mode _cmode = set.second;

        if(_cmode == create && check_f()) {
          return true;
        }
        if(_cmode == create && !check_f()) {
          bool res = new_f();
          if(res) {
            log << "Created new " << set_name << ": '" << set.first << "'\n";

          }
          else {
            err_log << "Could not create new " << set_name << ": '" << set.first << "'\n";
          }
          return res;
        }

        if(_cmode == do_not_create && !check_f()) {
          err_log << "Error: The " << set_name << " named '" << set.first << "' does not exist.\n";
          err_log << "  Check your input or use --new-" << set_name << " to create it first.\n";

          return false;
        }

        if(_cmode == do_not_create && check_f()) {
          return true;
        }

        // shouldn't be able to get here
        throw std::runtime_error("Unknown error updating CASM settings");
      }
    };
  }


  // ///////////////////////////////////////
  // 'settings' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int settings_command(const CommandArgs &args) {


    std::string single_input;
    std::vector<std::string> multi_input;
    po::variables_map vm;


    try {
      /// Set command line options using boost program_options
      Completer::SettingsOption settings_opt;
      const po::options_description &desc = settings_opt.desc(); //I'm tired of fixing merge conflicts, this is ugly and should not stay like this.

      try {
        po::store(po::parse_command_line(args.argc, args.argv, desc), vm); // can throw

        bool call_help = false;

        std::vector<std::string> all_opt = {"list", "desc",
                                            "new-property", "new-bset", "new-calctype", "new-ref", "new-eci",
                                            "new-clex", "set-formation-energy", "erase-clex",
                                            "set-default-clex", "set-property", "set-bset", "set-calctype", "set-ref", "set-eci", "set-all",
                                            "set-cxx", "set-cxxflags", "set-soflags",
                                            "set-casm-prefix", "set-casm-includedir", "set-casm-libdir",
                                            "set-boost-prefix", "set-boost-includedir", "set-boost-libdir",
                                            "set-view-command"
                                           };
        int option_count = 0;
        for(int i = 0; i < all_opt.size(); i++) {
          option_count += vm.count(all_opt[i]);
        }

        // must call one and only one option at a time:
        if(!vm.count("help")) {
          if(option_count == 0) {
            args.log << "Error in 'casm settings'. No option selected." << std::endl;
            call_help = true;
          }
          else if(option_count > 1) {
            args.log << "Error in 'casm settings'. Use one option (other than --clex) at a time." << std::endl;
            call_help = true;
          }
        }

        // --help option
        if(vm.count("help") || call_help) {
          args.log << "\n";
          args.log << desc << std::endl;

          return 0;
        }

        if(vm.count("desc")) {
          args.log << "\n";
          args.log << desc << std::endl;
          args.log << "DESCRIPTION" << std::endl;
          args.log << "\n";
          args.log << "    Often it is useful to try multiple different basis sets, \n" <<
                   "    calculation settings, references, or ECI fits in a single\n" <<
                   "    project. The 'casm settings' option helps to organize    \n" <<
                   "    these within a project and quickly switch between        \n" <<
                   "    different settings.                                      \n";
          args.log << "\n";
          args.log << "    Examples:\n";
          args.log << "                                                             \n" <<
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
                   "      casm settings --new-calctype 'my_new_calctype'         \n" <<
                   "      casm settings --new-ref 'my_new_ref'                   \n" <<
                   "      casm settings --new-eci 'my_new_eci'                   \n" <<
                   "      - Creates new settings directories with appropriate    \n" <<
                   "        names                                                \n" <<
                   "      - For --new-property, new 'default' calctype, ref, and \n" <<
                   "        eci are created.                                     \n" <<
                   "      - For --new-calctype, a new 'default' ref is created.  \n" <<
                   "      - For --new-ref, a new reference is created for the    \n" <<
                   "        current calctype                                     \n" <<
                   "      - For --new-eci, a new eci directory is created for the\n" <<
                   "         current clex, calctype and ref.                     \n\n" <<

                   "      casm settings --set-property 'other_property'          \n" <<
                   "      casm settings --set-bset 'other_bset'                  \n" <<
                   "      casm settings --set-calctype 'other_calctype'          \n" <<
                   "      casm settings --set-ref 'other_ref'                    \n" <<
                   "      casm settings --set-eci 'other_eci'                    \n" <<
                   "      casm settings --set-all 'property' 'calctype', 'ref', 'bset', 'eci'\n" <<
                   "      - Switch the current settings                          \n" <<
                   "      - For --set-property, standard options are:            \n" <<
                   "        - 'formation_energy' (the only standard option for now)\n" <<
                   "        After switching, the 'default' calctype, ref, bset,  \n" <<
                   "        and eci are used"
                   "      - For --set-calctype, the current property and bset are\n" <<
                   "        maintained, and the 'default' ref, and eci are used. \n" <<
                   "      - For --set-ref, the current property, calctype, and   \n" <<
                   "        bset are maintained, and the 'default' eci is used.  \n" <<
                   "      - For --set-bset, the current property, calctype, and  \n" <<
                   "        ref are maintained, and the 'default' eci is used.   \n" <<
                   "      - For --set-eci, the current property, calctype, ref,  \n" <<
                   "        and bset are maintained.                             \n" <<
                   "      - For --set-all, all settings are switched at once.    \n\n" <<

                   "      casm settings --set-cxx 'cxx'                          \n" <<
                   "      - Specifies compiler to use. In order of priority: \n"
                   "        1) User specified by 'casm settings --set-cxx' \n"
                   "           (use '' to clear) \n"
                   "        2) $CASM_CXX \n"
                   "        3) $CXX \n"
                   "        4) \"g++\" \n\n"

                   "      casm settings --set-cxxflags 'cxxflags'                \n"
                   "      - Specifies compiler options. In order of priority:    \n"
                   "        1) User specified by 'casm settings --set-cxxflags' \n"
                   "           (use '' to clear) \n"
                   "        2) $CASM_CXXFLAGS \n"
                   "        3) \"-O3 -Wall -fPIC --std=c++11\" \n\n"

                   "      casm settings --set-soflags 'soflags' \n"
                   "      - Specifies shared object construction flags. In order \n"
                   "        of priority: \n"
                   "        1) User specified by 'casm settings --set-soflags' \n"
                   "           (use '' to clear) \n"
                   "        2) $CASM_SOFLAGS \n"
                   "        3) \"-shared -lboost_system\" \n\n"

                   "      casm settings --set-casm-prefix 'casm_prefix' \n"
                   "      casm settings --set-casm-includedir 'casm_includedir' \n"
                   "      casm settings --set-casm-libdir 'casm_libdir' \n"
                   "      - Specifies location to find CASM header files and shared\n"
                   "        libraries for runtime compilation and linking.         \n"
                   "        In order of priority: \n"
                   "        1) User specified by via 'casm settings' (use '' to clear) \n"
                   "        2) $CASM_INCLUDEDIR and $CASM_LIBDIR \n"
                   "        3) $CASM_PREFIX/include and $CASM_PREFIX/lib \n"
                   "        4) (default search paths) \n\n"

                   "      casm settings --set-boost-prefix 'boost_prefix' \n"
                   "      casm settings --set-boost-includedir 'boost_includedir' \n"
                   "      casm settings --set-boost-libdir 'boost_libdir' \n"
                   "      - Specifies location to find boost header files and shared\n"
                   "        libraries for runtime compilation and linking.         \n"
                   "        In order of priority: \n"
                   "        1) User specified by via 'casm settings' (use '' to clear) \n"
                   "        2) $CASM_BOOST_INCLUDEDIR and $CASM_BOOST_LIBDIR \n"
                   "        3) $CASM_BOOST_PREFIX/include and $CASM_BOOST_PREFIX/lib \n"
                   "        4) (default search paths) \n\n"

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

        single_input = settings_opt.input_str();
        multi_input = settings_opt.input_vec();

      }
      catch(po::error &e) {
        args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
        args.err_log << desc << std::endl;
        return 1;
      }
    }
    catch(std::exception &e) {
      args.err_log << "Unhandled Exception reached the top of main: "
                   << e.what() << ", application will now exit" << std::endl;
      return 1;

    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << "current_path: " << fs::current_path() << std::endl;
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

    Data d(args.primclex, dir, set, clex_desc, args.log, args.err_log);

    typedef Data::pair_type pair_type;

    // disply all settings
    if(vm.count("list")) {
      set.print_summary(args.log);
      return 0;
    }

    // create new cluster_expansions/clex.property directory and sub-directories
    else if(vm.count("new-property")) {

      d.property = pair_type(single_input, create);
      d.calctype = pair_type("default", create);
      d.ref = pair_type("default", create);
      d.bset = pair_type("default", create);
      d.eci = pair_type("default", create);

      return d.update();
    }

    // create new bset, and default eci for current calctype and ref
    else if(vm.count("new-bset")) {

      d.bset = pair_type(single_input, create);
      d.eci = pair_type("default", create);

      return d.update();
    }

    // create new calctype, and ref (either as specified, or default)
    else if(vm.count("new-calctype")) {

      d.calctype = pair_type(single_input, create);
      d.ref = pair_type("default", create);
      d.eci = pair_type("default", create);

      return d.update();
    }

    // create new ref for current calctype
    else if(vm.count("new-ref")) {

      d.ref = pair_type(single_input, create);
      d.eci = pair_type("default", create);

      return d.update();
    }

    // create new eci directory
    else if(vm.count("new-eci")) {

      d.eci = pair_type(single_input, create);

      return d.update();
    }

    // create new cluster expansion settings group
    else if(vm.count("new-clex")) {
      clex_desc.name = single_input;
      if(set.new_clex(clex_desc)) {
        args.log << "Created new cluster expansion named '" << single_input << "'\n\n";
      }
      else {
        args.log << "Could not create new cluster expansion named '" << single_input << "'.\n\n";
        return 1;
      }

      if(clex_exists(dir, clex_desc) && set.set_default_clex(clex_desc)) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        args.log << "Set '" << clex_desc.name << "' as default cluster expansion.\n\n";
        return 0;
      }
      else {
        args.log << "Could not set '" << clex_desc.name << "' as default cluster expansion.\n\n";
        return 1;
      }
    }

    // erase cluster expansion settings group
    else if(vm.count("erase-clex")) {
      if(set.default_clex().name == single_input) {
        args.log << "Coult not erase the cluster expansion named '" << single_input << "' "
                 << "because it is the default cluster expansion.\n\n";
        return 1;
      }

      if(set.cluster_expansions().size() == 1) {
        args.log << "Could not erase the cluster expansion named '" << single_input << "' "
                 << "because it is the only cluster expansion.\n\n";
        return 1;
      }

      if(set.erase_clex(set.clex(single_input))) {
        set.commit();
        if(args.primclex) {
          args.primclex->refresh(true, false, true, true, true);
        }
        args.log << "Erased cluster expansion named '" << single_input << "'\n\n";
        return 0;
      }
      else {
        args.log << "Could not erase cluster expansion named '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set clex to use by default for formation energy
    else if(vm.count("set-formation-energy")) {

      if(clex_desc.property != "formation_energy") {
        args.err_log << "Attempting to use cluster expansion '" << clex_desc.name
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
        args.log << "Now using cluster expansion '" << single_input
                 << "' as the default for formation_energy.\n\n";
        return 0;
      }
      else {
        args.log << "Could not use cluster expansion '" << single_input
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
        args.log << "Switched to cluster expansion '" << single_input << "'.\n\n";
        return 0;
      }
      else {
        args.log << "Could not switch to cluster expansion '" << single_input << "'.\n\n";
        return 1;
      }
    }

    // set property
    else if(vm.count("set-property")) {

      d.property = pair_type(single_input, do_not_create);
      d.calctype = pair_type("default", do_not_create);
      d.ref = pair_type("default", do_not_create);
      d.bset = pair_type("default", do_not_create);
      d.eci = pair_type("default", create);

      return d.update();
    }

    // set bset
    else if(vm.count("set-bset")) {

      d.bset = pair_type(single_input, do_not_create);
      d.eci = pair_type("default", create);

      return d.update();
    }

    // set calctype
    else if(vm.count("set-calctype")) {

      d.calctype = pair_type(single_input, do_not_create);
      d.ref = pair_type("default", do_not_create);
      d.eci = pair_type("default", create);

      return d.update();
    }

    // set ref
    else if(vm.count("set-ref")) {

      d.ref = pair_type(single_input, do_not_create);
      d.eci = pair_type("default", create);

      return d.update();
    }

    // set eci
    else if(vm.count("set-eci")) {

      d.eci = pair_type(single_input, do_not_create);

      return d.update();
    }

    // set all
    else if(vm.count("set-all")) {

      if(multi_input.size() != 5) {
        args.log << "Error using --set-all: Expected 5 arguments: 'property' 'calctype' 'ref' 'bset' 'eci'" << std::endl;
        return ERR_INVALID_ARG;
      }

      d.property = pair_type(multi_input[0], do_not_create);
      d.calctype = pair_type(multi_input[1], do_not_create);
      d.ref = pair_type(multi_input[2], do_not_create);
      d.bset = pair_type(multi_input[3], do_not_create);
      d.eci = pair_type(multi_input[4], create);

      return d.update();
    }

    // set cxx
    else if(vm.count("set-cxx")) {
      set.set_cxx(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("cxx", set.cxx());
      args.log << "Compile command is now: '" << set.compile_options() << "'\n\n";
      args.log << "Shared object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set cxxflags
    else if(vm.count("set-cxxflags")) {
      set.set_cxxflags(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("cxxflags", set.cxxflags());
      args.log << "Compile command is now: '" << set.compile_options() << "'\n\n";

      return 0;
    }

    // set soflags
    else if(vm.count("set-soflags")) {
      set.set_soflags(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("soflags", set.soflags());
      args.log << "Shared object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set casm prefix
    else if(vm.count("set-casm-prefix")) {
      set.set_casm_prefix(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("casm_includedir", set.casm_includedir());
      args.log << "Set " << _wdefaultval("casm_libdir", set.casm_libdir());
      args.log << "Compile command is now: '" << set.compile_options() << "'\n\n";
      args.log << "Shared object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set casm includedir
    else if(vm.count("set-casm-includedir")) {
      set.set_casm_includedir(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("casm_includedir", set.casm_includedir());
      args.log << "Compile command is now: '" << set.compile_options() << "'\n\n";
      args.log << "Shared object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set casm libdir
    else if(vm.count("set-casm-libdir")) {
      set.set_casm_libdir(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("casm_libdir", set.casm_libdir());
      args.log << "Compile command is now: '" << set.compile_options() << "'\n\n";
      args.log << "Shared object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set boost prefix
    else if(vm.count("set-boost-prefix")) {
      set.set_boost_prefix(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("boost_includedir", set.boost_includedir());
      args.log << "Set " << _wdefaultval("boost_libdir", set.boost_libdir());
      args.log << "Compile command is now: '" << set.compile_options() << "'\n\n";
      args.log << "Shared object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set boost includedir
    else if(vm.count("set-boost-includedir")) {
      set.set_boost_includedir(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("boost_includedir", set.boost_includedir());
      args.log << "Compile command is now: '" << set.compile_options() << "'\n\n";
      args.log << "Shared object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set boost libdir
    else if(vm.count("set-boost-libdir")) {
      set.set_boost_libdir(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, true);
      }

      args.log << "Set " << _wdefaultval("boost_libdir", set.boost_libdir());
      args.log << "Compile command is now: '" << set.compile_options() << "'\n\n";
      args.log << "Shared object compile command is now: '" << set.so_options() << "'\n\n";

      return 0;
    }

    // set 'casm view' command
    else if(vm.count("set-view-command")) {
      set.set_view_command(single_input);
      set.commit();
      if(args.primclex) {
        args.primclex->refresh(true, false, false, false, false);
      }

      args.log << "Set view command to: '" << set.view_command() << "'\n\n";

      return 0;
    }

    args.log << std::endl;

    return 0;
  };


}
