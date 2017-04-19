#include "casm/app/casm_functions.hh"
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/Selection.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  void select_help(const DataFormatterDictionary<Configuration> &_dict, std::ostream &_stream, std::vector<std::string> help_opt) {
    _stream << "DESCRIPTION" << std::endl
            << "\n"
            << "    Use '[--set | --set-on | --set-off] [criteria]' for specifying or editing a selection.\n";

    for(const std::string &str : help_opt) {
      if(str.empty()) {
        continue;
      }

      if(str[0] == 'o') {
        _stream << "Available operators for use within selection criteria:" << std::endl;
        _dict.print_help(_stream, BaseDatumFormatter<Configuration>::Operator);
      }
      else if(str[0] == 'p') {
        _stream << "Available property tags are currently:" << std::endl;
        _dict.print_help(_stream, BaseDatumFormatter<Configuration>::Property);
      }
      _stream << std::endl;
    }
    _stream << std::endl;
  }

  namespace Completer {
    SelectOption::SelectOption(): OptionHandlerBase("select") {}

    const std::vector<std::string> &SelectOption::criteria_vec() const {
      return m_criteria_vec;
    }

    void SelectOption::initialize() {
      add_general_help_suboption();
      add_configlists_suboption();
      add_output_suboption("MASTER");

      m_desc.add_options()
      ("json", "Write JSON output (otherwise CSV, unless output extension is '.json' or '.JSON')")
      ("subset", "Only write selected configurations to output. Can be used by itself or in conjunction with other options")
      ("xor", "Performs logical XOR on two configuration selections")
      ("not", "Performs logical NOT on configuration selection")
      ("or", "Write configurations selected in any of the input lists. Equivalent to logical OR")
      ("and", "Write configurations selected in all of the input lists. Equivalent to logical AND")
      ("set-on", po::value<std::vector<std::string> >(&m_criteria_vec)->multitoken()->zero_tokens(), "Add configurations to selection if they meet specified criteria.  Call using 'casm select --set-on [\"criteria\"]'")
      ("set-off", po::value<std::vector<std::string> >(&m_criteria_vec)->multitoken()->zero_tokens(), "Remove configurations from selection if they meet specified criteria.  Call using 'casm select --set-off [\"criteria\"]'")
      ("set", po::value<std::vector<std::string> >(&m_criteria_vec)->multitoken(), "Create a selection of Configurations that meet specified criteria.  Call using 'casm select --set [\"criteria\"]'")
      ("force,f", "Overwrite output file");

      return;
    }

  }

  void write_selection_stats(Index Ntot, const DB::Selection<Configuration> &config_select, Log &log, bool only_selected) {

    auto Nselected = config_select.selected_size();
    auto Ninclude = only_selected ? Nselected : config_select.size();

    log << "# configurations in this project: " << Ntot << "\n";
    log << "# configurations included in this list: " << Ninclude << "\n";
    log << "# configurations selected in this list: " << Nselected << "\n";
  }


  // ///////////////////////////////////////
  // 'select' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int select_command(const CommandArgs &args) {

    //casm enum [—supercell min max] [—config supercell ] [—hopconfigs hop.background]
    //- enumerate supercells and configs and hop local configurations

    std::vector<std::string> criteria_vec, help_opt_vec;
    std::vector<fs::path> selection;

    fs::path out_path;
    COORD_TYPE coordtype;
    po::variables_map vm;

    /// Set command line options using boost program_options
    // NOTE: multitoken() is used instead of implicit_value() because implicit_value() is broken on some systems -- i.e., braid.cnsi.ucsb.edu
    //       (not sure if it's an issue with a particular shell, or boost version, or something else)
    Completer::SelectOption select_opt;

    std::string cmd;
    std::vector<std::string> allowed_cmd = {"and", "or", "xor", "not", "set-on", "set-off", "set"};

    try {
      po::store(po::parse_command_line(args.argc, args.argv, select_opt.desc()), vm); // can throw

      Index num_cmd(0);
      for(const std::string &cmd_str : allowed_cmd) {
        if(vm.count(cmd_str)) {
          num_cmd++;
          cmd = cmd_str;
        }
      }

      if(!vm.count("help")) {
        if(num_cmd > 1) {
          args.err_log << "Error in 'casm select'. Must use exactly one of --set-on, --set-off, --set, --and, --or, --xor, or --not." << std::endl;
          return ERR_INVALID_ARG;
        }
        else if(vm.count("subset") && vm.count("config") && selection.size() != 1) {
          args.err_log << "ERROR: 'casm select --subset' expects zero or one list as argument." << std::endl;
          return ERR_INVALID_ARG;
        }



        if(!vm.count("output") && (cmd == "or" || cmd == "and" || cmd == "xor" || cmd == "not")) {
          args.err_log << "ERROR: 'casm select --" << cmd << "' expects an --output file." << std::endl;
          return ERR_INVALID_ARG;
        }

      }

      // Start --help option
      if(vm.count("help")) {
        args.log << std::endl << select_opt.desc() << std::endl;
      }

      po::notify(vm); // throws on error, so do after help in case of problems

      criteria_vec = select_opt.criteria_vec();
      help_opt_vec = select_opt.help_opt_vec();
      selection = select_opt.selection_paths();
      out_path = select_opt.output_path();

      // Finish --help option
      if(vm.count("help")) {
        const fs::path &root = args.root;
        if(root.empty()) {
          auto dict = make_dictionary<Configuration>();
          select_help(dict, args.log, help_opt_vec);
        }
        else {
          // set status_stream: where query settings and PrimClex initialization messages are sent
          Log &status_log = (out_path.string() == "STDOUT") ? args.err_log : args.log;

          // If '_primclex', use that, else construct PrimClex in 'uniq_primclex'
          // Then whichever exists, store reference in 'primclex'
          std::unique_ptr<PrimClex> uniq_primclex;
          if(out_path.string() == "STDOUT") {
            args.log.set_verbosity(0);
          }
          PrimClex &primclex = make_primclex_if_not(args, uniq_primclex, status_log);

          select_help(primclex.settings().query_handler<Configuration>().dict(), args.log, help_opt_vec);
        }
        return 0;
      }

      if((vm.count("set-on") || vm.count("set-off") || vm.count("set")) && vm.count("config") && selection.size() != 1) {
        std::string cmd = "--set-on";
        if(vm.count("set-off")) {
          cmd = "--set-off";
        }
        if(vm.count("set")) {
          cmd = "--set";
        }

        args.err_log << "Error in 'casm select " << cmd << "'. " << selection.size() << " config selections were specified, but no more than one selection is allowed (MASTER list is used if no other is specified)." << std::endl;
        return ERR_INVALID_ARG;
      }

    }
    catch(po::error &e) {
      args.err_log << select_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log << select_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;
    }


    if(vm.count("output") && out_path != "MASTER") {

      //check now so we can exit early with an obvious error
      if(fs::exists(out_path) && !vm.count("force")) {
        args.err_log << "ERROR: File " << out_path << " already exists. Use --force to force overwrite." << std::endl;
        return ERR_EXISTING_FILE;
      }
    }

    bool only_selected(false);
    if(selection.empty()) {
      only_selected = true;
      selection.push_back("MASTER");
    }

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
    ProjectSettings &set = primclex.settings();

    // count total number of configurations in this project one time
    Index Ntot = primclex.db<Configuration>().size();

    // ---- load initial selection into config_select ----
    // this is also the selection that will be printed at end
    DB::Selection<Configuration> config_select(primclex.db<Configuration>(), selection[0]);

    std::stringstream ss;
    ss << selection[0];

    std::vector<DB::Selection<Configuration> > tselect;
    tselect.reserve(selection.size() - 1);
    for(int i = 1; i < selection.size(); ++i) {
      ss << ", " << selection[i];
      tselect.emplace_back(primclex.db<Configuration>(), selection[i]);
    }

    set.query_handler<Configuration>().set_selected(config_select);

    // ---- write starting stats ----
    args.log.custom("Input config list", selection[0].string());
    write_selection_stats(Ntot, config_select, args.log, false);
    args.log << std::endl;

    for(int i = 1; i < selection.size(); ++i) {
      args.log.custom("Input config list", selection[i].string());
      write_selection_stats(Ntot, tselect[i], args.log, false);
      args.log << std::endl;
    }

    // ---- perform requested logic ----
    if(vm.count("set-on") || vm.count("set-off") || vm.count("set")) {
      bool select_switch = vm.count("set-on");
      std::string criteria;
      if(criteria_vec.size() == 1) {
        criteria = criteria_vec[0];
      }
      else if(criteria_vec.size() > 1) {
        args.err_log << "ERROR: Selection criteria must be a single string.  You provided " << criteria_vec.size() << " strings:\n";
        for(const std::string &str : criteria_vec)
          args.err_log << "     - " << str << "\n";
        return ERR_INVALID_ARG;
      }

      if(vm.count("set-on")) {
        args.log.custom("set-on", criteria);
      }
      if(vm.count("set-off")) {
        args.log.custom("set-off", criteria);
      }
      if(vm.count("set")) {
        args.log.custom("set", criteria);
      }
      args.log.begin_lap();

      try {
        if(vm.count("set")) {
          config_select.set(set.query_handler<Configuration>().dict(), criteria);
        }
        else {
          config_select.set(set.query_handler<Configuration>().dict(), criteria, select_switch);
        }
      }
      catch(std::exception &e) {
        args.err_log << "ERROR: " << e.what() << "\n";
        return ERR_INVALID_ARG;
      }

      args.log << "selection time: " << args.log.lap_time() << " (s)\n" << std::endl;
    }

    if(vm.count("subset")) {
      args.log.custom("subset");
      args.log.begin_lap();
      args.log << "selection time: " << args.log.lap_time() << " (s)\n" << std::endl;
      only_selected = true;
    }

    if(vm.count("not")) {
      if(selection.size() != 1) {
        args.err_log << "ERROR: Option --not requires exactly 1 selection as argument\n";
        return ERR_INVALID_ARG;
      }

      args.log.custom(std::string("not ") + selection[0].string());
      args.log.begin_lap();

      // loop through other lists, keeping only configurations selected in the other lists
      for(auto it = config_select.all().begin(); it != config_select.all().end(); ++it) {
        it.is_selected() = !it.is_selected();
      }
      args.log << "selection time: " << args.log.lap_time() << " (s)\n" << std::endl;
    }

    if(vm.count("or")) {

      args.log.custom(std::string("or(") + ss.str() + ")");
      args.log.begin_lap();

      // loop through other lists, inserting all selected configurations
      for(int i = 1; i < selection.size(); i++) {
        for(const auto &val : tselect[i].data()) {
          if(!val.second) {
            continue;
          }
          auto res = config_select.data().insert({val.first, true});
          res.first->second = true;
        }
      }
      args.log << "selection time: " << args.log.lap_time() << " (s)\n" << std::endl;
      only_selected = true;
    }

    if(vm.count("and")) {
      args.log.custom(std::string("and(") + ss.str() + ")");
      args.log.begin_lap();

      // loop through other lists, keeping only configurations selected in the other lists
      for(int i = 1; i < selection.size(); i++) {
        for(const auto &val : config_select.data()) {
          if(!val.second) {
            continue;
          }
          auto res = config_select.data().insert({val.first, true});
          res.first->second = tselect[i].is_selected(val.first);
        }
      }

      args.log << "selection time: " << args.log.lap_time() << " (s)\n" << std::endl;
      only_selected = true;
    }

    if(vm.count("xor")) {
      if(selection.size() != 2) {
        args.err_log << "ERROR: Option --xor requires exactly 2 selections as argument\n";
        return 1;
      }

      args.log.custom(selection[0].string() + " xor " + selection[1].string());
      args.log.begin_lap();

      for(const auto &val : tselect[1].data()) {
        // if not selected in second, use 'config_select' 'is_selected' value
        if(!val.second) {
          continue;
        }
        // else, if selected in second:

        // if not in 'config_select' insert selected
        auto find_it = config_select.data().find(val.first);
        if(find_it == config_select.data().end()) {
          config_select.data().insert(val);
        }
        // else, use opposite of config_select 'is_selected' value
        else {
          find_it->second = !find_it->second;
        }

      }
      args.log << "selection time: " << args.log.lap_time() << " (s)\n" << std::endl;
      only_selected = true;
    }

    /// Only write selection to disk past this point

    args.log.write("Selection");
    int ret_code = config_select.write(
                     set.query_handler<Configuration>().dict(),
                     vm.count("force") || (out_path.string() == "MASTER"),
                     out_path,
                     vm.count("json"),
                     only_selected);
    args.log << "write: " << out_path << "\n" << std::endl;

    args.log.custom("Output config list", out_path.string());
    write_selection_stats(Ntot, config_select, args.log, only_selected);

    args.log << std::endl;
    return ret_code;

  };

}


