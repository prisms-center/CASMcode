#include "select.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"

namespace CASM {

  template<typename ConfigIterType>
  void set_selection(ConfigIterType begin, ConfigIterType end, const std::string &criteria, bool mk) {
    //boost::trim(criteria);

    if(criteria.size()) {
      DataFormatter<Configuration> tformat(ConfigIOParser::parse(criteria));
      for(; begin != end; ++begin) {
        ValueDataStream<bool> select_stream;
        select_stream << tformat(*begin);
        if(select_stream.value())
          begin.set_selected(mk);
      }
    }
    else {
      for(; begin != end; ++begin) {
        begin.set_selected(mk);
      }
    }

    return;
  }

  template<typename ConfigIterType>
  void set_selection(ConfigIterType begin, ConfigIterType end, const std::string &criteria) {
    //boost::trim(criteria);

    if(criteria.size()) {
      DataFormatter<Configuration> tformat(ConfigIOParser::parse(criteria));
      for(; begin != end; ++begin) {
        ValueDataStream<bool> select_stream;
        select_stream << tformat(*begin);
        begin.set_selected(select_stream.value());
      }
    }

    return;
  }

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

  void select_help(std::ostream &_stream, std::string help_opt = "") {
    _stream << "DESCRIPTION" << std::endl
            << "\n"
            << "    Use '[--set | --on | --off] [criteria]' for specifying or editing a selection.\n";

    if(help_opt == "operators" || help_opt == "operator") {
      _stream << "Available operators for use within selection criteria:" << std::endl;
      ConfigIOParser::print_help(_stream, BaseDatumFormatter<Configuration>::Operator);
    }
    else if(help_opt == "property" || help_opt == "properties") {
      _stream << "Available property tags are currently:" << std::endl;
      ConfigIOParser::print_help(_stream);
    }
    _stream << std::endl;
  }

  // ///////////////////////////////////////
  // 'select' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int select_command(int argc, char *argv[]) {

    //casm enum [—supercell min max] [—config supercell ] [—hopconfigs hop.background]
    //- enumerate supercells and configs and hop local configurations

    std::string criteria, help_opt;
    std::vector<std::string> selection;
    fs::path out_path;
    COORD_TYPE coordtype;
    po::variables_map vm;
    bool force;


    /// Set command line options using boost program_options
    po::options_description desc("'casm select' usage");
    desc.add_options()
    ("help,h", po::value<std::string>(&help_opt)->implicit_value(""), "Write help documentation. Use '--help properties' for a list of selectable properties or '--help operators' for a list of selection operators")
    ("config,c", po::value<std::vector<std::string> >(&selection)->multitoken(),
     "One or more configuration files to operate on. If not given, or if given the keyword \"MASTER\" the master list is used.")
    ("output,o", po::value<fs::path>(&out_path), "Name for output file")
    ("json", "Write JSON output (otherwise CSV, unless output extension in .json or .JSON)")
    ("subset", "Write selected configurations in input list")
    ("union", "Write configurations selected in any of the input lists")
    ("intersection", "Write configurations selected in all of the input lists")
    ("on", po::value<std::string>(&criteria), "Add configurations to selection if they meet specified criteria.  Call using 'casm select -on [criteria]'")
    ("off", po::value<std::string>(&criteria), "Remove configurations from selection if they meet specified criteria.  Call using 'casm select -off [criteria]'")
    ("set", po::value<std::string>(&criteria), "Create a selection of Configurations that meet specified criteria.  Call using 'casm select -set [criteria]'")
    ("force,f", po::value(&force)->zero_tokens(), "Overwrite output file");

    try {
      po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

      if(!vm.count("help")) {
        if(vm.count("on") + vm.count("off") + vm.count("set") + vm.count("subset") + vm.count("union") + vm.count("intersection") != 1) {
          std::cout << desc << std::endl;
          std::cout << "Error in 'casm select'. Must use exactly one of --on, --off, --set, --subset, --union, or --intersection." << std::endl;
          return 1;
        }

        if(vm.count("subset") && !vm.count("output")) {
          std::cout << "Error in 'casm select --subset'. Please specify an --output file." << std::endl;
          return 1;
        }

        if(vm.count("union") && !vm.count("output")) {
          std::cout << "Error in 'casm select --union'. Please specify an --output file." << std::endl;
          return 1;
        }

        if(vm.count("intersection") && !vm.count("output")) {
          std::cout << "Error in 'casm select --intersection'. Please specify an --output file." << std::endl;
          return 1;
        }
      }

      /** Start --help option
       */
      if(vm.count("help")) {
        std::cout << std::endl << desc << std::endl;
      }

      po::notify(vm); // throws on error, so do after help in case of problems

      /** Finish --help option
       */
      if(vm.count("help")) {
        fs::path root = find_casmroot(fs::current_path());
        fs::path alias_file = root / ".casm/query_alias.json";
        if(fs::exists(alias_file)) {
          ConfigIOParser::load_aliases(alias_file);
        }
        select_help(std::cout, help_opt);
        return 0;
      }


      if((vm.count("on") || vm.count("off") || vm.count("set")) && vm.count("config") && selection.size() != 1) {
        std::string cmd = "--on";
        if(vm.count("off"))
          cmd = "--off";
        if(vm.count("set"))
          cmd = "--set";

        std::cout << "Error in 'casm select " << cmd << "'. " << selection.size() << " config selections were specified, but no more than one selection is allowed (MASTER list is used if none is specified)." << std::endl;
        return 1;
      }

      if(vm.count("subset") && vm.count("config") && selection.size() != 1) {
        std::cout << "Error in 'casm select --subset'. Zero (implies using Master list) or one lists should be input." << std::endl;
        return 1;
      }

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

    if(vm.count("on") || vm.count("off") || vm.count("set")) {
      bool select_switch = vm.count("on");
      bool only_selected = false;

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
          if(vm.count("set"))
            set_selection(primclex.config_begin(), primclex.config_end(), criteria);
          else
            set_selection(primclex.config_begin(), primclex.config_end(), criteria, select_switch);

          std::cout << "  DONE." << std::endl << std::endl;

          std::cout << "Writing config_list..." << std::endl;
          primclex.write_config_list();
          std::cout << "  DONE." << std::endl;
          return 0;
        }
        else {
          ConfigSelection<true> config_select(primclex);
          if(vm.count("set"))
            set_selection(config_select.config_begin(), config_select.config_end(), criteria);
          else
            set_selection(config_select.config_begin(), config_select.config_end(), criteria, select_switch);

          std::cout << "  DONE." << std::endl << std::endl;
          return write_selection(config_select, vm.count("force"), out_path, vm.count("json"), only_selected);
        }
      }
      else {
        ConfigSelection<true> config_select(primclex, selection[0]);

        if(vm.count("set"))
          set_selection(config_select.config_begin(), config_select.config_end(), criteria);
        else
          set_selection(config_select.config_begin(), config_select.config_end(), criteria, select_switch);

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


