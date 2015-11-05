#include "select.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/CASM_global_definitions.hh"

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

  template<bool IsConst>
  bool write_selection(const ConfigSelection<IsConst> &config_select, bool force, const fs::path &out_path, bool write_json, bool only_selected) {
    if(fs::exists(out_path) && !force) {
      std::cerr << "File " << out_path << " already exists. Use --force to force overwrite." << std::endl;
      return 1;
    }

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
    return 0;
  }

  void select_help(std::ostream &_stream, std::vector<std::string> help_opt) {
    _stream << "DESCRIPTION" << std::endl
            << "\n"
            << "    Use '[--set | --set-on | --set-off] [criteria]' for specifying or editing a selection.\n";

    for(const std::string &str : help_opt){
      if(str == "operators" || str == "operator") {
        _stream << "Available operators for use within selection criteria:" << std::endl;
        ConfigIOParser::print_help(_stream, BaseDatumFormatter<Configuration>::Operator);
      }
      else if(str == "property" || str == "properties") {
        _stream << "Available property tags are currently:" << std::endl;
        ConfigIOParser::print_help(_stream);
      }
      _stream << std::endl;
    }
    _stream << std::endl;
  }

  // ///////////////////////////////////////
  // 'select' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int select_command(int argc, char *argv[]) {

    //casm enum [—supercell min max] [—config supercell ] [—hopconfigs hop.background]
    //- enumerate supercells and configs and hop local configurations

    std::vector<std::string> criteria_vec, help_opt_vec;
    std::vector<std::string> selection;

    fs::path out_path;
    COORD_TYPE coordtype;
    po::variables_map vm;

    /// Set command line options using boost program_options
    // NOTE: multitoken() is used instead of implicit_value() because implicit_value() is broken on some systems -- i.e., braid.cnsi.ucsb.edu
    //       (not sure if it's an issue with a particular shell, or boost version, or something else)
    po::options_description desc("'casm select' usage");
    desc.add_options()
      ("help,h", po::value<std::vector<std::string> >(&help_opt_vec)->multitoken()->zero_tokens(), "Write help documentation. Use '--help properties' for a list of selectable properties or '--help operators' for a list of selection operators")
      ("config,c", po::value<std::vector<std::string> >(&selection)->multitoken(),
       "One or more configuration files to operate on. If not given, or if given the keyword \"MASTER\" the master list is used.")
      ("output,o", po::value<fs::path>(&out_path), "Name for output file")
      ("json", "Write JSON output (otherwise CSV, unless output extension is '.json' or '.JSON')")
      ("subset", "Only write selected configurations to output. Can be used by itself or in conjunction with other options")
      ("xor", "Performs logical XOR on two configuration selections")
      ("not", "Performs logical NOT on configuration selection")
      ("or", "Write configurations selected in any of the input lists. Equivalent to logical OR")
      ("and", "Write configurations selected in all of the input lists. Equivalent to logical AND")
      ("set-on", po::value<std::vector<std::string> >(&criteria_vec)->multitoken()->zero_tokens(), "Add configurations to selection if they meet specified criteria.  Call using 'casm select --set-on [\"criteria\"]'")
      ("set-off", po::value<std::vector<std::string> >(&criteria_vec)->multitoken()->zero_tokens(), "Remove configurations from selection if they meet specified criteria.  Call using 'casm select --set-off [\"criteria\"]'")
      ("set", po::value<std::vector<std::string> >(&criteria_vec)->multitoken(), "Create a selection of Configurations that meet specified criteria.  Call using 'casm select --set [\"criteria\"]'")
      ("force,f", "Overwrite output file");

    std::string cmd;
    std::vector<std::string> allowed_cmd={"and","or","xor","not","set-on","set-off","set"};
    try {
      po::store(po::parse_command_line(argc, argv, desc), vm); // can throw
      Index num_cmd(0);
      for(const std::string &cmd_str : allowed_cmd){
        if(vm.count(cmd_str)){
          num_cmd++;
          cmd=cmd_str;
        }
      }

      if(!vm.count("help")) {
        if(num_cmd>1){
          std::cout << desc << std::endl;
          std::cout << "Error in 'casm select'. Must use exactly one of --set-on, --set-off, --set, --and, --or, --xor, or --not." << std::endl;
          return 1;
        }
        else if(vm.count("subset") && vm.count("config") && selection.size() != 1) {
          std::cout << "ERROR: 'casm select --subset' expects zero or one list as argument." << std::endl;
          return 1;
        }

        

        if(!vm.count("output") && (cmd=="or" || cmd=="and" || cmd=="xor" || cmd=="not")){
          std::cout << "ERROR: 'casm select --" << cmd << "' expects an --output file." << std::endl;
          return 1;
        }

      }

      // Start --help option
      if(vm.count("help")) {
        std::cout << std::endl << desc << std::endl;
      }

      po::notify(vm); // throws on error, so do after help in case of problems

      // Finish --help option
      if(vm.count("help")) {
        fs::path root = find_casmroot(fs::current_path());
        fs::path alias_file = root / ".casm/query_alias.json";
        if(fs::exists(alias_file)) {
          ConfigIOParser::load_aliases(alias_file);
        }
        select_help(std::cout, help_opt_vec);
        return 0;
      }

      if((vm.count("set-on") || vm.count("set-off") || vm.count("set")) && vm.count("config") && selection.size() != 1) {
        std::string cmd = "--set-on";
        if(vm.count("set-off"))
          cmd = "--set-off";
        if(vm.count("set"))
          cmd = "--set";

        std::cout << "Error in 'casm select " << cmd << "'. " << selection.size() << " config selections were specified, but no more than one selection is allowed (MASTER list is used if no other is specified)." << std::endl;
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


    if(vm.count("output") && out_path!="MASTER"){

      out_path = fs::absolute(out_path);
      //check now so we can exit early with an obvious error
      if(fs::exists(out_path) && !vm.count("force")) {
        std::cerr << desc << std::endl << std::endl;
        std::cerr << "ERROR: File " << out_path << " already exists. Use --force to force overwrite." << std::endl;
        return 1;
      }
    }

    bool only_selected(false);
    for(int i = 0; i < selection.size(); i++) {
      if(selection[i] != "MASTER") {
        selection[i] = fs::absolute(fs::path(selection[i])).string();
      }
    }
    if(selection.empty()) {
      only_selected=true;
      selection.push_back("MASTER");
    }
    
    // switch to root directory
    fs::path orig_path = fs::current_path();
    fs::path root = find_casmroot(orig_path);
    if(root.empty()) {
      std::cerr << "Error: No casm project found." << std::endl;
      return 1;
    }
    fs::current_path(root);


    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    // load initial selection into config_select -- this is also the selection that will be printed at end
    ConfigSelection<false> config_select(primclex,selection[0]);


    if(vm.count("set-on") || vm.count("set-off") || vm.count("set")) {
      bool select_switch = vm.count("set-on");
      std::string criteria;
      if(criteria_vec.size()==1){
        criteria=criteria_vec[0];
      }
      else if(criteria_vec.size()>1){
        std::cerr << "ERROR: Selection criteria must be a single string.  You provided " << criteria_vec.size() << " strings:\n";
        for(const std::string &str : criteria_vec)
          std::cerr << "     - " << str << "\n";
        return 1;
      }
      std::cout << "Set selection: " << criteria << std::endl << std::endl;
      
      /// Prepare for calculating correlations. Maybe this should get put into Clexulator.
      const DirectoryStructure &dir = primclex.dir();
      const ProjectSettings &set = primclex.settings();
      if(fs::exists(dir.clexulator_src(set.name(), set.bset()))) {
        primclex.read_global_orbitree(dir.clust(set.bset()));
        primclex.generate_full_nlist();
        primclex.generate_supercell_nlists();

      }

      
      if(vm.count("set"))
        set_selection(config_select.config_begin(), config_select.config_end(), criteria);
      else
        set_selection(config_select.config_begin(), config_select.config_end(), criteria, select_switch);
      
      std::cout << "  DONE." << std::endl << std::endl;
      
      //return write_selection(config_select, vm.count("force"), out_path, vm.count("json"), only_selected);
      
    }

    if(vm.count("subset")) {
      only_selected = true;
    }

    if(vm.count("not")) {
      if(selection.size()!=1){
        std::cerr << "ERROR: Option --not requires exactly 1 selection as argument\n";
        return 1;
      }
      // loop through other lists, keeping only configurations selected in the other lists
      auto it = config_select.config_begin();
      for(; it != config_select.config_end(); ++it) {
        it.set_selected(!it.selected());
      }

    }

    if(vm.count("or")) {
      
      // loop through other lists, inserting all selected configurations
      for(int i = 1; i < selection.size(); i++) {
        ConstConfigSelection tselect(primclex, selection[i]);
        for(auto it = tselect.selected_config_cbegin(); it != tselect.selected_config_cend(); ++it) {
          config_select.set_selected(it.name(), true);
        }
      }

      only_selected = true;
    }

    if(vm.count("and")) {
      // loop through other lists, keeping only configurations selected in the other lists
      for(int i = 1; i < selection.size(); i++) {
        ConstConfigSelection tselect(primclex, selection[i]);

        auto it = config_select.selected_config_begin();
        for(; it != config_select.selected_config_end(); ++it) {
          it.set_selected(tselect.selected(it.name()));
        }
      }

      only_selected = true;
      
    }

    if(vm.count("xor")) {
      if(selection.size()!=2){
        std::cerr << "ERROR: Option --xor requires exactly 2 selections as argument\n";
        return 1;
      }

      ConstConfigSelection tselect(primclex, selection[1]);
      for(auto it = tselect.selected_config_begin(); it != tselect.selected_config_end(); ++it) {
        //If selected in both lists, deselect it
        if(config_select.selected(it.name())){
          config_select.set_selected(it.name(), false);
        }
        else // If only selected in tselect, add it to config_select
          config_select.set_selected(it.name(), true);

      }
      only_selected = true;
    }

    /// Only write selection to disk past this point
    if(!vm.count("output") || out_path=="MASTER") {
      auto pc_it=primclex.config_begin(), pc_end=primclex.config_end();
      for(; pc_it!=pc_end; ++pc_it){
        pc_it->set_selected(false);
      }

      auto it=config_select.selected_config_begin(), it_end=config_select.selected_config_end();
      for(; it!=it_end; ++it){
        it->set_selected(true);
      }
      
      std::cout << "Writing config_list..." << std::endl;
      primclex.write_config_list();
      std::cout << "  DONE." << std::endl;

      std::cout << "\n***************************\n" << std::endl;
      std::cout << std::endl;

      return 0;    
    }
    else{
      
      std::cout << "Writing selection to " << out_path << std::endl;
      int ret_code=write_selection(config_select, vm.count("force"), out_path, vm.count("json"), only_selected);
      std::cout << "  DONE." << std::endl;
      std::cout << "\n***************************\n" << std::endl;
      
      std::cout << std::endl;
      return ret_code;
    }

    return 0;
  };

}


