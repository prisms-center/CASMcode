#include "query.hh"

#include <string>
#include <boost/algorithm/string.hpp>
#include "casm_functions.hh"
#include "casm/CASM_classes.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/clex/ConfigIOSelected.hh"

namespace CASM {
  void query_help(std::ostream &_stream, std::vector<std::string > help_opt_vec) {
    _stream << "Prints the properties for a set of configurations for the set of currently selected" << std::endl
            << "configurations or for a set of configurations specifed by a selection file." << std::endl
            << std::endl
            << "Property values are output in column-separated (default) or JSON format.  By default, " << std::endl
            << "entries for 'name' and 'selected' values are included in the output. " << std::endl
            << std::endl;
    for(std::string &help_opt : help_opt_vec){
      boost::trim(help_opt);
      if(help_opt == "operators" || help_opt == "operator") {
        _stream << "Available operators for use within queries:" << std::endl;
        ConfigIOParser::print_help(_stream, BaseDatumFormatter<Configuration>::Operator);
      }
      else if(help_opt == "property" || help_opt == "properties") {
        _stream << "Available property tags are currently:" << std::endl;
        ConfigIOParser::print_help(_stream);
      }
      _stream << std::endl;
    }
    _stream << std::endl;
  }

  int query_command(int argc, char *argv[]) {

    std::string new_alias;
    fs::path config_path, out_path;
    std::vector<std::string> columns,help_opt_vec;
    po::variables_map vm;
    bool json_flag(false), no_header(false), verbatim_flag(false);

    po::options_description desc("'casm query' usage");
    // Set command line options using boost program_options
    desc.add_options()
      ("help,h", po::value<std::vector<std::string> >(&help_opt_vec)->multitoken(), "Print general help. Use '--help properties' for a list of query-able properties or '--help operators' for a list of query operators")
    ("config,c", po::value<fs::path>(&config_path), "config_list files containing configurations for which to collect energies")
    ("columns,k", po::value<std::vector<std::string> >(&columns)->multitoken(), "List of values you want printed as columns")
    ("learn,l", po::value<std::string>(&new_alias), "Teach casm a new command that will persist within this project. Ex: 'casm query --learn is_Ni_dilute = lt(atom_frac(Ni),0.10001)'")
    ("json,j", po::value(&json_flag)->zero_tokens(), "Print in JSON format (CSV otherwise, unless output extension is .json or .JSON)")
    ("verbatim,v", po::value(&verbatim_flag)->zero_tokens(), "Print exact properties specified, without prepending 'name' and 'selected' entries")
    ("output,o", po::value<fs::path>(&out_path), "Name for output file")
    //("force,f", po::value(&force)->zero_tokens(), "Overrwrite output file")
    ("no-header,n", po::value(&no_header)->zero_tokens(), "Print without header (CSV only)");


    try {
      po::store(po::parse_command_line(argc, argv, desc),vm); // can throw
      
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
        query_help(std::cout, help_opt_vec);
        return 0;
      }


    }
    catch(po::error &e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      std::cerr << desc << std::endl;
      return 1;
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return 1;
    }

    if(!vm.count("learn")&&!vm.count("columns")){
      std::cout << std::endl << desc << std::endl;
    }

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cerr << "Error in 'casm query': No casm project found." << std::endl;
      return 1;
    }
    fs::current_path(root);

    fs::path alias_file = root / ".casm/query_alias.json";
    if(vm.count("learn")) {
      jsonParser mjson;
      if(fs::exists(alias_file)) {
        mjson.read(alias_file);
      }
      else {
        mjson.put_obj();
      }
      auto it = std::find(new_alias.cbegin(), new_alias.cend(), '=');
      std::string alias_name = boost::trim_copy(std::string(new_alias.cbegin(), it));
      std::string alias_command = boost::trim_copy(std::string(++it, new_alias.cend()));
      try {
        ConfigIOParser::add_custom_formatter(datum_formatter_alias<Configuration>(alias_name, alias_command));
      }
      catch(std::runtime_error &e) {
        std::cerr << "That's confusing. I can't learn malformed input:\n"
                  << e.what() << std::endl;
        return 1;
      }
      if(mjson.contains(alias_name)) {
        std::cerr << "WARNING: I already know '" << alias_name << "' as:\n"
                  << "             " << mjson[alias_name].get<std::string>() << "\n"
                  << "         I will forget it and learn '" << alias_name << "' as:\n"
                  << "             " << alias_command << std::endl;
      }
      mjson[alias_name] = alias_command;
      mjson.write(alias_file);
      std::cout << "  DONE" << std::endl;
      return 0;
    }
    if(!vm.count("columns")) {
      std::cerr << "ERROR: the option '--columns' is required but missing" << std::endl;
      return 1;
    }
    //else{ //option is "columns"
    if(fs::exists(alias_file)) {
      ConfigIOParser::load_aliases(alias_file);
    }

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    out_path = fs::absolute(out_path);

    std::cout << "Print:" << std::endl;
    for(int p = 0; p < columns.size(); p++) {
      std::cout << "   - " << columns[p] << std::endl;
    }

    if(vm.count("config"))
      std::cout << "to " << out_path << std::endl << std::endl;

    std::ofstream output_file;
    if(vm.count("output"))
      output_file.open(out_path.string().c_str());

    std::ostream &output_stream(vm.count("output") ? output_file : std::cout);
    output_stream << FormatFlag(output_stream).print_header(!no_header);
    DataFormatter<Configuration> formatter;
    ConstConfigSelection selection;

    /// Prepare for calculating correlations. Maybe this should get put into Clexulator.
    const DirectoryStructure &dir = primclex.dir();
    const ProjectSettings &set = primclex.settings();
    if(fs::exists(dir.clexulator_src(set.name(), set.bset()))) {
      primclex.read_global_orbitree(dir.clust(set.bset()));
      primclex.generate_full_nlist();
      primclex.generate_supercell_nlists();
    }

    try {
      if(vm.count("config"))
        selection = ConstConfigSelection(primclex, fs::absolute(config_path));
      else
        selection = ConstConfigSelection(primclex);

      auto it(columns.cbegin());
      std::vector<std::string> all_columns;
      if(!verbatim_flag) {
        formatter.push_back(ConfigIO::configname());
        formatter.push_back(datum_formatter_alias("selected", ConfigIO::selected_in(selection)));

        while(it != columns.cend() && ((*it) == "configname" || (*it) == "selected")) {
          ++it;
        }
      }
      all_columns.insert(all_columns.end(), it, columns.cend());
      formatter.append(ConfigIOParser::parse(all_columns));
    }
    catch(std::exception &e) {
      std::cerr << "Parsing error: " << e.what() << "\n\n";
      return 1;
    }

    try {
      // JSON output block
      if(json_flag || out_path.extension() == ".json" || out_path.extension() == ".JSON") {
        jsonParser json;

        //std::cout << "Read in config selection... it is:\n" << selection;
        json = formatter(selection.selected_config_begin(), selection.selected_config_end());

        output_stream << json;
      }
      // CSV output block
      else {
        //std::cout << "Read in config selection... it is:\n" << selection;
        output_stream << formatter(selection.selected_config_begin(), selection.selected_config_end());
      }
    }
    catch(std::exception &e) {
      std::cerr << "Initialization error: " << e.what() << "\n\n";
      return 1;
    }

    if(vm.count("output"))
      output_file.close();
    else {
      std::cerr << "\n   -Output printed to terminal, since no output file specified-\n";
    }

    std::cout << "  DONE." << std::endl << std::endl;

    return 0;
  };

}

