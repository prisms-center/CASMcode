#include "query.hh"

#include <string>
#include <boost/algorithm/string.hpp>
#include "casm_functions.hh"
#include "casm/CASM_classes.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIOSelected.hh"

namespace CASM {
  
  void query_help(std::ostream &_stream, std::vector<std::string > help_opt_vec) {
    _stream << "Prints the properties for a set of configurations for the set of currently selected" << std::endl
            << "configurations or for a set of configurations specifed by a selection file." << std::endl
            << std::endl
            << "Property values are output in column-separated (default) or JSON format.  By default, " << std::endl
            << "entries for 'name' and 'selected' values are included in the output. " << std::endl
            << std::endl;

    for(const std::string &help_opt : help_opt_vec) {
      if(help_opt.empty())
        continue;

      if(help_opt[0] == 'o') {
        _stream << "Available operators for use within queries:" << std::endl;
        ConfigIOParser::print_help(_stream, BaseDatumFormatter<Configuration>::Operator);
      }
      else if(help_opt[0] == 'p') {
        _stream << "Available property tags are currently:" << std::endl;
        ConfigIOParser::print_help(_stream);
      }
      _stream << std::endl;
    }
    _stream << std::endl;
  }

  int query_command(int argc, char *argv[], PrimClex* _primclex, std::ostream& sout, std::ostream& serr) {

    std::string new_alias, selection_str;
    fs::path config_path, out_path;
    std::vector<std::string> columns, help_opt_vec;
    po::variables_map vm;
    bool json_flag(false), no_header(false), verbatim_flag(false), gz_flag(false);

    po::options_description desc("'casm query' usage");
    // Set command line options using boost program_options
    desc.add_options()
    ("help,h", po::value<std::vector<std::string> >(&help_opt_vec)->multitoken()->zero_tokens(), "Print general help. Use '--help properties' for a list of query-able properties or '--help operators' for a list of query operators")
    ("config,c", po::value<std::string>(&selection_str)->default_value("MASTER"), "config_list files containing configurations for which to collect energies")
    ("columns,k", po::value<std::vector<std::string> >(&columns)->multitoken()->zero_tokens(), "List of values you want printed as columns")
    ("learn,l", po::value<std::string>(&new_alias), "Teach casm a new command that will persist within this project. Ex: 'casm query --learn is_Ni_dilute = lt(atom_frac(Ni),0.10001)'")
    ("json,j", po::value(&json_flag)->zero_tokens(), "Print in JSON format (CSV otherwise, unless output extension is .json/.JSON)")
    ("verbatim,v", po::value(&verbatim_flag)->zero_tokens(), "Print exact properties specified, without prepending 'name' and 'selected' entries")
    ("output,o", po::value<fs::path>(&out_path), "Name for output file. Use STDOUT to print results without extra messages. CSV format unless extension is .json/.JSON, or --json option used.")
    ("gzip,z", po::value(&gz_flag)->zero_tokens(), "Write gzipped output file.")
    ("all,a", "Print results all configurations in input selection, whether or not they are selected.")
    ("no-header,n", po::value(&no_header)->zero_tokens(), "Print without header (CSV only)");


    try {
      po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

      /** Start --help option
       */
      if(vm.count("help")) {
        sout << std::endl << desc << std::endl;
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
        ConfigIOParser::add_custom_formatter(
          datum_formatter_alias(
            "selected", 
            ConfigIO::selected_in(), 
            "Returns true if configuration is specified in the input selection"
          )
        );
        query_help(sout, help_opt_vec);
        return 0;
      }


    }
    catch(po::error &e) {
      serr << "ERROR: " << e.what() << std::endl << std::endl;
      serr << desc << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      serr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;
    }

    if(!vm.count("learn") && !vm.count("columns")) {
      sout << std::endl << desc << std::endl;
    }
    
    // set current path to project root
    fs::path root;
    if(!_primclex) {
      root = find_casmroot(fs::current_path());
      if(root.empty()) {
        serr << "Error in 'casm query': No casm project found." << std::endl;
        return ERR_NO_PROJ;
      }
    }
    else {
      root = _primclex->get_path();
    }
    
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
        serr << "That's confusing. I can't learn malformed input:\n"
                  << e.what() << std::endl;
        return 1;
      }
      if(mjson.contains(alias_name)) {
        serr << "WARNING: I already know '" << alias_name << "' as:\n"
                  << "             " << mjson[alias_name].get<std::string>() << "\n"
                  << "         I will forget it and learn '" << alias_name << "' as:\n"
                  << "             " << alias_command << std::endl;
      }
      mjson[alias_name] = alias_command;
      mjson.write(alias_file);
      sout << "  DONE" << std::endl;
      return 0;
    }
    if(!vm.count("columns")) {
      serr << "ERROR: the option '--columns' is required but missing" << std::endl;
      return ERR_INVALID_ARG;
    }
    
    
    // -------------------------------------------------------------------------
    // perform query operation
    
    if(fs::exists(alias_file)) {
      ConfigIOParser::load_aliases(alias_file);
    }
    
    
    auto check_gz = [=](fs::path p) {
      if(p.extension() == ".gz" || p.extension() == ".GZ") {
        return true;
      }
      return false;
    };
    
    auto check_json = [=](fs::path p) {
      if(p.extension() == ".json" || p.extension() == ".JSON") {
        return true;
      }
      return false;
    };
    
    
    // Checks for: X.json.gz / X.json / X.gz  (also accepts .JSON or .GZ)
    if(check_gz(out_path)) {
      gz_flag = true;
      json_flag = check_json(out_path.stem()) || json_flag;
    }
    else {
      json_flag = check_json(out_path) || json_flag;
    }
    
    // set output_stream: where the query results are written
    std::unique_ptr<std::ostream> uniq_fout;
    std::ostream& output_stream = make_ostream_if(vm.count("output"), sout, uniq_fout, out_path, gz_flag);
    output_stream << FormatFlag(output_stream).print_header(!no_header);
    
    // set status_stream: where query settings and PrimClex initialization messages are sent
    std::ostream &status_stream = (out_path.string() == "STDOUT") ? serr : sout;
    
    // If '_primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex& primclex = make_primclex_if_not(_primclex, uniq_primclex, root, status_stream);
    
    // Get configuration selection
    ConstConfigSelection selection(primclex, selection_str);
    
    // Print info
    status_stream << "Print:" << std::endl;
    for(int p = 0; p < columns.size(); p++) {
      status_stream << "   - " << columns[p] << std::endl;
    }
    if(vm.count("output"))
      status_stream << "to " << fs::absolute(out_path) << std::endl;
    status_stream << std::endl;
    
    // Construct DataFormatter
    DataFormatter<Configuration> formatter;
    ConfigIOParser::add_custom_formatter(
      datum_formatter_alias(
        "selected", 
        ConfigIO::selected_in(selection), 
        "Returns true if configuration is specified in the input selection"
      )
    );
        
    try {
      
      std::vector<std::string> all_columns;
      if(!verbatim_flag) {
        all_columns.push_back("configname");
        all_columns.push_back("selected");
      }
      all_columns.insert(all_columns.end(), columns.cbegin(), columns.cend());
      
      formatter.append(ConfigIOParser::parse(all_columns));
      
    }
    catch(std::exception &e) {
      serr << "Parsing error: " << e.what() << "\n\n";
      return ERR_INVALID_ARG;
    }

    try {
      
      auto begin = vm.count("all") ? selection.config_begin() : selection.selected_config_begin();
      auto end = vm.count("all") ? selection.config_end() : selection.selected_config_end();
      
      // JSON output block
      if(json_flag) {
        jsonParser json;

        //sout << "Read in config selection... it is:\n" << selection;
        json = formatter(begin, end);

        output_stream << json;
      }
      // CSV output block
      else {
        //sout << "Read in config selection... it is:\n" << selection;
        output_stream << formatter(begin, end);
      }

    }
    catch(std::exception &e) {
      serr << "Initialization error: " << e.what() << "\n\n";
      return ERR_UNKNOWN;
    }

    if(!uniq_fout) {
      status_stream << "\n   -Output printed to terminal, since no output file specified-\n";
    }

    status_stream << "  DONE." << std::endl << std::endl;

    return 0;
  };
}


