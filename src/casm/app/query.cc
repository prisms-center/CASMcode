#include <string>
#include <boost/algorithm/string.hpp>
#include "casm/app/casm_functions.hh"
#include "casm/CASM_classes.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIOSelected.hh"
#include "casm/completer/Complete.hh"


namespace CASM {

  void query_help(const DataFormatterDictionary<Configuration> &_dict, std::ostream &_stream, std::vector<std::string > help_opt_vec) {
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
        _dict.print_help(_stream, BaseDatumFormatter<Configuration>::Operator);
      }
      else if(help_opt[0] == 'p') {
        _stream << "Available property tags are currently:" << std::endl;
        _dict.print_help(_stream, BaseDatumFormatter<Configuration>::Property);
      }
      _stream << std::endl;
    }
    _stream << std::endl;
  }

  namespace Completer {
    ///Generate the option_description for the `query` command. You must pass all the required values to this function, and you must do it by reference.
    void add_query_options
    (po::options_description &desc,
     std::string &selection_str,
     fs::path &out_path,
     std::vector<std::string> &columns,
     std::vector<std::string> &help_opt_vec,
     std::vector<std::string> &new_alias,
     bool &json_flag,
     bool &no_header,
     bool &verbatim_flag,
     bool &gz_flag) {
      // Set command line options using boost program_options
      desc.add_options()
      ("help,h", po::value<std::vector<std::string> >(&help_opt_vec)->multitoken()->zero_tokens(), "Print general help. Use '--help properties' for a list of query-able properties or '--help operators' for a list of query operators")
      ("config,c", po::value<std::string>(&selection_str)->default_value("MASTER")->value_name(ArgHandler::path()), "config_list files containing configurations for which to collect energies")
      ("columns,k", po::value<std::vector<std::string> >(&columns)->multitoken()->zero_tokens()->value_name(ArgHandler::query()), "List of values you want printed as columns")
      ("json,j", po::value(&json_flag)->zero_tokens(), "Print in JSON format (CSV otherwise, unless output extension is .json/.JSON)")
      ("verbatim,v", po::value(&verbatim_flag)->zero_tokens(), "Print exact properties specified, without prepending 'name' and 'selected' entries")
      ("output,o", po::value<fs::path>(&out_path)->value_name(ArgHandler::path()), "Name for output file. Use STDOUT to print results without extra messages. CSV format unless extension is .json/.JSON, or --json option used.")
      ("gzip,z", po::value(&gz_flag)->zero_tokens(), "Write gzipped output file.")
      ("all,a", "Print results all configurations in input selection, whether or not they are selected.")
      ("no-header,n", po::value(&no_header)->zero_tokens(), "Print without header (CSV only)")
      ("alias", po::value<std::vector<std::string> >(&new_alias)->multitoken(),
       "Create an alias for a query that will persist within this project. "
       "Ex: 'casm query --alias is_Ni_dilute = lt(atom_frac(Ni),0.10001)'");

      return;
    }

  }

  int query_command(const CommandArgs &args) {

    std::string selection_str;
    fs::path out_path;
    std::vector<std::string> columns, help_opt_vec, new_alias;
    bool json_flag(false), no_header(false), verbatim_flag(false), gz_flag(false);

    po::options_description desc("'casm query' usage");
    // Set command line options using boost program_options
    Completer::add_query_options(desc, selection_str, out_path, columns, help_opt_vec, new_alias, json_flag, no_header, verbatim_flag, gz_flag);

    fs::path config_path;
    po::variables_map vm;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, desc), vm); // can throw

      /** Start --help option
       */
      if(vm.count("help")) {
        args.log << std::endl << desc << std::endl;
      }

      po::notify(vm); // throws on error, so do after help in case of problems

      /** Finish --help option
       */
      if(vm.count("help")) {
        if(args.root.empty()) {
          auto dict = make_dictionary<Configuration>();
          query_help(dict, std::cout, help_opt_vec);
        }
        else {
          ProjectSettings set(args.root);
          query_help(set.config_io(), args.log, help_opt_vec);
        }
        return 0;
      }


    }
    catch(po::error &e) {
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log << desc << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log << "Unhandled Exception reached the top of main: "
                   << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;
    }

    if(!vm.count("alias") && !vm.count("columns")) {
      args.log << std::endl << desc << std::endl;
    }

    // set current path to project root
    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    if(vm.count("alias")) {

      ProjectSettings set(root);

      // get user input
      std::string new_alias_str;
      for(auto const &substr : new_alias) {
        new_alias_str += substr;
      }

      // parse new_alias_str to create formatter
      auto it = std::find(new_alias_str.cbegin(), new_alias_str.cend(), '=');
      std::string alias_name = boost::trim_copy(std::string(new_alias_str.cbegin(), it));
      std::string alias_command = boost::trim_copy(std::string(++it, new_alias_str.cend()));

      try {
        set.add_alias(alias_name, alias_command, args.err_log);
        set.commit();
        return 0;
      }
      catch(std::runtime_error &e) {
        args.err_log << "Unable to learn alias\n"
                     << "   \"" << alias_name << " = " << alias_command << "\"\n"
                     << e.what() << std::endl;
        return ERR_UNKNOWN;
      }

    }
    if(!vm.count("columns")) {
      args.err_log << "ERROR: the option '--columns' is required but missing" << std::endl;
      return ERR_INVALID_ARG;
    }


    // -------------------------------------------------------------------------
    // perform query operation

    auto check_gz = [ = ](fs::path p) {
      if(p.extension() == ".gz" || p.extension() == ".GZ") {
        return true;
      }
      return false;
    };

    auto check_json = [ = ](fs::path p) {
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
    std::ostream &output_stream = make_ostream_if(vm.count("output"), args.log, uniq_fout, out_path, gz_flag);
    output_stream << FormatFlag(output_stream).print_header(!no_header);

    // If '_primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    // Get configuration selection
    ConstConfigSelection selection(primclex, selection_str);

    // set status_stream: where query settings and PrimClex initialization messages are sent
    Log &status_log = (out_path.string() == "STDOUT") ? args.err_log : args.log;

    // Print info
    status_log << "Print:" << std::endl;
    for(int p = 0; p < columns.size(); p++) {
      status_log << "   - " << columns[p] << std::endl;
    }
    if(vm.count("output"))
      status_log << "to " << fs::absolute(out_path) << std::endl;
    status_log << std::endl;

    // Construct DataFormatter
    primclex.settings().set_selected(selection);
    DataFormatter<Configuration> formatter;

    try {

      std::vector<std::string> all_columns;
      if(!verbatim_flag) {
        all_columns.push_back("configname");
        all_columns.push_back("selected");
      }
      all_columns.insert(all_columns.end(), columns.cbegin(), columns.cend());

      formatter.append(primclex.settings().config_io().parse(all_columns));

    }
    catch(std::exception &e) {
      args.err_log << "Parsing error: " << e.what() << "\n\n";
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
      args.err_log << "Initialization error: " << e.what() << "\n\n";
      return ERR_UNKNOWN;
    }

    if(!uniq_fout) {
      status_log << "\n   -Output printed to terminal, since no output file specified-\n";
    }

    status_log << "  DONE." << std::endl << std::endl;

    return 0;
  };
}


