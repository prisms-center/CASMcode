#include "query.hh"

#include <string>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"
#include "casm/app/ProjectSettings.hh"

namespace CASM {
  void query_help(std::ostream &_stream) {
    _stream << "Prints the properties for a set of configurations for the set of currently selected" << std::endl
            << "configurations or for a set of configurations specifed by a selection file." << std::endl
            << std::endl
            << "Property values are output in column-separated (default) or JSON format.  By default, " << std::endl
            << "entries for 'name' and 'selected' values are included in the output. " << std::endl
            << std::endl
            << "Available property tags are currently:" << std::endl;
    ConfigIOParser::print_help(_stream);
    _stream << std::endl;
  }

  int query_command(int argc, char *argv[]) {

    fs::path config_path, out_path;
    std::vector<std::string> columns;
    po::variables_map vm;
    bool json_flag(false), no_header(false), verbatim_flag(false);

    po::options_description desc("'casm query' usage");
    // Set command line options using boost program_options
    desc.add_options()
    ("help,h", "Print help message")
    ("config,c", po::value<fs::path>(&config_path), "config_list files containing configurations for which to collect energies")
    ("columns,k", po::value<std::vector<std::string> >(&columns)->multitoken()->required(), "List of values you want printed as columns")
    ("json,j", po::value(&json_flag)->zero_tokens(), "Print in JSON format (CSV otherwise, unless output extension is .json or .JSON)")
    ("verbatim,v", po::value(&verbatim_flag)->zero_tokens(), "Print exact properties specified, without prepending 'name' and 'selected' entries")
    ("output,o", po::value<fs::path>(&out_path), "Name for output file")
    //("force,f", po::value(&force)->zero_tokens(), "Overrwrite output file")
    ("no-header,n", po::value(&no_header)->zero_tokens(), "Print without header (CSV only)");


    try {

      po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        std::cout << std::endl << desc << std::endl;
        query_help(std::cout);
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
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return 1;
    }

    //If output file exists quit runnint
    /*if(fs::exists(out_path) && !force) {
      std::cerr << "File " << out_path << " already exists. I'd hate to erase your data." << std::endl;
      return 2;
    }
    */

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cerr << "Error in 'casm query': No casm project found." << std::endl;
      return 1;
    }
    fs::current_path(root);



    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    out_path = fs::absolute(out_path);

    std::cout << "Print:" << std::endl;
    for(int p = 0; p < columns.size(); p++) {
      std::cout << "-" << columns[p] << std::endl;
    }

    if(vm.count("config"))
      std::cout << "to " << out_path << std::endl << std::endl;

    std::ofstream output_file;
    if(vm.count("output"))
      output_file.open(out_path.string().c_str());

    const DirectoryStructure &dir = primclex.dir();
    ProjectSettings &set = primclex.settings();

    /// Prepare for calculating correlations. Maybe this should get put into Clexulator.
    if(fs::exists(dir.clexulator_src(set.name(), set.bset()))) {
      primclex.read_global_orbitree(dir.clust(set.bset()));
      primclex.generate_full_nlist();
      primclex.generate_supercell_nlists();

    }


    std::ostream &output_stream(vm.count("output") ? output_file : std::cout);
    output_stream << FormatFlag(output_stream).print_header(!no_header);

    auto it(columns.cbegin());
    std::vector<std::string> all_columns;
    if(!verbatim_flag) {
      all_columns = std::vector<std::string>({"configname", "selected"});

      while(it != columns.cend() && ((*it) == "configname" || (*it) == "selected")) {
        ++it;
      }
    }
    all_columns.insert(all_columns.end(), it, columns.cend());

    // JSON output block
    try {
      if(json_flag || out_path.extension() == ".json" || out_path.extension() == ".JSON") {
        jsonParser json;
        if(vm.count("config")) {
          ConstConfigSelection selection(primclex, fs::absolute(config_path));
          //std::cout << "Read in config selection... it is:\n" << selection;

          json = ConfigIOParser::parse(all_columns)(selection.selected_config_begin(), selection.selected_config_end());
        }
        else {
          json =  ConfigIOParser::parse(all_columns)(primclex.selected_config_begin(), primclex.selected_config_end());
        }
        output_stream << json;
      }
      // CSV output block
      else {
        if(vm.count("config")) {
          ConstConfigSelection selection(primclex, fs::absolute(config_path));
          //std::cout << "Read in config selection... it is:\n" << selection;
          output_stream << ConfigIOParser::parse(all_columns)(selection.selected_config_begin(), selection.selected_config_end());
        }
        else {
          output_stream << ConfigIOParser::parse(all_columns)(primclex.selected_config_begin(), primclex.selected_config_end());
        }
      }
    }
    catch(std::exception &e) {
      std::cerr << "Parsing error: " << e.what() << "\n\n";
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

