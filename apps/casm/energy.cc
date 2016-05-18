#include "energy.hh"

#include <string>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"

namespace CASM {

  int energy_command(int argc, char *argv[]) {

    fs::path abs_outpath, abs_outpath_hull, abs_outpath_hulljson;
    po::variables_map vm;
    bool force;


    try {

      // Set command line options using boost program_options
      po::options_description desc("'casm energy' usage");
      desc.add_options()
      ("help,h", "Print help message")
      ("force,f", po::value(&force)->zero_tokens(), "Overrwrite output file");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "Use this to print energies from your runs. Specify a list of        \n" <<
                    "configurations using 'casm select'. The following files are printed:\n" <<
                    "  $ROOT/properties/relaxed_energy/$CURR_CALCTYPE/$CURR_REF/energy   \n" <<
                    "  $ROOT/properties/relaxed_energy/$CURR_CALCTYPE/$CURR_REF/hull     \n" <<
                    "  $ROOT/properties/relaxed_energy/$CURR_CALCTYPE/$CURR_REF/hull.json\n\n";

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

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cout << "Error in 'casm energy': No casm project found." << std::endl;
      return 1;
    }

    // initialize primclex
    Log log(std::cout);
    PrimClex primclex(root, log);
    
    // Make sure directories exist:
    fs::path dir = root / "properties";
    if(!fs::exists(dir))
      fs::create_directory(dir);

    // relaxed_energy
    dir /= "relaxed_energy";
    if(!fs::exists(dir))
      fs::create_directory(dir);

    // calctype
    dir /= primclex.get_curr_calctype();
    if(!fs::exists(dir))
      fs::create_directory(dir);

    // ref
    dir /= primclex.get_curr_ref();
    if(!fs::exists(dir))
      fs::create_directory(dir);


    //Sorry, you no longer get to choose where to print stuff. Use casm awk for that
    abs_outpath = root / "properties" / "relaxed_energy" / primclex.get_curr_calctype() / primclex.get_curr_ref() / "energy";
    abs_outpath_hull = root / "properties" / "relaxed_energy" / primclex.get_curr_calctype() / primclex.get_curr_ref() / "hull";
    abs_outpath_hulljson = root / "properties" / "relaxed_energy" / primclex.get_curr_calctype() / primclex.get_curr_ref() / "hull.json";

    //If output files exists quit runnint
    if(fs::exists(abs_outpath) && !force) {
      std::cerr << "File " << abs_outpath.string() << " already exists. I'd hate to erase your data." << std::endl;
      return 2;
    }

    else if(fs::exists(abs_outpath_hull) && !force) {
      std::cerr << "File " << abs_outpath_hull.string() << " already exists. I'd hate to erase your data." << std::endl;
      return 2;
    }

    else if(fs::exists(abs_outpath_hulljson) && !force) {
      std::cerr << "File " << abs_outpath_hulljson.string() << " already exists. I'd hate to erase your data." << std::endl;
      return 2;
    }


    std::cout << "Make sure the hull gets generated and generate JSON object..." << std::endl << std::endl;
    jsonParser hulljson;
    ConfigSelection<false> config_select(primclex);
    hulljson = update_hull_props(primclex, config_select.begin(), config_select.end());
    std::cout << "  DONE." << std::endl << std::endl;

    std::cout << "Update Configuration files..." << std::endl << std::endl;
    primclex.write_config_list();
    std::cout << "  DONE." << std::endl << std::endl;

    std::string target = abs_outpath.string();
    std::cout << "Print energies to " << target << std::endl << std::endl;
    ConfigPrintStream cpstream(target);
    cpstream.add_printer("energy");
    cpstream.add_printer("dummy");
    cpstream.add_printer("parametric composition");
    cpstream.add_printer("distance from hull");
    cpstream.add_printer("name");
    cpstream.print_header();
    primclex.generic_print(cpstream);
    std::cout << "  DONE." << std::endl << std::endl;

    //You got the energy file, but you're not done yet, now we want a similar file that only has groundstates
    //Start by switching off everything that's not a groundstate
    //Under the current way the configuration selector works the provided config_list will be ignored here!!!
    target = abs_outpath_hull.string();
    std::cout << "Print groundstate energies to " << target << std::endl << std::endl;
    CASM::Array<std::string> justgroundstate;
    justgroundstate.push_back("off");
    justgroundstate.push_back("groundstate");
    justgroundstate.push_back("NOT");

    bool selection;
    for(auto it = primclex.config_begin(); it != primclex.config_end(); ++it) {
      it->set_selected(get_selection(justgroundstate, *it, it->selected()));
    }

    ConfigPrintStream cpstreamhull(target);
    cpstreamhull.add_printer("energy");
    cpstreamhull.add_printer("dummy");
    cpstreamhull.add_printer("parametric composition");
    cpstreamhull.add_printer("distance from hull");
    cpstreamhull.add_printer("name");
    cpstreamhull.print_header();
    primclex.generic_print(cpstreamhull);
    std::cout << "  DONE." << std::endl << std::endl;

    //Also write out the hull in json format
    target = abs_outpath_hulljson.string();
    std::cout << "Print groundstate energies in JSON format " << target << std::endl << std::endl;
    hulljson.write(abs_outpath_hulljson);
    std::cout << "  DONE." << std::endl << std::endl;

    return 0;
  };

}

