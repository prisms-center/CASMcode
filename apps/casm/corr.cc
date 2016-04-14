#include "corr.hh"

#include <string>

#include <casm/core>

#include "casm_functions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"

namespace CASM {

  int corr_command(int argc, char *argv[]) {

    fs::path abs_outpath;
    std::vector<fs::path> config_path, abs_config_path;
    po::variables_map vm;
    std::string outfile, cspecsfile;
    bool force;

    try {

      // Set command line options using boost program_options
      po::options_description desc("'casm corr' usage");
      desc.add_options()
      ("help,h", "Print help message")
      ("config,c", po::value<std::vector<fs::path> >(&config_path)->multitoken()->required(), "List of config_list files containing configurations for which to calculate correlations")
      ("output,o", po::value<std::string>(&outfile), "Name for output file")
      ("force,f", po::value(&force)->zero_tokens(), "Overwrite output file");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;
          std::cout << "Use this to make an old CASM style corr.in file to fit your cluster expansion. It will read in the FCLUST.json "
                    << "file in your $CASMROOT/clusters/clex.$CLEXTYPE directory. Pick the type of basis functions you want and an "
                    << "output file. Getting correlations has never been so easy!" << std::endl;
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
      std::cout << "Error in 'casm corr': No casm project found." << std::endl;
      return 1;
    }

    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    DirectoryStructure dir(root);
    ProjectSettings set(root);

    if(vm.count("output")) {
      abs_outpath = outfile;
      abs_outpath = fs::absolute(abs_outpath);
    }
    else {
      abs_outpath = dir.corr_in(set.bset());
    }

    // if(vm.count("cspecs")) {
    //   abs_cspecs_path = cspecsfile;
    //   abs_cspecs_path = fs::absolute(abs_cspecs_path);
    // }
    // else {
    //   abs_cspecs_path = root / "clusters" / primclex.get_curr_clex() / "CSPECS";
    // }

    //If output file exists quit runnint
    if(fs::exists(abs_outpath) && !force) {
      std::cerr << "File " << outfile << " already exists. I'd hate to erase your data." << std::endl;
      return 2;
    }

    //Test if FCLUST exists
    if(!fs::exists(dir.clust(set.bset()))) {
      std::cerr << "There is no clust.json in " << dir.clust(set.bset()) << std::endl;
      std::cerr << "Try re-running casm corr after you run casm clusters --enumerate" << std::endl;
      return 3;
    }

    // want absolute paths
    for(int i = 0; i < config_path.size(); i++)
      abs_config_path.push_back(fs::absolute(config_path[i]));

    jsonParser json_bspecs(dir.bspecs(set.bset()));
    std::string basis_type_str;
    basis_type_str = json_bspecs["orbitree_specs"]["basis_functions"].get<std::string>();

    //Fill up basis functions in primitive structure and orbitree
    char basis_type = basis_type_str[0];
    std::cout << "Fill basis tables and read in global orbitree." << std::endl << std::endl;
    primclex.populate_basis_tables(basis_type);
    //Reading in the global_orbitree from $CASMROOT/clusters/$CURRCLEX
    primclex.read_global_orbitree(dir.clust(set.bset()));
    // primclex.generate_global_orbitree();
    primclex.populate_cluster_basis_function_tables();
    primclex.generate_full_nlist();
    primclex.generate_supercell_nlists();
    std::cout << "  global orbitree basis set size: " << primclex.get_global_orbitree().basis_set_size() << std::endl;
    std::cout << "  DONE." << std::endl << std::endl;

    //Make the Clexulator
    //std::cout << "Print '" << primclex.get_curr_clexulator() << "' in " << dir.clexulator_src(set.name(), set.bset()) << std::endl;
    //primclex.print_global_clexulator();
    //std::cout << "  DONE." << std::endl << std::endl;


    Clexulator clexulator(set.global_clexulator(), dir.clexulator_dir(set.bset()), set.compile_options(), set.so_options());

    //Make the correlations
    std::cout << "Calculate global scalar correlations using the orbitree in  " << dir.FCLUST(set.bset()) << std::endl << std::endl;
    //std::ofstream outstream(outfile.c_str());
    //primclex.generate_global_scalar_correlations();
    PrimClex::config_iterator it = primclex.config_begin();
    for(; it != primclex.config_end(); ++it) {
      if(it->selected()) {
        it->set_correlations(clexulator);
      }
      // for(int i = 0; i < scel_index.size(); i++) {
      //   primclex.populate_global_correlations(scel_index[i], config_index[i]);
    }
    std::cout << "  DONE." << std::endl << std::endl;

    std::cout << "Update Configuration files..." << std::endl << std::endl;
    primclex.write_config_list();
    std::cout << "  DONE." << std::endl << std::endl;

    std::cout << "Count number of configurations and correlations for header..." << std::endl << std::endl;
    int num_configs = primclex.amount_selected();

    //Instead of just grabbing the first one, grab the first one that is turned on
    int num_corrs = -1;
    for(Index i = 0; i < primclex.get_supercell_list().size(); i++) {
      for(Index j = 0; j < primclex.get_supercell_list()[i].get_config_list().size(); j++) {
        if(primclex.get_supercell_list()[i].get_config_list()[j].selected()) {
          num_corrs = primclex.get_supercell(i).get_config(j).get_correlations().get_unrolled_correlations().size();
          break;
        }
      }
      if(num_corrs > 0) {
        break;
      }
    }

    if(num_corrs < 0) {
      std::cerr << "ERROR in casm corr" << std::endl;
      std::cerr << "No configurations selected!" << std::endl;
      std::cerr << "Exiting without writing anything" << std::endl;
      return 4;
    }

    else if(num_corrs == 0) {
      std::cerr << "ERROR in casm corr" << std::endl;
      std::cerr << "Shenanigans in your CSPECS files. You seem to have zero clusters" << std::endl;
      return 5;
    }

    std::stringstream numconfstream;
    numconfstream << num_configs;
    std::stringstream numcorrstream;
    numcorrstream << num_corrs;

    std::string file_header = numcorrstream.str() + " #clusters\n" + numconfstream.str() + " #configurations\nclusters";
    std::cout << "  DONE." << std::endl << std::endl;

    std::cout << "Print correlations to " << abs_outpath << std::endl << std::endl;
    ConfigPrintStream cpstream(abs_outpath);
    cpstream.print_header(file_header);
    cpstream.add_printer("correlations");
    primclex.generic_print(cpstream);
    std::cout << "  DONE." << std::endl << std::endl;

    return 0;
  };

}

