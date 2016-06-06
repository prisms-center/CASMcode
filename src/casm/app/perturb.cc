#include "casm/app/casm_functions.hh"
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/clex/ConfigEnumIterator.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'perturb' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int perturb_command(const CommandArgs &args) {

    throw std::runtime_error("Error: 'casm perturb' is being re-implemented");

    /*

    double tol = CASM::TOL;
    bool is_trans = false;
    fs::path cspecs_path, abs_cspecs_path;
    fs::path selection;
    COORD_TYPE coordtype = CASM::FRAC;
    po::variables_map vm;
    std::vector<Index> subgrids;
    std::vector<double> mags;
    double inc_value;
    std::string strain_mode;
    Index poly_order;

    /// Set command line options using boost program_options
    po::options_description desc("'casm perturb' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    ("occ", "Perturb occupations")
    ("cspecs", po::value<fs::path>(&cspecs_path), "Cluster specifications file defining perturbation")
    ("config,c", po::value<fs::path>(&selection),
     "Selected configurations are used reference for generating perturbations. If not specified, or 'MASTER' given, uses master list selection.");
    //("strain,s", "Generate strain perturbations")
    //("gridsize", po::value<std::vector<Index> > (&subgrids)->multitoken(), "Size of grid for each subspace")
    //("mag", po::value<std::vector<double> > (&mags)->multitoken(), "Magnitude of grid")
    //("poly-order", po::value<Index> (&poly_order)->default_value(6), "Max order of strain polynomials")
    //("strainmode", po::value<std::string> (&strain_mode)->default_value("GL"), "Strain mode name");


    try {
      po::store(po::parse_command_line(args.argc, args.argv, desc), vm); // can throw

      // --help option
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << desc << std::endl;

        std::cout << "DESCRIPTION" << std::endl;
        std::cout << "    Generate supercells that are perturbations of a reference\n";
        std::cout << "    configuration.                                           \n";
        std::cout << "    - using the --cspecs option, a bspecs.json type file is  \n";
        std::cout << "      required to determine the extent of the perturbations. \n";
        std::cout << "      Currently only 'orbit_branch_specs' are supported.     \n";
        std::cout << "    - perturbations are generated about selected reference   \n";
        std::cout << "      configurations                                         \n";
        std::cout << std::endl;

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

    COORD_MODE C(coordtype);

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

    DirectoryStructure dir(root);
    ProjectSettings set(root);

    if(vm.count("strain")) {
      StrainConverter sconvert(strain_mode);
      sconvert.set_symmetrized_sop(primclex.prim().point_group());

      Eigen::MatrixXd axes = sconvert.sop_transf_mat();
      std::vector<Index> mult;
      std::vector<Eigen::MatrixXd> wedges = sconvert.irreducible_wedges(primclex.prim().point_group(), mult);

      Index num_sub = wedges.size();

      BasisSet strain_vars;
      Array<ContinuousDoF> tvars;
      for(Index i = 1; i <= 6; i++) {
        tvars.push_back(ContinuousDoF("E", i, -1e+15, 1e+15));
        tvars.back().lock_ID();
      }
      strain_vars.set_variable_basis(tvars, sconvert.symrep_ID());
      Eigen::MatrixXd trans_mat(2, 6);
      trans_mat << 0, 1, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0;
      BasisSet sub_vars(strain_vars.transform_copy(trans_mat));
      BasisSet poly;
      for(Index i = 0; i <= poly_order; i++) {
        BasisSet tmono;
        tmono.construct_invariant_polynomials(Array<BasisSet const *>(1, &strain_vars), primclex.prim().point_group(), i);
        poly.append(tmono);
      }

      poly.accept(VariableLabeler("E(:,%n)"));
      fs::ofstream mfile(root / "poly.m");
      mfile << "corr =[\n";
      for(Index i = 0; i < poly.size(); i++) {
        mfile << "    " << poly[i]->formula();
        if(i + 1 < poly.size())
          mfile << ",...";
        mfile << "\n";
      }
      mfile << "];\n";
      mfile.close();
      if(num_sub != subgrids.size() || num_sub != mags.size()) {
        std::cout << "Option --strain selected.  Based on crystal symmetry, strains can be independently enumerated in the following subspaces:\n";

        std::cout.precision(8);
        std::cout.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
        Index nc = 0;
        for(Index i = 0; i < num_sub; i++) {
          std::cout << " Subspace " << i + 1 << ":\n";
          std::cout << wedges[i].transpose() << "\n\n";
        }
        std::cout << "To proceed, you must specify " << num_sub << " values for both '--mag' and '--subgrids'\n";
        return 1;
      }

      ConfigSelection<false> config_select;
      if(!vm.count("config") || selection == "MASTER") {
        config_select = ConfigSelection<false>(primclex);
      }
      else {
        config_select = ConfigSelection<false>(primclex, selection);
      }

      std::cout << "\n***************************\n" << std::endl;

      std::cout << "Generating perturbations about configurations " << std::endl << std::endl;

      bool verbose = false;
      bool print = true;
      for(auto it = config_select.selected_config_begin(); it != config_select.selected_config_end(); ++it) {
        Index num_before = (it->supercell()).config_list().size();
        ConfigEnumStrain<Configuration> enumerator(it->supercell(), *it, subgrids, mags, strain_mode);
        (it->supercell()).add_unique_canon_configs(enumerator.begin(), enumerator.end());
        std::cout << "Enumerated " << (it->supercell()).config_list().size() - num_before << " deformations.\n";
      }
    }
    else if(vm.count("occ")) {
      // want absolute paths
      abs_cspecs_path = fs::absolute(cspecs_path);



      ConfigSelection<false> config_select;
      if(!vm.count("config") || selection == "MASTER") {
        config_select = ConfigSelection<false>(primclex);
      }
      else {
        config_select = ConfigSelection<false>(primclex, selection);
      }

      std::cout << "\n***************************\n" << std::endl;

      std::cout << "Generating perturbations about configurations " << std::endl << std::endl;

      bool verbose = false;
      bool print = true;
      for(auto it = config_select.selected_config_begin(); it != config_select.selected_config_end(); ++it) {
        std::cout << "  " << it->supercell().name() << "/" << it->id() << std::endl;
        it->supercell().enumerate_perturb_configurations(*it, abs_cspecs_path, tol, verbose, print);
      }
    }
    else {
      std::cout << "\n";
      std::cout << desc << std::endl;

      std::cout << "DESCRIPTION" << std::endl;
      std::cout << "    Generate supercells that are perturbations of a reference\n";
      std::cout << "    configuration.                                           \n";
      std::cout << "    - using the --cspecs option, a bspecs.json type file is  \n";
      std::cout << "      required to determine the extent of the perturbations. \n";
      std::cout << "      Currently only 'orbit_branch_specs' are supported.     \n";
      std::cout << "    - perturbations are generated about selected reference   \n";
      std::cout << "      configurations                                         \n";
      std::cout << std::endl;

      return 1;
    }

    std::cout << std::endl << "  DONE." << std::endl << std::endl;

    std::cout << "Writing config_list..." << std::endl;
    primclex.write_config_list();
    std::cout << "  DONE" << std::endl;

    std::cout << std::endl;

    return 0;
    */
  };
}

