#include "casm/app/casm_functions.hh"
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

namespace Completer {
PerturbOption::PerturbOption() : OptionHandlerBase("perturb") {}

const fs::path &PerturbOption::cspecs_path() const { return m_cspecs_path; }

void PerturbOption::initialize() {
  add_help_suboption();
  add_configlist_suboption();
  m_desc.add_options()("occ", "Perturb occupations")(
      "cspecs",
      po::value<fs::path>(&m_cspecs_path)->value_name(ArgHandler::path()),
      "Cluster specifications file defining perturbation");
  //("strain,s", "Generate strain perturbations")
  //("gridsize", po::value<std::vector<Index> > (&subgrids)->multitoken(), "Size
  //of grid for each subspace")
  //("mag", po::value<std::vector<double> > (&mags)->multitoken(), "Magnitude of
  //grid")
  //("poly-order", po::value<Index> (&poly_order)->default_value(6), "Max order
  //of strain polynomials")
  //("strainmode", po::value<std::string> (&strain_mode)->default_value("GL"),
  //"Strain mode name");
  return;
}
}  // namespace Completer

// ///////////////////////////////////////
// 'perturb' function for casm
//    (add an 'if-else' statement in casm.cpp to call this)

int perturb_command(const CommandArgs &args) {
  throw std::runtime_error("Error: 'casm perturb' is being re-implemented");

  /*

  double tol = CASM::TOL;
  fs::path cspecs_path, abs_cspecs_path;
  fs::path selection;
  COORD_TYPE coordtype = CASM::FRAC;
  po::variables_map vm;
  std::vector<Index> subgrids;
  std::vector<double> mags;
  std::string strain_mode;
  //Index poly_order;

  /// Set command line options using boost program_options
  Completer::PerturbOption perturb_opt;

  try {
    po::store(po::parse_command_line(args.argc(), args.argv(),
  perturb_opt.desc()), vm); // can throw

    // --help option
    if(vm.count("help")) {
      log() << "\n";
      log() << perturb_opt.desc() << std::endl;

      return 0;
    }

    if(vm.count("desc")) {
      log() << "\n";
      log() << perturb_opt.desc() << std::endl;

      log() << "DESCRIPTION" << std::endl;
      log() << "    Generate supercells that are perturbations of a
  reference\n"; log() << "    configuration. \n"; log() << "    - using the
  --cspecs option, a bspecs.json type file is  \n"; log() << "      required to
  determine the extent of the perturbations. \n"; log() << "      Currently only
  'orbit_branch_specs' are supported.     \n"; log() << "    - perturbations are
  generated about selected reference   \n"; log() << "      configurations \n";
      log() << std::endl;

      return 0;
    }

    po::notify(vm); // throws on error, so do after help in case
    // there are any problems

    cspecs_path = perturb_opt.cspecs_path();
    selection = perturb_opt.selection_path();
  }
  catch(po::error &e) {
    err_log() << "ERROR: " << e.what() << std::endl << std::endl;
    err_log() << perturb_opt.desc() << std::endl;
    return 1;
  }
  catch(std::exception &e) {
    err_log() << "Unhandled Exception reached the top of main: "
                 << e.what() << ", application will now exit" << std::endl;
    return 1;

  }

  COORD_MODE C(coordtype);

  const fs::path &root = args.root;
  if(root.empty()) {
    err_log().error("No casm project found");
    err_log() << std::endl;
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
    std::vector<Eigen::MatrixXd> wedges =
  sconvert.irreducible_wedges(primclex.prim().point_group(), mult);

    Index num_sub = wedges.size();

    BasisSet strain_vars;
    Array<ContinuousDoF> tvars;
    for(Index i = 1; i <= 6; i++) {
      tvars.push_back(ContinuousDoF(traits, "E", i, -1e+15, 1e+15));
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
      tmono.construct_invariant_polynomials(Array<BasisSet const *>(1,
  &strain_vars), primclex.prim().point_group(), i); poly.append(tmono);
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
      log() << "Option --strain selected.  Based on crystal symmetry, strains
  can be independently enumerated in the following subspaces:\n";

      static_cast<std::ostream &>(log()).precision(8);
      static_cast<std::ostream &>(log()).flags(std::ios::showpoint |
  std::ios::fixed | std::ios::right); Index nc = 0; for(Index i = 0; i <
  num_sub; i++) { log() << " Subspace " << i + 1 << ":\n"; log() <<
  wedges[i].transpose() << "\n\n";
      }
      log() << "To proceed, you must specify " << num_sub << " values for both
  '--mag' and '--subgrids'\n"; return 1;
    }

    ConfigSelection<false> config_select;
    if(!vm.count("config") || selection == "MASTER") {
      config_select = ConfigSelection<false>(primclex);
    }
    else {
      config_select = ConfigSelection<false>(primclex, selection);
    }

    log() << "\n***************************\n" << std::endl;

    log() << "Generating perturbations about configurations " << std::endl <<
  std::endl;

    bool verbose = false;
    bool print = true;
    for(auto it = config_select.selected_config_begin(); it !=
  config_select.selected_config_end(); ++it) { Index num_before =
  (it->get_supercell()).get_config_list().size(); ConfigEnumStrain
  enumerator(it->get_supercell(), *it, subgrids, mags, strain_mode);
      (it->get_supercell()).add_unique_canon_configs(enumerator.begin(),
  enumerator.end()); log() << "Enumerated " <<
  (it->get_supercell()).get_config_list().size() - num_before << "
  deformations.\n";
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

    log() << "\n***************************\n" << std::endl;

    log() << "Generating perturbations about configurations " << std::endl <<
  std::endl;

    bool verbose = false;
    bool print = true;
    for(auto it = config_select.selected_config_begin(); it !=
  config_select.selected_config_end(); ++it) { log() << "  " <<
  it->get_supercell().get_name() << "/" << it->get_id() << std::endl;
      it->get_supercell().enumerate_perturb_configurations(*it, abs_cspecs_path,
  tol, verbose, print);
    }
  }
  else {
    log() << "\n";
    log() << perturb_opt.desc() << std::endl;

    log() << "DESCRIPTION" << std::endl;
    log() << "    Generate supercells that are perturbations of a reference\n";
    log() << "    configuration.                                           \n";
    log() << "    - using the --cspecs option, a bspecs.json type file is  \n";
    log() << "      required to determine the extent of the perturbations. \n";
    log() << "      Currently only 'orbit_branch_specs' are supported.     \n";
    log() << "    - perturbations are generated about selected reference   \n";
    log() << "      configurations                                         \n";
    log() << std::endl;

    return 1;
  }

  log() << std::endl << "  DONE." << std::endl << std::endl;

  log() << "Writing config_list..." << std::endl;
  primclex.write_config_list();
  log() << "  DONE" << std::endl;

  log() << std::endl;

  return 0;
  */
};
}  // namespace CASM
