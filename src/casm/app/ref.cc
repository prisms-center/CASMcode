#include "casm/app/casm_functions.hh"
#include "casm/casm_io/json_io/clex.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace ref_impl {

    int initialize_global(fs::path chem_ref_path, const PrimClex &primclex, const jsonParser &json_ref, double lin_alg_tol) {

      jsonParser json = jsonParser::object();
      json["chemical_reference"]["global"] = json_ref;

      ChemicalReference chem_ref = read_chemical_reference(json, primclex.get_prim(), lin_alg_tol);

      primclex.log() << "Initializing the chemical reference to: \n\n";
      ChemicalReferencePrinter p(primclex.log(), chem_ref);
      p.print_all();
      write_chemical_reference(chem_ref, chem_ref_path);
      return 0;
    }

    int update_global(fs::path chem_ref_path, const PrimClex &primclex, const jsonParser &json_ref, double lin_alg_tol) {

      ChemicalReference chem_ref =
        read_chemical_reference(chem_ref_path, primclex.get_prim(), lin_alg_tol);

      auto input = one_chemical_reference_from_json(primclex.get_prim(), json_ref);

      if(input.second.empty()) {
        chem_ref.set_global(input.first);
      }
      else {
        chem_ref.set_global(input.second.begin(), input.second.end(), lin_alg_tol);
      }

      primclex.log() << "Updating the project-wide chemical reference to: \n";
      ChemicalReferencePrinter p(primclex.log(), chem_ref);
      p.print_global();
      write_chemical_reference(chem_ref, chem_ref_path);
      return 0;
    }

    int update_config(std::string configname,
                      fs::path chem_ref_path,
                      const PrimClex &primclex,
                      const jsonParser &json_ref,
                      double lin_alg_tol) {

      if(!fs::exists(chem_ref_path)) {
        primclex.err_log() << "Error using 'casm ref --set --configname': No reference found.\n";
        primclex.err_log() << "  Expected file at: " << chem_ref_path << "\n";
        primclex.err_log() << "Use 'casm ref --set' or 'casm ref --set-auto' to set a project-wide reference first.\n";
        return ERR_MISSING_INPUT_FILE;
      }

      ChemicalReference chem_ref = read_chemical_reference(chem_ref_path, primclex.get_prim(), lin_alg_tol);

      try {
        const Configuration &config = primclex.configuration(configname);
        (void) config;
      }
      catch(...) {
        primclex.err_log() << "Error using 'casm ref --set --configname': \n"
                           "  Could not find configuration with name: " << configname << "\n";
        return ERR_INVALID_ARG;
      }

      auto input = one_chemical_reference_from_json(primclex.get_prim(), json_ref);
      if(input.second.empty()) {
        chem_ref.set_config(configname, input.first);
      }
      else {
        chem_ref.set_config(configname, input.second.begin(), input.second.end(), lin_alg_tol);
      }

      primclex.log() << "Updating the " << configname << " specialized reference to: \n";
      ChemicalReferencePrinter p(primclex.log(), chem_ref);
      p.print_config(configname);
      write_chemical_reference(chem_ref, chem_ref_path);
      return 0;

    }

    int update_supercell(std::string scelname,
                         fs::path chem_ref_path,
                         const PrimClex &primclex,
                         const jsonParser &json_ref,
                         double lin_alg_tol) {

      if(!fs::exists(chem_ref_path)) {
        primclex.err_log() << "Error using 'casm ref --set --scelname': No reference found.\n";
        primclex.err_log() << "  Expected file at: " << chem_ref_path << "\n";
        primclex.err_log() << "Use 'casm ref --set' or 'casm ref --set-auto' to set a project-wide reference first.\n";
        return ERR_MISSING_INPUT_FILE;
      }

      ChemicalReference chem_ref = read_chemical_reference(chem_ref_path, primclex.get_prim(), lin_alg_tol);

      try {
        const Supercell &scel = primclex.get_supercell(scelname);
        (void) scel;
      }
      catch(...) {
        primclex.err_log() << "Error using 'casm ref --set --scelname': \n"
                           "  Could not find supercell with name: " << scelname << "\n";
        return ERR_INVALID_ARG;
      }

      auto input = one_chemical_reference_from_json(primclex.get_prim(), json_ref);
      if(input.second.empty()) {
        chem_ref.set_supercell(scelname, input.first);
      }
      else {
        chem_ref.set_supercell(scelname, input.second.begin(), input.second.end(), lin_alg_tol);
      }

      primclex.log() << "Updating the " << scelname << " specialized reference to: \n";
      ChemicalReferencePrinter p(primclex.log(), chem_ref);
      p.print_supercell(scelname);
      write_chemical_reference(chem_ref, chem_ref_path);
      return 0;

    }

  }

  namespace Completer {
    RefOption::RefOption(): OptionHandlerBase("ref") {}

    const std::string &RefOption::set_str() const {
      return m_set_str;
    }

    void RefOption::initialize() {
      add_help_suboption();
      add_configname_suboption();
      add_scelname_suboption();

      m_desc.add_options()
      //("composition-space", "Display information on current composition space")
      ("display,d", "Display current reference states")
      ("set-auto", "Automatically set project level reference states using DFT results")
      ("set", po::value<std::string>(&m_set_str),
       "Set reference states using user specified compositions and energies "
       "(Default: set project-wide references). \n"
       "See examples below for the form of expected input.")
      ("erase", "Erase reference states (Default: clear project-wide references).")
      ("clex", po::value<std::string>(), "Name of the cluster expansion using the reference");

      return;
    }
  }

  // ///////////////////////////////////////
  // 'ref' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int ref_command(const CommandArgs &args) {

    po::variables_map vm;
    std::string scelname, configname, set_str;
    std::string species_order_string = "\n\n";

    const fs::path &root = args.root;
    if(!root.empty()) {
      std::stringstream ss;
      DirectoryStructure dir(root);
      Structure prim(read_prim(dir.prim()));

      ss << "       For this project, the expected order is:\n"
         << "        '[";
      auto names = prim.get_struc_molecule_name();
      for(int i = 0; i < names.size(); i++) {
        ss << names[i];
        if(i != names.size() - 1) {
          ss << ", ";
        }
      }
      ss << "]'\n\n";

      species_order_string = ss.str();
    }


    Completer::RefOption ref_opt;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, ref_opt.desc()), vm);

      bool call_help = false;

      //quit out if there are no arguments
      if(!vm.count("help") && !vm.count("desc")) {
        if(vm.count("set") + vm.count("display") + vm.count("erase") + vm.count("set-auto") != 1) {
          args.log << "Error in 'casm ref'. Please select one of --display, \n";
          args.log << "--set, --set-auto, or --erase to use this option." << std::endl;

          call_help = true;
        }

        if(vm.count("set")) {
          if(vm.count("scelname") + vm.count("configname") > 1) {
            args.err_log << "Error in 'casm ref --set'. Please select only one of --scelname, --configname \n";

            call_help = true;
          }
        }

        if(vm.count("erase")) {
          if(vm.count("scelname") + vm.count("configname") > 1) {
            args.err_log << "Error in 'casm ref --erase'. Please select only one of --scelname, --configname \n";

            call_help = true;
          }
        }
      }

      /** --help option
       */
      if(vm.count("help") || call_help) {
        args.log << std::endl;
        args.log << ref_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << ref_opt.desc() << std::endl;

        args.log << "DESCRIPTION" << std::endl;
        args.log << "    The chemical reference determines the value of the formation energy  \n"
                 "    and chemical potentials calculated by CASM.                          \n\n"

                 "    Chemical references states are set by specifying a hyperplane in     \n"
                 "    energy/atom - composition (as atom_frac) space. This may be done by  \n"
                 "    specifying the hyperplane explicitly, or by specifying several       \n"
                 "    reference states with energy/atom and composition (as atom_frac) for \n"
                 "    enough states to span the composition space of the allowed occupants \n"
                 "    specified in the prim. For consistency with other CASM projects,     \n"
                 "    additional reference states extending to other compositional         \n"
                 "    dimensions may be included also. The pure Va reference is always 0.  \n\n";

        args.log << "    The input to '--set' can be one of three forms:                      \n\n"

                 "    1) Input the energy_per_species for pure states:                     \n" <<
                 R"(       '{"A": X, "B": X, "C": X}')" << "\n\n" <<

                 "    2) Input reference state composition and energy_per_species:         \n" <<
                 R"(       '[)" << "\n" <<
                 R"(          {"A": 3.4, "C": 2.0, "energy_per_species": 2.0},)" << "\n" <<
                 R"(          {"B": 2.0, "energy_per_species": 4.0}, )" << "\n" <<
                 R"(          {"C": 1.0, "energy_per_species": 3.0}  )" << "\n" <<
                 R"(        ]')" << "\n\n" <<

                 "    3) Input an array of energy_per_species, for each species in prim,   \n"
                 "       including 0.0 for vacancy:                                        \n"
                 "        '[X, X, X]'                                                      \n"
                 << species_order_string;

        args.log << "    When using '--set' it is also possible to specialize the chemical    \n"
                 "    reference at the supercell or configuration level by adding the      \n"
                 "    --scelname or --configname option.                                   \n\n";



        args.log << "    Examples:\n";
        //args.log << "      casm ref --composition-space \n";
        //args.log << "      - Print composition space column matrix of the primitive\n";
        //args.log << "      - Print null space column matrix\n";
        //args.log << "\n";
        args.log << "      casm ref --display \n";
        args.log << "      - Print chemical reference\n";
        args.log << "\n";
        args.log << "      casm ref --set-auto\n";
        args.log << "      - set all reference states using DFT results for configurations with\n";
        args.log << "        extreme compositions.\n";
        args.log << "      - set reference for compositions outside range of this project to 0.0\n";
        args.log << "\n";
        args.log << "      casm ref --set \n"
                 "        '[{\"Zr\":1, \"energy_per_species\":-8.546979385}, \n"
                 "          {\"Zr\":1, \"O\":1, \"energy_per_species\":-9.090697345}]'\n"
                 "      - set Zr and ZrO, with given energy per species, as reference states\n\n";

        args.log << "      casm ref --scelname SCEL3_3_1_1_0_2_2 --set \n"
                 "        '[{\"Zr\":1, \"energy_per_species\":-8.546979385}, \n"
                 "          {\"Zr\":1, \"O\":1, \"energy_per_species\":-9.090697345}]'\n"
                 "      - set reference states as specified for configurations in supercell SCEL3_3_1_1_0_2_2\n\n";

        args.log << "      casm ref --configname SCEL3_3_1_1_0_2_2/2 --set \n"
                 "        '[{\"Zr\":1, \"energy_per_species\":-8.546979385}, \n"
                 "          {\"Zr\":1, \"O\":1, \"energy_per_species\":-9.090697345}]'\n"
                 "      - set reference states as specified for configuration SCEL3_3_1_1_0_2_2/2\n\n";

        args.log << "      casm ref --scelname SCEL3_3_1_1_0_2_2 --erase \n"
                 "      - erase specialized reference states for configurations in supercell SCEL3_3_1_1_0_2_2\n\n";

        args.log << "      casm ref --configname SCEL3_3_1_1_0_2_2/2 --erase \n"
                 "      - erase specialized reference states for configuration SCEL3_3_1_1_0_2_2/2\n\n";


        if(call_help)
          return ERR_INVALID_ARG;

        return 0;
      }

      po::notify(vm);

      scelname = ref_opt.supercell_str();
      configname = ref_opt.config_str();
      set_str = ref_opt.set_str();
    }
    catch(po::error &e) {
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log << ref_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log << "Unhandled Exception reached the top of main: "
                   << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }


    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);
    const ProjectSettings &set = primclex.settings();
    double lin_alg_tol = primclex.settings().lin_alg_tol();

    ClexDescription clex_desc;
    if(!vm.count("clex")) {
      clex_desc = set.default_clex();
    }
    else {
      auto it = set.cluster_expansions().find(vm["clex"].as<std::string>());
      if(it == set.cluster_expansions().end()) {
        args.err_log.error("Invalid --clex value");
        args.err_log << vm["clex"].as<std::string>() << " not found.";
        return ERR_INVALID_ARG;
      }
      clex_desc = it->second;
    }

    std::string calctype = clex_desc.calctype;
    std::string ref = clex_desc.ref;
    fs::path chem_ref_path = primclex.dir().chemical_reference(calctype, ref);
    int result_code = 0;

    if(vm.count("display")) {
      if(!primclex.has_chemical_reference()) {
        args.err_log << "Error using 'casm ref --display': No reference found.\n";
        args.err_log << "  Expected file at: " << chem_ref_path << "\n";
        args.err_log << "Use 'casm ref --set' or 'casm ref --set-auto' to set a reference\n";
        return ERR_MISSING_INPUT_FILE;
      }

      ChemicalReferencePrinter p(args.log, primclex.chemical_reference());
      p.print_all();

      result_code = 0;
    }
    else if(vm.count("set-auto")) {
      try {
        args.log << "  Set reference states automatically.\n\n" << std::endl;
        ChemicalReference chem_ref = auto_chemical_reference(primclex, lin_alg_tol);
        ChemicalReferencePrinter p(args.log, chem_ref);
        p.print_all();
        write_chemical_reference(chem_ref, chem_ref_path);
        result_code = 0;
      }
      catch(std::exception &e) {
        args.err_log << "Error setting reference states automatically.\n\n";
        args.err_log << e.what() << std::endl;
        return ERR_UNKNOWN;
      }
    }
    else if(vm.count("set")) {

      using namespace ref_impl;

      jsonParser json_ref;
      try {
        json_ref = jsonParser::parse(set_str);
      }
      catch(std::exception &e) {
        args.err_log << "Error parsing JSON input for 'casm ref --set ' with: \n"
                     << set_str << std::endl;
        args.err_log << e.what() << std::endl;
        return ERR_INVALID_ARG;
      }

      // --- Set project-wide ref

      if(!(vm.count("scelname") + vm.count("configname"))) {

        if(!fs::exists(chem_ref_path)) {
          // -- Initializing ref
          result_code = initialize_global(chem_ref_path, primclex, json_ref, lin_alg_tol);
        }
        else {
          // -- Updating project-wide ref
          result_code = update_global(chem_ref_path, primclex, json_ref, lin_alg_tol);
        }
      }

      // --- Set config specific ref

      else if(vm.count("configname")) {

        result_code = update_config(configname, chem_ref_path, primclex, json_ref, lin_alg_tol);
      }

      // --- Set supercell specific ref

      else {

        result_code = update_supercell(scelname, chem_ref_path, primclex, json_ref, lin_alg_tol);

      }

    }
    else if(vm.count("erase")) {

      // --- Erase project-wide ref

      if(!(vm.count("scelname") + vm.count("configname"))) {

        if(!fs::exists(chem_ref_path)) {
          args.err_log << "No chemical reference found. \n";
          return ERR_INVALID_ARG;
        }
        else {
          fs::remove(chem_ref_path);
          args.log << "Erased chemical reference" << std::endl;
        }
      }

      ChemicalReference chem_ref = primclex.chemical_reference();

      // --- Erase configuration specific ref

      if(vm.count("configname")) {
        if(!chem_ref.erase_config(configname)) {
          args.err_log << "No " << configname << " specialized reference found. \n";
          return ERR_INVALID_ARG;
        }
        else {
          args.log << "Erased specialized reference for " << configname << std::endl;
          write_chemical_reference(chem_ref, chem_ref_path);
        }
      }

      // --- Erase supercell specific ref

      else {
        if(!chem_ref.erase_supercell(scelname)) {
          args.err_log << "No " << scelname << " specialized reference found. \n";
          return ERR_INVALID_ARG;
        }
        else {
          args.log << "Erased specialized reference for " << scelname << std::endl;
          write_chemical_reference(chem_ref, chem_ref_path);
        }
      }

      result_code = 0;

    }

    if(!result_code) {
      if(args.primclex) {
        args.primclex->refresh(false, false, true, false);
      }
    }

    return result_code;
  }

}
