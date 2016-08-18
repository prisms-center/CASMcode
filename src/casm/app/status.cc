#include "casm/CASM_global_definitions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/casm_functions.hh"
#include "casm/clex/PrimClex.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {
  void status_unitialized() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Initialize a CASM project\n\
- Create and cd to the directory where you want the project to be located.\n\
  This will be called the 'project root directory' or project's 'location'.\n\
- Add a 'prim.json' file to the directory describing the primitive cell.  \n\
  See 'casm format --prim' for the format of the 'prim.json' file.        \n\
- Execute: 'casm init'                                                    \n\
- Several directories are created: 'symmetry', 'basis_sets',              \n\
  'training_data', and 'cluster_expansions'  \n\
- If necessary, set configuration options for runtime compilation and     \n\
  linking by using the 'casm settings' command or by setting environment  \n\
  variables. \n\
                                                                          \n\
    'cxx': \n\
      Specifies compiler to use. In order of priority: \n\
        1) User specified by 'casm settings --set-cxx' (use '' to clear) \n\
        2) $CASM_CXX \n\
        3) $CXX \n\
        4) \"g++\" \n\
\n\
    'cxxflags': \n\
      Compiler flags. In order of priority: \n\
        1) User specified by 'casm settings --set-cxxflags' \n\
        2) $CASM_CXXFLAGS \n\
        3) \"-O3 -Wall -fPIC --std=c++11\" \n\
\n\
    'soflags': \n\
      Shared object construction flags. In order of priority: \n\
        1) User specified by 'casm settings --set-soflags' \n\
        2) $CASM_SOFLAGS \n\
        3) \"-shared -lboost_system\" \n\
\n\
    'casm_prefix': \n\
      If not in a standard search path, CASM header files are expected in \n\
      '$CASM_PREFIX/include', and shared libraries in '$CASM_PREFIX/lib'. \n\
      In order of priority: \n\
        1) User specified by 'casm settings --set-casm-prefix' \n\
        2) $CASM_PREFIX \n\
        3) (default search paths) \n\
\n\
    Note: For the 'casm' Python package, $LIBCASM and $LIBCCASM, have \n\
    highest priority for locating libcasm and libccasm, respectively. \n\
\n\
    'boost_prefix': \n\
      If not in a standard search path, boost libraries are expected in \n\
      '$CASM_BOOST_PREFIX/lib'. \n\
      In order of priority: \n\
        1) User specified by 'casm settings --set-boost-prefix' \n\
        2) $CASM_BOOST_PREFIX \n\
        3) (default search paths) \n\
\n\
    Note: If shared libraries are installed in non-standard locations, you \n\
    may need to set: \n\
      (Linux) export LD_LIBRARY_PATH=$CASM_PREFIX/lib:$CASM_BOOST_PREFIX/lib:$LD_LIBRARY_PATH \n\
      (Mac)   export DYLD_FALLBACK_LIBRARY_PATH=$CASM_PREFIX/lib:$CASM_BOOST_PREFIX/lib:$DYLD_FALLBACK_LIBRARY_PATH \n\
\n\
- Subsequently, work on the CASM project can be done by executing 'casm'  \n\
  from the project's root directory or any subdirectory.                  \n\
- See 'casm format --prim' for description and location of the 'prim.json' file.\n";
  }

  void composition_unselected() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Select composition axes\n\
- Execute: 'casm composition -d' to display standard composition axes.    \n\
- Then execute 'casm composition -s <#>' to select one of the listed axes.\n\
- If no standard composition axis is satisfactory, edit the file          \n\
  'composition_axes.json' to add your own custom composition axes to the  \n\
  'custom_axes' JSON object.\n\
- See 'casm format --comp' for description and the location of  \n\
   the 'composition_axes.json' file.\n\n";

  }

  void supercells_ungenerated() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Enumerate supercells\n\
- Execute: 'casm enum --supercells --max V' to enumerate supercells up to \n\
  volume V (units: number of primitive cells).                            \n\
- Supercells are listed in the SCEL file.\n\
- See 'casm enum --desc' for extended help documentation on how to use the\n\
  '--matrix' and '--lattice-directions' options to perform restricted     \n\
  supercell enumeration (i.e. 2d, 1d, multiples of other supercells).     \n\
- See 'casm format' for a description and location of the  \n\
   'SCEL' file.\n\n";

  }

  void configs_ungenerated() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Enumerate configurations\n\
- Several options are possible:                                        \n\
- Execute: 'casm enum --configs --all' to enumerate configurations for \n\
  for all supercells.                                                  \n\
- Execute: 'casm enum --configs --min MINV --max MAXV' to enumerate    \n\
  configurations for supercells ranging in volume from MINV to MAXV    \n\
  (units: number of primitive cells).                                  \n\
- Execute: 'casm enum --configs --scellname NAME' to enumerate         \n\
  configurations for a particular supercell.                           \n\
- Generated configurations are listed in the 'config_list.json' file.  \n\
  This file should not usually be edited manually.                     \n\
- Use the 'casm view' command to quickly view configurations in your   \n\
  favorite visualization program. See 'casm view -h' for help.         \n\
- See 'casm enum --desc' for extended help documentation on how to use \n\
  '--filter' command to perform restricted enumeration of              \n\
  configurations.                                                      \n\
- Once you have a cluster expansion, see 'casm format --monte' for     \n\
  a description of how to save configurations enumerated during Monte  \n\
  Carlo calculations.                                                  \n\
- See 'casm format --config' for a description and location of         \n\
   the 'config_list.json' file.                                        \n\n";
  }

  void configs_uncalculated() {
    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Calculate configuration properties\n\
                                                                       \n\
Instructions for volume relaxed VASP energies:                         \n\n\
- Create INCAR, KPOINTS, POSCAR, SPECIES, and 'relax.json' files for   \n\
  VASP in: '$ROOT/training_data/settings/$CURR_CALCTYPE'. See          \n\
  'casm format --vasp' for a description and location of the VASP      \n\
  settings files.                                                      \n\
- Select which configurations to calculate properties for using the    \n\
  'casm select' command. Use 'casm select --set on' to select all      \n\
  configurations. By default, the 'selected' state of each             \n\
  configuration is stored by CASM in the master config_list.json file, \n\
  located in the hidden '.casm' directory. The standard selections     \n\
  'MASTER', 'CALCULATED', 'ALL', or 'NONE' may always be used.         \n\
- You can also save additional selection using the 'casm select -o'    \n\
  option to write a selection to a file.                               \n\
- Selections may be operated on to create new selections that are      \n\
  subsets, unions, or intersections of existing selections.            \n\
- Selection files may also be edited manually or via programs for more \n\
  complex selections than currently supported by 'casm select'. For all\n\
  options related to selection configurations, see 'casm select -h'.   \n\
- Selections may be used to query the properties of particular         \n\
  configurations using the 'casm query' command. See 'casm query -h'   \n\
  for the complete list of options.                                    \n\
- Execute: 'casm-calc --setup' to setup VASP input files for all       \n\
  selected configuration, but not submit the jobs. This is often a     \n\
  useful first step to check that input files have been prepared       \n\
  correctly.                                                           \n\
- Execute: 'casm-calc --submit' to submit VASP jobs for all selected   \n\
  configurations. Only configurations which have not yet been          \n\
  calculated will run.                                                 \n\
- See 'casm-calc -h' for help and other options.                       \n\
- VASP results will be stored at:                                      \n\
    '$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/properties.calc.json'\n\
  Results in 'properties.calc.json' are expected to be ordered to match\n\
  the 'POS' file at '$ROOT/training_data/$SCELNAME/$CONFIGID/POS'      \n\
- Execute 'casm update' to read in calculation results from the        \n\
  'properties.calc.json' files once completed. \n\n";

  }

  void references_unset() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Set chemical reference\n"
              "                                                                       \n"
              "- The chemical reference determines the value of the formation energy  \n"
              "  and chemical potentials calculated by CASM.                          \n\n"

              "- Chemical references states are set by specifying a hyperplane in     \n"
              "  energy/atom - composition (as atom_frac) space. This may be done by  \n"
              "  specifying the hyperplane explicitly, or by specifying several       \n"
              "  reference states with energy/atom and composition (as atom_frac) for \n"
              "  enough states to span the composition space of the allowed occupants \n"
              "  specified in the prim. For consistency with other CASM projects,     \n"
              "  additional reference states extending to other compositional         \n"
              "  dimensions may be included also.                                     \n\n"

              "- Execute 'casm ref --set-auto' to automatically set project level     \n"
              "  references using DFT calculated energies from configurations with    \n"
              "  extreme parametric compositions.\n\n"

              "- Execute 'casm ref --set '...JSON...'' to manually set the project    \n"
              "  level reference energies. See 'casm ref --help' for more information.\n\n"

              "- It is also possible to specialize the chemical reference at the      \n"
              "  supercell or configuration level.                                    \n\n"

              "- See 'casm format' for a description and location of the              \n"
              "  'chemical_reference.json' file.                                      \n\n";

  }

  void bset_uncalculated() {
    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Generate basis functions \n\
                                                                       \n\
Instructions for generating basis functions:                           \n\n\
- Write a '$ROOT/basis_sets/$CURR_BSET/bspecs.json' file containing    \n\
  specifications for how which orbits to include and which site basis  \n\
  functions to use in generating cluster expansion basis functions.    \n\
  See 'casm format --bspecs' for an example file.                      \n\
- Execute 'casm bset -u' to generate basis functions. If you edit the  \n\
  'bspecs.json' file, execute 'casm bset -u' again to update basis     \n\
  functions.                                                           \n\
- See 'casm format --bspecs' for description and location of the       \n\
  'bspecs.json' file.\n\n";
  }

  void eci_uncalculated() {
    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Fit effective cluster interactions (ECI)\n\
                                                                       \n\
Instructions for fitting ECI:                                          \n\n\
- Create a new directory within the CASM project, for example:         \n\
    mkdir fit_1 && cd fit_1                                            \n\
- Select which configurations to use as the training data with the     \n\
  'casm select' command. To select all calculated configurations:      \n\
    casm select --set 'is_calculated' -o train                         \n\
- See 'casm select -h' for more options.                               \n\
- Create a 'casm-learn' input file. Several example input files can be \n\
  generated from 'casm-learn --exMethodName'. For example:             \n\
    casm-learn --exGeneticAlgorithm > fit_1_ga.json                    \n\
  This file can be edited to adjust the problem being solved (training \n\
  data, weighting scheme, cross validation sets and scoring, linear    \n\
  estimator method, feature selection method, etc.)                    \n\
- See 'casm-learn --settings-format' for description and help with the \n\
  input file.                                                          \n\
- Execute: 'casm-learn -s fit_1_ga.json'                               \n\
- Results are stored in a Hall Of Fame file containing the best        \n\
  solutions as determined from cross validation scores.                \n\
- Different estimator methods (LinearRegression, Lasso, etc.) and      \n\
  different feature selection methods (GeneticAlgorithm, RFE, etc.) can\n\
  be used with the same problem specs (training data, weighting scheme,\n\
  cross validation sets and scoring) and compared in a single Hall Of  \n\
  Fame.                                                                \n\
- When some candidate ECI have been stored in a Hall Of Fame, use the  \n\
  'casm-learn --checkhull' option to check if ground state configurations \n\
  are accurately predicted by the cluster expansion.                   \n\
- When ready, use 'casm-learn --select' to write an 'eci.json' file to \n\
  use for Monte Carlo. \n\
- See 'casm format --eci' for a description and location of the        \n\
  'eci.json' files.\n\n";
  }

  void montecarlo() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Monte Carlo calculations\n\
                                                                       \n\
- Use 'casm monte' to run Monte Carlo calculations.                    \n\
- See 'casm monte --format' and 'casm monte -h' for help.              \n\n";
  }

  int update_eci_format(fs::path root) {

    DirectoryStructure dir(root);
    ProjectSettings set(root);

    // convert eci.out to eci.json (if eci.json does not exist)
    for(auto property : dir.all_property()) {
      for(auto bset : dir.all_bset()) {
        for(auto calctype : dir.all_calctype()) {
          for(auto ref : dir.all_ref(calctype)) {
            for(auto eci : dir.all_eci(property, calctype, ref, bset)) {

              auto eci_path = dir.eci(property, calctype, ref, bset, eci);
              auto eci_out = dir.eci_out(property, calctype, ref, bset, eci);
              fs::path basis_json_path = dir.basis(bset);

              if(!fs::exists(eci_path) && fs::exists(eci_out) && fs::exists(basis_json_path)) {
                // create an eci.json from eci.out and basis.json

                std::cout << "Converting: \n  " << eci_out << "\nto:\n  " << eci_path << std::endl;

                auto eci = read_eci_out(eci_out);

                fs::path basis_json_path = dir.basis(bset);
                // read basis.json
                jsonParser basis_json = jsonParser::parse(basis_json_path);

                // set eci
                for(int i = 0; i < eci.index().size(); ++i) {
                  auto index = eci.index()[i];
                  auto value = eci.value()[i];
                  basis_json["cluster_functions"][index]["eci"] = value;
                }

                //write eci.json
                basis_json.write(eci_path);
                std::cerr << "DONE" << std::endl << std::endl;

              }
            }
          }
        }
      }
    }

    return 0;
  }

  int update_format(fs::path root) {
    return update_eci_format(root);
  }

  namespace Completer {
    StatusOption::StatusOption(): OptionHandlerBase("status") {}

    void StatusOption::initialize() {
      add_help_suboption();

      m_desc.add_options()
      ("next,n", "Write next steps")
      ("warning,w", "Suppress warnings")
      ("details,d", "Print detailed information")
      ("all,a", "Print all 'casm status -n' help messages");

      return;
    }
  }

  int status_command(const CommandArgs &args) {

    po::variables_map vm;

    /// Set command line options using boost program_options
    Completer::StatusOption status_opt;
    try {
      po::store(po::parse_command_line(args.argc, args.argv, status_opt.desc()), vm); // can throw

      /** --help option
      */
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << status_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        std::cout << "\n";
        std::cout << status_opt.desc() << std::endl;
        std::cout << "DESCRIPTION" << std::endl;
        std::cout << "    Get status information for the current CASM project.\n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems
    }
    catch(po::error &e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      std::cerr << status_opt.desc() << std::endl;
      return 1;
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return 1;

    }

    /// 1) Check if a project exists

    std::cout << "\n#################################\n\n";

    std::cout << "CASM status:\n\n";

    if(vm.count("all")) {
      std::cout << "\n#################################\n\n";
      status_unitialized();
    }


    const fs::path &root = args.root;

    if(root.empty()) {
      std::cout << "1) Project initialized: FALSE\n\n";

      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        status_unitialized();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }

    {
      DirectoryStructure dir(root);

      if(!fs::exists(dir.prim())) {
        std::cout << " ERROR\n\n";

        std::cout << "- Found a CASM project, but no '" << dir.prim() << "' file." << std::endl;
        std::cout << "- CASMP project location: " << root << std::endl;
        std::cout << "Please add a prim.json file, or rm the '.casm' directory." << std::endl << std::endl;

        return 1;
      }
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    const DirectoryStructure &dir = primclex.dir();
    const ProjectSettings &settings = primclex.settings();
    const ClexDescription &desc = settings.default_clex();

    std::string property = desc.property;
    std::string calctype = desc.calctype;
    std::string ref = desc.ref;
    std::string bset = desc.bset;
    std::string eci = desc.eci;


    std::cout << "TRUE\n";
    std::cout << "- Project name: " << primclex.settings().name() << std::endl;
    std::cout << "- Project location: " << primclex.dir().root_dir().string() << std::endl;

    // it'd be nice to just read this...
    SymGroup prim_pg;
    primclex.prim().lattice().generate_point_group(prim_pg);
    prim_pg.character_table();
    std::cout << "- Lattice point group size: " << prim_pg.size() << std::endl;
    std::cout << "- Lattice point group is " << prim_pg.get_name() << std::endl;
    std::cout << "- Factor group size: " << primclex.prim().factor_group().size() << std::endl;
    std::cout << "- Crystal point group is: " << primclex.prim().point_group().get_name() << std::endl;
    if(!vm.count("warning")) {
      if(primclex.prim().factor_group().size() > prim_pg.size()) {
        std::cout << "*** Warning: Finding a factor group that is larger than the lattice \n"
                  << "             point group implies that your structure is not primitive." << std::endl;
      }
    }
    std::cout << std::endl << std::endl;

    /// 2) Composition axes

    if(vm.count("all")) {
      std::cout << "\n#################################\n\n";
      std::cout << "\n2) Composition axes \n\n";
      composition_unselected();

    }

    std::cout << "2) Composition axes \n";

    std::cout << "- Composition axes selected: ";

    if(!primclex.has_composition_axes()) {
      std::cout << "FALSE\n\n";

      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        composition_unselected();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }

    std::cout << "TRUE\n\n\n";
    /*
        // It'd be nice to note standard vs. custom axes, and just '*' the current composition axes
        std::cout << "- Standard & custom composition axes: " << std::endl << std::endl;
        primclex.get_param_comp().print_composition_axes(std::cout);
        std::cout << std::endl;
        std::cout << "- Current composition axes: " << std::endl << std::endl;
        primclex.get_param_comp().print_curr_composition_axes(std::cout);
        std::cout << std::endl << std::endl;
    */

    /// 3) Configuration generation


    if(vm.count("all")) {
      std::cout << "\n#################################\n\n";
      std::cout << "3) Generate configurations \n\n";

      supercells_ungenerated();

      configs_ungenerated();

    }

    std::cout << "\n3) Generate configurations \n";

    int tot_gen = 0;
    int tot_calc = 0;
    int tot_sel = 0;

    for(int i = 0; i < primclex.supercell_list().size(); i++) {
      int gen = primclex.supercell(i).config_list().size();
      int calc = 0, sel = 0;
      const Supercell &scel = primclex.supercell(i);
      for(int j = 0; j < scel.config_list().size(); j++) {
        if(scel.config(j).selected()) {
          sel++;
        }
        if(scel.config(j).calc_properties().contains("relaxed_energy")) {
          calc++;
        }
      }
      tot_gen += gen;
      tot_calc += calc;
      tot_sel += sel;
    }

    std::cout << "- Number of supercells generated: " << primclex.supercell_list().size() << "\n";
    std::cout << "- Number of configurations generated: " << tot_gen << "\n";
    std::cout << "- Number of configurations currently selected: " << tot_sel << "\n";

    if(primclex.supercell_list().size() == 0) {

      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        supercells_ungenerated();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }

    if(tot_gen == 0) {

      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        configs_ungenerated();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }

    std::cout << std::endl << std::endl;

    /// 4) Calculate configuration properties


    if(vm.count("all")) {
      std::cout << "\n#################################\n\n";
      std::cout << "4) Calculate configuration properties\n\n";
      configs_uncalculated();

    }

    std::cout << "4) Calculate configuration properties\n";
    std::cout << "- Current calctype: " << primclex.settings().calctype() << "\n";
    std::cout << "- Current cluster expansion: " << primclex.settings().clex() << "\n";
    std::cout << "- Number of configurations calculated: " << tot_calc << " / " << tot_gen << " generated (Update with 'casm update')\n\n";

    if(vm.count("details")) {
      //std::cout << std::setw(6) << " " << " " << std::setw(30) << " " << "     " << "#CONFIGS" << std::endl;
      std::cout << std::setw(6) << "INDEX" << " " << std::setw(30) << "SUPERCELL" << "     " << "#CONFIGS G / C / S" << std::endl;
      std::cout << "---------------------------------------------------------------------------" << std::endl;
      for(int i = 0; i < primclex.supercell_list().size(); i++) {
        int tot_gen = 0;
        int tot_calc = 0;
        int tot_sel = 0;

        int gen = primclex.supercell(i).config_list().size();
        int calc = 0, sel = 0;
        const Supercell &scel = primclex.supercell(i);
        for(int j = 0; j < scel.config_list().size(); j++) {
          if(scel.config(j).selected()) {
            sel++;
          }
          if(scel.config(j).calc_properties().contains("relaxed_energy")) {
            calc++;
          }
        }
        tot_gen += gen;
        tot_calc += calc;
        tot_sel += sel;
        std::cout << std::setw(6) << i << " " << std::setw(30) << primclex.supercell(i).name() << "     " << gen << " / " << calc << " / " << sel << std::endl;
      }
      std::cout << "---------------------------------------------------------------------------" << std::endl;
      std::cout << std::setw(6) << " " << " " << std::setw(30) << "TOTAL" << "     " << tot_gen << " / " << tot_calc << " / " << tot_sel << std::endl;
      std::cout << "\nG:Generated, C:Calculated, S:Selected" << std::endl << std::endl;
    }
    else {
      std::cout << "For the number of configurations generated, calculated,\n and selected by supercell, run 'casm status -d'\n\n";
    }

    if(tot_calc == 0) {
      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        configs_uncalculated();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }

    std::cout << std::endl;


    /// 5) Choose chemical reference


    if(vm.count("all")) {
      std::cout << "\n#################################\n\n";
      std::cout << "5) Choose chemical reference\n\n";
      references_unset();
    }

    std::cout << "5) Choose chemical reference\n";

    std::cout << "- Chemical reference set: ";
    if(primclex.has_chemical_reference()) {
      std::cout << "TRUE" << "\n";
    }
    else {
      std::cout << "FALSE" << "\n";
    }
    std::cout << "\n";

    if(primclex.has_chemical_reference()) {
      std::cout << "To show the chemical reference, run 'casm ref -d'\n\n";
    }
    else {

      std::cout << "No chemical reference set." << std::endl << std::endl;

      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        references_unset();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }

    std::cout << std::endl;


    /// 6) Generate basis functions:


    if(vm.count("all")) {
      std::cout << "\n#################################\n\n";
      std::cout << "6) Generate basis functions: \n\n";
      bset_uncalculated();
    }

    std::cout << "6) Generate basis functions: ";

    if(!fs::exists(dir.clexulator_src(settings.name(), bset))) {
      std::cout << "FALSE\n\n";

      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        bset_uncalculated();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }
    std::cout << "TRUE\n\n\n";

    /// 7) Fit effective cluster interactions (ECI):


    if(vm.count("all")) {
      std::cout << "\n#################################\n\n";
      std::cout << "7) Fit effective cluster interactions (ECI): \n\n";
      eci_uncalculated();
    }

    std::cout << "7) Fit effective cluster interactions (ECI): ";

    if(!fs::exists(dir.eci(property, calctype, ref, bset, eci))) {
      std::cout << "FALSE\n\n";

      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        eci_uncalculated();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }
    std::cout << "TRUE\n\n\n";

    /// 7) Monte Carlo

    std::cout << std::endl;

    if(vm.count("next") || vm.count("all")) {
      std::cout << "\n#################################\n\n";
      std::cout << "8) Monte Carlo Calculations: \n\n";
      montecarlo();
    }
    else {
      std::cout << "For next steps, run 'casm status -n'\n\n";
    }

    return 0;

  };

}


