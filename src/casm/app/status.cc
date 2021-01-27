#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/casm_functions.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/database/Selection.hh"
#include "casm/global/definitions.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
void status_unitialized(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Initialize a CASM project\n\
- Create and cd to the directory where you want the project to be located.\n\
  This will be called the 'project root directory' or project's 'location'.\n\
- Add a 'prim.json' file to the directory describing the primitive cell.  \n\
  See 'casm format --prim' for the format of the 'prim.json' file.        \n\
- Execute: 'casm init'                                                    \n\
- Several directories are created: 'basis_sets', 'cluster_expansions',    \n\
  'reports', 'symmetry', and 'training_data'                              \n\
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
        3) \"-O3 -Wall -fPIC --std=c++17\" \n\
\n\
    'soflags': \n\
      Shared object construction flags. In order of priority: \n\
        1) User specified by 'casm settings --set-soflags' \n\
        2) $CASM_SOFLAGS \n\
        3) \"-shared -lboost_system\" \n\
\n\
    'casm headers and libraries': \n\
      CASM header files and shared libraries are expected in the following\n\
      locations.                                                          \n\
      In order of priority: \n\
        1) User specified by 'casm settings --set-casm-includedir' and \n\
           'casm settings --set-casm-libdir' \n\
        2) $CASM_INCLUDEDIR and $CASM_LIBDIR \n\
        3) $CASM_PREFIX/include and $CASM_PREFIX/lib \n\
        3) (default search paths) \n\
\n\
    Note: For the 'casm' Python package, $LIBCASM and $LIBCCASM, have \n\
    highest priority for locating libcasm and libccasm, respectively. \n\
\n\
    'boost headers and libraries': \n\
      The boost libraries are expected in the following locations.        \n\
      In order of priority: \n\
        1) User specified by 'casm settings --set-boost-includedir' and \n\
           'casm settings --set-boost-libdir' and \n\
        2) $CASM_BOOST_INCLUDEDIR and $CASM_BOOST_LIBDIR \n\
        3) $CASM_BOOST_PREFIX/include $CASM_BOOST_PREFIX/lib \n\
        4) (default search paths) \n\
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

void composition_unselected(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Select composition axes\n\
- Execute: 'casm composition -d' to display standard composition axes.    \n\
- Then execute 'casm composition -s <#>' to select one of the listed axes.\n\
- If no standard composition axis is satisfactory, edit the file          \n\
  'composition_axes.json' to add your own custom composition axes to the  \n\
  'custom_axes' JSON object.\n\
- See 'casm format --comp' for description and the location of  \n\
   the 'composition_axes.json' file.\n\n";
}

void supercells_ungenerated(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Enumerate supercells\n\
- Execute: 'casm enum --method ScelEnum --max V' to enumerate supercells up to \n\
  volume V (units: number of primitive cells).                            \n\
- Supercells are listed in the SCEL file as well as scel_list.json        \n\
  This file should not usually be edited manually.\n\
- See 'casm enum --desc ScelEnum' for extended help documentation on how to use the\n\
  '--matrix' and '--lattice-directions' options to perform restricted     \n\
  supercell enumeration (i.e. 2d, 1d, multiples of other supercells).     \n\
- See 'casm format' for a description and location of the  \n\
   'SCEL' file.\n\n";
}

void configs_ungenerated(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Enumerate configurations\n\
- Several options are possible:                                        \n\
- Execute: 'casm enum --method ConfigEnumAllOccupations --all' to      \n\
  enumerate configurations for all supercells.                         \n\
- Execute: 'casm enum --method ConfigEnumAllOccupations --min MINV --max MAXV' \n\
  to enumerate configurations for supercells ranging in volume from    \n\
  MINV to MAXV (units: number of primitive cells).                     \n\
- Execute: 'casm enum --method ConfigEnumAllOccupations --scelname NAME' \n\
  to enumerate configurations for a particular supercell.              \n\
- Generated configurations are listed in the 'config_list.json' file.  \n\
  This file should not usually be edited manually.                     \n\
- Use the 'casm view' command to quickly view configurations in your   \n\
  favorite visualization program. See 'casm view -h' for help.         \n\
- See 'casm enum --desc ConfigEnumAllOccupations' for extended help documentation on how to use \n\
  '--filter' command to perform restricted enumeration of              \n\
  configurations.                                                      \n\
- Once you have a cluster expansion, see 'casm format --monte' for     \n\
  a description of how to save configurations enumerated during Monte  \n\
  Carlo calculations.                                                  \n\
                                                                      \n\n";
}

void configs_uncalculated(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Calculate configuration properties\n\
                                                                       \n\
Instructions for volume relaxed VASP energies:                         \n\n\
- Create INCAR, KPOINTS, POSCAR, SPECIES, and 'relax.json' files for   \n\
  VASP in: '$ROOT/training_data/settings/$CURR_CALCTYPE'. See          \n\
  'casm format --vasp' for a description and location of the VASP      \n\
  settings files.                                                      \n\
Instructions for volume relaxed Quantum Espresso energies:             \n\n\
- Create $inputfile, SPECIES, and 'relax.json' (with calculator tag set) files for\n\
  Quantum Espresso in: '$ROOT/training_data/settings/$CURR_CALCTYPE'. See\n\
  'casm format --qe' for a description and location of the Quantum Espresso\n\
  settings files.                                                      \n\n\
For either of the choices above do the following:             \n\n\
- Select which configurations to calculate properties for using the    \n\
  'casm select' command. Use 'casm select --set-on -c ALL' to select all\n\
  configurations. By default, the 'selected' state of each             \n\
  configuration is stored by CASM in the master_selection file,       \n\
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
    '$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/properties.calc.json'\n\
- Execute: 'casm-calc --setup' to setup VASP/QuantumEspresso input files for all\n\
  selected configurations, but not submit the jobs. This is often a     \n\
  useful first step to check that input files have been prepared       \n\
  correctly.                                                           \n\
- Execute: 'casm-calc --submit' to submit VASP/QuantumEspresso jobs for all selected   \n\
  configurations. Only configurations which have not yet been          \n\
  calculated will run.                                                 \n\
- See 'casm-calc -h' for help and other options.                       \n\
- VASP/QuantumEspresso results will be stored at:                      \n\
    '$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/properties.calc.json'\n\
  Results in 'properties.calc.json' are expected to be ordered to match\n\
  the 'POS' file at '$ROOT/training_data/$SCELNAME/$CONFIGID/POS'      \n\
- Execute 'casm update' to read in calculation results from the        \n\
  'properties.calc.json' files once completed. \n\n";
}

void references_unset(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Set chemical reference\n"
           "                                                                   "
           "    \n"
           "- The chemical reference determines the value of the formation "
           "energy  \n"
           "  and chemical potentials calculated by CASM.                      "
           "    \n\n"

           "- Chemical references states are set by specifying a hyperplane in "
           "    \n"
           "  energy/atom - composition (as atom_frac) space. This may be done "
           "by  \n"
           "  specifying the hyperplane explicitly, or by specifying several   "
           "    \n"
           "  reference states with energy/atom and composition (as atom_frac) "
           "for \n"
           "  enough states to span the composition space of the allowed "
           "occupants \n"
           "  specified in the prim. For consistency with other CASM projects, "
           "    \n"
           "  additional reference states extending to other compositional     "
           "    \n"
           "  dimensions may be included also.                                 "
           "    \n\n"

           "- Execute 'casm ref --set-auto' to automatically set project level "
           "    \n"
           "  references using DFT calculated energies from configurations "
           "with    \n"
           "  extreme parametric compositions.\n\n"

           "- Execute 'casm ref --set '...JSON...'' to manually set the "
           "project    \n"
           "  level reference energies. See 'casm ref --help' for more "
           "information.\n\n"

           "- It is also possible to specialize the chemical reference at the  "
           "    \n"
           "  supercell or configuration level.                                "
           "    \n\n"

           "- See 'casm format' for a description and location of the          "
           "    \n"
           "  'chemical_reference.json' file.                                  "
           "    \n\n";
}

void bset_uncalculated(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Generate basis functions \n\
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

void eci_uncalculated(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Fit effective cluster interactions (ECI)\n\
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

void montecarlo(const CommandArgs &args) {
  log() << "NEXT STEPS:\n\n";

  log() << "Monte Carlo calculations\n\
                                                                       \n\
- Use 'casm monte' to run Monte Carlo calculations.                    \n\
- See 'casm monte --format' and 'casm monte -h' for help.              \n\n";
}

namespace Completer {
StatusOption::StatusOption() : OptionHandlerBase("status") {}

void StatusOption::initialize() {
  add_help_suboption();

  m_desc.add_options()("next,n", "Write next steps")(
      "warning,w", "Suppress warnings")("details,d",
                                        "Print detailed information")(
      "all,a", "Print all 'casm status -n' help messages");

  return;
}
}  // namespace Completer

struct PrintDetails {
  PrintDetails(const CommandArgs &_args, const PrimClex &_primclex)
      : args(_args), primclex(_primclex) {}

  const CommandArgs &args;
  const PrimClex &primclex;

  template <typename ConfigType>
  void eval() const {
    int tot_gen = 0;
    int tot_calc = 0;
    int tot_sel = 0;

    log() << traits<ConfigType>::name << ":" << std::endl;
    log() << std::setw(6) << "INDEX"
          << " " << std::setw(30) << "SUPERCELL"
          << "     "
          << "#CONFIGS G / C / S" << std::endl;
    log() << "-----------------------------------------------------------------"
             "----------"
          << std::endl;
    Index i = 0;
    DB::Selection<ConfigType> master_selection(primclex);
    for (const Supercell &scel : primclex.db<Supercell>()) {
      int gen = primclex.db<ConfigType>().scel_range_size(scel.name());
      int calc = 0, sel = 0;
      for (const auto &config :
           primclex.db<ConfigType>().scel_range(scel.name())) {
        if (master_selection.data()[config.name()]) {
          sel++;
        }
        if (is_calculated(config)) {
          calc++;
        }
      }
      tot_gen += gen;
      tot_calc += calc;
      tot_sel += sel;
      log() << std::setw(6) << i << " " << std::setw(30) << scel.name()
            << "     " << gen << " / " << calc << " / " << sel << std::endl;
      ++i;
    }
    log() << "-----------------------------------------------------------------"
             "----------"
          << std::endl;
    log() << std::setw(6) << " "
          << " " << std::setw(30) << "TOTAL"
          << "     " << tot_gen << " / " << tot_calc << " / " << tot_sel
          << std::endl;
    log() << "\nG:Generated, C:Calculated, S:Selected" << std::endl
          << std::endl;
  }
};

int status_command(const CommandArgs &args) {
  po::variables_map vm;

  /// Set command line options using boost program_options
  Completer::StatusOption status_opt;
  try {
    po::store(
        po::parse_command_line(args.argc(), args.argv(), status_opt.desc()),
        vm);  // can throw

    /** --help option
     */
    if (vm.count("help")) {
      log() << "\n";
      log() << status_opt.desc() << std::endl;

      return 0;
    }

    if (vm.count("desc")) {
      log() << "\n";
      log() << status_opt.desc() << std::endl;
      log() << "DESCRIPTION" << std::endl;
      log() << "    Get status information for the current CASM project.\n\n";

      return 0;
    }

    po::notify(vm);  // throws on error, so do after help in case
    // there are any problems
  } catch (po::error &e) {
    err_log() << "ERROR: " << e.what() << std::endl << std::endl;
    err_log() << status_opt.desc() << std::endl;
    return 1;
  } catch (std::exception &e) {
    err_log() << "Unhandled Exception reached the top of main: " << e.what()
              << ", application will now exit" << std::endl;
    return 1;
  }

  /// 1) Check if a project exists

  log() << "\n#################################\n\n";

  log() << "CASM status:\n\n";

  if (vm.count("all")) {
    log() << "\n#################################\n\n";
    status_unitialized(args);
  }

  const fs::path &root = args.root;

  if (root.empty()) {
    log() << "1) Project initialized: FALSE\n\n";

    if (vm.count("next")) {
      log() << "\n#################################\n\n";

      status_unitialized(args);
    } else {
      log() << "For next steps, run 'casm status -n'\n\n";
    }

    return 0;
  }

  {
    DirectoryStructure dir(root);

    if (!fs::exists(dir.prim())) {
      log() << " ERROR\n\n";

      log() << "- Found a CASM project, but no " << dir.prim() << " file."
            << std::endl;
      log() << "- CASM project location: " << root << std::endl;
      log() << "Please add a prim.json file, or rm the '.casm' directory."
            << std::endl
            << std::endl;

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

  log() << "1) Project initialized: TRUE\n\n";
  log() << "- Project name: " << primclex.settings().project_name()
        << std::endl;
  log() << "- Project location: " << primclex.dir().root_dir().string()
        << std::endl;

  // it'd be nice to just read this...
  SymGroup prim_pg(SymGroup::lattice_point_group(primclex.prim().lattice()));
  log() << "- Lattice point group size: " << prim_pg.size() << std::endl;
  log() << "- Lattice point group is " << prim_pg.get_name() << std::endl;
  log() << "- Factor group size: " << primclex.prim().factor_group().size()
        << std::endl;
  log() << "- Crystal point group is: "
        << primclex.prim().point_group().get_name() << std::endl;
  if (!vm.count("warning")) {
    if (primclex.prim().factor_group().size() > prim_pg.size()) {
      log() << "*** Warning: Finding a factor group that is larger than the "
               "lattice \n"
            << "             point group implies that your structure is not "
               "primitive."
            << std::endl;
    }
  }
  log() << std::endl << std::endl;

  /// 2) Composition axes

  if (vm.count("all")) {
    log() << "\n#################################\n\n";
    log() << "\n2) Composition axes \n\n";
    composition_unselected(args);
  }

  log() << "2) Composition axes \n";

  log() << "- Composition axes selected: ";

  if (!primclex.has_composition_axes()) {
    log() << "FALSE\n\n";

    if (vm.count("next")) {
      log() << "\n#################################\n\n";

      composition_unselected(args);
    } else {
      log() << "For next steps, run 'casm status -n'\n\n";
    }

    return 0;
  }

  log() << "TRUE\n\n\n";
  /*
      // It'd be nice to note standard vs. custom axes, and just '*' the current
     composition axes log() << "- Standard & custom composition axes: " <<
     std::endl << std::endl;
      primclex.get_param_comp().print_composition_axes(log());
      log() << std::endl;
      log() << "- Current composition axes: " << std::endl << std::endl;
      primclex.get_param_comp().print_curr_composition_axes(log());
      log() << std::endl << std::endl;
  */

  /// 3) Configuration generation

  if (vm.count("all")) {
    log() << "\n#################################\n\n";
    log() << "3) Generate configurations \n\n";

    supercells_ungenerated(args);

    configs_ungenerated(args);
  }

  log() << "\n3) Generate configurations \n";

  DB::Selection<Configuration> master_selection(primclex.db<Configuration>());
  int tot_gen = master_selection.size();
  int tot_sel = master_selection.selected_size();
  int tot_calc = 0;
  for (const auto &config : master_selection.all()) {
    tot_calc += is_calculated(config);
  }

  log() << "- Number of supercells generated: "
        << primclex.db<Supercell>().size() << "\n";
  log() << "- Number of configurations generated: " << tot_gen << "\n";
  log() << "- Number of configurations currently selected: " << tot_sel << "\n";

  if (primclex.db<Supercell>().size() == 0) {
    if (vm.count("next")) {
      log() << "\n#################################\n\n";

      supercells_ungenerated(args);
    } else {
      log() << "For next steps, run 'casm status -n'\n\n";
    }

    return 0;
  }

  if (tot_gen == 0) {
    if (vm.count("next")) {
      log() << "\n#################################\n\n";

      configs_ungenerated(args);
    } else {
      log() << "For next steps, run 'casm status -n'\n\n";
    }

    if (vm.count("details")) {
      DB::for_each_config_type(PrintDetails(args, primclex));
    } else {
      log() << "For the number of configurations generated, calculated,\n and "
               "selected by supercell, run 'casm status -d'\n\n";
    }

    return 0;
  }

  log() << std::endl << std::endl;

  /// 4) Calculate configuration properties

  if (vm.count("all")) {
    log() << "\n#################################\n\n";
    log() << "4) Calculate configuration properties\n\n";
    configs_uncalculated(args);
  }

  log() << "4) Calculate configuration properties\n";
  log() << "- Current calctype: " << calctype << "\n";
  log() << "- Current cluster expansion: " << desc.name << "\n";
  log() << "- Number of configurations calculated: " << tot_calc << " / "
        << tot_gen << " generated (Update with 'casm update')\n\n";

  if (vm.count("details")) {
    DB::for_each_config_type(PrintDetails(args, primclex));
  } else {
    log() << "For the number of configurations generated, calculated,\n and "
             "selected by supercell, run 'casm status -d'\n\n";
  }

  if (tot_calc == 0) {
    if (vm.count("next")) {
      log() << "\n#################################\n\n";

      configs_uncalculated(args);
    } else {
      log() << "For next steps, run 'casm status -n'\n\n";
    }

    return 0;
  }

  log() << std::endl;

  /// 5) Choose chemical reference

  if (vm.count("all")) {
    log() << "\n#################################\n\n";
    log() << "5) Choose chemical reference\n\n";
    references_unset(args);
  }

  log() << "5) Choose chemical reference\n";

  log() << "- Chemical reference set: ";
  if (primclex.has_chemical_reference()) {
    log() << "TRUE"
          << "\n";
  } else {
    log() << "FALSE"
          << "\n";
  }
  log() << "\n";

  if (primclex.has_chemical_reference()) {
    log() << "To show the chemical reference, run 'casm ref -d'\n\n";
  } else {
    log() << "No chemical reference set." << std::endl << std::endl;

    if (vm.count("next")) {
      log() << "\n#################################\n\n";

      references_unset(args);
    } else {
      log() << "For next steps, run 'casm status -n'\n\n";
    }

    return 0;
  }

  log() << std::endl;

  /// 6) Generate basis functions:

  if (vm.count("all")) {
    log() << "\n#################################\n\n";
    log() << "6) Generate basis functions: \n\n";
    bset_uncalculated(args);
  }

  log() << "6) Generate basis functions: ";

  if (!fs::exists(dir.clexulator_src(settings.project_name(), bset))) {
    log() << "FALSE\n\n";

    if (vm.count("next")) {
      log() << "\n#################################\n\n";

      bset_uncalculated(args);
    } else {
      log() << "For next steps, run 'casm status -n'\n\n";
    }

    return 0;
  }
  log() << "TRUE\n\n\n";

  /// 7) Fit effective cluster interactions (ECI):

  if (vm.count("all")) {
    log() << "\n#################################\n\n";
    log() << "7) Fit effective cluster interactions (ECI): \n\n";
    eci_uncalculated(args);
  }

  log() << "7) Fit effective cluster interactions (ECI): ";

  if (!fs::exists(dir.eci(property, calctype, ref, bset, eci))) {
    log() << "FALSE\n\n";

    if (vm.count("next")) {
      log() << "\n#################################\n\n";

      eci_uncalculated(args);
    } else {
      log() << "For next steps, run 'casm status -n'\n\n";
    }

    return 0;
  }
  log() << "TRUE\n\n\n";

  /// 7) Monte Carlo

  log() << std::endl;

  if (vm.count("next") || vm.count("all")) {
    log() << "\n#################################\n\n";
    log() << "8) Monte Carlo Calculations: \n\n";
    montecarlo(args);
  } else {
    log() << "For next steps, run 'casm status -n'\n\n";
  }

  return 0;
};

}  // namespace CASM
