#include "status.hh"

#include<cstring>

#include "casm/CASM_classes.hh"
#include "casm/app/casm_functions.hh"

namespace CASM {


  void status_unitialized() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Initialize a CASM project\n\
- Create and cd to the directory where you want the project to be located.\n\
  This will be called the 'project root directory' or project's 'location'.\n\
- Add a 'prim.json' file to the directory describing the primitive cell.  \n\
  See 'casm format --prim' for the format of the 'prim.json' file.        \n\
- Execute: 'casm init --name myproject'                                   \n\
- The 'basis_sets' and 'cluster_expansions' directories are created with a\n\
  default format.                                                         \n\
- If necessary, set compilation options using                             \n\
    'casm settings --set-compile-options' and                             \n\
    'casm settings --set-so-options'.                                     \n\
  This may be necessary if, for instance, the CASM header files are       \n\
  installed in a location that is not in your default compiler search path.\n\
- Subsequently, work on 'myproject' can be done by executing 'casm' from  \n\
  the project's root directory or any subdirectory.                       \n\
\n\
- See 'casm format' for descriptions and locations of the 'prim.json' file.\n";
  }

  void standard_composition_uncalculated() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Calculate standard composition axes\n\
- Execute: 'casm composition --calc'                                      \n\
- If successful, standard composition axes will be printed to screen and  \n\
  saved in the 'composition_axes.json' file.                              \n\
- Then execute 'casm composition -s key' to select one of the listed axes.\n\
- If none of the composition axes are satisfactory, edit the file         \n\
  'composition_axes.json' to add your own custom composition axes to the  \n\
  'custom_axes' JSON object.\n\n";

    std::cout <<
              "- See 'casm format' for a description and the location of  \n\
   the 'composition_axes.json' file.\n\n";

  }

  void composition_unselected() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Select composition axes\n\
- Execute: 'casm composition -d' to display composition axes.             \n\
- Then execute 'casm composition -s N' to select one of the listed axes.  \n\
- If none of the composition axes are satisfactory, edit the file         \n\
  'composition_axes.json' to add your own custom composition axes to the  \n\
  'custom_axes' JSON object.\n\n";

    std::cout <<
              "- See 'casm format' for a description and the location of  \n\
   the 'composition_axes.json' file.\n\n";

  }

  void supercells_ungenerated() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Enumerate supercells\n\
- Execute: 'casm enum --supercells --max V' to enumerate supercells up to \n\
  volume V (units: number of primitive cells).                            \n\
- Supercells are listed in the SCEL file.\n\n";

    std::cout <<
              "- See 'casm format' for a description and location of the  \n\
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
  This file should not usually be edited manually.                     \n\n";

    std::cout <<
              "- See 'casm format' for a description and location of   \n\
   the 'config_list.json' file.                                        \n\
 - See 'casm format' for a description and location of                 \n\
   the data files related to a particular configuration.\n\n";

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
  configurations. By default, the 'is selected?' state of each         \n\
  configuration is stored by CASM in the master config_list.json file, \n\
  located in the hidden '.casm' directory. You can also save additional\n\
  selection using the 'casm select -o' option to write a selection to a\n\
  file. Selections may be operated on to create new selections that    \n\
  are subsets, unions, or intersections of existing selections.        \n\
  Selection files may also be edited manually or via programs for more \n\
  complex selections than currently supported by 'casm select'. For all\n\
  options related to selection configurations, see 'casm select -h'.   \n\
- Selections may be used to query the properties of particular         \n\
  configurations using the 'casm query' command. See 'casm query -h'   \n\
  for the complete list of options.                                    \n\
- Execute 'casm run -e \"vasp.relax\" --write-pos' to submit  \n\
  VASP jobs for all selected configurations. This depends on the python\n\
  modules 'pbs', 'casm', and 'vasp' being installed and the script     \n\
  'vasp.relax' being in the PATH. Only configurations which have not   \n\
  yet been calculated will run.                                        \n\
  *Note: You can also use 'casm run -e \"vasp.setup\" --write-pos to   \n\
  setup VASP input files for all selected configuration, but not submit\n\
  the jobs. This is often a useful first step to check that input files\n\
  have been prepared correctly.*                                       \n\
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
  functions.                                                           \n\n";

    std::cout <<
              "- See 'casm format --bspecs' for a description and location of   \n\
   the 'bspecs.json' files.\n\n";
  }

  void eci_uncalculated() {
    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Fit effective cluster interactions (ECI)\n\
                                                                       \n\
Instructions for fitting ECI:                                          \n\n\
- Select which configurations to use as the training data with the     \n\
  'casm select' command. Use 'casm select --set on' to select all      \n\
  configurations. See 'casm select -h' for options.                    \n\
- Execute 'casm fit' to generate input files for fitting eci with the  \n\
  program 'eci_search'.                                                \n\
- Execute 'eci_search -h' for descriptions of the available fitting    \n\
  options                                                              \n\
- Results will be stored at:                                           \n\
    'root/cluster_expansions/clex.formation_energy/SCELNAME/CURR_CALCTYPE/CURR_REF/CURR_ECI/energy\n\
    'root/cluster_expansions/clex.formation_energy/SCELNAME/CURR_CALCTYPE/CURR_REF/CURR_ECI/eci.in\n\
    'root/cluster_expansions/clex.formation_energy/SCELNAME/CURR_CALCTYPE/CURR_REF/CURR_ECI/corr.in\n\n";

    std::cout <<
              "- See 'casm format --fit' for a description and location of   \n\
   the 'energy', 'eci.in', and 'corr.in' files.\n\n";
  }

  void advanced_steps() {

    std::cout << "NEXT STEPS:\n\n";

    std::cout <<
              "Advanced steps\n\
                                                                       \n\
- Alternative calculation settings, composition axes, or reference     \n\
  states can be explored within a single CASM project.                 \n\
                                                                       \n\
- Use 'casm settings' to add a new calculation type. This will create  \n\
  directories with alternative VASP calculation settings, such as a    \n\
  different psuedopotential or spin polarization setting. Then use     \n\
  'casm run' as usual to calculate configuration properties using the  \n\
  alternative settings.                                                \n\
                                                                       \n\
- Use 'casm settings' to add an alternative reference states, then use \n\
  'casm composition' and 'casm ref' as usual to set alternative        \n\
  composition axes or reference states.\n\n";

  }

  int update_eci_format(fs::path root) {

    DirectoryStructure dir(root);
    ProjectSettings set(root);

    // convert eci.out to eci.json (if eci.json does not exist)
    for(auto clex : dir.all_clex()) {
      for(auto bset : dir.all_bset()) {
        for(auto calctype : dir.all_calctype()) {
          for(auto ref : dir.all_ref(calctype)) {
            for(auto eci : dir.all_eci(clex, calctype, ref, bset)) {

              auto eci_path = dir.eci(clex, calctype, ref, bset, eci);
              auto eci_out = dir.eci_out(clex, calctype, ref, bset, eci);
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

  int status_command(int argc, char *argv[]) {

    po::variables_map vm;

    try {

      /// Set command line options using boost program_options
      po::options_description desc("'casm status' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("next,n", "Write next steps")
      ("warning,w", "Suppress warnings")
      ("details,d", "Print detailed information")
      ("update,u", "Update file formats for current version");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION" << std::endl;
          std::cout << "    Get status information for the current CASM project.\n\n";

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

    if(vm.count("update")) {
      return update_format(find_casmroot(fs::current_path()));
    }


    /// 1) Check if a project exists

    std::cout << "\n#################################\n\n";

    std::cout << "CASM status:\n\n";
    std::cout << "1) Project initialized: ";

    fs::path root = find_casmroot(fs::current_path());

    if(root.empty()) {
      std::cout << " FALSE\n\n";

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


    std::ostringstream tmp;
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, tmp);
    std::cout << "  DONE." << std::endl << std::endl;

    const DirectoryStructure &dir = primclex.dir();
    const ProjectSettings &settings = primclex.settings();

    std::string calctype = settings.calctype();
    std::string ref = settings.ref();
    std::string bset = settings.bset();
    std::string clex = settings.clex();
    std::string eci = settings.eci();


    std::cout << "TRUE\n";
    std::cout << "- Project name: " << primclex.name() << std::endl;
    std::cout << "- Project location: " << primclex.get_path().string() << std::endl;

    // it'd be nice to just read this...
    SymGroup prim_pg;
    primclex.get_prim().lattice().generate_point_group(prim_pg);
    prim_pg.character_table();
    std::cout << "- Lattice point group size: " << prim_pg.size() << std::endl;
    std::cout << "- Lattice point group is " << prim_pg.get_name() << std::endl;
    std::cout << "- Factor group size: " << primclex.get_prim().factor_group().size() << std::endl;
    std::cout << "- Crystal point group is: " << primclex.get_prim().point_group().get_name() << std::endl;
    if(!vm.count("warning")) {
      if(primclex.get_prim().factor_group().size() > prim_pg.size()) {
        std::cout << "*** Warning: Finding a factor group that is larger than the lattice \n"
                  << "             point group implies that your structure is not primitive." << std::endl;
      }
    }
    std::cout << std::endl << std::endl;


    /// 2) Composition axes

    std::cout << "2) Composition axes \n";

    std::cout << "- Standard composition axes calculated: ";

    if(!fs::is_regular_file(primclex.dir().composition_axes(calctype, ref))) {
      std::cout << "FALSE\n\n";

      if(vm.count("next")) {
        std::cout << "\n#################################\n\n";

        standard_composition_uncalculated();
      }
      else {
        std::cout << "For next steps, run 'casm status -n'\n\n";
      }

      return 0;
    }

    std::cout << "TRUE\n";

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

    std::cout << "3) Generate configurations \n";

    int tot_gen = 0;
    int tot_calc = 0;
    int tot_sel = 0;

    for(int i = 0; i < primclex.get_supercell_list().size(); i++) {
      int gen = primclex.get_supercell(i).get_config_list().size();
      int calc = 0, sel = 0;
      const Supercell &scel = primclex.get_supercell(i);
      for(int j = 0; j < scel.get_config_list().size(); j++) {
        if(scel.get_config(j).selected()) {
          sel++;
        }
        if(scel.get_config(j).calc_properties().contains("relaxed_energy")) {
          calc++;
        }
      }
      tot_gen += gen;
      tot_calc += calc;
      tot_sel += sel;
    }

    std::cout << "- Number of supercells generated: " << primclex.get_supercell_list().size() << "\n";
    std::cout << "- Number of configurations generated: " << tot_gen << "\n";
    std::cout << "- Number of configurations currently selected: " << tot_sel << "\n";

    if(primclex.get_supercell_list().size() == 0) {

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

    std::cout << "4) Calculate configuration properties\n";
    std::cout << "- Current calctype: " << primclex.get_curr_calctype() << "\n";
    std::cout << "- Current cluster expansion: " << primclex.get_curr_clex() << "\n";
    std::cout << "- Number of configurations calculated: " << tot_calc << " / " << tot_gen << " generated (Update with 'casm update')\n\n";

    if(vm.count("details")) {
      //std::cout << std::setw(6) << " " << " " << std::setw(30) << " " << "     " << "#CONFIGS" << std::endl;
      std::cout << std::setw(6) << "INDEX" << " " << std::setw(30) << "SUPERCELL" << "     " << "#CONFIGS G / C / S" << std::endl;
      std::cout << "---------------------------------------------------------------------------" << std::endl;
      for(int i = 0; i < primclex.get_supercell_list().size(); i++) {
        int tot_gen = 0;
        int tot_calc = 0;
        int tot_sel = 0;

        int gen = primclex.get_supercell(i).get_config_list().size();
        int calc = 0, sel = 0;
        const Supercell &scel = primclex.get_supercell(i);
        for(int j = 0; j < scel.get_config_list().size(); j++) {
          if(scel.get_config(j).selected()) {
            sel++;
          }
          if(scel.get_config(j).calc_properties().contains("relaxed_energy")) {
            calc++;
          }
        }
        tot_gen += gen;
        tot_calc += calc;
        tot_sel += sel;
        std::cout << std::setw(6) << i << " " << std::setw(30) << primclex.get_supercell(i).get_name() << "     " << gen << " / " << calc << " / " << sel << std::endl;
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

    std::cout << "7) Fit effective cluster interactions (ECI): ";

    if(!fs::exists(dir.eci(clex, calctype, ref, bset, eci))) {
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


    /// 7) Advanced steps

    std::cout << std::endl;

    if(vm.count("next")) {
      std::cout << "\n#################################\n\n";

      advanced_steps();
    }
    else {
      std::cout << "For next steps, run 'casm status -n'\n\n";
    }




    return 0;

  };

}


