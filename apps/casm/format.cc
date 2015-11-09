#include "format.hh"

#include <cstring>

#include "casm/CASM_classes.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'format' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int format_command(int argc, char *argv[]) {

    po::variables_map vm;

    try {
      po::options_description desc("'casm format' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("dir,d", "CASM project directory structure summary")
      ("project_settings", "Description and location of 'project_settings' file")
      ("prim", "Description and location of 'prim.json' and 'PRIM' files")
      ("config_list", "Description and location of 'config_list.json' file")
      ("sym", "Description and location of 'lattice_point_group.json', 'factor_group.json' and 'crystal_point_group.json' files")
      ("vasp", "Description and location of VASP settings files")
      ("comp", "Description and location of 'composition_axes.json' file")
      ("bspecs", "Description and location of 'bspecs.json' file")
      ("ref_state", "Description and location of 'properties.ref_state.i.json' files")
      ("scel", "Description and location of 'SCEL' file")
      ("lat", "Description and location of 'LAT' files")
      ("pos", "Description and location of 'POS' files")
      ("fit", "Description and location of the 'energy', 'corr.in', and 'eci.in' files")
      ("monte", "Description and location of the Monte Carlo input file");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        /** --help option
         */
        if(vm.count("help") || vm.size() == 0) {
          std::cout << std::endl;
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION" << std::endl;
          std::cout << "    This option describes the files contained within a CASM project \n";
          std::cout << "    and where to find them. For a summary of the directory structure\n";
          std::cout << "    of a CASM project using VASP for calculating configuration use  \n";
          std::cout << "    the --dir option. Not all files are always present.             \n";

          return 0;
        }

        po::notify(vm);

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

    if(vm.count("dir")) {
      std::cout << "\n### dir ##################\n\n";

      std::cout << "  The expected CASM project directory structure with VASP settings  \n";
      std::cout << "  files:                                                            \n";
      std::cout << "                                                                    \n";
      std::cout << "    $ROOT/.casm                                                     \n";
      std::cout << "      project_settings.json                                         \n";
      std::cout << "      config_list.json                                              \n";
      std::cout << "    $ROOT/                                                          \n";
      std::cout << "      prim.json                                                     \n";
      std::cout << "      (PRIM)                                                        \n";
      std::cout << "      LOG                                                           \n";
      std::cout << "    $ROOT/symmetry/                                                 \n";
      std::cout << "      lattice_point_group.json                                      \n";
      std::cout << "      factor_group.json                                             \n";
      std::cout << "      crystal_point_group.json                                      \n";
      std::cout << "    $ROOT/basis_sets/$CURR_BSET/                                    \n";
      std::cout << "      bspecs.json                                                   \n";
      std::cout << "      clust.json                                                    \n";
      std::cout << "      $NAME_Clexulator.cc                                           \n";
      std::cout << "    $ROOT/training_data/                                            \n";
      std::cout << "      SCEL                                                          \n";
      std::cout << "    $ROOT/training_data/settings/$CURR_CALCTYPE/                    \n";
      std::cout << "      relax.json                                                    \n";
      std::cout << "      INCAR                                                         \n";
      std::cout << "      SPECIES                                                       \n";
      std::cout << "      KPOINTS                                                       \n";
      std::cout << "      POSCAR                                                        \n";
      std::cout << "    $ROOT/training_data/settings/$CURR_CALCTYPE/$CURR_REF/          \n";
      std::cout << "      composition_axes.json                                         \n";
      std::cout << "      properties.ref_state.i.json                                   \n";
      std::cout << "    $ROOT/training_data/$SCELNAME/                                  \n";
      std::cout << "      LAT                                                           \n";
      std::cout << "    $ROOT/training_data/$SCELNAME/$CONFIGID                         \n";
      std::cout << "      POS                                                           \n";
      std::cout << "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE          \n";
      std::cout << "      (VASP results)                                                \n";
      std::cout << "      properties.calc.json                                          \n";
      std::cout << "    $ROOT/cluster_expansions/clex.formation_energy/$CURR_BSET/$CURR_CALCTYPE/$CURR_REF/$CURR_ECI\n";
      std::cout << "      energy                                                        \n";
      std::cout << "      corr.in                                                       \n";
      std::cout << "      eci.in                                                        \n";
      std::cout << "      eci.out                                                       \n";
      std::cout << " \n";
      std::cout << " \n";
      std::cout << "    Variable descriptions:                                          \n";
      std::cout << " \n";
      std::cout << "    $ROOT: root directory of the CASM project                       \n";
      std::cout << " \n";
      std::cout << "    $CURR_BSET: Current basis set, by default this is 'bset.default'.\n";
      std::cout << "    The current value can be inspected via 'casm settings -l'.      \n";
      std::cout << " \n";
      std::cout << "    $CURR_CALCTYPE: Current calctype, by default this is 'calctype.default'.\n";
      std::cout << "    The current value can be inspected via 'casm settings -l'.      \n";
      std::cout << " \n";
      std::cout << "    $CURR_REF: Current composition axes and reference states, by    \n";
      std::cout << "    default this is 'ref.default'. The current value can be inspected\n";
      std::cout << "    via 'casm settings -l'.                                         \n";
      std::cout << " \n";
      std::cout << "    $SCELNAME: Supercell name, in the form SCELV_A_B_C_D_E_F. 'V' is\n";
      std::cout << "    volume of the supercell in terms of the primitive cell, and     \n";
      std::cout << "    'A'-'F' are the values of the hermite normal form of the        \n";
      std::cout << "    transformation matrix.                                          \n";
      std::cout << " \n";
      std::cout << "    $CONFIGID: Configuration id, a unique integer.                  \n";
      std::cout << "\n";
      std::cout << "    Note: The 'settings' heirarchy can be located at the project    \n";
      std::cout << "    level as shown above, or at the supercell or configuration level\n";
      std::cout << "    in order to override calculation, composition, or reference     \n";
      std::cout << "    state settings at the supercell or configuration level.  The    \n";
      std::cout << "    most local settings are always used for a configuration.        \n";
      std::cout << " \n";

    }

    if(vm.count("project_settings")) {
      std::cout << "\n### project_settings.json ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/.casm/project_settings.json\n\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "Current CASM project settings.\n\n\n";

      std::cout << "EXAMPLE:\n";
      std::cout << "-------\n";
      std::cout <<
                "{\n  \"compile_options\" : \"g++ -O3 -Wall -fPIC --std=c++11\",\n  \"curr_bset\" : \"default\",\n  \"curr_calctype\" : \"default\",\n  \"curr_clex\" : \"formation_energy\",\n  \"curr_eci\" : \"default\",\n  \"curr_properties\" : [ \"relaxed_energy\" ],\n  \"curr_ref\" : \"default\",\n  \"name\" : \"ZrO\",\n  \"so_options\" : \"g++ -shared -lboost_system\",\n  \"tol\" : 0.000010000000\n}\n";
      std::cout << "-------\n";
      std::cout << std::endl << std::endl;
    }

    if(vm.count("prim")) {
      std::cout << "\n### prim.json ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/prim.json\n";
      std::cout << "$ROOT/PRIM\n\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "'prim.json' describes the primitive cell structure. It includes the lattice \n";
      std::cout << "vectors, crystal basis sites and a list of possible occupant molecules on each\n";
      std::cout << "basis site.\n\n";

      std::cout << "- Molecule names are case sensitive.\n";
      std::cout << "- 'Va' is reserved for vacancies.\n";
      std::cout << "- The default tolerance for checking symmetry is 1e-5, so basis site coordinates\n";
      std::cout << "  should include 6 significant digits or more.\n\n\n";

      std::cout << "EXAMPLE 1: An FCC ternary alloy of elements A, B, and C\n";
      std::cout << "-------\n";
      std::cout <<
                "{\n  \"basis\" : [\n    {\n      \"coordinate\" : [ 0.000000000000, 0.000000000000, 0.000000000000 ],\n      \"occupant_dof\" : [ \"A\", \"B\", \"C\" ]\n    }\n  ],\n  \"coordinate_mode\" : \"Fractional\",\n  \"description\" : \"Face-centered Cubic (FCC, cF)\",\n  \"lattice_vectors\" : [\n    [ 2.000000000000, 2.000000000000, 0.000000000000 ],\n    [ 0.000000000000, 2.000000000000, 2.000000000000 ],\n    [ 2.000000000000, 0.000000000000, 2.000000000000 ]\n  ],\n  \"title\" : \"ABC\"\n}\n";
      std::cout << "-------\n\n";

      std::cout << "EXAMPLE 2: HCP Zr with O in octahedral interstitial positions\n";
      std::cout << "-------\n";
      std::cout <<
                "\n{\n  \"basis\" : [\n    {\n      \"coordinate\" : [ 0.0, 0.0, 0.0 ],\n      \"occupant_dof\" : [ \"Zr\" ]\n    },\n    {\n      \"coordinate\" : [ 0.666666, 0.333333, 0.5 ],\n      \"occupant_dof\" : [ \"Zr\" ]\n    },\n    {\n      \"coordinate\" : [ 0.333333, 0.666666, 0.25 ],\n      \"occupant_dof\" : [ \"Va\", \"O\" ]\n    },\n    {\n      \"coordinate\" : [ 0.333333, 0.666666, 0.75 ],\n      \"occupant_dof\" : [ \"Va\", \"O\" ]\n    }\n  ],\n  \"coordinate_mode\" : \"Fractional\",\n  \"description\" : \"hcp Zr with oct (O) \",\n  \"lattice_vectors\" : [\n    [ 3.233986860000, 0.000000000000, 0.000000000000 ],\n    [ -1.616993430000, 2.800714770000, 0.000000000000 ],\n    [ -0.000000000000, 0.000000000000, 5.168678340000 ]\n  ],\n  \"title\" : \"ZrO\"\n}\n";
      std::cout << "-------\n\n";

      std::cout << "\n### PRIM ##################\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "PRIM is the input file used by previous version of casm. It can be read and        \n";
      std::cout << "converted to 'prim.json'. The format of PRIM is very similar to the VASP POSCAR    \n";
      std::cout << "except a list of possible occupant molecules is included with each basis site.     \n\n";

      std::cout << "- Molecule names are case sensitive.\n";
      std::cout << "- 'Va' is reserved for vacancies.\n";
      std::cout << "- The default tolerance for checking symmetry is 1e-5, so basis site coordinates\n";
      std::cout << "  should include 6 significant digits or more.\n\n\n";

      std::cout << "EXAMPLE 1: An FCC ternary alloy of elements A, B, and C\n";
      std::cout << "-------\n";
      std::cout <<
                "Face-centered Cubic (FCC, cF)\n\
1.0\n\
0 2.0 2.0\n\
2.0 0 2.0\n\
2.0 2.0 0\n\
1\n\
D\n\
0.00 0.00 0.00 A B C\n";
      std::cout << "-------\n\n";
      std::cout << "EXAMPLE 2: A structure with some alloying sites and some non-alloying sites\n";
      std::cout << "-------\n";
      std::cout <<
                "LiTiO2 - Bronze\n\
 1.00000000\n\
       1.91357600      0.00000000     -6.23799200\n\
       6.08935000     -1.87060000      0.00000000\n\
       0.00000000     -3.74120000      0.00000000\n\
 5 4 8\n\
Direct\n\
   0.0000000   0.0000000   0.0000000 Li Va\n\
   0.3800000   0.9000000   0.5500000 Li Va\n\
   0.6200000   0.1000000   0.4500000 Li Va\n\
   0.0000000   0.2600000   0.3700000 Li Va\n\
   0.0000000   0.7400000   0.6300000 Li Va\n\
   0.7080000   0.3940000   0.8030000 Ti\n\
   0.2920000   0.6060000   0.1970000 Ti\n\
   0.2950000   0.1980000   0.9010000 Ti\n\
   0.7050000   0.8020000   0.0990000 Ti\n\
   0.9960000   0.2640000   0.8680000 O\n\
   0.0040000   0.7360000   0.1320000 O\n\
   0.3470000   0.5280000   0.7360000 O\n\
   0.6530000   0.4720000   0.2640000 O\n\
   0.6290000   0.1200000   0.9400000 O\n\
   0.3710000   0.8800000   0.0600000 O\n\
   0.7070000   0.7240000   0.6380000 O\n\
   0.2930000   0.2760000   0.3620000 O\n";
      std::cout << "-------\n";
      std::cout << std::endl << std::endl;
    }

    if(vm.count("config_list")) {
      std::cout << "\n### config_list.json ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/.casm/config_list.json\n\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "A list of generated configurations. This file is generated at the   \n";
      std::cout << "project level once 'casm enum' has been used to generate            \n";
      std::cout << "configurations.                                                     \n";
      std::cout << "                                                                    \n";
      std::cout << "Contains basic information describing the configuration:            \n\n" <<

                "supercells:supercell_name:configid:                                   \n" <<
                "  Configurations are organized in the JSON file first by the SCELNAME \n" <<
                "  and then listed by the configuration's CONFIGID.                    \n\n" <<

                "source:                                                             \n" <<
                "  Describes the possibly mutiple ways in which this configuration   \n" <<
                "  was generated. Possibilities are \"enumerated\": the configuration\n" <<
                "  was generated when all symmetrically unique configurations of the \n" <<
                "  supercell were enumerated, or \"perturbation\": the configuration \n" <<
                "  was generated as a perturbation with respect to some other        \n" <<
                "  configuration.                                                    \n\n" <<

                "occupation:                                                         \n" <<
                "  An array of int that describes the decoration for a particular    \n" <<
                "  configuration in a given supercell. The length of the array is the\n" <<
                "  same length as your basis, with each value corresponding to a     \n" <<
                "  particular basis site on the supercell. The value for a particular\n" <<
                "  site determines what type of molecule is occupying said site, with\n" <<
                "  values going from zero to the number of allowed occupants for that\n" <<
                "  particular site. The array is ordered in blocks of basis sites of \n" <<
                "  your primitive structure.                                         \n\n" <<

                "  For example, for an \"occupation\" corresponding to a volume 5    \n" <<
                "  supercell of a PRIM with 2 basis sites, the first 5 values in the \n" <<
                "  array represent basis sites for the first site in the PRIM, and   \n" <<
                "  the next 5 represent basis sites for the second site in the PRIM. \n" <<
                "  The occupation values can also be read from the \'config\' column \n" <<
                "  of the config_list file.                                          \n\n" <<

                "... other degrees of freedom will also be included in this file ... \n\n" <<

                "properties:calc:                                                    \n" <<
                "  This contains the values of properties calculated using a         \n" <<
                "  particular calctype. Properties are normalized per unit cell where\n" <<
                "  appropriate (relaxed_energy, volume).                             \n\n" <<

                "properties:ref:                                                     \n" <<
                "  This contains the reference values of properties based on the     \n" <<
                "  values of those properties in the current reference states for    \n" <<
                "  this configuration. The reference values are determined by a      \n" <<
                "  linear interpolation of their values in the reference states to   \n" <<
                "  the composition of this configuration. The most local reference   \n" <<
                "  states are used; see 'casm format --ref_state' for details.       \n\n" <<

                "properties:delta:                                                   \n" <<
                "  This contains the difference between the properties listed in     \n" <<
                "  \'properties:calc\' from those in \'properties:ref\'. These are   \n" <<
                "  the values to be cluster expanded.                                \n\n";
      std::cout << "\n\n\n";

      std::cout << "EXAMPLE:\n";
      std::cout << "-------\n";
      std::cout <<
                "{\n  \"supercells\" : {\n    \"SCEL1_1_1_1_0_0_0\" : {\n      \"0\" : {\n        \"calctype.default\" : {\n          \"ref.default\" : {\n            \"properties\" : {\n              \"calc\" : {\n                \"basis_deformation\" : 0.000000000000,\n                \"data_timestamp\" : 1441172550,\n                \"lattice_deformation\" : 0.000000676576,\n                \"relaxation_strain\" : [ 0.001443293898, 0.001443293305, 0.002332246990, 0.000000000000, 0.000000000000, -0.000000001264 ],\n                \"relaxed_energy\" : -17.093958770000,\n                \"rms_force\" : 0.000000000000,\n                \"volume_relaxation\" : 1.005222845232\n              },\n              \"delta\" : {\n                \"relaxed_energy\" : 0.000000000000\n              },\n              \"ref\" : {\n                \"relaxed_energy\" : -17.093958770000\n              }\n            }\n          }\n        },\n        \"dof\" : {\n          \"occupation\" : [ 0, 0, 0, 0 ]\n        },\n        \"selected\" : false,\n        \"source\" : [ \"occupation_enumeration\" ]\n      },\n\
      ... other configurations ...\n\
    },\n\
    ... other supercells ... \n\
  }\n\
}\n";
      std::cout << "-------\n";
      std::cout << std::endl << std::endl;
    }

    if(vm.count("sym")) {
      std::cout << "\n### sym ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/symmetry/lattice_point_group.json\n";
      std::cout << "$ROOT/symmetry/factor_group.json\n";
      std::cout << "$ROOT/symmetry/crystal_point_group.json\n\n\n";



      std::cout << "DESCRIPTION:\n" <<
                "Symgroup files report each element of a group as representative     \n" <<
                "linear transformations of 3-dimensional space (i.e., the symmetry   \n" <<
                "operation). The symmetery operation transforms a spatial coordinate \n" <<
                "'x' according to x' = A*x+b, where 'A' is the 3x3 'operation matrix'\n" <<
                "and 'b' is the 'shift' vector. Operations are printed either in     \n" <<
                "direct (fractional) coordinates or Cartesian coordinates, depending \n" <<
                "on command-line execution options.\n\n" <<

                "For a crystal or PRIM file, CASM reports the following groups:      \n\n" <<

                "lattice_point_group.json:                                           \n" <<
                "  This is the point group of the Bravais lattice, and is the list of\n" <<
                "  operations that map the lattice vectors onto themselves. The      \n" <<
                "  'shift' vectors will always be zero.                              \n\n" <<

                "factor_group.json:                                                  \n" <<
                "  This is a finite description of the crystal spacegroup, in which  \n" <<
                "  all redundant operations that differ only by a 'shift' are        \n" <<
                "  represented by a single operation, whose 'shift' lies within the  \n" <<
                "  primitive cell. Formally, this is a group formed by the cosets of \n" <<
                "  'T' in 'S', where 'T' is the translation group of the Bravais     \n" <<
                "  lattice and 'S' is the crystal space group.                       \n\n" <<

                "crystal_point_group.json:                                           \n" <<
                "  This is a group of point operations formed by taking the          \n" <<
                "  factor group operations and setting their 'shifts' to zero.       \n" <<
                "  Macroscopic properties of the crystal must exhibit the symmetries \n" <<
                "  of the crystal point group. It is, by definition a subgroup of    \n" <<
                "  the lattice point group.\n";

      std::cout << "\n\n";
    }

    if(vm.count("vasp")) {
      std::cout << "\n### vasp ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n\n";

      std::cout << "INPUT SETTINGS:\n";
      std::cout << "$CALC_SETTINGS/relax.json\n";
      std::cout << "$CALC_SETTINGS/INCAR\n";
      std::cout << "$CALC_SETTINGS/SPECIES\n";
      std::cout << "$CALC_SETTINGS/KPOINTS\n";
      std::cout << "$CALC_SETTINGS/POSCAR\n\n";

      std::cout << "For global settings:\n";
      std::cout << "  CALC_SETTINGS = $ROOT/training_data/settings/$CURR_CALCTYPE\n";
      std::cout << "For supercell specific settings:\n";
      std::cout << "  CALC_SETTINGS = $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE\n";
      std::cout << "For configuration specific settings:\n";
      std::cout << "  CALC_SETTINGS = $ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE\n\n";

      std::cout << "RESULTS:\n";
      std::cout << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/(VASP results)\n";
      std::cout << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/properties.calc.json (read)\n";

      std::cout << "\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "CASM comes with wrappers for using VASP to calculate the properties \n" <<
                "of configurations, but is designed so that any type of calculation  \n" <<
                "software or method could be used if an appropriate set of wrapper   \n" <<
                "scripts are available. By convention, input settings for software   \n" <<
                "used to calculate the properties of a particular configuration      \n" <<
                "should be checked for in the following directories:                 \n" <<
                "  1) $ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE\n" <<
                "  2) $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE          \n" <<
                "  3) $ROOT/training_data/settings/$CURR_CALCTYPE                    \n\n" <<

                "The VASP wrappers included with CASM check for input settings files \n" <<
                "in the above directories, using the most local settings for a       \n" <<
                "particular configuration. In most cases, the global settings files  \n" <<
                "are stored in $ROOT/training_data/settings/$CURR_CALCTYPE and used  \n" <<
                "for all configurations. Settings files are searched for on a file-by-file\n" <<
                "basis, so to set supercell or configuration specific settings it is \n" <<
                "sufficient to only include the particular files necessary in the    \n" <<
                "supercell or configuration level settings folder.                   \n\n" <<

                "PBS job submission using the VASP wrappers depends on using the pbs \n" <<
                "python module available here: https://github.com/prisms-center/pbs  \n\n" <<

                "Included with CASM, the 'vasp.relax' script can be executed by the  \n" <<
                "'casm run' command to submit a batch of VASP jobs that for selected \n" <<
                "configurations. For each selected configuration, VASP is re-run     \n" <<
                "using the output of the previous calculation until full convergence \n" <<
                "is achieved. The convergence criteria is: if the cell shape and     \n" <<
                "volume remain constant (ISIF=0, 1, or 2) then a single calculation  \n" <<
                "is performed; else the calculation is converged if at least 2 jobs  \n" <<
                "are complete, and: 1) the last job completed with <= 3 ionic steps  \n" <<
                " or, if \"nrg_convergence\" is set in the 'relax.json' file, 2) the \n" <<
                "last two calculations had final E0 differ by less than the value of \n" <<
                " \"nrg_convergence\". Once converged, a final constant volume       \n" <<
                "calculation is performed with the following INCAR settings: (ISIF=2,\n" <<
                "ISMEAR=-5, NSW=0, IBRION=-1).                                       \n\n" <<

                "relax.json:                                                         \n" <<
                "  This JSON file contains a single JSON object which contains       \n" <<
                "  parameters used to control PBS job submission settings.           \n" <<
                "  Required keys are:                                                \n" <<
                "    \"queue\": queue to submit job in                               \n" <<
                "    \"ppn\": processors (cores) per node to request                 \n" <<
                "    \"atom_per_proc\": max number of atoms per processor (core)     \n" <<
                "    \"walltime\": walltime to request (ex. \"48:00:00\")            \n\n" <<

                " Optional keys are:                                                 \n" <<
                "    \"account\": account to submit job under (default None)         \n" <<
                "    \"pmem\": string for requested memory (default None)            \n" <<
                "    \"priority\": requested job priority (default \"0\")            \n" <<
                "    \"message\": when to send messages about jobs (ex. \"abe\",     \n" <<
                "                 default \"a\")                                     \n" <<
                "    \"email\": where to send messages (ex. \"me@fake.com\", default \n" <<
                "               None)                                                \n" <<
                "    \"npar\": vasp incar setting (default None)                     \n" <<
                "    \"ncore\": vasp incar setting (default None)                    \n" <<
                "    \"vasp_cmd\": vasp execution command (default is \"vasp\" if    \n" <<
                "                  ncpus=1, else \"mpirun -np {NCPUS} vasp\"         \n" <<
                "    \"ncpus\": number of cpus (cores) to run on (default $PBS_NP)   \n" <<
                "    \"run_limit\": number of vasp runs until \"not_converging\"     \n" <<
                "                   (default 10)                                     \n" <<
                "    \"nrg_convergence\": converged if last two runs complete and    \n" <<
                "                         differ in energy by less than this amount  \n" <<
                "                         (default None)                             \n\n";

      std::cout << "EXAMPLE: relax.json \n";
      std::cout << "-------\n";
      std::cout <<
                "{\n\
  \"account\":\"prismsprojectdebug_flux\",\n\
  \"queue\":\"flux\",\n\
  \"priority\":\"-200\",\n\
  \"walltime\":\"1:00:00\",\n\
  \"pmem\":\"3800mb\",\n\
  \"email\":\"username@univ.edu\",\n\
  \"ppn\":\"16\",\n\
  \"atom_per_proc\":\"2\",\n\
  \"run_limit\":10,\n\
  \"nrg_convergence\":0.002\n\
}\n";
      std::cout << "-------\n\n\n";


      std::cout << "SPECIES:                                                            \n" <<
                "  This file contains information for selecting POTCARs and specifing\n" <<
                "  parameters that must be set on an atom-by-atom basis in the INCAR,\n" <<
                "  such as MAGMOM. The first line in the file specifies the value of \n" <<
                "  'POTCAR_DIR_PATH', which is the base path used to find POTCAR     \n" <<
                "  files. The second line contains column headings (at least 4), and \n" <<
                "  then there are lines for each distinct species. The first column  \n" <<
                "  specifies the 'SPECIES' and must match a species names in the PRIM\n" <<
                "  file. The second column gives an 'ALIAS' name for the species which\n" <<
                "  is used for ordering like atoms in the generated POSCAR files. The\n" <<
                "  third column should be either '0' or '1', such that only one      \n" <<
                "  species with a given ALIAS has a '1'. For that species the fourth \n" <<
                "  column must contain the path that should be appended to the       \n" <<
                "  POTCAR_DIR_PATH to specify the POTCAR file for that species.      \n\n" <<

                "  Additional columns, such as 'MAGMOM' in the example below are     \n\n" <<
                "  and used to specify the value used for a particular species in the\n" <<
                "  INCAR. The column heading must match a valid VASP INCAR setting.  \n\n";

      std::cout << "EXAMPLE: SPECIES \n";
      std::cout << "-------\n" <<
                "POTCAR_DIR_PATH = /absolute/path/to/vasp_potentials\n" <<
                "SPECIES    ALIAS    POTCAR  POTCAR_location    MAGMOM\n" <<
                "Mn3        Mn       0       -                  3\n" <<
                "Mn4        Mn       1       PAW_PBE/Mn         4\n";
      std::cout << "-------\n\n\n";

      std::cout << "INCAR:                                                              \n" <<
                "  This is a template INCAR used for VASP calculations. The settings \n" <<
                "  are generally used as given though some may be automatically set  \n" <<
                "  based on settings in the 'relax.json' or 'SPECIES' files. Also,   \n" <<
                "  some settings might be added or changed if certain errors are     \n" <<
                "  during calculation. The actual INCAR used for each calculation is \n" <<
                "  saved.                                                            \n\n";

      std::cout << "EXAMPLE: INCAR \n";
      std::cout << "-------\n" <<
                "System = Test of VASP submission\n\
ISPIN = 1 # non-spin polarized\n\
PREC = Accurate \n\
IBRION = 2 # conjugate gradient ionic minimization\n\
NSW = 61\n\
ISIF= 3 # relax ions and volume\n\
ENMAX = 400 \n\
ISMEAR = 1 # for metals\n\
SIGMA = 0.2 \n\
LWAVE = .FALSE.\n\
LCHARG = .FALSE.\n";
      std::cout << "-------\n\n\n";

      std::cout << "KPOINTS, POSCAR:                                                    \n" <<
                "  This is a template KPOINTS file used for VASP calculations. If the\n" <<
                "  mode (third line) is set to 'Auto', this file is used as is for all\n" <<
                "  VASP calculations.  Otherwise, if the mode is 'Gamma' or 'M', a   \n" <<
                "  reference POSCAR must also exist and a scaling method is used to  \n" <<
                "  calculate the kpoint mesh for a supercell, such that it has an    \n" <<
                "  equal or greater kpoint density than in the reference POSCAR.     \n\n";

      std::cout << "properties.calc.json:                                               \n" <<
                "  Results of calculations for a particular configuration should be  \n" <<
                "  stored in the directory                                           \n" <<
                "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE,         \n" <<
                "  and calculated properties summarized in the file                  \n" <<
                "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/properties.calc.json \n" <<
                "  The 'properties.calc.json' file is read by CASM to extract the    \n" <<
                "  first-principles calculted properties of interest. If the         \n" <<
                "  'properties.calc.json' file does not exist in the                 \n" <<
                "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE directory\n" <<
                "  CASM assumes that no data is available for that configuration.    \n\n";

      std::cout << "EXAMPLE:\n";
      std::cout << "-------\n";
      std::cout << "{\n\
    \"atom_type\": [\n\
        \"A\", \n\
        \"B\"\n\
    ], \n\
    \"atoms_per_type\": [\n\
        1, \n\
        2\n\
    ], \n\
    \"coord_mode\": \"direct\", \n\
    \"is_complete\": true, \n\
    \"relaxed_basis\": [\n\
        [0.6666667, 0.6666667, 0.6666667],\n\
        [0.00255632, 0.99488736, 0.00255632],\n\
        [0.33077698, 0.33844594, 0.33077698]\n\
    ], \n\
    \"relaxed_energy\": -16.27773537, \n\
    \"relaxed_forces\": [\n\
        [0.0, 0.0, 0.0], \n\
        [0.0, 0.00987362, -0.00987362], \n\
        [0.0, -0.00987362, 0.00987362]\n\
    ], \n\
    \"relaxed_lattice\": [\n\
        [0.0, 1.9174843, 1.9174843], \n\
        [1.61158655, -1.88219884, 3.79968315], \n\
        [3.22317311, 0.0, 0.0]\n\
    ]\n\
}\n";
      std::cout << "-------\n";

    }

    if(vm.count("comp")) {
      std::cout << "\n### composition_axes.json ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/training_data/settings/$CURR_CALCTYPE/$CURR_REF/composition_axes.json\n";
      std::cout << "\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "This JSON file contains the currently selected composition axes, and \n" <<
                "a list of possible standard or custom composition axes.              \n\n" <<

                "standard_axes:                                                      \n" <<
                "  A JSON object containing each possible standard composition axes  \n" <<
                "  as an attribute with its index as the key.                        \n\n" <<

                "custom_axes:                                                        \n" <<
                "  A JSON object containing each custom composition axes as an       \n" <<
                "  attribute with its index as the key. The keys should not be       \n" <<
                "  repeats of any of the standard_axes.                              \n\n" <<

                "standard_axes/composition_axes:components                           \n" <<
                "  A JSON array containing the names of possible species.            \n\n" <<

                "standard_axes/composition_axes:independent_compositions             \n" <<
                "  The number of independent composition axes.                       \n\n" <<

                "standard_axes/composition_axes:origin                               \n" <<
                "  The composition of origin the of composition axes in terms of     \n" <<
                "  number of each component species per primitive cell, ordered as in\n" <<
                "  the 'components' array.                                           \n\n" <<

                "standard_axes/composition_axes:a, b, c, ...                         \n" <<
                "  The composition of end members a, b, c, etc. in terms of number of\n" <<
                "  each component species per primitive cell, ordered as in the      \n" <<
                "  'components' array.                                               \n\n" <<

                "standard_axes/composition_axes:param_formula:                       \n" <<
                "  The formula that converts 'comp_n' (# of each component per       \n" <<
                "  primitive cell) to 'comp' (composition relative the selected      \n" <<
                "  composition axes).                                                \n\n" <<

                "standard_axes/composition_axes:mol_formula:                         \n" <<
                "  The formula that converts 'comp' (composition relative the        \n" <<
                "  selected composition axes) to 'comp_n' (# of each component per   \n" <<
                "  primitive cell).                                                  \n\n\n";


      std::cout << "EXAMPLE:\n";
      std::cout << "-------\n";
      std::cout <<
                "{\n  \"current_axes\" : \"0\",\n  \"standard_axes\" : {\n    \"0\" : {\n      \"a\" : [\n        [ 2.000000000000 ],\n        [ 0.000000000000 ],\n        [ 2.000000000000 ]\n      ],\n      \"components\" : [ \"Zr\", \"Va\", \"O\" ],\n      \"independent_compositions\" : 1,\n      \"mol_formula\" : \"Zr(2)Va(2-2a)O(2a)\",\n      \"origin\" : [\n        [ 2.000000000000 ],\n        [ 2.000000000000 ],\n        [ 0.000000000000 ]\n      ],\n      \"param_formula\" : \"a(0.5-0.25Va+0.25O)\"\n    },\n    \"1\" : {\n      \"a\" : [\n        [ 2.000000000000 ],\n        [ 2.000000000000 ],\n        [ 0.000000000000 ]\n      ],\n      \"components\" : [ \"Zr\", \"Va\", \"O\" ],\n      \"independent_compositions\" : 1,\n      \"mol_formula\" : \"Zr(2)Va(2a)O(2-2a)\",\n      \"origin\" : [\n        [ 2.000000000000 ],\n        [ 0.000000000000 ],\n        [ 2.000000000000 ]\n      ],\n      \"param_formula\" : \"a(0.5+0.25Va-0.25O)\"\n    }\n  }\n}\n";

      std::cout << "-------\n";

    }

    if(vm.count("bspecs")) {
      std::cout << "\n### bspecs.json ##################\n\n";

      std::cout << "LOCATION:\n";
      std::cout << "$ROOT/basis_sets/$CURR_BSET/bspecs.json\n";
      std::cout << "\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "This JSON file contains specifications for generating the cluster\n" <<
                "basis functions.                                                    \n\n";

      std::cout << "The 'site_basis_functions' may be 'occupation' or 'chebychev'.      \n\n";

      std::cout << "The JSON object 'orbit_branch_specs' specifies the maximum size of pair,   \n" <<
                "triplet, quadruplet, etc. clusters in terms of the maximum distance \n" <<
                "between any two sites in the cluster.\n\n";

      std::cout << "The JSON array 'orbit_specs' allows specifying particular custom orbits \n" <<
                "by providing the prototype cluster coordinates. The 'include_subclusters'\n" <<
                "option allows including all orbits of subclusters of the specified cluster\n\n\n";


      std::cout << "EXAMPLE:\n";
      std::cout << "-------\n";
      std::cout <<
                "{\n  \"basis_functions\" : {\n    \"site_basis_functions\" : \"occupation\"\n  },\n  \"orbit_branch_specs\" : {\n    \"2\" : {\"max_length\" : 4.01},\n    \"3\" : {\"max_length\" : 3.01}\n  },\n  \"orbit_specs\" : [\n    {\n      \"coordinate_mode\" : \"Direct\",\n      \"prototype\" : [\n        [ 0.000000000000, 0.000000000000, 0.000000000000 ],\n        [ 1.000000000000, 0.000000000000, 0.000000000000 ],\n        [ 2.000000000000, 0.000000000000, 0.000000000000 ],\n        [ 3.000000000000, 0.000000000000, 0.000000000000 ]\n      ],\n      \"include_subclusters\" : true  \n    },\n    {\n      \"coordinate_mode\" : \"Direct\",\n      \"prototype\" : [\n        [ 0.000000000000, 0.000000000000, 0.000000000000 ],\n        [ 0.000000000000, 1.000000000000, 0.000000000000 ],\n        [ 0.000000000000, 0.000000000000, 1.000000000000 ],\n        [ 1.000000000000, 1.000000000000, 1.000000000000 ]\n      ],\n      \"include_subclusters\" : true\n    }\n  ]\n}\n";
      std::cout << "-------\n";

    }

    if(vm.count("ref_state")) {
      std::cout << "\n### ref_state ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n\n";
      std::cout << "For global reference states:\n";
      std::cout << "  $ROOT/training_data/settings/$CURR_CALCTYPE/$CURR_REF/properties.ref_state.i.json\n";
      std::cout << "For supercell specific reference states:\n";
      std::cout << "  $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE/$CURR_REF/properties.ref_state.i.json\n";
      std::cout << "For configuration specific reference states:\n";
      std::cout << "  $ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE/$CURR_REF/properties.ref_state.i.json\n";
      std::cout << "\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "For a particular configuration (SCELNAME, CONFIGID), the following \n" <<
                "directories are checked for reference state files:                 \n" <<
                "    1) $ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE/$CURR_REF\n" <<
                "    2) $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE/$CURR_REF        \n" <<
                "    3) $ROOT/training_data/$CURR_CALCTYPE/$CURR_REF                        \n\n" <<

                "The most local references found are always used. There must be N     \n" <<
                "reference states set, where N = param_composition.size()+1. See      \n" <<
                "'casm ref' for how to set reference states.                        \n\n" <<

                "The 'properties.ref_state.i.json' JSON files contain 'supercell_name',\n" <<
                "'configid', 'param_composition' identifying the configuration used \n" <<
                "as a reference, and 'ref_state', a JSON object that contains the   \n" <<
                "values of each property in the reference state.                    \n\n" <<

                "If you create custom 'properties.ref_state.i.json' files, only the \n" <<
                "'param_composition' and 'ref_state' members are required.          \n\n\n";


      std::cout << "EXAMPLE: properties.ref_state.0.json\n";
      std::cout << "-------\n";
      std::cout <<
                "{\n\
  \"configid\" : \"0\",\n\
  \"param_composition\" : [ 0.0, 0.0 ],\n\
  \"ref_state\" : {\n\
    \"relaxed_energy\" : -3.699287870000\n\
  },\n\
  \"supercell_name\" : \"SCEL1_1_1_1_0_0_0\"\n\
}\n";
      std::cout << "-------\n";

    }

    if(vm.count("scel")) {
      std::cout << "\n### scel ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/training_data/SCEL\n\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "Contains a list of all the generated and imported supercells. Each  \n" <<
                "entry gives the name of a supercell, its volume in number of        \n" <<
                "primitive volumes, and the transformation matrix to go from the     \n" <<
                "primitive cell to the supercell. The convention is                  \n" <<
                "            LAT.scel = LAT.prim*transf_matrix,                      \n" <<
                "where the columns of the LAT matrices are the lattice vectors.      \n\n\n";

      std::cout << "EXAMPLE:\n";
      std::cout << "-------\n";
      std::cout <<
                "Supercell Name: 'SCEL1_1_1_1_0_0_0' Number: 0 Volume: 1\n\
Supercell Transformation Matrix: \n\
1 0 0\n\
0 1 0\n\
0 0 1\n\
\n\
Supercell Name: 'SCEL2_1_1_2_0_0_0' Number: 1 Volume: 2\n\
Supercell Transformation Matrix: \n\
1 0 -1\n\
0 1 0\n\
0 0 2\n\
\n\
Supercell Name: 'SCEL2_1_2_1_0_0_1' Number: 2 Volume: 2\n\
Supercell Transformation Matrix: \n\
1 -1 0\n\
0 1 -1\n\
0 1 1\n\
...\n";
      std::cout << "-------\n";
    }

    if(vm.count("lat")) {
      std::cout << "\n### lat ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/training_data/$SCELNAME/LAT\n\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "Contains the lattice vectors of a particular supercell of your CASM \n" <<
                "project. The format is the same as the first lines in a standard    \n" <<
                "vasp POSCAR file, excluding the title (scaling followed by the three\n" <<
                "lattice vectors as rows).\n\n\n";

      std::cout << "EXAMPLE:\n";
      std::cout << "-------\n";
      std::cout <<
                " 1.00000000\n\
      10.52850134      0.00000000      0.00000000\n\
      0.00000000     10.52850134      0.00000000\n\
      0.00000000      0.00000000     10.52850134\n";
      std::cout << "-------\n";

    }

    if(vm.count("pos")) {
      std::cout << "\n### pos ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/training_data/$SCELNAME/$CONFIGID/POS\n\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "This file is generated using the '--write-pos' option for 'casm run'.\n";
      std::cout << "Decorated configuration for a particular supercell. It is a         \n" <<
                "supercell of your primitive structure after the enumeration on the  \n" <<
                "alloying sites. The format is standard vasp 5.x format, and the     \n" <<
                "coordinates of the sites for the configuration are the ideal sites  \n" <<
                "specified in the PRIM file.\n\n\n";

      std::cout << "EXAMPLE:\n";
      std::cout << "-------\n";
      std::cout <<
                "SCEL3_1_1_3_0_0_0\n\
1.00000000\n\
      0.00000000      1.75475022      1.75475022\n\
      1.75475022      0.00000000      1.75475022\n\
      3.50950045      3.50950045     -3.50950045\n\
Al Ni\n\
2 1\n\
Direct\n\
   0.3333333   0.3333333   0.3333333\n\
   0.6666667   0.6666667   0.6666667\n\
   0.0000000   0.0000000   0.0000000\n\
\n";
      std::cout << "-------\n";

    }

    if(vm.count("fit")) {
      std::cout << "\n### fit ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "$ROOT/cluster_expansions/clex.formation_energy/$CURR_BSET/$CURR_CALCTYPE/$CURR_REF/$CURR_ECI/energy\n";
      std::cout << "$ROOT/cluster_expansions/clex.formation_energy/$CURR_BSET/$CURR_CALCTYPE/$CURR_REF/$CURR_ECI/corr.in\n";
      std::cout << "$ROOT/cluster_expansions/clex.formation_energy/$CURR_BSET/$CURR_CALCTYPE/$CURR_REF/$CURR_ECI/eci.in\n\n\n";

      std::cout << "DESCRIPTION:\n";
      std::cout << "The 'energy' file contains information about every selected    \n" <<
                "configuration that will be included as training data for fitting ECI.\n\n" <<

                "1st column:                                                        \n" <<
                "  Formation energy determined from the reference states.           \n" <<
                "  (See 'casm ref' and 'casm format --ref_state' for details)       \n\n" <<

                "2nd column:                                                        \n" <<
                "  Weight to be placed on configuration when fitting energies with  \n" <<
                "  eci_search.                                                      \n\n" <<

                "3rd and following columns:                                         \n" <<
                "  Composition of configuration. For a system with N independent    \n" <<
                "  occupants there will be N-1 columns (see 'casm comp')            \n\n" <<

                "2nd column from back:                                              \n" <<
                "  Distance to convex hull. Groundstates will have a value of 0.0000.\n\n" <<

                "Last column:                                                       \n" <<
                "  Path to configuration.                                           \n\n" <<

                "The energy file is to be used together with the corr.in and eci.in \n" <<
                "files to fit the cluster expansion using the eci_search program.   \n\n";

      std::cout << "The 'corr.in' file contains a matrix of correlations for each   \n" <<
                "selected configuration.    \n\n";

      std::cout << "The 'eci.in' file contains a list of calculated correlations and \n" <<
                "can be used to control with correlations are fit by 'eci_search'. \n\n";

      std::cout << "The 'eci.out' file contains the fitted ECI as calculated by 'eci_search'.\n\n";

      std::cout << "EXAMPLE: energy\n";
      std::cout << "-------\n";
      std::cout <<
                "#formation_energy    n/a    n/a    n/a    path    \n\
0.0000000000000  1.000000000000  1.000000000000  0.000000000000  /home/user/science/supercells/SCEL1_1_1_1_0_0_0/0  \n\
0.0000000000000  1.000000000000  0.000000000000  0.000000000000  /home/user/science/supercells/SCEL1_1_1_1_0_0_0/1  \n\
-0.415501770000  1.000000000000  0.500000000000  0.243052905000  /home/user/science/supercells/SCEL2_1_1_2_0_0_0/0  \n\
-0.658554675000  1.000000000000  0.500000000000  0.000000000000  /home/user/science/supercells/SCEL2_1_2_1_0_0_1/0  \n\
-0.307639756667  1.000000000000  0.666666666667  0.131644213333  /home/user/science/supercells/SCEL3_1_1_3_0_0_0/0  \n\
-0.243993753333  1.000000000000  0.333333333333  0.277377812022  /home/user/science/supercells/SCEL3_1_1_3_0_0_0/1  \n\
-0.388569660000  1.000000000000  0.666666666667  0.050714310000  /home/user/science/supercells/SCEL3_1_3_1_0_0_1/0  \n\
-0.444539536667  1.000000000000  0.333333333333  0.076832028688  /home/user/science/supercells/SCEL3_1_3_1_0_0_1/1  \n\
-0.377047050000  1.000000000000  0.666666666667  0.062236920000  /home/user/science/supercells/SCEL3_1_3_1_0_0_2/0  \n";
      std::cout << "-------\n";

    }
    
    if(vm.count("monte")) {
      std::cout << "\n### monte ##################\n\n";

      std::cout << "LOCATION WHEN GENERATED:\n";
      std::cout << "  User determined\n\n\n";
      
      std::cout << "DESCRIPTION:\n";
      std::cout << "  The Monte Carlo input file does not need to be in any particular \n" <<
                   "  location, as long as it is somewhere inside the CASM project     \n" <<
                   "  directory. The input file contains a JSON object with a \"type\" \n" <<
                   "  specifying the type of Monte Carlo calculation to be performed   \n" <<
                   "  and three main subcategories: \"initialization\", \"data\", and  \n" <<
                   "  \"driver\".                                                      \n\n" <<
                   
                   "  Restarts: Monte Carlo calculations that are run in \"incremental\"\n" <<
                   "   drive mode and are stopped before the entire path has been      \n" <<
                   "   calculated can be restarted as long as the input settings do not\n" <<
                   "   change. Upon restart, the settings file is checked for changes, \n" <<
                   "   and the results summary file is checked for the last finished   \n" <<
                   "   conditions. Then the path is resumed from the next set of       \n" <<
                   "   conditions.                                                     \n\n" <<
                   
                   "Input file parameters:                                             \n\n" <<
                   
                   "\"type\" (string):                                                 \n\n" <<
                   
                   "  \"grand_canonical\": Currently the only option.                  \n" <<
                   "    Semi-grand canonical Monte Carlo calculation in which the total\n" <<
                   "    number of sites is fixed, but the occupants on each site may   \n" <<
                   "    vary. One occupant change at a time is attempted.              \n\n\n" <<
                   
                   
                   "\"initialization\": (JSON object) Supercell and parameterization   \n"
                   "    initialization options.                                        \n\n" <<
                   
                   "  /\"clex\", /\"bset\", /\"calctype\", /\"ref\", /\"eci\": (string)\n" <<
                   "    The CASM project settings that should be used for the monte    \n" <<
                   "    carlo calculation.                                             \n\n" <<
                   
                   "  /\"matrix\": (3x3 JSON arrays of integers) The supercell        \n" <<
                   "    transformation matrix.                                         \n\n" <<
                   
                   "  /\"motif\": (JSON object) Specifies the initial occupation of   \n" <<
                   "      the supercell.                                               \n\n" <<
                   
                   "    /\"configname\": (string) The configuration that is tiled to  \n" <<
                   "      fill the supercell (ex. \"SCEL3_3_1_1_0_2_2/0\").            \n" <<
                   "      An error is thrown if the \"matrix\" given is not a supercell\n" <<
                   "      of the specified configuration.                              \n\n\n" <<

                   
                   "\"data\": (JSON object) Data collection options                    \n\n" <<
                   
                   "  /\"sample_by\": (string) Specify unit for the period between    \n" <<
                   "    samples.  May be either \"step\" (sample after every \"sample_period\"\n" <<
                   "    proposed Monte Carlo events), or \"pass\" (sample after the    \n" <<
                   "    \"sample_period\" number of passes), where 1 pass is a number of\n" <<
                   "    steps equal to the number of sites in the supercell that have  \n" <<
                   "    variable occupation).                                          \n\n" <<
                   
                   "  /\"sample_period\": (integer) Specify how many steps or passes  \n" <<
                   "    to wait between data samples.                                  \n\n" <<
                   
                   "  /\"measurements\": (JSON array containing JSON objects)         \n" <<
                   "    Specifies which properties to sample. Each JSON object should  \n" <<
                   "    include \"quantity\" (string) specifying a property to be      \n" <<
                   "    sampled. Optionally, it may also include \"precision\" (number),\n" <<
                   "    indicating the required (absolute) precision in the average of \n" <<
                   "    the quantity for the calculation to be considered converged. If\n" <<
                   "    a precision is given for any quantity, then the Monte Carlo    \n" <<
                   "    calculations run in automatic convergence mode and continue    \n" <<
                   "    until all quantities with a specified precision are converged  \n" <<
                   "    to level requested.                                            \n\n" <<
                   
                   "    Possible options for \"quantity\" are:                         \n" <<
                   "      \"comp\": composition, relative the composition axes         \n" <<
                   "      \"comp_n\": composition, number of atoms per unit cell       \n" <<
                   "      \"site_frac\": composition, normalized per basis site        \n" <<
                   "      \"atom_frac\": composition, normalized per total number of atoms\n" <<
                   "      \"formation_energy\": formation energy (per unit cell)       \n" <<
                   "      \"potential_energy\": potential energy (per unit cell),      \n" <<
                   "        (= formation_energy - sum_i(mu_i*comp_i))                  \n" <<
                   "      \"non_zero_eci_correlations\": correlations (per unit cell)  \n" <<
                   "        which have non-zero eci values.                            \n" <<
                   "      \"all_correlations\": correlations (per unit cell)           \n\n" <<
                   
                   "  /\"confidence\": (number, range (0.0, 1.0), default 0.95) The   \n" <<
                   "    confidence level used for calculating the precision in the     \n" <<
                   "    average value of sampled quantities.                           \n\n" <<
                   
                   "  /\"min_pass\", /\"min_step\", /\"min_sample\": (integer) If in\n" <<
                   "    automatic convergence mode, prevents the calculation from a    \n" <<
                   "    minimum number of passes, steps, or samples have occurred.     \n\n" <<
                   
                   "  /\"max_pass\", /\"max_step\", /\"max_sample\": (integer) If in\n" <<
                   "    automatic convergence mode, stops the calculation if the       \n" <<
                   "    specified number of passes, steps, or samples have occurred.   \n\n" <<
                   
                   "  /\"N_pass\", /\"N_step\", /\"N_sample\": (integer) When not in\n" <<
                   "    automatic convergence mode (no precision has been specified for\n" <<
                   "    any quantities being sampled), stops the calculation when the  \n" <<
                   "    specified number of passes, steps, or samples have occurred.   \n\n" <<
                   
                   "  /\"equilibration_passes_first_run\": (integer) If included, the \n" <<
                   "    requested number of passes will be performed at the initial    \n" <<
                   "    conditions as a preliminary step before the actual run begins. \n" <<
                   "    This may be useful when not running in automatic convergence   \n" <<
                   "    mode.                                                          \n\n" <<
                   
                   "  /\"equilibration_passes_each_run\": (integer) If included, the \n" <<
                   "    requested number of passes will be performed at each condition \n" <<
                   "    as a preliminary step before the actual run begins. This may be\n" <<
                   "    useful when not running in automatic convergence mode.         \n\n" <<
                   
                   "  /\"storage\": (JSON object) Options for writing results.        \n\n" <<
                   
                   "    /\"output_format\": (string or JSON array of string) Specifies\n" <<
                   "      the type or types of output files. Current options are \"csv\"\n" <<
                   "      or \"json\". Type names with either all lower case or all    \n" <<
                   "      upper case are accepted.                                     \n\n" <<
                   
                   "    /\"write_observations\": (boolean, default false) If true,    \n" <<
                   "      all individual observations of the quantities requested to be\n" <<
                   "      sampled will be written to compressed files:                 \n" <<
                   "        \"output_directory\"/conditions.i/observations.ext.gz      \n" <<
                   "      where 'i' is the condition index and 'ext' is the output     \n" <<
                   "      format.                                                      \n\n" <<
                   
                   "    /\"write_trajectory\": (boolean, default false) If true,      \n" <<
                   "      the value of all degrees of freedom at the time of each      \n" <<
                   "      sample will be written to compressed files:                  \n" <<
                   "        \"output_directory\"/conditions.i/trajectory.ext.gz        \n" <<
                   "      where 'i' is the condition index and 'ext' is the output     \n" <<
                   "      format.                                                      \n\n\n" <<
                   
                   
                   "\"driver\": (JSON object) Contains options controlling the         \n" <<
                   "    conditions at which calculations are performed.                \n\n" <<
                   
                   "  /\"mode\": (string) Specify the drive mode.                     \n\n" <<
                   
                   "    Possible options for \"mode\" are:                             \n" <<
                   "      \"single\": perform one calculation at the initial conditions\n" <<
                   "      \"incremental\": perform several calculations, starting at the\n" <<
                   "        initial conditions and incrementing by the incremental     \n" <<
                   "        conditions up to (and including) the final conditions.     \n\n" <<
                   
                   "  /\"initial_conditions\",\n" <<
                   "  /\"incremental_conditions\", \n" <<
                   "  /\"final_conditions\": \n" <<
                   "    (JSON object) Specifies the applied conditions for the         \n" <<
                   "    calculation. For \"incremental_conditions\", specifies the     \n" <<
                   "    change in conditions between individual calculations. Each JSON\n" <<
                   "    object should include:                                         \n\n" <<
                   
                   "    /\"temperature\": (number) The temperature in K.              \n\n" <<
                   
                   "    /\"mu\" (JSON object): The chemical potential(s)              \n\n" <<
                   
                   "      /\"a\", /\"b\", ...: (number) Each chemical potential. For \n" <<
                   "      example, \"a\" specifies dG/dcomp(a) and \"b\" specifies     \n" <<
                   "      dG/dcomp(b). The number of chemical potentials provided must \n" <<
                   "      match the number of independent compositions.                \n\n" <<
                   
                   "    /\"tolerance\": (number) For \"incremental\" drive mode,      \n" <<
                   "      specifies a numerical tolerance for comparing conditions.    \n\n\n";

      std::cout << "EXAMPLE: Settings for an incremental calculation with increasing temperature in manual convergence mode.\n";
      std::cout << "-------\n";
      std::cout << "{\n  \"comment\" : \"This is a sample input file. Unrecognized attributes (like this one) are ignored.\",\n  \"type\" : \"grand_canonical\",\n  \"initialization\" : {\n    \"clex\" : \"formation_energy\",\n    \"bset\" : \"default\",\n    \"calctype\" : \"default\",\n    \"ref\" : \"default\",\n    \"eci\" : \"default\",\n    \"matrix\" : [\n      [9, 0, 0],\n      [0, 9, 0],\n      [0, 0, 9]\n    ],\n    \"motif\" : {\n      \"configname\" : \"SCEL3_3_1_1_0_2_2/0\"\n    }\n  },\n  \"data\" : {\n    \"sample_by\" : \"pass\",\n    \"sample_period\" : 10,\n    \"equilibration_passes_each_run\" : 1000,\n    \"N_pass\" : 10000,\n    \"measurements\" : [ \n      { \n        \"quantity\" : \"formation_energy\"\n      },\n      { \n        \"quantity\" : \"potential_energy\",\n        \"precision\" : 1e-3\n      },\n      { \n        \"quantity\" : \"comp\",\n        \"precision\" : 1e-3\n      },\n      { \n        \"quantity\" : \"comp_n\"\n      },\n      { \n        \"quantity\" : \"all_correlations\"\n      }\n    ],\n    \"storage\" : {\n      \"output_format\" : \"json\"\n    }\n  },\n  \"driver\" : {\n    \"mode\" : \"incremental\",\n    \"initial_conditions\" : {\n      \"mu\" : {\n        \"a\" : -1.75\n      },\n      \"temperature\" : 200.0,\n      \"tolerance\" : 0.001\n    },\n    \"final_conditions\" : {\n      \"mu\" : {\n        \"a\" : -1.75\n      },\n      \"temperature\" : 800.0,\n      \"tolerance\" : 0.001\n    },\n    \"incremental_conditions\" : {\n      \"mu\" : {\n        \"a\" : 0.00\n      },\n      \"temperature\" : 10.0,\n      \"tolerance\" : 0.001\n    }\n  }\n}\n";
      std::cout << "-------\n\n";
      
      std::cout << "EXAMPLE: Settings for an incremental calculation with increasing mu in automatic convergence mode.\n";
      std::cout << "-------\n";
      std::cout << "{\n  \"comment\" : \"This is a sample input file. Unrecognized attributes (like this one) are ignored.\",\n  \"type\" : \"grand_canonical\",\n  \"initialization\" : {\n    \"clex\" : \"formation_energy\",\n    \"bset\" : \"default\",\n    \"calctype\" : \"default\",\n    \"ref\" : \"default\",\n    \"eci\" : \"default\",\n    \"matrix\" : [\n      [9, 0, 0],\n      [0, 9, 0],\n      [0, 0, 9]\n    ],\n    \"motif\" : {\n      \"configname\" : \"SCEL3_3_1_1_0_2_2/0\"\n    }\n  },\n  \"data\" : {\n    \"sample_by\" : \"pass\",\n    \"sample_period\" : 1,\n    \"max_pass\" : 1000000,\n    \"measurements\" : [ \n      { \n        \"quantity\" : \"formation_energy\"\n      },\n      { \n        \"quantity\" : \"potential_energy\",\n        \"precision\" : 1e-3\n      },\n      { \n        \"quantity\" : \"atom_frac\"\n      },\n      { \n        \"quantity\" : \"site_frac\"\n      },\n      { \n        \"quantity\" : \"comp\",\n        \"precision\" : 1e-3\n      },\n      { \n        \"quantity\" : \"comp_n\"\n      },\n      { \n        \"quantity\" : \"non_zero_eci_correlations\"\n      }\n    ],\n    \"storage\" : {\n      \"write_observations\" : true,\n      \"write_trajectory\" : true,\n      \"output_format\" : \"csv\"\n    }\n  },\n  \"driver\" : {\n    \"mode\" : \"incremental\",\n    \"initial_conditions\" : {\n      \"mu\" : {\n        \"a\" : -1.75\n      },\n      \"temperature\" : 800.0,\n      \"tolerance\" : 0.001\n    },\n    \"final_conditions\" : {\n      \"mu\" : {\n        \"a\" : -1.00\n      },\n      \"temperature\" : 800.0,\n      \"tolerance\" : 0.001\n    },\n    \"incremental_conditions\" : {\n      \"mu\" : {\n        \"a\" : 0.01\n      },\n      \"temperature\" : 0.0,\n      \"tolerance\" : 0.001\n    }\n  }\n}\n";
      std::cout << "-------\n";

    }

    return 0;
  }

}
