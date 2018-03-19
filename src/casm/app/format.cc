#include <cstring>

#include "casm/CASM_global_definitions.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {


  namespace Completer {
    FormatOption::FormatOption(): OptionHandlerBase("format") {}

    void FormatOption::initialize() {
      add_help_suboption();

      m_desc.add_options()
      ("dir,d", "CASM project directory structure summary")
      ("project_settings", "Description and location of 'project_settings' file")
      ("prim", "Description and location of 'prim.json' and 'PRIM' files")
      ("config_list", "Description and location of 'config_list.json' file")
      ("sym", "Description and location of 'lattice_point_group.json', 'factor_group.json' and 'crystal_point_group.json' files")
      ("vasp", "Description and location of VASP settings files")
      ("properties", "Description and location of properties.calc.json files")
      ("qe", "Description and location of Quantum Espresso settings files")
      ("comp", "Description and location of 'composition_axes.json' file")
      ("bspecs", "Description and location of 'bspecs.json' file")
      ("clust", "Description and location of 'clust.json' file")
      ("basis", "Description and location of 'basis.json' file")
      ("clex", "Description and location of '$TITLE_Clexulator.*' files")
      ("ref", "Description and location of 'chemical_reference.json' files")
      ("scel", "Description and location of 'SCEL' file")
      ("lat", "Description and location of 'LAT' files")
      ("pos", "Description and location of 'POS' files")
      ("eci", "Description and location of 'eci.json' file")
      ("monte", "Description and location of the Monte Carlo input file");
      return;
    }
  }

  // ///////////////////////////////////////
  // 'format' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int format_command(const CommandArgs &args) {

    po::variables_map vm;

    Completer::FormatOption format_opt;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, format_opt.desc()), vm);

      /** --help option
       */
      if(vm.count("help") || vm.size() == 0) {
        args.log << std::endl;
        args.log << format_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << format_opt.desc() << std::endl;

        args.log << "DESCRIPTION" << std::endl;
        args.log << "    This option describes the files contained within a CASM project \n";
        args.log << "    and where to find them. For a summary of the directory structure\n";
        args.log << "    of a CASM project using VASP for calculating configuration use  \n";
        args.log << "    the --dir option. Not all files are always present.             \n";

        return 0;
      }

      po::notify(vm);

    }
    catch(po::error &e) {
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log << format_opt.desc() << std::endl;
      return 1;
    }
    catch(std::exception &e) {
      args.err_log << "Unhandled Exception reached the top of main: "
                   << e.what() << ", application will now exit" << std::endl;
      return 1;

    }

    if(vm.count("dir")) {
      args.log << "\n### dir ##################\n\n";

      args.log << "  The expected CASM project directory structure with VASP settings  \n";
      args.log << "  files:                                                            \n";
      args.log << "                                                                    \n";
      args.log << "    $ROOT/                                                          \n";
      args.log << "      prim.json                                                     \n";
      args.log << "      (PRIM)                                                        \n";
      args.log << "      LOG                                                           \n";
      args.log << "    $ROOT/.casm                                                     \n";
      args.log << "      project_settings.json                                         \n";
      args.log << "      config_list.json                                              \n";
      args.log << "      composition_axes.json                                         \n";
      args.log << "    $ROOT/symmetry/                                                 \n";
      args.log << "      lattice_point_group.json                                      \n";
      args.log << "      factor_group.json                                             \n";
      args.log << "      crystal_point_group.json                                      \n";
      args.log << "    $ROOT/basis_sets/$CURR_BSET/                                    \n";
      args.log << "      bspecs.json                                                   \n";
      args.log << "      basis.json                                                    \n";
      args.log << "      clust.json                                                    \n";
      args.log << "      $TITLE_Clexulator.*                                           \n";
      args.log << "    $ROOT/training_data/                                            \n";
      args.log << "      SCEL                                                          \n";
      args.log << "    $ROOT/training_data/settings/$CURR_CALCTYPE/                    \n";
      args.log << "      relax.json                                                    \n";
      args.log << "      INCAR                                                         \n";
      args.log << "      SPECIES                                                       \n";
      args.log << "      KPOINTS                                                       \n";
      args.log << "      POSCAR                                                        \n";
      args.log << "    $ROOT/training_data/settings/$CURR_CALCTYPE/$CURR_REF/          \n";
      args.log << "      chemical_reference.json                                       \n";
      args.log << "    $ROOT/training_data/$SCELNAME/                                  \n";
      args.log << "      LAT                                                           \n";
      args.log << "    $ROOT/training_data/$SCELNAME/$CONFIGID                         \n";
      args.log << "      POS                                                           \n";
      args.log << "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE          \n";
      args.log << "      (VASP/QE results)                                             \n";
      args.log << "      properties.calc.json                                          \n";
      args.log << "    $ROOT/cluster_expansions/clex.formation_energy/$CURR_BSET/$CURR_CALCTYPE/$CURR_REF/$CURR_ECI\n";
      args.log << "      eci.json                                                      \n";
      args.log << " \n";
      args.log << " \n";
      args.log << "    Variable descriptions:                                          \n";
      args.log << " \n";
      args.log << "    $ROOT: root directory of the CASM project                       \n";
      args.log << " \n";
      args.log << "    $CURR_BSET: Current basis set, by default this is 'bset.default'.\n";
      args.log << "    The current value can be inspected via 'casm settings -l'.      \n";
      args.log << " \n";
      args.log << "    $CURR_CALCTYPE: Current calctype, by default this is 'calctype.default'.\n";
      args.log << "    The current value can be inspected via 'casm settings -l'.      \n";
      args.log << " \n";
      args.log << "    $CURR_REF: Current composition axes and reference states, by    \n";
      args.log << "    default this is 'ref.default'. The current value can be inspected\n";
      args.log << "    via 'casm settings -l'.                                         \n";
      args.log << " \n";
      args.log << "    $SCELNAME: Supercell name, in the form SCELV_A_B_C_D_E_F. 'V' is\n";
      args.log << "    volume of the supercell in terms of the primitive cell, and     \n";
      args.log << "    'A'-'F' are the values of the hermite normal form of the        \n";
      args.log << "    transformation matrix.                                          \n";
      args.log << " \n";
      args.log << "    $CONFIGID: Configuration id, a unique integer.                  \n";
      args.log << " \n";
      args.log << "    $TITLE: Title of the CASM project                               \n";
      args.log << "\n";
      args.log << "    Note: The 'settings' heirarchy can be located at the project    \n";
      args.log << "    level as shown above, or at the supercell or configuration level\n";
      args.log << "    in order to override calculation, composition, or reference     \n";
      args.log << "    state settings at the supercell or configuration level.  The    \n";
      args.log << "    most local settings are always used for a configuration.        \n";
      args.log << " \n";

    }

    if(vm.count("project_settings")) {
      args.log << "\n### project_settings.json ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "$ROOT/.casm/project_settings.json\n\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "Current CASM project settings.\n\n\n";

      args.log << "EXAMPLE:\n";
      args.log << "-------\n";
      args.log <<
               "{\n  \"cluster_expansions\" : {\n    \"formation_energy\" : {\n      \"bset\" : \"default\",\n      \"calctype\" : \"default\",\n      \"eci\" : \"default\",\n      \"name\" : \"formation_energy\",\n      \"property\" : \"formation_energy\",\n      \"ref\" : \"default\"\n    }\n  },\n  \"crystallography_tol\" : 1.000000000000000082e-05,\n  \"curr_properties\" : [ \"relaxed_energy\" ],\n  \"default_clex\" : \"formation_energy\",\n  \"lin_alg_tol\" : 1.000000000000000036e-10,\n  \"name\" : \"ZrO\",\n  \"nlist_sublat_indices\" : [ 2, 3 ],\n  \"nlist_weight_matrix\" : [\n    [ 2, -1, 0 ],\n    [ -1, 2, 0 ],\n    [ 0, 0, 5 ]\n  ],\n  \"query_alias\" : {\n  },\n  \"view_command\" : \"casm.view \\\"open -a /Applications/VESTA/VESTA.app\\\"\"\n}" << std::endl;
      args.log << "-------\n";
      args.log << std::endl << std::endl;
    }

    if(vm.count("prim")) {
      args.log << "\n### prim.json ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "$ROOT/prim.json\n";
      args.log << "$ROOT/PRIM (legacy)\n\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "'prim.json' describes the primitive cell structure. It includes the lattice \n";
      args.log << "vectors, crystal basis sites and a list of possible occupant molecules on each\n";
      args.log << "basis site.\n\n";

      args.log <<  "'prim.json' parameters:                                            \n\n"

               "\"title\" (string):                                                \n"
               "  A title for the project. Must consist of alphanumeric characters \n"
               "  and underscores only. The first character may not be a number.   \n\n"

               "\"lattice_vectors\" (JSON array of 3 JSON arrays of 3 numbers):    \n"
               "  Lattice vectors for the primitive structure, in Angstroms.       \n\n"

               "\"coordinate_mode\" (string):                                      \n"
               "  Coordinate mode for basis sites. One of:                         \n"
               "    \"Fractional\" or \"Direct\",                                  \n"
               "    \"Cartesian\"                                                  \n\n"

               "\"basis\" (JSON array of JSON objects):                            \n\n"

               "  /\"coordinate\" (JSON array of 3 numbers):                       \n"
               "    Coordinate of the basis site with units as specified by the    \n"
               "    the \"coordinate_mode\" parameter. The default tolerance for   \n"
               "    checking symmetry is 1e-5, so basis site coordinates should    \n"
               "    include 6 significant digits or more.                          \n"

               "  /\"occupant_dof\" (JSON array of string):                        \n"
               "    A list of the possible occupant atoms (and in future versions  \n"
               "    CASM, molecules) that on each site. The names are case         \n"
               "    sensitive, and \"Va\" is reserved for vacancies.               \n\n\n";


      args.log << "EXAMPLE 1: An FCC ternary alloy of elements A, B, and C\n";
      args.log << "-------\n";
      args.log <<
               "{\n  \"basis\" : [\n    {\n      \"coordinate\" : [ 0.000000000000, 0.000000000000, 0.000000000000 ],\n      \"occupant_dof\" : [ \"A\", \"B\", \"C\" ]\n    }\n  ],\n  \"coordinate_mode\" : \"Fractional\",\n  \"description\" : \"Face-centered Cubic (FCC, cF)\",\n  \"lattice_vectors\" : [\n    [ 2.000000000000, 2.000000000000, 0.000000000000 ],\n    [ 0.000000000000, 2.000000000000, 2.000000000000 ],\n    [ 2.000000000000, 0.000000000000, 2.000000000000 ]\n  ],\n  \"title\" : \"ABC\"\n}\n";
      args.log << "-------\n\n";

      args.log << "EXAMPLE 2: HCP Zr with O in octahedral interstitial positions\n";
      args.log << "-------\n";
      args.log <<
               "\n{\n  \"basis\" : [\n    {\n      \"coordinate\" : [ 0.0, 0.0, 0.0 ],\n      \"occupant_dof\" : [ \"Zr\" ]\n    },\n    {\n      \"coordinate\" : [ 0.666666, 0.333333, 0.5 ],\n      \"occupant_dof\" : [ \"Zr\" ]\n    },\n    {\n      \"coordinate\" : [ 0.333333, 0.666666, 0.25 ],\n      \"occupant_dof\" : [ \"Va\", \"O\" ]\n    },\n    {\n      \"coordinate\" : [ 0.333333, 0.666666, 0.75 ],\n      \"occupant_dof\" : [ \"Va\", \"O\" ]\n    }\n  ],\n  \"coordinate_mode\" : \"Fractional\",\n  \"description\" : \"hcp Zr with oct (O) \",\n  \"lattice_vectors\" : [\n    [ 3.233986860000, 0.000000000000, 0.000000000000 ],\n    [ -1.616993430000, 2.800714770000, 0.000000000000 ],\n    [ -0.000000000000, 0.000000000000, 5.168678340000 ]\n  ],\n  \"title\" : \"ZrO\"\n}\n";
      args.log << "-------\n\n";

      args.log << "\n### PRIM ##################\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "PRIM is the input file used by previous version of casm. It can be read and        \n";
      args.log << "converted to 'prim.json'. The format of PRIM is very similar to the VASP POSCAR    \n";
      args.log << "except a list of possible occupant molecules is included with each basis site.     \n\n";

      args.log << "- Molecule names are case sensitive.\n";
      args.log << "- 'Va' is reserved for vacancies.\n";
      args.log << "- The default tolerance for checking symmetry is 1e-5, so basis site coordinates\n";
      args.log << "  should include 6 significant digits or more.\n\n\n";

      args.log << "EXAMPLE 1: An FCC ternary alloy of elements A, B, and C\n";
      args.log << "-------\n";
      args.log <<
               "Face-centered Cubic (FCC, cF)\n\
1.0\n\
0 2.0 2.0\n\
2.0 0 2.0\n\
2.0 2.0 0\n\
1\n\
D\n\
0.00 0.00 0.00 A B C\n";
      args.log << "-------\n\n";
      args.log << "EXAMPLE 2: A structure with some alloying sites and some non-alloying sites\n";
      args.log << "-------\n";
      args.log <<
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
      args.log << "-------\n";
      args.log << std::endl << std::endl;
    }

    if(vm.count("config_list")) {
      args.log << "\n### config_list.json ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "$ROOT/.casm/config_list.json\n\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "A list of generated configurations. This file is generated at the   \n";
      args.log << "project level once 'casm enum' has been used to generate            \n";
      args.log << "configurations.                                                     \n";
      args.log << "                                                                    \n";
      args.log << "Contains basic information describing the configuration:            \n\n" <<

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
               "  states are used; see 'casm format --ref' for details.             \n\n" <<

               "properties:delta:                                                   \n" <<
               "  This contains the difference between the properties listed in     \n" <<
               "  \'properties:calc\' from those in \'properties:ref\'. These are   \n" <<
               "  the values to be cluster expanded.                                \n\n";
      args.log << "\n\n\n";

      args.log << "EXAMPLE:\n";
      args.log << "-------\n";
      args.log <<
               "{\n  \"supercells\" : {\n    \"SCEL1_1_1_1_0_0_0\" : {\n      \"0\" : {\n        \"calctype.default\" : {\n          \"ref.default\" : {\n            \"properties\" : {\n              \"calc\" : {\n                \"basis_deformation\" : 0.000000000000,\n                \"data_timestamp\" : 1441172550,\n                \"lattice_deformation\" : 0.000000676576,\n                \"relaxation_strain\" : [ 0.001443293898, 0.001443293305, 0.002332246990, 0.000000000000, 0.000000000000, -0.000000001264 ],\n                \"relaxed_energy\" : -17.093958770000,\n                \"rms_force\" : 0.000000000000,\n                \"relaxed_magmom\" : -4.125372200000,\n                \"volume_relaxation\" : 1.005222845232\n              },\n              \"delta\" : {\n                \"relaxed_energy\" : 0.000000000000\n              },\n              \"ref\" : {\n                \"relaxed_energy\" : -17.093958770000\n              }\n            }\n          }\n        },\n        \"dof\" : {\n          \"occupation\" : [ 0, 0, 0, 0 ]\n        },\n        \"selected\" : false,\n        \"source\" : [ \"occupation_enumeration\" ]\n      },\n\
      ... other configurations ...\n\
    },\n\
    ... other supercells ... \n\
  }\n\
}\n";
      args.log << "-------\n";
      args.log << std::endl << std::endl;
    }

    if(vm.count("sym")) {
      args.log << "\n### sym ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "$ROOT/symmetry/lattice_point_group.json\n";
      args.log << "$ROOT/symmetry/factor_group.json\n";
      args.log << "$ROOT/symmetry/crystal_point_group.json\n\n\n";



      args.log << "DESCRIPTION:\n" <<
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

      args.log << "\n\n";
    }

    if(vm.count("vasp")) {
      args.log << "\n### vasp ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n\n";

      args.log << "INPUT SETTINGS:\n";
      args.log << "$CALC_SETTINGS/relax.json\n";
      args.log << "$CALC_SETTINGS/INCAR\n";
      args.log << "$CALC_SETTINGS/SPECIES\n";
      args.log << "$CALC_SETTINGS/KPOINTS\n";
      args.log << "$CALC_SETTINGS/POSCAR\n\n";

      args.log << "For global settings:\n";
      args.log << "  CALC_SETTINGS = $ROOT/training_data/settings/$CURR_CALCTYPE\n";
      args.log << "For supercell specific settings:\n";
      args.log << "  CALC_SETTINGS = $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE\n";
      args.log << "For configuration specific settings:\n";
      args.log << "  CALC_SETTINGS = $ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE\n\n";

      args.log << "RESULTS:\n";
      args.log << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/(VASP results)\n";
      args.log << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/properties.calc.json (read)\n";

      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "CASM comes with wrappers for using VASP to calculate the properties \n" <<
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
               "               default \"a\")                                       \n" <<
               "    \"email\": where to send messages (ex. \"me@fake.com\", default \n" <<
               "             None)                                                  \n" <<
               "    \"qos\": quality of service, 'qos' option (ex. \"fluxoe\")      \n" <<
               "    \"npar\": vasp incar setting (default None)                     \n" <<
               "    \"ncore\": vasp incar setting (default None)                    \n" <<
               "    \"kpar\": vasp incar setting (default None)                     \n" <<
               "    \"vasp_cmd\": vasp execution command (default is \"vasp\" if    \n" <<
               "                ncpus=1, else \"mpirun -np {NCPUS} vasp\"           \n" <<
               "    \"ncpus\": number of cpus (cores) to run on (default $PBS_NP)   \n" <<
               "    \"run_limit\": number of vasp runs until \"not_converging\"     \n" <<
               "                 (default 10)                                       \n" <<
               "    \"nrg_convergence\": converged if last two runs complete and    \n" <<
               "                       differ in energy by less than this amount    \n" <<
               "                       (default None)                               \n" <<
               "    \"move\": files to move at the end of a run (ex. \"POTCAR\",    \n" <<
               "            \"WAVECAR\"], default [\"POTCAR\"])                     \n" <<
               "    \"copy\": files to copy from run to run (ex. [\"INCAR\",        \n" <<
               "            \"KPOINTS\"], default [\"INCAR, KPOINTS\"])             \n" <<
               "    \"remove\": files to remove at the end of a run (ex. [\"IBZKPT\",\n" <<
               "              \"CHGCAR\"], default [\"IBKZPT\", \"CHG\", \"CHGCAR\",\n" <<
               "              \"WAVECAR\", \"TMPCAR\", \"EIGENVAL\", \"DOSCAR\",    \n" <<
               "              \"PROCAR\", \"PCDAT\", \"XDATCAR\", \"LOCPOT\", \"ELFCAR\",\n" <<
               "              \"PROOUT\"]                                           \n" <<
               "    \"compress\": files to compress at the end of a run (ex.        \n" <<
               "                [\"OUTCAR\", \"vasprun.xml\"], default [])          \n" <<
               "    \"backup\": files to compress to backups at the end of a run,   \n" <<
               "              used in conjunction with move (ex. [\"WAVECAR\"])     \n" <<
               "    \"encut\": [START, STOP, STEP] values for converging ENCUT to   \n" <<
               "             within nrg_convergence (ex. [\"450\", \"Auto\",        \n" <<
               "             \"10\"], default [\"Auto\", \"Auto\", \"10\"] where    \n" <<
               "             \"Auto\" is either the largest ENMAX in all POTCARS    \n" <<
               "             called in SPECIES for START, or 2.0 * largest ENMAX    \n" <<
               "             for STOP)                                              \n" <<
               "    \"kpoints\": [start, stop, step] values for converging KPOINTS  \n" <<
               "               to within nrg_convergence (ex. [\"5\", \"50\", \"1\"],\n" <<
               "               default [\"5\", \"Auto\", \"1\"] where \"Auto\" can  \n" <<
               "               only be used for STOP and means to keep increasing   \n" <<
               "               the KPOINTS length by STEP until either              \n" <<
               "               nrg_convergence or walltime is reached). For meaning \n" <<
               "               of the KPOINTS length parameter see the VASP         \n" <<
               "               documentation at http://cms.mpi.univie.ac.at/vasp/   \n" <<
               "               vasp/Automatic_k_mesh_generation.html                \n" <<
               "    \"extra_input_files\": extra input files to be copied from the  \n" <<
               "                         settings directory, e.g., a vdW kernel     \n" <<
               "                         file.                                      \n" <<
               "    \"initial\": location of INCAR with tags for the initial run,   \n" <<
               "               if desired (e.g. to generate a PBE WAVECAR for use   \n" <<
               "               with M06-L)                                          \n" <<
               "    \"final\": location of INCAR with tags for the final run, if    \n" <<
               "             desired (e.g. \"ISMEAR = -5\", etc). Otherwise, the    \n" <<
               "             settings enforced are (\"ISMEAR = -5\", \"NSW = 0\",   \n" <<
               "              \"IBRION = -1\", \"ISIF = 2\")                        \n" <<
               "    \"err_types\": list of errors to check for. Allowed entries are \n" <<
               "                 \"IbzkptError\" and \"SubSpaceMatrixError\".       \n" <<
               "                 Default: [\"SubSpaceMatrixError\"]                 \n" <<
               "\n";

      args.log << "EXAMPLE: relax.json \n";
      args.log << "-------\n";
      args.log <<
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
      args.log << "-------\n\n\n";


      args.log << "SPECIES:                                                            \n" <<
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

      args.log << "EXAMPLE: SPECIES \n";
      args.log << "-------\n" <<
               "POTCAR_DIR_PATH = /absolute/path/to/vasp_potentials\n" <<
               "SPECIES    ALIAS    POTCAR  POTCAR_location    MAGMOM\n" <<
               "Mn3        Mn       0       -                  3\n" <<
               "Mn4        Mn       1       PAW_PBE/Mn         4\n";
      args.log << "-------\n\n\n";

      args.log << "INCAR:                                                              \n" <<
               "  This is a template INCAR used for VASP calculations. The settings \n" <<
               "  are generally used as given though some may be automatically set  \n" <<
               "  based on settings in the 'relax.json' or 'SPECIES' files. Also,   \n" <<
               "  some settings might be added or changed if certain errors are     \n" <<
               "  during calculation. The actual INCAR used for each calculation is \n" <<
               "  saved.                                                            \n\n";

      args.log << "EXAMPLE: INCAR \n";
      args.log << "-------\n" <<
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
      args.log << "-------\n\n\n";

      args.log << "KPOINTS, POSCAR:                                                    \n" <<
               "  This is a template KPOINTS file used for VASP calculations. If the\n" <<
               "  mode (third line) is set to 'Auto', this file is used as is for all\n" <<
               "  VASP calculations.  Otherwise, if the mode is 'Gamma' or 'M', a   \n" <<
               "  reference POSCAR must also exist and a scaling method is used to  \n" <<
               "  calculate the kpoint mesh for a supercell, such that it has an    \n" <<
               "  equal or greater kpoint density than in the reference POSCAR.     \n\n";

      args.log << "-------\n";

    }

    if(vm.count("qe")) {
      args.log << "\n### quantum espresso ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n\n";

      args.log << "INPUT SETTINGS:\n";
      args.log << "$CALC_SETTINGS/relax.json\n";
      args.log << "$CALC_SETTINGS/$CUSTOM_INFILE_NAME\n";
      args.log << "$CALC_SETTINGS/SPECIES\n";

      args.log << "For global settings:\n";
      args.log << "  CALC_SETTINGS = $ROOT/training_data/settings/$CURR_CALCTYPE\n";
      args.log << "For supercell specific settings:\n";
      args.log << "  CALC_SETTINGS = $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE\n";
      args.log << "For configuration specific settings:\n";
      args.log << "  CALC_SETTINGS = $ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE\n\n";

      args.log << "RESULTS:\n";
      args.log << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/(quantum espresso results)\n";
      args.log << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/properties.calc.json (read)\n";

      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "CASM comes with wrappers for using Quantum Espresso to calculate the properties \n" <<
               "of configurations, but is designed so that any type of calculation  \n" <<
               "software or method could be used if an appropriate set of wrapper   \n" <<
               "scripts are available. By convention, input settings for software   \n" <<
               "used to calculate the properties of a particular configuration      \n" <<
               "should be checked for in the following directories:                 \n" <<
               "  1) $ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE\n" <<
               "  2) $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE          \n" <<
               "  3) $ROOT/training_data/settings/$CURR_CALCTYPE                    \n\n" <<

               "The Quantum Espresso wrappers included with CASM check for input settings files \n" <<
               "in the above directories, using the most local settings for a       \n" <<
               "particular configuration. In most cases, the global settings files  \n" <<
               "are stored in $ROOT/training_data/settings/$CURR_CALCTYPE and used  \n" <<
               "for all configurations. Settings files are searched for on a file-by-file\n" <<
               "basis, so to set supercell or configuration specific settings it is \n" <<
               "sufficient to only include the particular files necessary in the    \n" <<
               "supercell or configuration level settings folder.                   \n\n" <<

               "PBS job submission using the Quantum Espresso wrappers depends on using the pbs \n" <<
               "python module available here: https://github.com/prisms-center/pbs  \n\n" <<

               "Included with CASM, the 'qe.relax' script can be executed by the  \n" <<
               "'casm run' command to submit a batch of Quantum Espresso jobs that for selected \n" <<
               "configurations. For each selected configuration, Quantum Espresso is re-run\n" <<
               "using the output of the previous calculation until full convergence \n" <<
               "is achieved. The convergence criteria is: if the cell shape and     \n" <<
               "volume remain constant (calculation != vc-relax) then a single calculation  \n" <<
               "is performed; else the calculation is converged if at least 2 jobs  \n" <<
               "are complete, and: 1) the last job completed with <= 3 ionic steps  \n" <<
               " or, if \"nrg_convergence\" is set in the 'relax.json' file, 2) the \n" <<
               "last two calculations had final energy differ by less than the value of \n" <<
               " \"nrg_convergence\". Once converged, a final constant volume       \n" <<
               "calculation is performed with the following setting: (calculation = 'relax')\n" <<

               "relax.json:                                                         \n" <<
               "  This JSON file contains a single JSON object which contains       \n" <<
               "  parameters used to control PBS job submission settings.           \n" <<
               "  Required keys are:                                                \n" <<
               "    \"queue\": queue to submit job in                               \n" <<
               "    \"ppn\": processors (cores) per node to request                 \n" <<
               "    \"atom_per_proc\": max number of atoms per processor (core)     \n" <<
               "    \"walltime\": walltime to request (ex. \"48:00:00\")            \n\n" <<
               "    \"software\": needs to be quantumespresso for quantum espresso to be used\n\n" <<

               " Optional keys are:                                                 \n" <<
               "    \"account\": account to submit job under (default None)         \n" <<
               "    \"pmem\": string for requested memory (default None)            \n" <<
               "    \"priority\": requested job priority (default \"0\")            \n" <<
               "    \"message\": when to send messages about jobs (ex. \"abe\",     \n" <<
               "               default \"a\")                                       \n" <<
               "    \"email\": where to send messages (ex. \"me@fake.com\", default \n" <<
               "             None)                                                  \n" <<
               "    \"qos\": quality of service, 'qos' option (ex. \"fluxoe\")      \n" <<
               "    \"qe_cmd\": quantum espresso execution command (default is \"pw.x < {INFILE} > {OUTFILE}\" if    \n" <<
               "                ncpus=1, else \"mpirun -np {NCPUS} pw.x < {INFILE} > {OUTFILE}\"           \n" <<
               "    \"infile\": quantum espresso input file name (default is \"std.in\"\n" <<
               "    \"outfile\": quantum espresso output file name (default is \"std.out\"\n" <<
               "    \"ncpus\": number of cpus (cores) to run on (default $PBS_NP)   \n" <<
               "    \"run_limit\": number of vasp runs until \"not_converging\"     \n" <<
               "                 (default 10)                                       \n" <<
               "    \"nrg_convergence\": converged if last two runs complete and    \n" <<
               "                       differ in energy by less than this amount    \n" <<
               "                       (default None)                               \n" <<
               "    \"move\": files to move at the end of a run (ex. \"\",    \n" <<
               "            \".wfc\"], default [])                     \n" <<
               "    \"copy\": files to copy from run to run  default [$infilename]) \n" <<
               "    \"remove\": files to remove at the end of a run                \n" <<
               "               default [\".wfc\", \".igk\", \".save\"]             \n" <<
               "    \"compress\": files to compress at the end of a run (ex.        \n" <<
               "                [$outfilename], default [])          \n" <<
               "    \"backup\": files to compress to backups at the end of a run,   \n" <<
               "              used in conjunction with move (ex. [\".wfc\"])     \n" <<
               "    \"extra_input_files\": extra input files to be copied from the  \n" <<
               "                         settings directory, e.g., an OCCUPATIONS     \n" <<
               "                         file.                                      \n" <<
               "    \"initial\": location of $infile with tags for the initial run,   \n" <<
               "               if desired                                            \n" <<
               "    \"final\": location of $infile with tags for the final run, if  \n" <<
               "             desired                                                \n" <<
               "    \"err_types\": list of errors to check for. Not Implemented yet  \n" <<
               "\n";

      args.log << "EXAMPLE: relax.json \n";
      args.log << "-------\n";
      args.log <<
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
  \"calculator\":\"quantumespresso\"\n\
  \"infilename\":\"LixCoO2.in\"\n\
  \"outfilename\":\"LixCoO2.out\"\n\
}\n";
      args.log << "-------\n\n\n";


      args.log << "SPECIES:                                                            \n" <<
               "  This file contains information for selecting pseudopotentials and specifing\n" <<
               "  parameters that must be set on an atom-by-atom basis in the infile,\n" <<
               "  such as magnetic moment (non currently implemented).\n" <<
               "  The first line in the file specifies the value of \n" <<
               "  'PSEUDO_DIR_PATH', which is the base path used to find UPF     \n" <<
               "  files. The second line contains column headings (at least 4), and \n" <<
               "  then there are lines for each distinct species. The first column  \n" <<
               "  specifies the 'SPECIES' and must match a species names in the PRIM\n" <<
               "  file. The second column gives an 'ALIAS' name for the species which\n" <<
               "  is used for ordering like atoms in the generated input files. The\n" <<
               "  third column should be either '0' or '1', such that only one      \n" <<
               "  species with a given ALIAS has a '1'. For that species the fourth \n" <<
               "  column must contain the path that should be appended to the       \n" <<
               "  PSEUDO_DIR_PATH to specify the UPF file for that species.      \n\n" <<
               "  Additional columns, such as 'if_pos' in the example below are     \n\n" <<
               "  and used to specify the value used for a particular species in the\n" <<
               "  infile. The column heading must match a valid quantum espresso input setting.\n"
               "  For now only supported additional tag is if_pos, a way to fixed certain lattice positions.\n\n";

      args.log << "EXAMPLE: SPECIES \n";
      args.log << "-------\n" <<
               "PSEUDO_DIR_PATH = /absolute/path/to/quantumespresso_potentials\n" <<
               "SPECIES    ALIAS    UPF  UPF_location     if_pos\n" <<
               "Ni         Ni       1       PAW_PBE/Ni.UPF     1,1,1\n" <<
               "Al        Al       1       PAW_PBE/Al.UPF      1,1,1\n";
      args.log << "-------\n\n\n";

      args.log << "$infilename:                                                              \n" <<
               "  This is a template input file used for Quantum Espresso calculations. The settings \n" <<
               "  are generally used as given though some may be automatically set  \n" <<
               "  based on settings in the 'relax.json' or 'SPECIES' files. Also,   \n" <<
               "  some settings might be added or changed if certain errors are     \n" <<
               "  during calculation. The actual input file used for each calculation is \n" <<
               "  saved.                                                            \n\n";
      args.log << "Note:                                                    \n" <<
               "  K_POINTS will be adjusted accordingly such that the density is maintained\n" <<
               "  over all configurations in the project for all Quantum Espresso calculations\n" <<
               "  this uses the CELL_PARAMETERS and the K_POINTS cards in the input file to calculate\n" <<
               "  a density and rescale configurations k-point mesh accordingly\n";

      args.log << "EXAMPLE: Mg2Ti4S8.in \n";
      args.log << "-------\n" <<
               "System = Test of Quantum Espresso submission\n\
&CONTROL\n\
 calculation = 'vc-relax',\n\
 pseudo_dir = '/home/skolli/quantum_espresso/pseudo/',\n\
 tprnfor = .true.,\n\
 prefix = 'Mg2Ti4S8',\n\
 restart_mode = 'from_scratch',\n\
 tstress = .true.,\n\
/\n\
&SYSTEM\n\
 ecutwfc = 45.0,\n\
 occupations = 'fixed',\n\
 celldm(1) = 7.3794,\n\
 ibrav = 0,\n\
 nat = 14,\n\
 ntyp = 3,\n\
 ecutrho = 200.0,\n\
/\n\
&ELECTRONS\n\
 diagonalization = 'cg',\n\
 mixing_mode = 'plain',\n\
 mixing_beta = 0.7,\n\
 conv_thr = 1e-08,\n\
/\n\
&IONS\n\
 ion_dynamics = 'bfgs',\n\
/\n\
&CELL\n\
 press = 0.1,\n\
 cell_factor = 1.6,\n\
 cell_dynamics = 'bfgs',\n\
/\n\
\n\
ATOMIC_SPECIES\n\
 Mg 24.31 Mg.pbe-nsp-bpaw.UPF\n\
 Ti 47.88 Ti.pbe-sp-hgh.UPF\n\
 S 32.07 S.pbe-n-kjpaw_psl.0.1.UPF\n\
\n\
CELL_PARAMETERS angstrom\n\
 0.0000000000000000 5.1756022947592379 5.1756022947592388\n\
 5.1756022947592388 0.0000000000000000 5.1756022947592388\n\
 5.1756022947592388 5.1756022947592379 0.0000000000000000\n\
\n\
ATOMIC_POSITIONS crystal\n\
Mg 0.000000000 0.000000000 0.000000000\n\
Mg 0.250000000 0.250000000 0.250000000\n\
Ti 0.625000000 0.625000000 0.625000000\n\
Ti 0.125000000 0.625000000 0.625000000\n\
Ti 0.625000000 0.125000000 0.625000000\n\
Ti 0.625000000 0.625000000 0.125000000\n\
S 0.3842989149764762 0.3842989149764762 0.3842989149764762\n\
S 0.8657010850235238 0.8657010850235238 0.8657010850235238\n\
S 0.3842989149764762 0.8471032550705786 0.3842989149764762\n\
S 0.3842989149764762 0.3842989149764762 0.8471032550705786\n\
S 0.8471032550705786 0.3842989149764762 0.3842989149764762\n\
S 0.8657010850235238 0.8657010850235238 0.4028967449294214\n\
S 0.8657010850235238 0.4028967449294214 0.8657010850235238\n\
S 0.4028967449294214 0.8657010850235238 0.8657010850235238\n\
\n\
K_POINTS automatic\n\
 6 6 6 0 0 0\n\
\n";
      args.log << "-------\n\n\n";

    }
    if(vm.count("properties")) {
      args.log << "\n### properties.calc.json ##################\n\n";

      args.log << "properties.calc.json:                                               \n" <<
               "  Results of calculations for a particular configuration should be  \n" <<
               "  stored in the directory                                           \n" <<
               "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE,         \n" <<
               "  and calculated properties summarized in the file                  \n" <<
               "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/properties.calc.json \n" <<
               "  The 'properties.calc.json' file is read by CASM to extract the    \n" <<
               "  first-principles calculted properties of interest. If the         \n" <<
               "  'properties.calc.json' file does not exist in the                 \n" <<
               "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE directory\n" <<
               "  CASM assumes that no data is available for that configuration.    \n" <<
               "  The 'properties.calc.json' uses CASM standard units eV and Angstroms\n\n" ;

      args.log << "EXAMPLE:\n";
      args.log << "-------\n";
      args.log << "{\n\
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
          \"relaxed_mag_basis\": [\n\
              -3.93,\n\
               3.82,\n\
               1.198\n\
          ], \n\
          \"relaxed_magmom\": -1.3086, \n\
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
      args.log << "-------\n";

    }
    if(vm.count("comp")) {
      args.log << "\n### composition_axes.json ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "$ROOT/.casm/composition_axes.json\n";
      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "This JSON file contains the currently selected composition axes, and \n" <<
               "a list of possible standard or custom composition axes.              \n\n" <<

               "standard_axes:                                                      \n" <<
               "  A JSON object containing each possible standard composition axes  \n" <<
               "  as an attribute with its index as the key.                        \n\n" <<

               "custom_axes:                                                        \n" <<
               "  A JSON object containing each custom composition axes as an       \n" <<
               "  attribute with its index as the key. The keys should not be       \n" <<
               "  repeats of any of the standard_axes.                              \n\n" <<

               "standard_axes/custom_axes:components                                \n" <<
               "  A JSON array containing the names of possible species.            \n\n" <<

               "standard_axes/custom_axes:independent_compositions                  \n" <<
               "  The number of independent composition axes.                       \n\n" <<

               "standard_axes/custom_axes:origin                                    \n" <<
               "  The composition of origin the of composition axes in terms of     \n" <<
               "  number of each component species per primitive cell, ordered as in\n" <<
               "  the 'components' array.                                           \n\n" <<

               "standard_axes/custom_axes:a, b, c, ...                              \n" <<
               "  The composition of end members a, b, c, etc. in terms of number of\n" <<
               "  each component species per primitive cell, ordered as in the      \n" <<
               "  'components' array.                                               \n\n" <<

               "standard_axes/custom_axes:param_formula:                            \n" <<
               "  The formula that converts 'comp_n' (# of each component per       \n" <<
               "  primitive cell) to 'comp' (composition relative the selected      \n" <<
               "  composition axes).                                                \n\n" <<

               "standard_axes/custom_axes:mol_formula:                              \n" <<
               "  The formula that converts 'comp' (composition relative the        \n" <<
               "  selected composition axes) to 'comp_n' (# of each component per   \n" <<
               "  primitive cell).                                                  \n\n\n";


      args.log << "EXAMPLE:\n";
      args.log << "-------\n";
      args.log <<
               "{\n  \"current_axes\" : \"0\",\n  \"standard_axes\" : {\n    \"0\" : {\n      \"a\" : [\n        [ 2.000000000000 ],\n        [ 0.000000000000 ],\n        [ 2.000000000000 ]\n      ],\n      \"components\" : [ \"Zr\", \"Va\", \"O\" ],\n      \"independent_compositions\" : 1,\n      \"mol_formula\" : \"Zr(2)Va(2-2a)O(2a)\",\n      \"origin\" : [\n        [ 2.000000000000 ],\n        [ 2.000000000000 ],\n        [ 0.000000000000 ]\n      ],\n      \"param_formula\" : \"a(0.5-0.25Va+0.25O)\"\n    },\n    \"1\" : {\n      \"a\" : [\n        [ 2.000000000000 ],\n        [ 2.000000000000 ],\n        [ 0.000000000000 ]\n      ],\n      \"components\" : [ \"Zr\", \"Va\", \"O\" ],\n      \"independent_compositions\" : 1,\n      \"mol_formula\" : \"Zr(2)Va(2a)O(2-2a)\",\n      \"origin\" : [\n        [ 2.000000000000 ],\n        [ 0.000000000000 ],\n        [ 2.000000000000 ]\n      ],\n      \"param_formula\" : \"a(0.5+0.25Va-0.25O)\"\n    }\n  }\n}\n";

      args.log << "-------\n";

    }

    if(vm.count("bspecs")) {
      args.log << "\n### bspecs.json ##################\n\n";

      args.log << "LOCATION:\n";
      args.log << "$ROOT/basis_sets/$CURR_BSET/bspecs.json\n";
      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "This JSON file contains specifications for generating the cluster\n" <<
               "basis functions.                                                    \n\n";

      std::cout << "'site_basis_functions' may specify a string, which can be either 'occupation' or \n"
                << "'chebychev'. Otherwise, specifies a JSON object containing a composition vector or\n"
                << "a JSON array containing multiple composition vectors. A single composition vector\n"
                << "is formatted as, e.g.\n"
                << "   \"composition\" : {\"Au\" : 0.25, \"Cu\" : 0.75} \n"
                << "The site basis functions will then be constructed as to be optimized for that composition.\n\n"

                << "To specify different compositions on multiple sublattices, an array can be used. \n"
                << "As an example, the following specifies a different composition on sublattice 0 than\n"
                << "on sublattices 1 and 3: \n\n"

                << "   \"site_basis_functions\" : [\n"
                << "                                {\n"
                << "                                  \"composition\" : {\"Ga\" : 0.3, \"In\" : 0.7},\n"
                << "                                  \"sublat_indices\" : [0]\n"
                << "                                },\n"
                << "                                {\n"
                << "                                  \"composition\" : {\"Ga\" : 1.0, \"In\" : 0.0},\n"
                << "                                  \"sublat_indices\" : [1,2]\n"
                << "                                }\n"
                << "                             ]\n\n"

                << "Sublattices are specified in the same order as in prim.json. Sublattice compositions\n"
                << "are not allowed to break the symmetry of the crystal. If equivalent sublattices are\n"
                << "assigned inequivalent compositions, one will be chosen arbitrarily and propagated to\n"
                << "all equivalent sublattices.  The resulting site basis functions can be reviewed using\n"
                << "'casm bset --functions'\n\n";


      args.log << "The JSON object 'orbit_branch_specs' specifies the maximum size of pair,   \n" <<
               "triplet, quadruplet, etc. clusters in terms of the maximum distance \n" <<
               "between any two sites in the cluster.\n\n";

      args.log << "The JSON array 'orbit_specs' allows specifying particular custom orbits    \n" <<
               "by providing the prototype cluster coordinates. The 'include_subclusters'  \n" <<
               "option allows including all orbits of subclusters of the specified cluster.\n" <<
               "The cluster coordinates may be in \"Direct\"/\"Fractional\" coordinates,   \n"
               "\"Cartesian\" coordinates, or \"Integral\" coordinates. \"Integral\"       \n"
               "coordinates are 4-element integer arrays indicating sublattice index, b,   \n"
               "followed by unit cell indices, i, j, k.                                    \n\n\n";


      args.log << "EXAMPLE:\n";
      args.log << "-------\n";
      args.log <<
               "{\n  \"basis_functions\" : {\n    \"site_basis_functions\" : \"occupation\"\n  },\n  \"orbit_branch_specs\" : {\n    \"2\" : {\"max_length\" : 4.01},\n    \"3\" : {\"max_length\" : 3.01}\n  },\n  \"orbit_specs\" : [\n    {\n      \"coordinate_mode\" : \"Direct\",\n      \"prototype\" : [\n        [ 0.000000000000, 0.000000000000, 0.000000000000 ],\n        [ 1.000000000000, 0.000000000000, 0.000000000000 ],\n        [ 2.000000000000, 0.000000000000, 0.000000000000 ],\n        [ 3.000000000000, 0.000000000000, 0.000000000000 ]\n      ],\n      \"include_subclusters\" : true  \n    },\n    {\n      \"coordinate_mode\" : \"Integral\",\n      \"prototype\" : [\n        [ 0, 0, 0, 0 ],\n        [ 1, 0, 0, 0 ],\n        [ 0, 0, 0, 3 ]\n      ],\n      \"include_subclusters\" : true\n    }\n  ]\n}\n";
      args.log << "-------\n";

    }

    if(vm.count("clust")) {
      args.log << "\n### clust.json ##################\n\n";

      args.log << "LOCATION:\n";
      args.log << "$ROOT/basis_sets/$CURR_BSET/clust.json\n";
      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "This JSON file contains the coordinates of sites in the prototype   \n" <<
               "clusters generated using the 'bspecs.json' specifications.          \n\n";


      args.log << "Prototype clusters can be accessed via:                            \n"
               "  [\"branches\"][branch_index][\"orbits\"][orbit_index][\"prototype\"]\n\n"

               "\"prototype\": (JSON object)                                       \n"

               "  /\"max_length\": (number)                                        \n"
               "     Maximum pair distance between sites in the cluster            \n\n"

               "  /\"min_length\": (number)                                        \n"
               "     Minimum pair distance between sites in the cluster            \n\n"

               "  /\"sites\": (JSON array of Integral coordinates)                 \n"
               "     An array listing sites in the prototype cluster using Integral\n"
               "     coordinates. Integral coordinates are 4-element integer arrays\n"
               "     indicating sublattice index, b, followed by unit cell indices,\n"
               "     i, j, k.                                                      \n\n"

               "\"bspecs\": (JSON object)                                          \n"
               "  For reference, the contents of the 'bspecs.json' file used to    \n"
               "  generate these clusters is reproduced here.                      \n\n"

               "\"lattice\": (JSON object)                                         \n"
               "  For reference, so that the Integral coordinates can be converted \n"
               "  into Fractional or Cartesian coordinates, the lattice vectors    \n"
               "  of the primitive structure are reproduced here.                  \n\n" << std::endl;
    }

    if(vm.count("basis")) {
      args.log << "\n### basis.json ##################\n\n";

      args.log << "LOCATION:\n";
      args.log << "$ROOT/basis_sets/$CURR_BSET/basis.json\n";
      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "This JSON file contains the basis functions generated using the    \n"
               "'bspecs.json' specifications.                                      \n\n";


      args.log << "\"site_functions\": (JSON array of JSON object)                    \n"
               "  Gives the site basis functions. One JSON object for each basis   \n"
               "  site. \n\n"

               "  /\"sublat\": (int)                                               \n"
               "    Basis site index.                                              \n\n"

               "  /\"asym_unit\": (int)                                            \n"
               "    Index of the asymmetric unit this basis site belongs to.       \n\n"

               "  /\"basis\": (JSON object)                                        \n"
               "     Gives the value of each site basis function for each possible \n"
               "     occupant. Of the form:                                        \n\n"
               "       { \n"
               "         \"\\\\phi_b_i\": { \n"
               "           \"A\": val, \n"
               "           \"B\": val, \n"
               "           ... \n"
               "         }, \n"
               "         ... \n"
               "       } \n"

               "\"cluster_functions\": (JSON array of JSON object)                 \n"
               "  Gives the cluster basis functions. One JSON object for each      \n"
               "  cluster basis function.                                          \n\n"

               "  /\"linear_function_index\": (int)                                \n"
               "    Linear function index. This corresponds to ECI indices.        \n\n"

               "  /\"mult\": (int)                                                 \n"
               "    Multiplicity of symmetrically equivalent cluter functions.     \n\n"

               "  /\"orbit\": (JSON array of 3 int)                                \n"
               "     Gives the cluster branch index, cluster orbit index, and index\n"
               "     of this basis function in the cluster basis.                  \n\n"

               "  /\"prototype\": (JSON object)                                    \n"
               "     Specifies the prototype cluster, as in the 'clust.json' file. \n\n"

               "  /\"prototype_function\": (string)                                \n"
               "     Latex-style function for the prototype cluster.               \n\n" << std::endl;
    }

    if(vm.count("clex")) {
      args.log << "\n### $TITLE_Clexulator.* ##################\n\n";

      args.log << "LOCATION:\n";
      args.log << "$ROOT/basis_sets/$CURR_BSET/$TITLE_Clexulator.*\n";
      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "$TITLE_Clexulator.cc contains C++ code generated by CASM for       \n"
               "the cluster basis functions. It is automatically compiled into     \n"
               "$TITLE_Clexulator.o and $TITLE_Clexulator.so for use by CASM.      \n\n" << std::endl;
    }

    if(vm.count("ref")) {
      args.log << "\n### ref ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n\n";
      args.log << "$ROOT/training_data/settings/$CURR_CALCTYPE/$CURR_REF/chemical_reference.json\n";
      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
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

      args.log << "    The reference states are stored in the 'chemical_reference.json' file\n"
               "    in one of two formats:                                               \n\n"

               "    1) Reference state composition and energy_per_species.               \n"
               "       In this format each reference state is represented by a JSON      \n"
               "       object storing the number of each species present in the reference\n"
               "       state and the energy_per_species for that reference state. Species\n"
               "       that are not in the primitive structure may also be included in   \n"
               "       the reference states as long as the composition space of the      \n"
               "       primitive structure is spanned by the hyperplane connecting the   \n"
               "       provided reference states.                                        \n"
               R"(       '[)" << "\n" <<
               R"(          {"A": 3.4, "C": 2.0, "energy_per_species": 2.0},)" << "\n" <<
               R"(          {"B": 2.0, "energy_per_species": 4.0}, )" << "\n" <<
               R"(          {"C": 1.0, "energy_per_species": 3.0}  )" << "\n" <<
               R"(        ]')" << "\n\n" <<

               "    2) Input an array of energy_per_species, for each species in prim,   \n"
               "       including 0.0 for vacancy:                                        \n"
               "        '[X, X, X]'                                                      \n\n";

      args.log << "    When using '--set' it is also possible to specialize the chemical    \n"
               "    reference at the supercell or configuration level by adding the      \n"
               "    --scelname or --configname option.                                   \n\n";


      args.log << "EXAMPLE: chemical_reference.json\n";
      args.log << "-------\n";
      args.log <<
               R"({
"chemical_reference" : {
"config" : {
"SCEL4_2_2_1_1_1_0 / 0" : [
{
"A" : 1.000000000000,
"energy_per_species" : -1.500000000000
},
{
"B" : 1.000000000000,
"energy_per_species" : -2.000100000000
},
{
"C" : 1.000000000000,
"energy_per_species" : -8.030000000000
}
],
"SCEL4_2_2_1_1_1_0 / 2" : [ -1.520000000000, -2.000100000000, -8.030000000000 ]
},
"global" : [
{
"A" : 0.500000000000,
"B" : 0.500000000000,
"energy_per_species" : -1.500000000000
},
{
"B" : 1.000000000000,
"energy_per_species" : -2.000000000000
},
{
"C" : 1.000000000000,
"energy_per_species" : -8.000000000000
},
{
"D" : 1.000000000000,
"energy_per_species" : -4.000000000000
}
],
"species_order" : [ "A", "B", "C" ],
"supercell" : {
"SCEL3_1_3_1_1_0_0" : [
{
"A" : 1.000000000000,
"energy_per_species" : -1.500000000000
},
{
"B" : 1.000000000000,
"energy_per_species" : -2.000000000000
},
{
"C" : 1.000000000000,
"energy_per_species" : -8.001000000000
}
],
"SCEL4_2_2_1_1_1_0" : [
{
"A" : 1.000000000000,
"energy_per_species" : -1.500000000000
},
{
"B" : 1.000000000000,
"energy_per_species" : -2.000000000000
},
{
"C" : 1.000000000000,
"energy_per_species" : -8.030000000000
}
]
}
}
})";
      args.log << "-------\n";

    }

    if(vm.count("scel")) {
      args.log << "\n### scel ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "$ROOT/training_data/SCEL\n\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "Contains a list of all the generated and imported supercells. Each  \n" <<
               "entry gives the name of a supercell, its volume in number of        \n" <<
               "primitive volumes, and the transformation matrix to go from the     \n" <<
               "primitive cell to the supercell. The convention is                  \n" <<
               "            LAT.scel = LAT.prim*transf_matrix,                      \n" <<
               "where the columns of the LAT matrices are the lattice vectors.      \n\n\n";

      args.log << "EXAMPLE:\n";
      args.log << "-------\n";
      args.log <<
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
      args.log << "-------\n";
    }

    if(vm.count("lat")) {
      args.log << "\n### lat ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "$ROOT/training_data/$SCELNAME/LAT\n\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "Contains the lattice vectors of a particular supercell of your CASM \n" <<
               "project. The format is the same as the first lines in a standard    \n" <<
               "vasp POSCAR file, excluding the title (scaling followed by the three\n" <<
               "lattice vectors as rows).\n\n\n";

      args.log << "EXAMPLE:\n";
      args.log << "-------\n";
      args.log <<
               " 1.00000000\n\
      10.52850134      0.00000000      0.00000000\n\
      0.00000000     10.52850134      0.00000000\n\
      0.00000000      0.00000000     10.52850134\n";
      args.log << "-------\n";

    }

    if(vm.count("pos")) {
      args.log << "\n### pos ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "$ROOT/training_data/$SCELNAME/$CONFIGID/POS\n\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "This file is generated using the '--write-pos' option for 'casm run'.\n";
      args.log << "Decorated configuration for a particular supercell. It is a         \n" <<
               "supercell of your primitive structure after the enumeration on the  \n" <<
               "alloying sites. The format is standard vasp 5.x format, and the     \n" <<
               "coordinates of the sites for the configuration are the ideal sites  \n" <<
               "specified in the PRIM file.\n\n\n";

      args.log << "EXAMPLE:\n";
      args.log << "-------\n";
      args.log <<
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
      args.log << "-------\n";

    }

    if(vm.count("eci")) {
      args.log << "\n### eci.json ##################\n\n";

      args.log << "LOCATION:\n";
      args.log << "$ROOT/cluster_expansions/clex.formation_energy/$CURR_BSET/$CURR_CALCTYPE/$CURR_REF/$CURR_ECI/eci.json\n";
      args.log << "\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "This is a copy of the $ROOT/basis_sets/$CURR_BSET/'basis.json' file \n"
               "with the following additions:                                      \n\n"

               "\"cluster_functions\": (JSON array of JSON object)                 \n\n"

               "  /\"eci\": (number, optional, default=0.0)                        \n"
               "     The value of the ECI for the cluster basis function. If not   \n"
               "     given, use 0.0.                                               \n\n"

               "\"fit\": (JSON object)                                             \n"
               "  Data from 'casm-learn' specifying how the ECI where generated and\n"
               "  some goodness of fit measures.                                   \n\n" << std::endl;

    }

    if(vm.count("monte")) {
      args.log << "\n### monte ##################\n\n";

      args.log << "LOCATION WHEN GENERATED:\n";
      args.log << "  User determined\n\n\n";

      args.log << "DESCRIPTION:\n";
      args.log << "  The Monte Carlo input file does not need to be in any particular \n" <<
               "  location, as long as it is somewhere inside the CASM project     \n" <<
               "  directory or subdirectories. The input file contains a JSON      \n" <<
               "  object with \"ensemble\", \"method\", \"model\", \"supercell\",  \n" <<
               "  \"data\", and \"driver\" attributes, as described below. An      \n" <<
               "  optional attribute \"debug\" may also be included to print       \n" <<
               "  information that may be useful for debugging an input file.      \n\n" <<

               "Input file parameters:                                             \n\n" <<

               "\"ensemble\" (string):                                             \n\n" <<

               "  Possible options for \"ensemble\" are:                           \n\n" <<

               "    \"GrandCanonical\" or \"grand_canonical\": Semi-grand canonical\n" <<
               "    Monte Carlo calculation in which the total number of sites is  \n" <<
               "    fixed, but the occupants on each site may vary. One occupant   \n" <<
               "    change at a time is attempted.                                 \n\n" <<

               "    \"Canonical\" or \"canonical\": Canonical Monte Carlo \n" <<
               "    calculation in which the total number of each type of occupant \n"
               "    is fixed. Each Monte Carlo step attempts to swap a pair of     \n"
               "    occupants.                                                     \n\n\n" <<


               "\"method\" (string):                                               \n\n" <<

               "  Possible options for \"method\" are:                             \n\n" <<

               "    \"Metropolis\" or \"metropolis\": Run Monte Carlo calculations \n" <<
               "    using the Metropolis algorithm.                                \n\n" <<

               "    \"LTE1\" or \"lte1\": Single spin flip low temperature         \n" <<
               "    expansion calculations.                                        \n\n\n" <<


               "\"model\": (JSON object)                                           \n\n" <<

               "  /\"formation_energy\": (string, optional, default=\"formation_energy\")\n" <<
               "    Specifies the cluster expansion to use to calculated formation \n"
               "    energy. Should be one of the ones listed by 'casm settings -l'.\n\n\n" <<


               "\"supercell\": (3x3 JSON arrays of integers)                      \n" <<
               "    The supercell transformation matrix.                           \n\n" <<


               "\"data\": (JSON object)                                            \n\n" <<

               "  /\"sample_by\": (string)                                         \n" <<
               "    Specify unit for the period between samples.  May be either    \n" <<
               "    \"step\" (sample after every \"sample_period\" proposed Monte  \n" <<
               "    Carlo events), or \"pass\" (sample after the \"sample_period\" \n" <<
               "    number of passes), where 1 pass is a number of steps equal to  \n" <<
               "    the number of sites in the supercell that have variable        \n" <<
               "    occupation).                                                   \n\n" <<

               "  /\"sample_period\": (integer)                                    \n" <<
               "    Specify how many steps or passes to wait between data samples. \n\n" <<

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
               "      \"all_correlations\": correlations (per unit cell)           \n" <<
               "      \"<anything else>\": is interpreted as a 'casm query' query  \n\n" <<

               "  /\"confidence\": (number, range (0.0, 1.0), default 0.95)        \n" <<
               "    The confidence level used for calculating the precision in the \n" <<
               "    average value of sampled quantities.                           \n\n" <<

               "  /\"min_pass\", /\"min_step\", /\"min_sample\": (integer)         \n" <<
               "    If in automatic convergence mode, prevents the calculation from\n" <<
               "    a minimum number of passes, steps, or samples have occurred.   \n\n" <<

               "  /\"max_pass\", /\"max_step\", /\"max_sample\": (integer)         \n" <<
               "    If in automatic convergence mode, stops the calculation if the \n" <<
               "    specified number of passes, steps, or samples have occurred.   \n\n" <<

               "  /\"N_pass\", /\"N_step\", /\"N_sample\": (integer)               \n" <<
               "    When not in automatic convergence mode (no precision has been  \n" <<
               "    specified for any quantities being sampled), stops the         \n" <<
               "    calculation when the specified number of passes, steps, or     \n" <<
               "    samples have occurred.                                         \n\n" <<

               "  /\"equilibration_passes_first_run\": (integer)                   \n" <<
               "    If included, the requested number of passes will be performed  \n" <<
               "    at the initial conditions as a preliminary step before the     \n" <<
               "    actual run begins. This may be useful when not running in      \n" <<
               "    automatic convergence mode.                                    \n\n" <<

               "  /\"equilibration_passes_each_run\": (integer)                    \n" <<
               "    If included, the requested number of passes will be performed  \n" <<
               "    at each condition as a preliminary step before the actual run  \n" <<
               "    begins. This may be useful when not running in automatic       \n" <<
               "    convergence mode.                                              \n\n" <<

               "  /\"storage\": (JSON object) Options for writing results.         \n\n" <<

               "    /\"output_format\": (string or JSON array of string)           \n" <<
               "      Specifies the type or types of output files. Current options \n" <<
               "      are \"csv\" or \"json\". Type names with either all lower    \n" <<
               "      case or all   upper case are accepted.                       \n\n" <<

               "    /\"write_observations\": (boolean, default false)              \n" <<
               "      If true, all individual observations of the quantities       \n" <<
               "      requested to be sampled will be written to compressed files: \n" <<
               "        \"output_directory\"/conditions.i/observations.ext.gz      \n" <<
               "      where 'i' is the condition index and 'ext' is the output     \n" <<
               "      format.                                                      \n\n" <<

               "    /\"write_trajectory\": (boolean, default false)                \n" <<
               "      If true, the value of all degrees of freedom at the time of  \n" <<
               "      each sample will be written to compressed files:             \n" <<
               "        \"output_directory\"/conditions.i/trajectory.ext.gz        \n" <<
               "      where 'i' is the condition index and 'ext' is the output     \n" <<
               "      format.                                                      \n\n" <<

               "  /\"enumeration\": (JSON object, optional)                        \n" <<
               "    If included, save configurations encountered during Monte      \n" <<
               "    Carlo calculations by keeping a 'hall of fame' of best scoring \n" <<
               "    configurations. After the calculation at a particular set of   \n" <<
               "    thermodynamic conditions completes, the configurations in the  \n" <<
               "    hall of fame are saved to the project configuration list.      \n\n" <<

               "    /\"check\": (string, default=\"eq(1,1)\")                      \n" <<
               "      A 'casm query'-like string that returns a boolean value      \n" <<
               "      indicating if (true) a configuration should be considered for\n" <<
               "      for the enumeration hall of fame. The default always returns \n" <<
               "      true.                                                        \n\n" <<

               "    /\"metric\": (string, default=\"clex_hull_dist(ALL)\")         \n" <<
               "      A 'casm query'-like string that provides a metric for ranking\n" <<
               "      ranking configurations as they are encountered during a Monte\n" <<
               "      Carlo calculation. The resulting value is used to create a   \n" <<
               "      hall of fame of 'best' configurations encountered during the \n" <<
               "      calculation. When the calculation is complete configurations \n" <<
               "      in the hall of fame are added to the CASM project config     \n" <<
               "      list. The 'casm query'-like command should evaluate to a     \n" <<
               "      number.                                                      \n\n" <<

               "      Besides the properties listed via 'casm query -h properties',\n" <<
               "      and 'casm query -h operators', both \"check\" and \"metric\" \n" <<
               "      can also use the property \"potential_energy\".              \n\n" <<

               "    /\"sample_mode\": (string, optional, default=\"on_sample\")    \n" <<
               "      Indicate when to attempt to insert configurations into the   \n" <<
               "      enumeration hall of fame. Options are:                       \n" <<
               "        \"on_accept\": after each accepted Monte Carlo event       \n" <<
               "        \"on_sample\": after each data sample                      \n\n" <<

               "    /\"check_existence\": (bool, optional, default=true)           \n" <<
               "      If true, only configurations that do not already exist in the\n" <<
               "      config list are inserted into the enumeration hall of fame.  \n\n" <<

               "    /\"insert_canonical\": (bool, optional, default=true)          \n" <<
               "      If true, configurations are inserted into the enumeration    \n" <<
               "      hall of fame in their canonical form. If 'check_existence' is\n" <<
               "      true, this must be set to true.                              \n\n" <<

               "    /\"N_halloffame\": (integer, optional, default=100)            \n" <<
               "      The number of configurations that are allowed in the         \n" <<
               "      enumeration hall of fame.                                    \n\n" <<

               "    /\"tolerance\": (number, optional, default=1e-8)               \n" <<
               "      Tolerance used for floating point comparison of configuration\n" <<
               "      scores in the enumeration hall of fame.                      \n\n\n" <<


               "\"driver\": (JSON object)                                          \n\n" <<

               "  /\"motif\": (JSON object)                                        \n" <<
               "      Specifies the initial occupation of the supercell.           \n\n" <<

               "      For canonical ensemble Monte Carlo calculations an additional\n" <<
               "      step changes the occupants on random sites to make the actual\n" <<
               "      composition as close as possible to the requested composition.\n\n" <<

               "    /\"configname\": (string, optional)                            \n" <<
               "      A configuration name, \"auto\", \"restricted_auto\", or      \n" <<
               "      \"default\".                                                 \n\n" <<

               "      Specifies the configuration that is tiled to fill the        \n" <<
               "      supercell. If necessary, symmetry operations may be applied  \n" <<
               "      An error is thrown if the specified configuration can not be \n" <<
               "      used to fill the \"supercell\".                              \n\n" <<

               "      Possible options for \"configname\" are:                     \n" <<
               "        A configuration name (ex. \"SCEL3_3_1_1_0_2_2/0\")         \n" <<
               "        \"auto\": (\"grand_canonical\" ensemble only) Enumerated   \n" <<
               "        configurations will be searched for the configuration with \n" <<
               "        the lowest potential energy to use as the motif.           \n" <<
               "        \"default\": If the value \"default\" is used, the initial \n" <<
               "        motif occupation is determined from the occupation order in\n" <<
               "        the PRIM.                                                  \n" <<
               "        \"restricted_auto\": (\"grand_canonical\" ensemble only)   \n" <<
               "        Same as \"auto\", but only configurations that can tile the\n" <<
               "        supercell are considered. As a last resort, \"default\" is \n" <<
               "        used.                     \n\n" <<

               "    /\"configdof\": (string, optional)                             \n" <<
               "      Specifies the path to a configdof JSON file, such as         \n" <<
               "      \"initial_state.json\" or \"final_state.json\", containing   \n" <<
               "      the degrees of freedom to initialize the supercell with      \n\n" <<

               "  /\"mode\": (string)                                              \n" <<
               "    Specify the drive mode.                                        \n\n" <<

               "    Possible options for \"mode\" are:                             \n" <<
               "      \"incremental\": perform one or more calculations, starting  \n" <<
               "        at the initial conditions and incrementing by the          \n" <<
               "        incremental conditions up to (and including) the final     \n" <<
               "        conditions.                                                \n\n" <<
               "      \"custom\": perform one or more calculations, as specified by\n" <<
               "        the \"custom_conditions\".                                 \n\n" <<

               "  /\"dependent_runs\": (boolean, default true)                     \n\n" <<

               "    If true, begin the next calculation with the final DoF from the\n" <<
               "    previous calculation. If false, begin each calculation with the\n" <<
               "    DoF specified for the \"motif\".\n\n" <<


               "  /\"initial_conditions\",\n" <<
               "  /\"incremental_conditions\", \n" <<
               "  /\"final_conditions\": (JSON object, optional)                    \n" <<
               "    Specifies the applied conditions for the calculation. For       \n" <<
               "    \"incremental_conditions\", specifies the change in conditions  \n" <<
               "    between individual calculations. Each JSON object should        \n" <<
               "    include:                                                        \n\n" <<

               "    /\"temperature\": (number)                                      \n" <<
               "      The temperature in K.                                         \n\n" <<

               "    /\"param_chem_pot\" (JSON object, \"grand_canonical\" ensemble only)\n" <<
               "      The parametric chemical potential(s)                          \n\n" <<

               "      /\"a\", /\"b\", ...: (number)                                 \n" <<
               "      The parametric chemical potentials conjugate to the parametric\n" <<
               "      compositions. The number of parametric chemical potentials    \n" <<
               "      provided must match the number of independent compositions.   \n\n" <<

               "    /\"comp\" (JSON object, \"canonical\" ensemble only, option 1)\n" <<
               "      The parametric composition(s)                          \n\n" <<

               "      /\"a\", /\"b\", ...: (number)                                 \n" <<
               "      The parametric composition. The number of entries provided    \n" <<
               "      must match the number of independent compositions.            \n\n" <<

               "    /\"comp\" (JSON array, \"canonical\" ensemble only, option 2)\n" <<
               "      The parametric composition(s)                          \n\n" <<

               "      [number, ...]                                                 \n" <<
               "      An array containing the parametric composition. The number of \n" <<
               "      entries provided must match the number of independent         \n" <<
               "      compositions.                                                 \n\n" <<

               "    /\"comp_n\" (JSON object, \"canonical\" ensemble only, option 3)\n" <<
               "      The mol composition per unitcell                          \n\n" <<

               "      /\"A\", /\"B\", ...: (number)                                 \n" <<
               "      The mol composition per unitcell. The entries provided must   \n" <<
               "      match occupant names in the 'prim.json' file. The values are  \n" <<
               "      summed, normalized, and then converted to parametric composition.\n\n" <<

               "    /\"comp_n\" (JSON array, \"canonical\" ensemble only, option 4)\n" <<
               "      The mol composition per unitcell                          \n\n" <<

               "      [number, ...]                                                 \n" <<
               "      An array containing the mol composition per unitcell. The     \n" <<
               "      number of entries provided must match the number components.  \n" <<
               "      The values are summed, normalized, and converted to parametric\n" <<
               "      composition.  \n\n" <<

               "    /\"tolerance\": (number)                                        \n" <<
               "      Specifies a numerical tolerance for comparing conditions.     \n\n" <<

               "  /\"custom_conditions\":\n" <<
               "    (JSON array of JSON objects) An array specifying a custom     \n" <<
               "    path of conditions.                                           \n\n" <<

               "  Restarts: Metropolis Monte Carlo calculations that are stopped   \n" <<
               "  before the entire path has been calculated can be restarted as   \n" <<
               "  long as the conditions of the existing calculations agree with   \n" <<
               "  the conditions specified in the input settings. This means that  \n" <<
               "  the \"final_conditions\" might be changed to increase the length \n" <<
               "  of a path, or additional \"custom_conditions\" might be added,   \n" <<
               "  but the \"incremental_conditions\" may not be changed. Upon      \n" <<
               "  restart, the results summary file is checked for the last        \n" <<
               "  finished conditions. Then the path is resumed from the next set  \n" <<
               "  of conditions. It is the responsibility of the user to ensure    \n" <<
               "  that other important settings, such as the \"model\" are not     \n" <<
               "  changed inappropriately.                                         \n\n\n" <<


               "\"debug\" (bool, default false):                                   \n\n" <<

               "  If true, will print as much information as possible to assist in \n" <<
               "  debugging input file settings.                                   \n\n\n";


      args.log << "EXAMPLE: Settings for an incremental Metropolis calculation     \n" <<
               "with increasing temperature in automatic convergence mode.\n";
      args.log << "-------\n";
      args.log << "{\n  \"comment\" : \"This is a sample input file. Unrecognized attributes (like the ones prepended with '_' are ignored.\",\n  \"debug\" : false,\n  \"ensemble\" : \"grand_canonical\",\n  \"method\" : \"metropolis\",\n  \"model\" : {\n    \"formation_energy\" : \"formation_energy\"\n  },\n  \"supercell\" : [\n    [10, 0, 0],\n    [0, 10, 0],\n    [0, 0, 10]\n  ],\n  \"data\" : {\n    \"sample_by\" : \"pass\",\n    \"sample_period\" : 1,\n    \"_N_sample\" : 1000, \n    \"_N_pass\" : 1000,\n    \"_N_step\" : 1000,\n    \"_max_pass\" : 10000,\n    \"min_pass\" : 1000,\n    \"_max_step\" : 10000,\n    \"_max_sample\" : 500,\n    \"_min_sample\" : 100,\n    \"confidence\" : 0.95,\n    \"measurements\" : [ \n      { \n        \"quantity\" : \"formation_energy\"\n      },\n      { \n        \"quantity\" : \"potential_energy\"\n      },\n      { \n        \"quantity\" : \"atom_frac\"\n      },\n      { \n        \"quantity\" : \"site_frac\"\n      },\n      { \n        \"quantity\" : \"comp\",\n        \"precision\" : 1e-3\n      },\n      { \n        \"quantity\" : \"comp_n\"\n      }\n    ],\n    \"storage\" : {\n      \"write_observations\" : false,\n      \"write_trajectory\" : false,\n      \"output_format\" : [\"csv\", \"json\"]\n    }\n  },\n  \"driver\" : {\n    \"mode\" : \"incremental\", \n    \"motif\" : {\n      \"configname\" : \"auto\",\n      \"_configname\" : \"SCEL3_3_1_1_0_2_2/0\",\n      \"_configdof\" : \"path/to/final_state.json\"\n    },\n    \"initial_conditions\" : {\n      \"param_chem_pot\" : {\n        \"a\" : -1.75\n      },\n      \"temperature\" : 100.0,\n      \"tolerance\" : 0.001\n    },\n    \"final_conditions\" : {\n      \"param_chem_pot\" : {\n        \"a\" : -1.75\n      },\n      \"temperature\" : 1000.0,\n      \"tolerance\" : 0.001\n    },\n    \"incremental_conditions\" : {\n      \"param_chem_pot\" : {\n        \"a\" : 0.0\n      },\n      \"temperature\" : 10.0,\n      \"tolerance\" : 0.001\n    }\n  }\n}\n";
      args.log << "-------\n\n";

      args.log << "EXAMPLE: Settings for an custom drive mode LTE1 calculation with\n" <<
               "increasing temperature.\n";
      args.log << "-------\n";
      args.log << "{\n  \"comment\" : \"This is a sample input file. Unrecognized attributes (like the ones prepended with '_' are ignored.\",\n  \"debug\" : false,\n  \"ensemble\" : \"grand_canonical\",\n  \"method\" : \"lte1\",\n  \"model\" : {\n    \"formation_energy\" : \"formation_energy\"\n  },\n  \"supercell\" : [\n    [9, 0, 0],\n    [0, 9, 0],\n    [0, 0, 9]\n  ],\n  \"data\" : {\n    \"storage\" : {\n      \"write_observations\" : false,\n      \"write_trajectory\" : false,\n      \"output_format\" : [\"csv\", \"json\"]\n    }\n  },\n  \"driver\" : {\n    \"mode\" : \"incremental\", \n    \"motif\" : {\n      \"configname\" : \"auto\",\n      \"_configname\" : \"SCEL3_3_1_1_0_2_2/0\",\n      \"_configdof\" : \"path/to/final_state.json\"\n    },\n    \"custom_conditions\" : [\n      {\n        \"param_chem_pot\" : {\n          \"a\" : 0.0\n        },\n        \"temperature\" : 100.0,\n        \"tolerance\" : 0.001\n      },\n      {\n        \"param_chem_pot\" : {\n          \"a\" : 0.0\n        },\n        \"temperature\" : 200.0,\n        \"tolerance\" : 0.001\n      },\n      {\n        \"param_chem_pot\" : {\n          \"a\" : 0.0\n        },\n        \"temperature\" : 400.0,\n        \"tolerance\" : 0.001\n      },\n      {\n        \"param_chem_pot\" : {\n          \"a\" : 0.0\n        },\n        \"temperature\" : 800.0,\n        \"tolerance\" : 0.001\n      }\n    ]\n  }\n}\n";
      args.log << "-------\n";

    }

    return 0;
  }

}
