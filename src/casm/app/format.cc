#include <cstring>

#include "casm/app/casm_functions.hh"
#include "casm/casm_io/Log.hh"
#include "casm/completer/Handlers.hh"
#include "casm/global/definitions.hh"

namespace CASM {

namespace Completer {
FormatOption::FormatOption() : OptionHandlerBase("format") {}

void FormatOption::initialize() {
  add_help_suboption();

  m_desc.add_options()("dir", "CASM project directory structure summary")(
      "project_settings",
      "Description and location of 'project_settings' file")(
      "prim", "Description and location of 'prim.json' and 'PRIM' files")(
      "sym",
      "Description and location of 'lattice_point_group.json', "
      "'factor_group.json' and 'crystal_point_group.json' files")(
      "vasp", "Description and location of VASP settings files")(
      "properties", "Description and location of properties.calc.json files")(
      "qe", "Description and location of Quantum Espresso settings files")(
      "comp", "Description and location of 'composition_axes.json' file")(
      "bspecs", "Description and location of 'bspecs.json' file")(
      "clust", "Description and location of 'clust.json' file")(
      "basis", "Description and location of 'basis.json' file")(
      "clex", "Description and location of '$TITLE_Clexulator.*' files")(
      "ref", "Description and location of 'chemical_reference.json' files")(
      "lat", "Description and location of 'LAT' files")(
      "pos", "Description and location of 'POS' files")(
      "eci", "Description and location of 'eci.json' file")(
      "monte", "Description and location of the Monte Carlo input file");
  return;
}
}  // namespace Completer

// ///////////////////////////////////////
// 'format' function for casm
//    (add an 'if-else' statement in casm.cpp to call this)

int format_command(const CommandArgs &args) {
  po::variables_map vm;

  Completer::FormatOption format_opt;

  try {
    po::store(
        po::parse_command_line(args.argc(), args.argv(), format_opt.desc()),
        vm);

    /** --help option
     */
    if (vm.count("help") || vm.size() == 0) {
      log() << std::endl;
      log() << format_opt.desc() << std::endl;

      return 0;
    }

    if (vm.count("desc")) {
      log() << "\n";
      log() << format_opt.desc() << std::endl;

      log() << "DESCRIPTION\n"
               "    This option describes the files contained within a CASM "
               "project \n"
               "    and where to find them. For a summary of the directory "
               "structure\n"
               "    of a CASM project using VASP for calculating configuration "
               "use  \n"
               "    the --dir option. Not all files are always present.        "
               "     \n";

      return 0;
    }

    po::notify(vm);

  } catch (po::error &e) {
    err_log() << "ERROR: " << e.what() << std::endl << std::endl;
    err_log() << format_opt.desc() << std::endl;
    return 1;
  } catch (std::exception &e) {
    err_log() << "Unhandled Exception reached the top of main: " << e.what()
              << ", application will now exit" << std::endl;
    return 1;
  }

  if (vm.count("dir")) {
    log() << "\n### dir ##################\n\n";

    log() << "  The expected CASM project directory structure with VASP "
             "settings  \n"
             "  files:                                                    \n"
             "                                                            \n"
             "    $ROOT/                                                  \n"
             "      prim.json                                             \n"
             "      (PRIM)                                                \n"
             "      LOG                                                   \n"
             "      .casm/                                                \n"
             "        composition_axes.json                               \n"
             "        project_settings.json                               \n"
             "        jsonDB/                                             \n"
             "          config_list.json                                  \n"
             "          scel_list.json                                    \n"
             "        query/                                              \n"
             "          Configuration/                                    \n"
             "          Supercell/                                        \n"
             "      symmetry/                                             \n"
             "        lattice_point_group.json                            \n"
             "        factor_group.json                                   \n"
             "        crystal_point_group.json                            \n"
             "      basis_sets/                                           \n"
             "        $CURR_BSET/                                         \n"
             "          bspecs.json                                       \n"
             "          basis.json                                        \n"
             "          clust.json                                        \n"
             "          $TITLE_Clexulator.*                               \n"
             "      training_data/                                        \n"
             "        SCEL                                                \n"
             "        settings/$CURR_CALCTYPE/                            \n"
             "          relax.json                                        \n"
             "          INCAR                                             \n"
             "          SPECIES                                           \n"
             "          KPOINTS                                           \n"
             "          POSCAR                                            \n"
             "          $CURR_CALCTYPE/                                   \n"
             "            $CURR_REF/                                      \n"
             "              chemical_reference.json                       \n"
             "        $SCELNAME/                                          \n"
             "          LAT                                               \n"
             "          $CONFIGID/                                        \n"
             "            POS                                             \n"
             "            config.json                                     \n"
             "            structure.json                                  \n"
             "            $CURR_CALCTYPE/                                 \n"
             "              (DFT files)                                   \n"
             "              properties.calc.json                          \n"
             "      cluster_expansions/                                   \n"
             "        clex.formation_energy/                              \n"
             "          $CURR_BSET/                                       \n"
             "            $CURR_CALCTYPE/                                 \n"
             "              $CURR_REF/                                    \n"
             "                $CURR_ECI/                                  \n"
             "                  eci.json                                  \n\n"

             "    Variable descriptions:                                  \n\n"

             "    $ROOT: root directory of the CASM project               \n\n"
             "    $CURR_BSET: Current basis set, by default this is "
             "'bset.default'.\n"
             "    The current value can be inspected via 'casm settings -l'.   "
             "   \n"
             " \n"
             "    $CURR_CALCTYPE: Current calctype, by default this is "
             "'calctype.default'.\n"
             "    The current value can be inspected via 'casm settings -l'.   "
             "   \n"
             " \n"
             "    $CURR_REF: Current composition axes and reference states, by "
             "   \n"
             "    default this is 'ref.default'. The current value can be "
             "inspected\n"
             "    via 'casm settings -l'.                                 \n"
             " \n"
             "    $SCELNAME: Supercell name, in the form SCELV_A_B_C_D_E_F. "
             "'V' is\n"
             "    volume of the supercell in terms of the primitive cell, and  "
             "   \n"
             "    'A'-'F' are the values of the hermite normal form of the\n"
             "    transformation matrix.                                  \n"
             " \n"
             "    $CONFIGID: Configuration id, a unique integer.          \n"
             " \n"
             "    $TITLE: Title of the CASM project                       \n"
             " \n"
             "    Note: The 'settings' heirarchy can be located at the project "
             "   \n"
             "    level as shown above, or at the supercell or configuration "
             "level\n"
             "    in order to override calculation, composition, or reference  "
             "   \n"
             "    state settings at the supercell or configuration level.  The "
             "   \n"
             "    most local settings are always used for a configuration.\n"
             " \n";
  }

  if (vm.count("project_settings")) {
    log() << "\n### project_settings.json ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n";
    log() << "$ROOT/.casm/project_settings.json\n\n\n";

    log() << "DESCRIPTION:\n";
    log() << "Current CASM project settings.\n\n\n";

    log() << "EXAMPLE:\n";
    log() << "-------\n";
    log() << "{\n"
             "  \"cluster_expansions\" : {\n"
             "    \"formation_energy\" : {\n"
             "      \"bset\" : \"default\",\n"
             "      \"calctype\" : \"default\",\n"
             "      \"eci\" : \"default\",\n"
             "      \"name\" : \"formation_energy\",\n"
             "      \"property\" : \"formation_energy\",\n"
             "      \"ref\" : \"default\"\n"
             "    }\n"
             "  },\n"
             "  \"crystallography_tol\" : 1.0000000000000000e-05,\n"
             "  \"curr_properties\" : [ \"energy\" ],\n"
             "  \"default_clex\" : \"formation_energy\",\n"
             "  \"lin_alg_tol\" : 1.0000000000000000e-10,\n"
             "  \"name\" : \"ZrO\",\n"
             "  \"nlist_sublat_indices\" : [ 2, 3 ],\n"
             "  \"nlist_weight_matrix\" : [\n"
             "    [ 2, -1, 0 ],\n"
             "    [ -1, 2, 0 ],\n"
             "    [ 0, 0, 5 ]\n"
             "  ],\n"
             "  \"query_alias\" : {\n"
             "  },\n"
             "  \"view_command\" : \"casm.view \\\"open -a "
             "/Applications/VESTA/VESTA.app\\\"\"\n"
             "}"
          << std::endl;
    log() << "-------\n";
    log() << std::endl << std::endl;
  }

  if (vm.count("prim")) {
    log()
        << "\n### prim.json ##################\n\n"

           "LOCATION WHEN GENERATED:\n"
           "$ROOT/prim.json\n"
           "$ROOT/PRIM (legacy)\n\n\n"

           "DESCRIPTION:\n"
           "'prim.json' describes the primitive cell structure and defines the "
           " \n"
           "allowed degrees of freedom. It lists vectors, crystal basis sites, "
           " \n"
           "global degrees of freedom, and site degrees of freedom, including  "
           " \n"
           "allowed occupant species on each basis site.\n\n"

           "'prim.json' parameters:                                            "
           "\n\n"

           "\"title\" (string):                                                "
           "\n"
           "  A title for the project. Must consist of alphanumeric characters "
           "\n"
           "  and underscores only. The first character may not be a number.   "
           "\n\n"

           "\"lattice_vectors\" (JSON array of 3 JSON arrays of 3 numbers):    "
           "\n"
           "  Lattice vectors for the primitive structure, in Angstroms.       "
           "\n\n"

           "\"coordinate_mode\" (string):                                      "
           "\n"
           "  Coordinate mode for basis sites. One of:                         "
           "\n"
           "    \"Fractional\" or \"Direct\",                                  "
           "\n"
           "    \"Cartesian\"                                                  "
           "\n\n"

           "\"dofs\" [OPTIONAL] (JSON dictionary of DoF JSON objects):      \n"
           "  For each allowed type of global site DoF (typically strain), an  "
           "\n"
           "  object specifying how it is defined.                             "
           "\n\n\n"

           "\"basis\" (JSON array of JSON objects):                            "
           "\n\n"

           "  /\"coordinate\" (JSON array of 3 numbers):                       "
           "\n"
           "    Coordinate of the basis site with units as specified by the    "
           "\n"
           "    the \"coordinate_mode\" parameter. The default tolerance for   "
           "\n"
           "    checking symmetry is 1e-5, so basis site coordinates should    "
           "\n"
           "    include 6 significant digits or more.                          "
           "\n\n"

           "  /\"occupants\" (JSON array of string):                           "
           "\n"
           "    A list of the possible occupant species that may reside at "
           "each\n"
           "    site. The names are case sensitive, and \"Va\" is reserved for "
           "\n"
           "    vacancies.                                                     "
           "\n\n"

           "  /\"dofs\" [OPTIONAL] (JSON dictionary of DoF JSON objects):   \n"
           "    For each allowed type of site DoF (e.g., displacement, "
           "magnetic\n"
           "    spin, etc.) an object specifying how it is defined. Within "
           "\"dofs\"\n"
           "    object, each DoF is given by the key/object pair\n"
           "      \"$DOFNAME\" : {DOFOBJECT}\n"
           "    where $DOFNAME is the name specifier of a particular DoF type  "
           "\n"
           "    and the associated object specifies non-default options.\n\n\n"

           "\"species\" [OPTIONAL] (JSON dictionary of Molecule JSON "
           "objects):\n"
           "  A dictionary used to define extended properties of any species \n"
           "  listed as an allowed occupant in \"basis\"/\"occupants\". \n\n\n"

           "DoF JSON object:\n"
           "  DoFs are continuous vectors having a standard basis that   "
           "\n"
           "  is related to the fixed reference frame of the crystal. The DoF "
           "\n"
           "  object encodes a user-specified basis in terms of the standard "
           "basis.\n"
           "  User-specified basis may fully span the standard basis or only a "
           "  "
           "\n"
           "  subspace.\n"
           "  EXAMPLE: Site displacement DoF:\n"
           "    \"disp\" : {\n"
           "      \"axis_names\" : [\"dxy\", \"dz\"],\n"
           "      \"basis\" : [[1.0, 1.0, 0.0],\n"
           "                   [0.0, 0.0, 1.0]]\n"
           "    }\n\n\n"

           "Molecule JSON object:\n"
           "  Used to define species the comprise multiple atoms, "
           "off-centered\n"
           "  atoms, or species with properties such as a magnetic spin or or\n"
           "  charge state.\n\n"

           "  Allowed fields:\n"
           "  \"atoms\" (JSON Array of JSON objects):\n"
           "    /\"name\" (string):\n"
           "      Name of atomic species.\n\n"

           "    /\"coordinate\" (3D JSON Array of doubles):\n"
           "      Position of the atom, relative to the basis site at which "
           "it\n"
           "      is placed. Coordinate mode is same as rest of prim.json.\n\n"

           "    /\"properties\" [OPTIONAL] (JSON dictionary of "
           "SpeciesProperty)\n"
           "      Additonal fixed properties of the atom, such as magnetic "
           "moment,\n"
           "      charge state, or selective dynamics flags. The name of each\n"
           "      attribute must correspond to a CASM-supported "
           "SpeciesProperty.\n\n"

           "  /\"properties\" [OPTIONAL] (JSON dictionary of "
           "SpeciesProperty)\n"
           "    Additonal fixed properties of the molecule as a whole, such as "
           "\n"
           "    magnetic spin, charge state, or selective dynamics flags. The "
           "\n"
           "    name of each attribute must correspond to a CASM-supported \n"
           "    SpeciesProperty.\n\n\n"

           "  /\"name\" [OPTIONAL] (string)\n"
           "    Chemical name of molecule, used to override its name in the\n"
           "    enclosing dictionary. This name is used for chemical "
           "comparisons\n"
           "    of molecules. Crystallographic and spatial comparisons are "
           "not\n"
           "    dependent on molecule names.\n\n\n"

           "SpeciesProperty JSON object:\n"
           "  Associates the discrete value of a vector property to an Atom or "
           "Moleule.\n\n"
           "  Allowed fields:\n"
           "  \"value\" (JSON Array of doubles):\n"
           "    Dimension of array must match the dimension of the specified \n"
           "    SpeciesProperty.\n\n";

    log() << "prim.json EXAMPLE 1: FCC ternary alloy of elements A, B, and C\n"
             "-------\n"
             "{\n"
             "  \"basis\" : [\n"
             "    {\n"
             "      \"coordinate\" : [ 0.000000000000, 0.000000000000, "
             "0.000000000000 ],\n"
             "      \"occupants\" : [ \"A\", \"B\", \"C\" ]\n"
             "    }\n"
             "  ],\n"
             "  \"coordinate_mode\" : \"Fractional\",\n"
             "  \"description\" : \"Face-centered Cubic (FCC, cF)\",\n"
             "  \"lattice_vectors\" : [\n"
             "    [ 2.000000000000, 2.000000000000, 0.000000000000 ],\n"
             "    [ 0.000000000000, 2.000000000000, 2.000000000000 ],\n"
             "    [ 2.000000000000, 0.000000000000, 2.000000000000 ]\n"
             "  ],\n"
             "  \"title\" : \"ABC\"\n"
             "}\n";
    log() << "-------\n\n";

    log() << "prim.json EXAMPLE 2: FCC binary alloy of elements A and B with "
             "Green-Lagrange strain DoF\n"
             "-------\n"
             "{\n"
             "  \"basis\" : [\n"
             "    {\n"
             "      \"coordinate\" : [ 0.000000000000, 0.000000000000, "
             "0.000000000000 ],\n"
             "      \"occupants\" : [ \"A\", \"B\" ]\n"
             "    }\n"
             "  ],\n"
             "  \"dofs\" : {\n"
             "    \"strainGL\" : {},\n"
             "  }\n"
             "  \"coordinate_mode\" : \"Fractional\",\n"
             "  \"description\" : \"Face-centered Cubic (FCC, cF)\",\n"
             "  \"lattice_vectors\" : [\n"
             "    [ 2.000000000000, 2.000000000000, 0.000000000000 ],\n"
             "    [ 0.000000000000, 2.000000000000, 2.000000000000 ],\n"
             "    [ 2.000000000000, 0.000000000000, 2.000000000000 ]\n"
             "  ],\n"
             "  \"title\" : \"AB_with_Strain\"\n"
             "}\n";
    log() << "-------\n\n";

    log() << "prim.json EXAMPLE 3: Zincblende GaAs with Hencky strain DoF "
             "along (x,x), (y,y), and (z,z) components\n"
             "-------\n"
             "{\n"
             "  \"basis\" : [\n"
             "    {\n"
             "      \"coordinate\" : [ 0.000000000000, 0.000000000000, "
             "0.000000000000 ],\n"
             "      \"occupants\" : [ \"Ga\" ]\n"
             "    },\n"
             "    {\n"
             "      \"coordinate\" : [ 0.250000000000, 0.250000000000, "
             "0.250000000000 ],\n"
             "      \"occupants\" : [ \"As\" ]\n"
             "    }\n"
             "  ],\n"
             "  \"dofs\" : {\n"
             "    \"strainH\" : {\n"
             "      \"axis_names\" : [ \"Exx\", \"Eyy\", \"Ezz\" ],\n"
             "      \"basis\" : [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],\n"
             "                   [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],\n"
             "                   [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]]\n"
             "    },\n"
             "  }\n"
             "  \"coordinate_mode\" : \"Fractional\",\n"
             "  \"description\" : \"Zincblende GaAs with deviatoric strain\",\n"
             "  \"lattice_vectors\" : [\n"
             "    [ 0.00000, 2.82663, 2.82663 ],\n"
             "    [ 2.82663, 0.00000, 2.82663 ],\n"
             "    [ 2.82663, 2.82663, 0.00000 ]\n"
             "  ],\n"
             "  \"title\" : \"GaAs\"\n"
             "}\n";
    log() << "-------\n\n";

    log() << "prim.json EXAMPLE 4: BCC lithium with atomic displacement DoF\n"
             "-------\n"
             "{\n"
             "  \"basis\" : [\n"
             "    {\n"
             "      \"coordinate\" : [ 0.000000000000, 0.000000000000, "
             "0.000000000000 ],\n"
             "      \"dofs\" : {\n"
             "        \"disp\" : {},\n"
             "      }\n"
             "      \"occupants\" : [ \"Li\" ]\n"
             "    }\n"
             "  ],\n"
             "  \"coordinate_mode\" : \"Fractional\",\n"
             "  \"description\" : \"Body-centered Cubic (BCC, cI)\",\n"
             "  \"lattice_vectors\" : [\n"
             "    [ -1.75445,  1.75445,  1.75445 ],\n"
             "    [  1.75445, -1.75445,  1.75445 ],\n"
             "    [  1.75445,  1.75445, -1.75445 ]\n"
             "  ],\n"
             "  \"title\" : \"Li_with_displacement\"\n"
             "}\n";
    log() << "-------\n\n";

    log() << "prim.json EXAMPLE 5: HCP Zr with O insertion at octahedral "
             "interstitial positions\n"
             "-------\n\n"
             "{\n"
             "  \"basis\" : [\n"
             "    {\n"
             "      \"coordinate\" : [ 0.0, 0.0, 0.0 ],\n"
             "      \"occupants\" : [ \"Zr\" ]\n"
             "    },\n"
             "    {\n"
             "      \"coordinate\" : [ 0.666666, 0.333333, 0.5 ],\n"
             "      \"occupants\" : [ \"Zr\" ]\n"
             "    },\n"
             "    {\n"
             "      \"coordinate\" : [ 0.333333, 0.666666, 0.25 ],\n"
             "      \"occupants\" : [ \"Va\", \"O\" ]\n"
             "    },\n"
             "    {\n"
             "      \"coordinate\" : [ 0.333333, 0.666666, 0.75 ],\n"
             "      \"occupants\" : [ \"Va\", \"O\" ]\n"
             "    }\n"
             "  ],\n"
             "  \"coordinate_mode\" : \"Fractional\",\n"
             "  \"description\" : \"hcp Zr with oct (O) \",\n"
             "  \"lattice_vectors\" : [\n"
             "    [  3.233986860000, 0.000000000000, 0.000000000000 ],\n"
             "    [ -1.616993430000, 2.800714770000, 0.000000000000 ],\n"
             "    [ -0.000000000000, 0.000000000000, 5.168678340000 ]\n"
             "  ],\n"
             "  \"title\" : \"ZrO\"\n"
             "}\n";
    log() << "-------\n\n";

    log() << "\n### PRIM ##################\n\n"

             "DESCRIPTION:\n"
             "PRIM is the input file used by previous version of casm. It can "
             "be read and        \n"
             "converted to 'prim.json'. The format of PRIM is very similar to "
             "the VASP POSCAR    \n"
             "except a list of possible occupant molecules is included with "
             "each basis site.     \n\n"

             "- Molecule names are case sensitive.\n"
             "- 'Va' is reserved for vacancies.\n"
             "- The default tolerance for checking symmetry is 1e-5, so basis "
             "site coordinates\n"
             "  should include 6 significant digits or more.\n\n\n"

             "EXAMPLE 1: An FCC ternary alloy of elements A, B, and C\n"
             "-------\n"

             "  Face-centered Cubic (FCC, cF)\n"
             "   1.0\n"
             "     0.00 2.00 2.00\n"
             "     2.00 0.00 2.00\n"
             "     2.00 2.00 0.00\n"
             "   1\n"
             "  Direct\n"
             "     0.00 0.00 0.00 A B C\n"
             "-------\n\n";
    log() << "EXAMPLE 2: A structure with some alloying sites and some "
             "non-alloying sites\n"
             "-------\n"
             "  LiTiO2 - Bronze\n"
             "   1.00000000\n"
             "     1.91357600      0.00000000     -6.23799200\n"
             "     6.08935000     -1.87060000      0.00000000\n"
             "     0.00000000     -3.74120000      0.00000000\n"
             "   5 4 8\n"
             "  Direct\n"
             "     0.0000000   0.0000000   0.0000000 Li Va\n"
             "     0.3800000   0.9000000   0.5500000 Li Va\n"
             "     0.6200000   0.1000000   0.4500000 Li Va\n"
             "     0.0000000   0.2600000   0.3700000 Li Va\n"
             "     0.0000000   0.7400000   0.6300000 Li Va\n"
             "     0.7080000   0.3940000   0.8030000 Ti\n"
             "     0.2920000   0.6060000   0.1970000 Ti\n"
             "     0.2950000   0.1980000   0.9010000 Ti\n"
             "     0.7050000   0.8020000   0.0990000 Ti\n"
             "     0.9960000   0.2640000   0.8680000 O\n"
             "     0.0040000   0.7360000   0.1320000 O\n"
             "     0.3470000   0.5280000   0.7360000 O\n"
             "     0.6530000   0.4720000   0.2640000 O\n"
             "     0.6290000   0.1200000   0.9400000 O\n"
             "     0.3710000   0.8800000   0.0600000 O\n"
             "     0.7070000   0.7240000   0.6380000 O\n"
             "     0.2930000   0.2760000   0.3620000 O\n"
             "-------\n";
    log() << std::endl << std::endl;
  }

  if (vm.count("sym")) {
    log() << "\n### sym ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n";
    log() << "$ROOT/symmetry/lattice_point_group.json\n";
    log() << "$ROOT/symmetry/factor_group.json\n";
    log() << "$ROOT/symmetry/crystal_point_group.json\n\n\n";

    log() << "DESCRIPTION:\n"
             "Symgroup files report each element of a group as            \n"
             "representative linear transformations of 3-dimensional space\n"
             "(i.e., the symmetry operation). The symmetery operation     \n"
             "transforms a spatial coordinate 'x' according to x' = A*x+b,\n"
             "where 'A' is the 3x3 'operation matrix' and 'b' is the      \n"
             "'shift' vector.                                             \n\n";

    log() << "For a crystal or PRIM file, CASM reports the following      \n"
             "groups:\n\n"

             "lattice_point_group.json:                                   \n"
             "  This is the point group of the Bravais lattice, and is the\n"
             "  list of operations that map the lattice vectors onto      \n"
             "  themselves. The 'shift' vectors will always be zero.      \n\n"

             "factor_group.json:                                          \n"
             "  This is a finite description of the crystal spacegroup, in\n"
             "  which all redundant operations that differ only by a 'shift'\n"
             "  are represented by a single operation, whose 'shift' lies \n"
             "  within the primitive cell. Formally, this is a group      \n"
             "  formed by the cosets of 'T' in 'S', where 'T' is the      \n"
             "  translation group of the Bravais lattice and 'S' is the   \n"
             "  crystal space group.                                      \n\n"

             "crystal_point_group.json:                                   \n"
             "  This is a group of point operations formed by taking the  \n"
             "  factor group operations and setting their 'shifts' to     \n"
             "  zero. Macroscopic properties of the crystal must exhibit  \n"
             "  the symmetries of the crystal point group. It is, by      \n"
             "  definition a subgroup of the lattice point group.         \n";

    log() << "\n\n";
  }

  if (vm.count("vasp")) {
    log() << "\n### vasp ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n\n";

    log() << "INPUT SETTINGS:\n";
    log() << "$CALC_SETTINGS/relax.json\n";
    log() << "$CALC_SETTINGS/INCAR\n";
    log() << "$CALC_SETTINGS/SPECIES\n";
    log() << "$CALC_SETTINGS/KPOINTS\n";
    log() << "$CALC_SETTINGS/POSCAR\n\n";

    log() << "For global settings:\n";
    log() << "  CALC_SETTINGS = $ROOT/training_data/settings/$CURR_CALCTYPE\n";
    log() << "For supercell specific settings:\n";
    log() << "  CALC_SETTINGS = "
             "$ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE\n";
    log() << "For configuration specific settings:\n";
    log() << "  CALC_SETTINGS = "
             "$ROOT/training_data/$SCELNAME/$CONFIGID/settings/"
             "$CURR_CALCTYPE\n\n";

    log() << "RESULTS:\n";
    log() << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/(VASP "
             "results)\n";
    log() << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/"
             "properties.calc.json (read)\n";

    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log() << "CASM comes with wrappers for using VASP to calculate the "
             "properties \n"
             "of configurations, but is designed so that any type of "
             "calculation  \n"
             "software or method could be used if an appropriate set of "
             "wrapper   \n"
             "scripts are available. By convention, input settings for "
             "software   \n"
             "used to calculate the properties of a particular configuration   "
             "   \n"
             "should be checked for in the following directories:         \n"
             "  1) "
             "$ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE\n"
             "  2) $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE  \n"
             "  3) $ROOT/training_data/settings/$CURR_CALCTYPE            \n\n"

             "The VASP wrappers included with CASM check for input settings "
             "files \n"
             "in the above directories, using the most local settings for a    "
             "   \n"
             "particular configuration. In most cases, the global settings "
             "files  \n"
             "are stored in $ROOT/training_data/settings/$CURR_CALCTYPE and "
             "used  \n"
             "for all configurations. Settings files are searched for on a "
             "file-by-file\n"
             "basis, so to set supercell or configuration specific settings it "
             "is \n"
             "sufficient to only include the particular files necessary in the "
             "   \n"
             "supercell or configuration level settings folder.           \n\n"

             "PBS job submission using the VASP wrappers depends on using the "
             "pbs \n"
             "python module available here: "
             "https://github.com/prisms-center/pbs  \n\n"

             "Included with CASM, the 'vasp.relax' script can be executed by "
             "the  \n"
             "'casm run' command to submit a batch of VASP jobs that for "
             "selected \n"
             "configurations. For each selected configuration, VASP is re-run  "
             "   \n"
             "using the output of the previous calculation until full "
             "convergence \n"
             "is achieved. The convergence criteria is: if the cell shape and  "
             "   \n"
             "volume remain constant (ISIF=0, 1, or 2) then a single "
             "calculation  \n"
             "is performed; else the calculation is converged if at least 2 "
             "jobs  \n"
             "are complete, and: 1) the last job completed with <= 3 ionic "
             "steps  \n"
             " or, if \"nrg_convergence\" is set in the 'relax.json' file, 2) "
             "the \n"
             "last two calculations had final E0 differ by less than the value "
             "of \n"
             " \"nrg_convergence\". Once converged, a final constant volume    "
             "   \n"
             "calculation is performed with the following INCAR settings: "
             "(ISIF=2,\n"
             "ISMEAR=-5, NSW=0, IBRION=-1).                               \n\n"

             "relax.json:                                                 \n"
             "  This JSON file contains a single JSON object which contains    "
             "   \n"
             "  parameters used to control PBS job submission settings.   \n"
             "  Required keys are:                                        \n"
             "    \"queue\": queue to submit job in                       \n"
             "    \"ppn\": processors (cores) per node to request         \n"
             "    \"atom_per_proc\": max number of atoms per processor (core)  "
             "   \n"
             "    \"walltime\": walltime to request (ex. \"48:00:00\")    \n\n"

             " Optional keys are:                                         \n"
             "    \"account\": account to submit job under (default None) \n"
             "    \"pmem\": string for requested memory (default None)    \n"
             "    \"priority\": requested job priority (default \"0\")    \n"
             "    \"message\": when to send messages about jobs (ex. \"abe\",  "
             "   \n"
             "               default \"a\")                               \n"
             "    \"email\": where to send messages (ex. \"me@fake.com\", "
             "default \n"
             "             None)                                          \n"
             "    \"qos\": quality of service, 'qos' option (ex. \"fluxoe\")   "
             "   \n"
             "    \"npar\": vasp incar setting (default None)             \n"
             "    \"ncore\": vasp incar setting (default None)            \n"
             "    \"kpar\": vasp incar setting (default None)             \n"
             "    \"vasp_cmd\": vasp execution command (default is \"vasp\" if "
             "   \n"
             "                ncpus=1, else \"mpirun -np {NCPUS} vasp\"   \n"
             "    \"ncpus\": number of cpus (cores) to run on (default "
             "$PBS_NP)   \n"
             "    \"run_limit\": number of vasp runs until \"not_converging\"  "
             "   \n"
             "                 (default 10)                               \n"
             "    \"nrg_convergence\": converged if last two runs complete and "
             "   \n"
             "                       differ in energy by less than this amount "
             "   \n"
             "                       (default None)                       \n"
             "    \"move\": files to move at the end of a run (ex. \"POTCAR\", "
             "   \n"
             "            \"WAVECAR\"], default [\"POTCAR\"])             \n"
             "    \"copy\": files to copy from run to run (ex. [\"INCAR\",\n"
             "            \"KPOINTS\"], default [\"INCAR, KPOINTS\"])     \n"
             "    \"remove\": files to remove at the end of a run (ex. "
             "[\"IBZKPT\",\n"
             "              \"CHGCAR\"], default [\"IBKZPT\", \"CHG\", "
             "\"CHGCAR\",\n"
             "              \"WAVECAR\", \"TMPCAR\", \"EIGENVAL\", \"DOSCAR\", "
             "   \n"
             "              \"PROCAR\", \"PCDAT\", \"XDATCAR\", \"LOCPOT\", "
             "\"ELFCAR\",\n"
             "              \"PROOUT\"]                                   \n"
             "    \"compress\": files to compress at the end of a run (ex.\n"
             "                [\"OUTCAR\", \"vasprun.xml\"], default [])  \n"
             "    \"backup\": files to compress to backups at the end of a "
             "run,   \n"
             "              used in conjunction with move (ex. [\"WAVECAR\"])  "
             "   \n"
             "    \"encut\": [START, STOP, STEP] values for converging ENCUT "
             "to   \n"
             "             within nrg_convergence (ex. [\"450\", \"Auto\",\n"
             "             \"10\"], default [\"Auto\", \"Auto\", \"10\"] where "
             "   \n"
             "             \"Auto\" is either the largest ENMAX in all POTCARS "
             "   \n"
             "             called in SPECIES for START, or 2.0 * largest ENMAX "
             "   \n"
             "             for STOP)                                      \n"
             "    \"kpoints\": [start, stop, step] values for converging "
             "KPOINTS  \n"
             "               to within nrg_convergence (ex. [\"5\", \"50\", "
             "\"1\"],\n"
             "               default [\"5\", \"Auto\", \"1\"] where \"Auto\" "
             "can  \n"
             "               only be used for STOP and means to keep "
             "increasing   \n"
             "               the KPOINTS length by STEP until either      \n"
             "               nrg_convergence or walltime is reached). For "
             "meaning \n"
             "               of the KPOINTS length parameter see the VASP \n"
             "               documentation at "
             "http://cms.mpi.univie.ac.at/vasp/   \n"
             "               vasp/Automatic_k_mesh_generation.html        \n"
             "    \"extra_input_files\": extra input files to be copied from "
             "the  \n"
             "                         settings directory, e.g., a vdW kernel  "
             "   \n"
             "                         file.                              \n"
             "    \"initial\": location of INCAR with tags for the initial "
             "run,   \n"
             "               if desired (e.g. to generate a PBE WAVECAR for "
             "use   \n"
             "               with M06-L)                                  \n"
             "    \"final\": location of INCAR with tags for the final run, if "
             "   \n"
             "             desired (e.g. \"ISMEAR = -5\", etc). Otherwise, the "
             "   \n"
             "             settings enforced are (\"ISMEAR = -5\", \"NSW = "
             "0\",   \n"
             "              \"IBRION = -1\", \"ISIF = 2\")                \n"
             "    \"err_types\": list of errors to check for. Allowed entries "
             "are \n"
             "                 \"IbzkptError\" and \"SubSpaceMatrixError\".    "
             "   \n"
             "                 Default: [\"SubSpaceMatrixError\"]         \n"
             "\n";

    log() << "EXAMPLE: relax.json \n";
    log() << "-------\n";
    log() << "{\n"
             "  \"account\":\"prismsprojectdebug_flux\",\n"
             "  \"queue\":\"flux\",\n"
             "  \"priority\":\"-200\",\n"
             "  \"walltime\":\"1:00:00\",\n"
             "  \"pmem\":\"3800mb\",\n"
             "  \"email\":\"username@univ.edu\",\n"
             "  \"ppn\":\"16\",\n"
             "  \"atom_per_proc\":\"2\",\n"
             "  \"run_limit\":10,\n"
             "  \"nrg_convergence\":0.002\n"
             "}\n";
    log() << "-------\n\n\n";

    log() << "SPECIES:                                                    \n"
             "  This file contains information for selecting POTCARs and "
             "specifing\n"
             "  parameters that must be set on an atom-by-atom basis in the "
             "INCAR,\n"
             "  such as MAGMOM. The first line in the file specifies the value "
             "of \n"
             "  'POTCAR_DIR_PATH', which is the base path used to find POTCAR  "
             "   \n"
             "  files. The second line contains column headings (at least 4), "
             "and \n"
             "  then there are lines for each distinct species. The first "
             "column  \n"
             "  specifies the 'SPECIES' and must match a species names in the "
             "PRIM\n"
             "  file. The second column gives an 'ALIAS' name for the species "
             "which\n"
             "  is used for ordering like atoms in the generated POSCAR files. "
             "The\n"
             "  third column should be either '0' or '1', such that only one   "
             "   \n"
             "  species with a given ALIAS has a '1'. For that species the "
             "fourth \n"
             "  column must contain the path that should be appended to the    "
             "   \n"
             "  POTCAR_DIR_PATH to specify the POTCAR file for that species.   "
             "   \n\n"

             "  Additional columns, such as 'MAGMOM' in the example below are  "
             "   \n\n"
             "  and used to specify the value used for a particular species in "
             "the\n"
             "  INCAR. The column heading must match a valid VASP INCAR "
             "setting.  \n\n";

    log() << "EXAMPLE: SPECIES \n";
    log() << "-------\n"
             "POTCAR_DIR_PATH = /absolute/path/to/vasp_potentials\n"
             "SPECIES    ALIAS    POTCAR  POTCAR_location    MAGMOM\n"
             "Mn3        Mn       0       -                  3\n"
             "Mn4        Mn       1       PAW_PBE/Mn         4\n";
    log() << "-------\n\n\n";

    log() << "INCAR:                                                      \n"
             "  This is a template INCAR used for VASP calculations. The "
             "settings \n"
             "  are generally used as given though some may be automatically "
             "set  \n"
             "  based on settings in the 'relax.json' or 'SPECIES' files. "
             "Also,   \n"
             "  some settings might be added or changed if certain errors are  "
             "   \n"
             "  during calculation. The actual INCAR used for each calculation "
             "is \n"
             "  saved.                                                    \n\n";

    log() << "EXAMPLE: INCAR \n";
    log() << "-------\n"
             "System = Test of VASP submission\n"
             "ISPIN = 1 # non-spin polarized\n"
             "PREC = Accurate \n"
             "IBRION = 2 # conjugate gradient ionic minimization\n"
             "NSW = 61\n"
             "ISIF= 3 # relax ions and volume\n"
             "ENMAX = 400 \n"
             "ISMEAR = 1 # for metals\n"
             "SIGMA = 0.2 \n"
             "LWAVE = .FALSE.\n"
             "LCHARG = .FALSE.\n";
    log() << "-------\n\n\n";

    log() << "KPOINTS, POSCAR:                                            \n"
             "  This is a template KPOINTS file used for VASP calculations. If "
             "the\n"
             "  mode (third line) is set to 'Auto', this file is used as is "
             "for all\n"
             "  VASP calculations.  Otherwise, if the mode is 'Gamma' or 'M', "
             "a   \n"
             "  reference POSCAR must also exist and a scaling method is used "
             "to  \n"
             "  calculate the kpoint mesh for a supercell, such that it has an "
             "   \n"
             "  equal or greater kpoint density than in the reference POSCAR.  "
             "   \n\n";

    log() << "-------\n";
  }

  if (vm.count("qe")) {
    log() << "\n### quantum espresso ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n\n";

    log() << "INPUT SETTINGS:\n";
    log() << "$CALC_SETTINGS/relax.json\n";
    log() << "$CALC_SETTINGS/$CUSTOM_INFILE_NAME\n";
    log() << "$CALC_SETTINGS/SPECIES\n";

    log() << "For global settings:\n";
    log() << "  CALC_SETTINGS = $ROOT/training_data/settings/$CURR_CALCTYPE\n";
    log() << "For supercell specific settings:\n";
    log() << "  CALC_SETTINGS = "
             "$ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE\n";
    log() << "For configuration specific settings:\n";
    log() << "  CALC_SETTINGS = "
             "$ROOT/training_data/$SCELNAME/$CONFIGID/settings/"
             "$CURR_CALCTYPE\n\n";

    log() << "RESULTS:\n";
    log() << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/(quantum "
             "espresso results)\n";
    log() << "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/"
             "properties.calc.json (read)\n";

    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log()
        << "CASM comes with wrappers for using Quantum Espresso to calculate "
           "the properties \n"
           "of configurations, but is designed so that any type of calculation "
           " \n"
           "software or method could be used if an appropriate set of wrapper  "
           " \n"
           "scripts are available. By convention, input settings for software  "
           " \n"
           "used to calculate the properties of a particular configuration     "
           " \n"
           "should be checked for in the following directories:                "
           " \n"
           "  1) "
           "$ROOT/training_data/$SCELNAME/$CONFIGID/settings/$CURR_CALCTYPE\n"
           "  2) $ROOT/training_data/$SCELNAME/settings/$CURR_CALCTYPE         "
           " \n"
           "  3) $ROOT/training_data/settings/$CURR_CALCTYPE                   "
           " \n\n"

           "The Quantum Espresso wrappers included with CASM check for input "
           "settings files \n"
           "in the above directories, using the most local settings for a      "
           " \n"
           "particular configuration. In most cases, the global settings files "
           " \n"
           "are stored in $ROOT/training_data/settings/$CURR_CALCTYPE and used "
           " \n"
           "for all configurations. Settings files are searched for on a "
           "file-by-file\n"
           "basis, so to set supercell or configuration specific settings it "
           "is \n"
           "sufficient to only include the particular files necessary in the   "
           " \n"
           "supercell or configuration level settings folder.                  "
           " \n\n"

           "PBS job submission using the Quantum Espresso wrappers depends on "
           "using the pbs \n"
           "python module available here: https://github.com/prisms-center/pbs "
           " \n\n"

           "Included with CASM, the 'qe.relax' script can be executed by the  "
           "\n"
           "'casm run' command to submit a batch of Quantum Espresso jobs that "
           "for selected \n"
           "configurations. For each selected configuration, Quantum Espresso "
           "is re-run\n"
           "using the output of the previous calculation until full "
           "convergence \n"
           "is achieved. The convergence criteria is: if the cell shape and    "
           " \n"
           "volume remain constant (calculation != vc-relax) then a single "
           "calculation  \n"
           "is performed; else the calculation is converged if at least 2 jobs "
           " \n"
           "are complete, and: 1) the last job completed with <= 3 ionic steps "
           " \n"
           " or, if \"nrg_convergence\" is set in the 'relax.json' file, 2) "
           "the \n"
           "last two calculations had final energy differ by less than the "
           "value of \n"
           " \"nrg_convergence\". Once converged, a final constant volume      "
           " \n"
           "calculation is performed with the following setting: (calculation "
           "= 'relax')\n"

           "relax.json:                                                        "
           " \n"
           "  This JSON file contains a single JSON object which contains      "
           " \n"
           "  parameters used to control PBS job submission settings.          "
           " \n"
           "  Required keys are:                                               "
           " \n"
           "    \"queue\": queue to submit job in                              "
           " \n"
           "    \"ppn\": processors (cores) per node to request                "
           " \n"
           "    \"atom_per_proc\": max number of atoms per processor (core)    "
           " \n"
           "    \"walltime\": walltime to request (ex. \"48:00:00\")           "
           " \n\n"
           "    \"software\": needs to be quantumespresso for quantum espresso "
           "to be used\n\n"

           " Optional keys are:                                                "
           " \n"
           "    \"account\": account to submit job under (default None)        "
           " \n"
           "    \"pmem\": string for requested memory (default None)           "
           " \n"
           "    \"priority\": requested job priority (default \"0\")           "
           " \n"
           "    \"message\": when to send messages about jobs (ex. \"abe\",    "
           " \n"
           "               default \"a\")                                      "
           " \n"
           "    \"email\": where to send messages (ex. \"me@fake.com\", "
           "default \n"
           "             None)                                                 "
           " \n"
           "    \"qos\": quality of service, 'qos' option (ex. \"fluxoe\")     "
           " \n"
           "    \"qe_cmd\": quantum espresso execution command (default is "
           "\"pw.x < {INFILE} > {OUTFILE}\" if    \n"
           "                ncpus=1, else \"mpirun -np {NCPUS} pw.x < {INFILE} "
           "> {OUTFILE}\"           \n"
           "    \"infile\": quantum espresso input file name (default is "
           "\"std.in\"\n"
           "    \"outfile\": quantum espresso output file name (default is "
           "\"std.out\"\n"
           "    \"ncpus\": number of cpus (cores) to run on (default $PBS_NP)  "
           " \n"
           "    \"run_limit\": number of vasp runs until \"not_converging\"    "
           " \n"
           "                 (default 10)                                      "
           " \n"
           "    \"nrg_convergence\": converged if last two runs complete and   "
           " \n"
           "                       differ in energy by less than this amount   "
           " \n"
           "                       (default None)                              "
           " \n"
           "    \"move\": files to move at the end of a run (ex. \"\",    \n"
           "            \".wfc\"], default [])                     \n"
           "    \"copy\": files to copy from run to run  default "
           "[$infilename]) \n"
           "    \"remove\": files to remove at the end of a run                "
           "\n"
           "               default [\".wfc\", \".igk\", \".save\"]             "
           "\n"
           "    \"compress\": files to compress at the end of a run (ex.       "
           " \n"
           "                [$outfilename], default [])          \n"
           "    \"backup\": files to compress to backups at the end of a run,  "
           " \n"
           "              used in conjunction with move (ex. [\".wfc\"])     \n"
           "    \"extra_input_files\": extra input files to be copied from the "
           " \n"
           "                         settings directory, e.g., an OCCUPATIONS  "
           "   \n"
           "                         file.                                     "
           " \n"
           "    \"initial\": location of $infile with tags for the initial "
           "run,   \n"
           "               if desired                                          "
           "  \n"
           "    \"final\": location of $infile with tags for the final run, if "
           " \n"
           "             desired                                               "
           " \n"
           "    \"err_types\": list of errors to check for. Not Implemented "
           "yet  \n"
           "\n";

    log() << "EXAMPLE: relax.json \n";
    log() << "-------\n";
    log() << "{\n"
             "  \"account\":\"prismsprojectdebug_flux\",\n"
             "  \"queue\":\"flux\",\n"
             "  \"priority\":\"-200\",\n"
             "  \"walltime\":\"1:00:00\",\n"
             "  \"pmem\":\"3800mb\",\n"
             "  \"email\":\"username@univ.edu\",\n"
             "  \"ppn\":\"16\",\n"
             "  \"atom_per_proc\":\"2\",\n"
             "  \"run_limit\":10,\n"
             "  \"nrg_convergence\":0.002\n"
             "  \"calculator\":\"quantumespresso\"\n"
             "  \"infilename\":\"LixCoO2.in\"\n"
             "  \"outfilename\":\"LixCoO2.out\"\n"
             "}\n";
    log() << "-------\n\n\n";

    log()
        << "SPECIES:                                                           "
           " \n"
           "  This file contains information for selecting pseudopotentials "
           "and specifing\n"
           "  parameters that must be set on an atom-by-atom basis in the "
           "infile,\n"
           "  such as magnetic moment (non currently implemented).\n"
           "  The first line in the file specifies the value of \n"
           "  'PSEUDO_DIR_PATH', which is the base path used to find UPF     \n"
           "  files. The second line contains column headings (at least 4), "
           "and \n"
           "  then there are lines for each distinct species. The first column "
           " \n"
           "  specifies the 'SPECIES' and must match a species names in the "
           "PRIM\n"
           "  file. The second column gives an 'ALIAS' name for the species "
           "which\n"
           "  is used for ordering like atoms in the generated input files. "
           "The\n"
           "  third column should be either '0' or '1', such that only one     "
           " \n"
           "  species with a given ALIAS has a '1'. For that species the "
           "fourth \n"
           "  column must contain the path that should be appended to the      "
           " \n"
           "  PSEUDO_DIR_PATH to specify the UPF file for that species.      "
           "\n\n"

           "  Additional columns, such as 'if_pos' in the example below are    "
           " \n\n"
           "  and used to specify the value used for a particular species in "
           "the\n"
           "  infile. The column heading must match a valid quantum espresso "
           "input setting.\n"
           "  For now only supported additional tag is if_pos, a way to fixed "
           "certain lattice positions.\n\n";

    log() << "EXAMPLE: SPECIES \n";
    log() << "-------\n"
             "PSEUDO_DIR_PATH = /absolute/path/to/quantumespresso_potentials\n"
             "SPECIES    ALIAS    UPF  UPF_location     if_pos\n"
             "Ni         Ni       1       PAW_PBE/Ni.UPF     1,1,1\n"
             "Al        Al       1       PAW_PBE/Al.UPF      1,1,1\n";
    log() << "-------\n\n\n";

    log() << "$infilename:                                                     "
             " \n"
             "  This is a template input file used for Quantum Espresso "
             "calculations. The settings \n"
             "  are generally used as given though some may be automatically "
             "set  \n"
             "  based on settings in the 'relax.json' or 'SPECIES' files. "
             "Also,   \n"
             "  some settings might be added or changed if certain errors are  "
             "   \n"
             "  during calculation. The actual input file used for each "
             "calculation is \n"
             "  saved.                                                    \n\n";
    log()
        << "Note:                                                    \n"
           "  K_POINTS will be adjusted accordingly such that the density is "
           "maintained\n"
           "  over all configurations in the project for all Quantum Espresso "
           "calculations\n"
           "  this uses the CELL_PARAMETERS and the K_POINTS cards in the "
           "input file to calculate\n"
           "  a density and rescale configurations k-point mesh accordingly\n";

    log() << "EXAMPLE: Mg2Ti4S8.in \n";
    log() << "-------\n"
             "System = Test of Quantum Espresso submission\n"
             "&CONTROL\n"
             " calculation = 'vc-relax',\n"
             " pseudo_dir = '/home/skolli/quantum_espresso/pseudo/',\n"
             " tprnfor = .true.,\n"
             " prefix = 'Mg2Ti4S8',\n"
             " restart_mode = 'from_scratch',\n"
             " tstress = .true.,\n"
             "/\n"
             "&SYSTEM\n"
             " ecutwfc = 45.0,\n"
             " occupations = 'fixed',\n"
             " celldm(1) = 7.3794,\n"
             " ibrav = 0,\n"
             " nat = 14,\n"
             " ntyp = 3,\n"
             " ecutrho = 200.0,\n"
             "/\n"
             "&ELECTRONS\n"
             " diagonalization = 'cg',\n"
             " mixing_mode = 'plain',\n"
             " mixing_beta = 0.7,\n"
             " conv_thr = 1e-08,\n"
             "/\n"
             "&IONS\n"
             " ion_dynamics = 'bfgs',\n"
             "/\n"
             "&CELL\n"
             " press = 0.1,\n"
             " cell_factor = 1.6,\n"
             " cell_dynamics = 'bfgs',\n"
             "/\n"
             "\n"
             "ATOMIC_SPECIES\n"
             " Mg 24.31 Mg.pbe-nsp-bpaw.UPF\n"
             " Ti 47.88 Ti.pbe-sp-hgh.UPF\n"
             " S 32.07 S.pbe-n-kjpaw_psl.0.1.UPF\n"
             "\n"
             "CELL_PARAMETERS angstrom\n"
             " 0.0000000000000000 5.1756022947592379 5.1756022947592388\n"
             " 5.1756022947592388 0.0000000000000000 5.1756022947592388\n"
             " 5.1756022947592388 5.1756022947592379 0.0000000000000000\n"
             "\n"
             "ATOMIC_POSITIONS crystal\n"
             "Mg 0.000000000 0.000000000 0.000000000\n"
             "Mg 0.250000000 0.250000000 0.250000000\n"
             "Ti 0.625000000 0.625000000 0.625000000\n"
             "Ti 0.125000000 0.625000000 0.625000000\n"
             "Ti 0.625000000 0.125000000 0.625000000\n"
             "Ti 0.625000000 0.625000000 0.125000000\n"
             "S 0.3842989149764762 0.3842989149764762 0.3842989149764762\n"
             "S 0.8657010850235238 0.8657010850235238 0.8657010850235238\n"
             "S 0.3842989149764762 0.8471032550705786 0.3842989149764762\n"
             "S 0.3842989149764762 0.3842989149764762 0.8471032550705786\n"
             "S 0.8471032550705786 0.3842989149764762 0.3842989149764762\n"
             "S 0.8657010850235238 0.8657010850235238 0.4028967449294214\n"
             "S 0.8657010850235238 0.4028967449294214 0.8657010850235238\n"
             "S 0.4028967449294214 0.8657010850235238 0.8657010850235238\n"
             "\n"
             "K_POINTS automatic\n"
             " 6 6 6 0 0 0\n"
             "\n";
    log() << "-------\n\n\n";
  }
  if (vm.count("properties")) {
    log() << "\n### properties.calc.json ##################\n\n";

    log() << "properties.calc.json:                                       \n"
             "  Results of calculations for a particular configuration    \n"
             "  should be stored in the directory                         \n"
             "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE, \n"
             "  and the final calculated structure and properties         \n"
             "  summarized in the file                                    \n"
             "    "
             "$ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE/"
             "properties.calc.json \n"
             "  The 'properties.calc.json' file is read by CASM to extract the "
             "  first-principles calculted properties of interest. If the \n"
             "  'properties.calc.json' file does not exist in the         \n"
             "    $ROOT/training_data/$SCELNAME/$CONFIGID/$CURR_CALCTYPE  \n"
             "  directory CASM assumes that no data is available for that \n"
             "  configuration.                                            \n\n"

             "  The 'properties.calc.json' uses CASM standard units eV and\n"
             "  Angstroms.\n\n";

    log() << "EXAMPLE: (to be updated)\n\n";
  }
  if (vm.count("comp")) {
    log() << "\n### composition_axes.json ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n";
    log() << "$ROOT/.casm/composition_axes.json\n";
    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log()
        << "This JSON file contains the currently selected composition axes, "
           "and \n"
           "a list of possible standard or custom composition axes.      \n\n"

           "possible_axes:                                              \n"
           "  A JSON object containing each possible set of composition axes "
           "  \n"
           "  as an attribute with a unique string as the key.               "
           "  \n\n"

           "current_axes:\n"
           "  A string denoting the key of currently selected composition "
           "axes\n"
           "  from among the list of \"possible_axes\".\n\n"

           "enumerated:\n"
           "  A list of strings denoting the keys of composition axes sets "
           "that\n"
           "  from among the list of \"possible_axes\" that were "
           "automaticatlly\n"
           "  enumerated by `casm composition --calc`.\n\n"

           "possible_axes:components                                    \n"
           "  A JSON array containing the names of possible species.    \n\n"

           "possible_axes:independent_compositions                  \n"
           "  The number of independent composition axes.               \n\n"

           "possible_axes:origin                                    \n"
           "  The composition of origin the of composition axes in terms of  "
           "   \n"
           "  number of each component species per primitive cell, ordered "
           "as in\n"
           "  the 'components' array.                                   \n\n"

           "possible_axes:a, b, c, ...                              \n"
           "  The composition of end members a, b, c, etc. in terms of "
           "number of\n"
           "  each component species per primitive cell, ordered as in the   "
           "   \n"
           "  'components' array.                                       \n\n"

           "possible_axes:param_formula:                            \n"
           "  The formula that converts 'comp_n' (# of each component per    "
           "   \n"
           "  primitive cell) to 'comp' (composition relative the selected   "
           "   \n"
           "  composition axes).                                        \n\n"

           "possible_axes:mol_formula:                              \n"
           "  The formula that converts 'comp' (composition relative the\n"
           "  selected composition axes) to 'comp_n' (# of each component "
           "per   \n"
           "  primitive cell).                                          \n\n\n";

    log() << "EXAMPLE:\n";
    log() << "-------\n";
    log() << "{\n"
             "  \"current_axes\" : \"0\",\n"
             "  \"enumerated\" : [\"0\", \"1\"],\n"
             "  \"possible_axes\" : {\n"
             "    \"0\" : {\n"
             "      \"a\" : [\n"
             "        [ 2.000000000000 ],\n"
             "        [ 0.000000000000 ],\n"
             "        [ 2.000000000000 ]\n"
             "      ],\n"
             "      \"components\" : [ \"Zr\", \"Va\", \"O\" ],\n"
             "      \"independent_compositions\" : 1,\n"
             "      \"mol_formula\" : \"Zr(2)Va(2-2a)O(2a)\",\n"
             "      \"origin\" : [\n"
             "        [ 2.000000000000 ],\n"
             "        [ 2.000000000000 ],\n"
             "        [ 0.000000000000 ]\n"
             "      ],\n"
             "      \"param_formula\" : \"a(0.5-0.25Va+0.25O)\"\n"
             "    },\n"
             "    \"1\" : {\n"
             "      \"a\" : [\n"
             "        [ 2.000000000000 ],\n"
             "        [ 2.000000000000 ],\n"
             "        [ 0.000000000000 ]\n"
             "      ],\n"
             "      \"components\" : [ \"Zr\", \"Va\", \"O\" ],\n"
             "      \"independent_compositions\" : 1,\n"
             "      \"mol_formula\" : \"Zr(2)Va(2a)O(2-2a)\",\n"
             "      \"origin\" : [\n"
             "        [ 2.000000000000 ],\n"
             "        [ 0.000000000000 ],\n"
             "        [ 2.000000000000 ]\n"
             "      ],\n"
             "      \"param_formula\" : \"a(0.5+0.25Va-0.25O)\"\n"
             "    }\n"
             "  }\n"
             "}\n";

    log() << "-------\n";
  }

  if (vm.count("bspecs")) {
    log() << "\n### bspecs.json ##################\n\n";

    log() << "LOCATION:\n";
    log() << "$ROOT/basis_sets/$CURR_BSET/bspecs.json\n";
    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log()
        << "This JSON file contains specifications for generating the     \n"
           "cluster basis functions. For a complete description of the    \n"
           "bspecs.json file, see the online documentation at:            \n\n"

           "  "
           "https://prisms-center.github.io/CASMcode_docs/pages/formats/casm/"
           "clex/ClexBasisSpecs.html \n\n";

    log() << "EXAMPLE: occupation DoF, periodic cluster expansion \n";
    log() << "-------\n";
    log() << "{\n"
             "  \"basis_function_specs\" : {\n"
             "    \"dof_specs\": {\n"
             "      \"occ\": {\n"
             "        \"site_basis_functions\" : \"occupation\"\n"
             "      }\n"
             "    }\n"
             "  },\n"
             "  \"cluster_specs\": {\n"
             "    \"method\": \"periodic_max_length\",\n"
             "    \"params\": {\n"
             "      \"orbit_branch_specs\" : {\n"
             "        \"2\" : {\"max_length\" : 4.01},\n"
             "        \"3\" : {\"max_length\" : 3.01}\n"
             "      },\n"
             "      \"orbit_specs\" : [\n"
             "        {\n"
             "          \"coordinate_mode\" : \"Direct\",\n"
             "          \"prototype\" : [\n"
             "            [ 0.0, 0.0, 0.0 ],\n"
             "            [ 1.0, 0.0, 0.0 ],\n"
             "            [ 2.0, 0.0, 0.0 ],\n"
             "            [ 3.0, 0.0, 0.0 ]\n"
             "          ],\n"
             "          \"include_subclusters\" : true  \n"
             "        },\n"
             "        {\n"
             "          \"coordinate_mode\" : \"Integral\",\n"
             "          \"prototype\" : [\n"
             "            [ 0, 0, 0, 0 ],\n"
             "            [ 1, 0, 0, 0 ],\n"
             "            [ 0, 0, 0, 3 ]\n"
             "          ],\n"
             "          \"include_subclusters\" : true\n"
             "        }\n"
             "      ]\n"
             "    }\n"
             "  }\n"
             "}\n";
    log() << "-------\n";
  }

  if (vm.count("clust")) {
    log() << "\n### clust.json ##################\n\n";

    log() << "LOCATION:\n";
    log() << "$ROOT/basis_sets/$CURR_BSET/clust.json\n";
    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log() << "This JSON file contains a description of all the cluster    \n"
             "orbits generated using the \"cluster_specs\" parameters in  \n"
             "the 'bspecs.json' specifications.  \n\n";
  }

  if (vm.count("basis")) {
    log() << "\n### basis.json ##################\n\n";

    log() << "LOCATION:\n";
    log() << "$ROOT/basis_sets/$CURR_BSET/basis.json\n";
    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log() << "This JSON file contains a description of the basis functions\n\n"
             "generated using the 'bspecs.json' specifications. \n\n";
  }

  if (vm.count("clex")) {
    log() << "\n### $TITLE_Clexulator.* ##################\n\n";

    log() << "LOCATION:\n";
    log() << "$ROOT/basis_sets/$CURR_BSET/$TITLE_Clexulator.*\n";
    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log() << "$TITLE_Clexulator.cc contains C++ code generated by CASM for     "
             "  \n"
             "the cluster basis functions. It is automatically compiled into   "
             "  \n"
             "$TITLE_Clexulator.o and $TITLE_Clexulator.so for use by CASM.    "
             "  \n\n"
          << std::endl;
  }

  if (vm.count("ref")) {
    log() << "\n### ref ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n\n";
    log() << "$ROOT/training_data/settings/$CURR_CALCTYPE/$CURR_REF/"
             "chemical_reference.json\n";
    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log() << "    The chemical reference determines the value of the formation "
             "energy  \n"
             "    and chemical potentials calculated by CASM.                  "
             "\n\n"

             "    Chemical references states are set by specifying a "
             "hyperplane in     \n"
             "    energy/atom - composition (as atom_frac) space. This may be "
             "done by  \n"
             "    specifying the hyperplane explicitly, or by specifying "
             "several       \n"
             "    reference states with energy/atom and composition (as "
             "atom_frac) for \n"
             "    enough states to span the composition space of the allowed "
             "occupants \n"
             "    specified in the prim. For consistency with other CASM "
             "projects,     \n"
             "    additional reference states extending to other compositional "
             "        \n"
             "    dimensions may be included also. The pure Va reference is "
             "always 0.  \n\n";

    log()
        << "    The reference states are stored in the "
           "'chemical_reference.json' file\n"
           "    in one of two formats:                                       "
           "\n\n"

           "    1) Reference state composition and energy_per_species.       \n"
           "       In this format each reference state is represented by a "
           "JSON      \n"
           "       object storing the number of each species present in the "
           "reference\n"
           "       state and the energy_per_species for that reference "
           "state. Species\n"
           "       that are not in the primitive structure may also be "
           "included in   \n"
           "       the reference states as long as the composition space of "
           "the      \n"
           "       primitive structure is spanned by the hyperplane "
           "connecting the   \n"
           "       provided reference states.                                \n"
           R"(       '[)"
        << "\n"
           R"(          {"A": 3.4, "C": 2.0, "energy_per_species": 2.0},)"
        << "\n"
           R"(          {"B": 2.0, "energy_per_species": 4.0}, )"
        << "\n"
           R"(          {"C": 1.0, "energy_per_species": 3.0}  )"
        << "\n"
           R"(        ]')"
        << "\n\n"

           "    2) Input an array of energy_per_species, for each species in "
           "prim,   \n"
           "       including 0.0 for vacancy:                                \n"
           "        '[X, X, X]'                                              "
           "\n\n";

    log() << "    When using '--set' it is also possible to specialize the "
             "chemical    \n"
             "    reference at the supercell or configuration level by adding "
             "the      \n"
             "    --scelname or --configname option.                           "
             "\n\n";

    log() << "EXAMPLE: chemical_reference.json\n";
    log() << "-------\n";
    log() << "{\n"
             "  \"chemical_reference\" : {\n"
             "    \"config\" : {\n"
             "      \"SCEL4_2_2_1_1_1_0 / 0\" : [\n"
             "        {\n"
             "          \"A\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -1.500000000000\n"
             "        },\n"
             "        {\n"
             "          \"B\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -2.000100000000\n"
             "        },\n"
             "        {\n"
             "          \"C\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -8.030000000000\n"
             "        }\n"
             "      ],\n"
             "      \"SCEL4_2_2_1_1_1_0 / 2\" : [ -1.520000000000, "
             "-2.000100000000, -8.030000000000 ]\n"
             "    },\n"
             "    \"global\" : [\n"
             "      {\n"
             "        \"A\" : 0.500000000000,\n"
             "        \"B\" : 0.500000000000,\n"
             "        \"energy_per_species\" : -1.500000000000\n"
             "      },\n"
             "      {\n"
             "        \"B\" : 1.000000000000,\n"
             "        \"energy_per_species\" : -2.000000000000\n"
             "      },\n"
             "      {\n"
             "        \"C\" : 1.000000000000,\n"
             "        \"energy_per_species\" : -8.000000000000\n"
             "      },\n"
             "      {\n"
             "        \"D\" : 1.000000000000,\n"
             "        \"energy_per_species\" : -4.000000000000\n"
             "      }\n"
             "    ],\n"
             "    \"species_order\" : [ \"A\", \"B\", \"C\" ],\n"
             "    \"supercell\" : {\n"
             "      \"SCEL3_1_3_1_1_0_0\" : [\n"
             "        {\n"
             "          \"A\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -1.500000000000\n"
             "        },\n"
             "        {\n"
             "          \"B\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -2.000000000000\n"
             "        },\n"
             "        {\n"
             "          \"C\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -8.001000000000\n"
             "        }\n"
             "      ],\n"
             "      \"SCEL4_2_2_1_1_1_0\" : [\n"
             "        {\n"
             "          \"A\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -1.500000000000\n"
             "        },\n"
             "        {\n"
             "          \"B\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -2.000000000000\n"
             "        },\n"
             "        {\n"
             "          \"C\" : 1.000000000000,\n"
             "          \"energy_per_species\" : -8.030000000000\n"
             "        }\n"
             "      ]\n"
             "    }\n"
             "  }\n"
             "}\n";
    log() << "-------\n";
  }

  if (vm.count("lat")) {
    log() << "\n### lat ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n";
    log() << "$ROOT/training_data/$SCELNAME/LAT\n\n\n";

    log() << "DESCRIPTION:\n";
    log() << "This file is generated using the '--write-pos' option for   \n"
             "'casm query' with '--type scel'. It contains the lattice    \n"
             "vectors of a particular supercell of as a column vector     \n"
             "matrix.  \n\n";

    log() << "EXAMPLE:\n";
    log() << "-------\n";
    log() << "      10.52850134      0.00000000      0.00000000\n"
             "      0.00000000     10.52850134      0.00000000\n"
             "      0.00000000      0.00000000     10.52850134\n";
    log() << "-------\n";
  }

  if (vm.count("pos")) {
    log() << "\n### pos ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n";
    log() << "$ROOT/training_data/$SCELNAME/$CONFIGID/POS\n\n\n";

    log() << "DESCRIPTION:\n";
    log() << "This file is generated using the '--write-pos' option for   \n"
             "'casm query' with '--type config'. It is VASP POSCAR type   \n"
             "file constructed from the DoF values of the particular      \n"
             "configuration. The format is standard vasp 5.x format.      \n\n";

    log() << "EXAMPLE:\n";
    log() << "-------\n";
    log() << "SCEL3_1_1_3_0_0_0\n"
             "1.00000000\n"
             "      0.00000000      1.75475022      1.75475022\n"
             "      1.75475022      0.00000000      1.75475022\n"
             "      3.50950045      3.50950045     -3.50950045\n"
             "Al Ni\n"
             "2 1\n"
             "Direct\n"
             "   0.3333333   0.3333333   0.3333333\n"
             "   0.6666667   0.6666667   0.6666667\n"
             "   0.0000000   0.0000000   0.0000000\n"
             "\n";
    log() << "-------\n";
  }

  if (vm.count("eci")) {
    log() << "\n### eci.json ##################\n\n";

    log() << "LOCATION:\n";
    log() << "$ROOT/cluster_expansions/clex.formation_energy/$CURR_BSET/"
             "$CURR_CALCTYPE/$CURR_REF/$CURR_ECI/eci.json\n";
    log() << "\n\n";

    log() << "DESCRIPTION:\n";
    log() << "This is a copy of the 'basis.json' file with the addition of\n"
             "\"eci: <value>\" for each basis function with non-zero      \n"
             "coefficient value. By convention, an additional attribute   \n"
             "\"fit\" may be added at the top level with a description of \n"
             "the fitting method used to fit the coefficients for         \n"
             "documentation purposes.                                     \n\n";
  }

  if (vm.count("monte")) {
    log() << "\n### monte ##################\n\n";

    log() << "LOCATION WHEN GENERATED:\n";
    log() << "  User determined\n\n\n";

    log() << "DESCRIPTION:\n";
    log()
        << "  The Monte Carlo input file does not need to be in any particular "
           "\n"
           "  location, as long as it is somewhere inside the CASM project     "
           "\n"
           "  directory or subdirectories. The input file contains a JSON      "
           "\n"
           "  object with \"ensemble\", \"method\", \"model\", \"supercell\",  "
           "\n"
           "  \"data\", and \"driver\" attributes, as described below. An      "
           "\n"
           "  optional attribute \"debug\" may also be included to print       "
           "\n"
           "  information that may be useful for debugging an input file.      "
           "\n\n"

           "Input file parameters:                                             "
           "\n\n"

           "\"ensemble\" (string):                                             "
           "\n\n"

           "  Possible options for \"ensemble\" are:                           "
           "\n\n"

           "    \"GrandCanonical\" or \"grand_canonical\": Semi-grand "
           "canonical\n"
           "    Monte Carlo calculation in which the total number of sites is  "
           "\n"
           "    fixed, but the occupants on each site may vary. One occupant   "
           "\n"
           "    change at a time is attempted.                                 "
           "\n\n"

           "    \"Canonical\" or \"canonical\": Canonical Monte Carlo \n"
           "    calculation in which the total number of each type of occupant "
           "\n"
           "    is fixed. Each Monte Carlo step attempts to swap a pair of     "
           "\n"
           "    occupants.                                                     "
           "\n\n\n"

           "\"method\" (string):                                               "
           "\n\n"

           "  Possible options for \"method\" are:                             "
           "\n\n"

           "    \"Metropolis\" or \"metropolis\": Run Monte Carlo calculations "
           "\n"
           "    using the Metropolis algorithm.                                "
           "\n\n"

           "    \"LTE1\" or \"lte1\": Single spin flip low temperature         "
           "\n"
           "    expansion calculations.                                        "
           "\n\n\n"

           "\"model\": (JSON object)                                           "
           "\n\n"

           "  /\"formation_energy\": (string, optional, "
           "default=\"formation_energy\")\n"
           "    Specifies the cluster expansion to use to calculated         \n"
           "    formation energy. Should be one of the ones listed by        \n"
           "    'casm settings -l'.                                        \n\n"

           "  /\"order_parameter\": (JSON object, optional)                  \n"
           "    A JSON object representing a DoF space defining the order    \n"
           "    parameter space and basis. This can be copied from the output\n"
           "    of:                                                          \n"
           "    - `casm sym --dof-space-analysis`: for a symmetry-adapted    \n"
           "      order parameter basis spanning an entire DoF vector        \n"
           "      subspace.                                                  \n"
           "    - `casm sym --config-space-analysis`: for a symmetry-adapted \n"
           "      order parameter basis spanning the DoF subspace spanned by \n"
           "      a choice of significant configurations. This method is     \n"
           "      faster than `--dof-space-analysis`, and often the resulting\n"
           "      axes do lie along high-symmetry directions, but they may   \n"
           "      not lie within a single irreducible subspace, and they are \n"
           "      not explicitly rotated to the highest symmetry directions."
           "\n\n"

           "    /\"subspaces\": (JSON list of list, optional)                \n"
           "      If added inside the \"order_parameter\" DoF space JSON     \n"
           "      object, this specifies order parameter subspaces for which \n"
           "      the l2-norm of the order parameter should be sampled and   \n"
           "      output as \"order_parameter_subspace(i)\". Example:        \n"
           "                                                                 \n"
           "         \"subspaces\": [                                        \n"
           "           [0],                  # a 1-d subspace                \n"
           "           [1, 2],               # a 2-d subspace                \n"
           "           [3, 4, 5]             # a 3-d subspace                \n"
           "           [6, 7, 8, 9, 10, 11]  # a 6-d subspace                \n"
           "         ]                                                       \n"
           "                                                                 \n"
           "      Only subspaces that should be sampled need to be included.   "
           "\n\n\n"

           "\"supercell\": (3x3 JSON arrays of integers)                      "
           "\n"
           "    The supercell transformation matrix.                           "
           "\n\n"

           "\"data\": (JSON object)                                            "
           "\n\n"

           "  /\"sample_by\": (string)                                         "
           "\n"
           "    Specify unit for the period between samples.  May be either    "
           "\n"
           "    \"step\" (sample after every \"sample_period\" proposed Monte  "
           "\n"
           "    Carlo events), or \"pass\" (sample after the \"sample_period\" "
           "\n"
           "    number of passes), where 1 pass is a number of steps equal to  "
           "\n"
           "    the number of sites in the supercell that have variable        "
           "\n"
           "    occupation).                                                   "
           "\n\n"

           "  /\"sample_period\": (integer)                                    "
           "\n"
           "    Specify how many steps or passes to wait between data samples. "
           "\n\n"

           "  /\"measurements\": (JSON array containing JSON objects)         "
           "\n"
           "    Specifies which properties to sample. Each JSON object should  "
           "\n"
           "    include \"quantity\" (string) specifying a property to be      "
           "\n"
           "    sampled. Optionally, it may also include \"precision\" "
           "(number),\n"
           "    indicating the required (absolute) precision in the average of "
           "\n"
           "    the quantity for the calculation to be considered converged. "
           "If\n"
           "    a precision is given for any quantity, then the Monte Carlo    "
           "\n"
           "    calculations run in automatic convergence mode and continue    "
           "\n"
           "    until all quantities with a specified precision are converged  "
           "\n"
           "    to level requested.                                            "
           "\n\n"

           "    Possible options for \"quantity\" are:                         "
           "\n"
           "      \"comp\": composition, relative the composition axes         "
           "\n"
           "      \"comp_n\": composition, number of atoms per unit cell       "
           "\n"
           "      \"site_frac\": composition, normalized per basis site        "
           "\n"
           "      \"atom_frac\": composition, normalized per total number of "
           "atoms\n"
           "      \"formation_energy\": formation energy (per unit cell)       "
           "\n"
           "      \"potential_energy\": potential energy (per unit cell),      "
           "\n"
           "        (= formation_energy - sum_i(mu_i*comp_i))                  "
           "\n"
           "      \"non_zero_eci_correlations\": correlations (per unit cell)  "
           "\n"
           "        which have non-zero eci values.                            "
           "\n"
           "      \"<anything else>\": is interpreted as a 'casm query' query  "
           "\n\n"

           "  /\"confidence\": (number, range (0.0, 1.0), default 0.95)        "
           "\n"
           "    The confidence level used for calculating the precision in the "
           "\n"
           "    average value of sampled quantities.                           "
           "\n\n"

           "  /\"min_pass\", /\"min_step\", /\"min_sample\": (integer)         "
           "\n"
           "    If in automatic convergence mode, prevents the calculation "
           "from\n"
           "    a minimum number of passes, steps, or samples have occurred.   "
           "\n\n"

           "  /\"max_pass\", /\"max_step\", /\"max_sample\": (integer)         "
           "\n"
           "    If in automatic convergence mode, stops the calculation if the "
           "\n"
           "    specified number of passes, steps, or samples have occurred.   "
           "\n\n"

           "  /\"N_pass\", /\"N_step\", /\"N_sample\": (integer)               "
           "\n"
           "    When not in automatic convergence mode (no precision has been  "
           "\n"
           "    specified for any quantities being sampled), stops the         "
           "\n"
           "    calculation when the specified number of passes, steps, or     "
           "\n"
           "    samples have occurred.                                         "
           "\n\n"

           "  /\"equilibration_passes_first_run\": (integer)                   "
           "\n"
           "    If included, the requested number of passes will be performed  "
           "\n"
           "    at the initial conditions as a preliminary step before the     "
           "\n"
           "    actual run begins. This may be useful when not running in      "
           "\n"
           "    automatic convergence mode.                                    "
           "\n\n"

           "  /\"equilibration_passes_each_run\": (integer)                    "
           "\n"
           "    If included, the requested number of passes will be performed  "
           "\n"
           "    at each condition as a preliminary step before the actual run  "
           "\n"
           "    begins. This may be useful when not running in automatic       "
           "\n"
           "    convergence mode.                                              "
           "\n\n"

           "  /\"storage\": (JSON object) Options for writing results.         "
           "\n\n"

           "    /\"output_format\": (string or JSON array of string)           "
           "\n"
           "      Specifies the type or types of output files. Current options "
           "\n"
           "      are \"csv\" or \"json\". Type names with either all lower    "
           "\n"
           "      case or all   upper case are accepted.                       "
           "\n\n"

           "    /\"save_state_details\": (bool, optional, default=false)     \n"
           "      If true, the initial and final state are written to a file \n"
           "      in the conditions directory, for each condition.         \n\n"

           "    /\"write_observations\": (boolean, default false)              "
           "\n"
           "      If true, all individual observations of the quantities       "
           "\n"
           "      requested to be sampled will be written to compressed files: "
           "\n"
           "        \"output_directory\"/conditions.i/observations.ext.gz      "
           "\n"
           "      where 'i' is the condition index and 'ext' is the output     "
           "\n"
           "      format.                                                      "
           "\n\n"

           "    /\"write_trajectory\": (boolean, default false)                "
           "\n"
           "      If true, the value of all degrees of freedom at the time of  "
           "\n"
           "      each sample will be written to compressed files:             "
           "\n"
           "        \"output_directory\"/conditions.i/trajectory.ext.gz        "
           "\n"
           "      where 'i' is the condition index and 'ext' is the output     "
           "\n"
           "      format.                                                      "
           "\n\n"

           "  /\"enumeration\": (JSON object, optional)                      \n"
           "    If included, save configurations encountered during Monte    \n"
           "    Carlo calculations by keeping a 'hall of fame' of best       \n"
           "    scoring configurations. After the calculation at a particular\n"
           "    set of thermodynamic conditions completes, the configurations\n"
           "    in the hall of fame are saved to the project configuration   \n"
           "    list.                                                      \n\n"

           "    /\"check\": (string, default=\"eq(1,1)\")                    \n"
           "      A 'casm query'-like string that returns a boolean value    \n"
           "      indicating if (true) a configuration should be considered  \n"
           "      for the enumeration hall of fame. The default always       \n"
           "      returns true.                                            \n\n"

           "    /\"metric\": (string, default=\"clex_hull_dist(ALL)\")       \n"
           "      A 'casm query'-like string that provides a metric for      \n"
           "      ranking configurations as they are encountered during a    \n"
           "      Monte Carlo calculation. The resulting value is used to    \n"
           "      create a hall of fame of 'best' configurations encountered \n"
           "      during the calculation. When the calculation is complete   \n"
           "      configurations in the hall of fame are added to the CASM   \n"
           "      project config list. The 'casm query'-like command should  \n"
           "      evaluate to a number.                                    \n\n"

           "      Besides the properties listed via 'casm query -h "
           "properties',\n"
           "      and 'casm query -h operators', both \"check\" and \"metric\" "
           "\n"
           "      can also use the property \"potential_energy\".          \n\n"

           "    /\"tolerance\": (number, optional, default=1e-8)             \n"
           "      Tolerance used for floating point comparison of            \n"
           "      configuration scores in the enumeration hall of fame.    \n\n"

           "    /\"sample_mode\": (string, optional, default=\"on_sample\")  \n"
           "      Indicate when to attempt to insert configurations into the \n"
           "      enumeration hall of fame. Options are:                     \n"
           "        \"on_accept\": after each accepted Monte Carlo event     \n"
           "        \"on_sample\": after each data sample                  \n\n"

           "    /\"check_existence\": (bool, optional, default=true)         \n"
           "      If true, only configurations that do not already exist in  \n"
           "      the config list are inserted into the enumeration hall of  \n"
           "      fame. If this is true, one of 'insert_canonical' or        \n"
           "      'insert_primitive_only' must also be true.               \n\n"

           "    /\"insert_canonical\": (bool, optional, default=true)        \n"
           "      If true, configurations are put into canonical form before \n"
           "      being inserted into the enumeration hall of fame.        \n\n"

           "    /\"insert_primitive_only\": (bool, optional, default=true)   \n"
           "      If true, configurations are put into primitive and         \n"
           "      canonical form before being inserted into the enumeration  \n"
           "      hall of fame.                                            \n\n"

           "    /\"save_primitive_only\": (bool, optional, default=true)     \n"
           "      If true, only the primitive form of configurations in the  \n"
           "      hall of fame are saved into the project database.        \n\n"

           "    /\"N_halloffame\": (integer, optional, default=100)          \n"
           "      The number of configurations that are allowed in the       \n"
           "      enumeration hall of fame.                                \n\n"

           "    /\"save_configs\": bool (optional, default=false) \n"
           "      If true, save hall of fame configurations in the project   \n"
           "      database after the run completes.\n\n"

           "    /\"save_configs_periodically\": bool (optional, "
           "default=false) \n"
           "      When `save_configs==true` and this is also true, save hall \n"
           "      of fame configurations in the project database periodically\n"
           "      during the run.\n\n"

           "    /\"save_configs_period\": bool (optional, default=10000)     \n"
           "      When `save_configs_periodically==true`, this is how often  \n"
           "      the configurations in the hall of fame are saved in the    \n"
           "      project database, in terms of the number of Monte Carlo    \n"
           "      samples taken.                                           \n\n"

           "    /\"output_configurations\": bool (optional, default=false)   \n"
           "      If true, output hall of fame as a data file in the         \n"
           "      conditions directory. Requires "
           "`data/storage/save_state_details==true`.\n"
           "      The \"insert\"/\"enum_output_file\" option of the method   \n"
           "      `casm info -m ConfigurationInfo` can be used to insert the \n"
           "      configurations written to this data file into the project  \n"
           "      database.                                                \n\n"

           "    /\"output_period\": int (optional, default=10000) \n"
           "      When `output_configurations==true`, this is how often the  \n"
           "      configurations in the hall of fame are output to a file, in\n"
           "      terms of the number of Monte Carlo samples taken.        \n\n"

           "    /\"output_configurations_options\": object (optional) \n"
           "      Output options for when `output_configurations==true`. \n\n"

           "      /\"path\": string (optional, default=\"enum.out\") \n"
           "        Output file name.\n"
           "      /\"json\": bool (optional, default=false) \n"
           "        If true, write JSON output files. Else CSV style.\n"
           "      /\"json_arrays\": bool (optional, default=false) \n"
           "        If true, write data in JSON arrays. \n"
           "      /\"compress\": bool (optional, default=false) \n"
           "        If true, compress data using gz. If `path` does not end \n"
           "        in '.gz' it will be appended. \n"
           "      /\"properties\": array of string (optional, default=[])    \n"
           "        Specify additional configuration properties to output,\n"
           "        using properties available from `casm query`.         "
           "\n\n\n"

           "\"driver\": (JSON object)                                          "
           "\n\n"

           "  /\"motif\": (JSON object)                                      \n"
           "      Specifies the initial occupation of the supercell.       \n\n"

           "      For canonical ensemble Monte Carlo calculations an         \n"
           "      additional step changes the occupants on random sites to   \n"
           "      make the actual composition as close as possible to the    \n"
           "      requested composition.                                   \n\n"

           "    /\"configname\": (string, optional)                          \n"
           "      A configuration name, \"auto\", \"restricted_auto\", or    \n"
           "      \"default\".                                             \n\n"

           "      Specifies the configuration that is tiled to fill the      \n"
           "      supercell. If necessary, symmetry operations may be applied\n"
           "      An error is thrown if the specified configuration can not  \n"
           "      be used to fill the \"supercell\".                       \n\n"

           "      Possible options for \"configname\" are:                   \n"
           "        A configuration name (ex. \"SCEL3_3_1_1_0_2_2/0\")       \n"
           "        \"auto\": (\"grand_canonical\" ensemble only) Enumerated \n"
           "        configurations will be searched for the configuration    \n"
           "        with the lowest potential energy to use as the motif.    \n"
           "        \"default\": If the value \"default\" is used, the       \n"
           "          initial motif occupation is determined from the        \n"
           "          occupation order in the PRIM. \n"
           "        \"restricted_auto\": (\"grand_canonical\" ensemble only) \n"
           "          Same as \"auto\", but only configurations that can tile\n"
           "          the supercell are considered. As a last resort,        \n"
           "          \"default\" is used.                                 \n\n"

           "    /\"configdof\": (string, optional)                           \n"
           "      Specifies the path to a configdof JSON file, such as       \n"
           "      \"initial_state.json\" or \"final_state.json\", containing \n"
           "      the degrees of freedom to initialize the supercell with  \n\n"

           "  /\"mode\": (string)                                            \n"
           "    Specify the drive mode.                                    \n\n"

           "    Possible options for \"mode\" are:                           \n"
           "      \"incremental\": perform one or more calculations, starting\n"
           "        at the initial conditions and incrementing by the        \n"
           "        incremental conditions up to (and including) the final   \n"
           "        conditions.                                            \n\n"
           "      \"custom\": perform one or more calculations, as specified \n"
           "        by the \"custom_conditions\".                          \n\n"

           "  /\"dependent_runs\": (boolean, default true)                   \n"
           "    If true, begin the next calculation with the final DoF from  \n"
           "    the previous calculation. If false, begin each calculation   \n"
           "    with the DoF specified for the \"motif\".                  \n\n"

           "  /\"initial_conditions\",\n"
           "  /\"incremental_conditions\", \n"
           "  /\"final_conditions\": (JSON object, optional)                 \n"
           "    Specifies the applied conditions for the calculation. For    \n"
           "    \"incremental_conditions\", specifies the change in          \n"
           "    conditions between individual calculations. Each JSON object \n"
           "    should include: \n\n"

           "    /\"temperature\": (number)                                   \n"
           "      The temperature in K.                                    \n\n"

           "    /\"param_chem_pot\" (JSON object, \"grand_canonical\" ensemble "
           "only)\n"
           "      The parametric chemical potential(s)                     \n\n"

           "      /\"a\", /\"b\", ...: (number)                              \n"
           "      The parametric chemical potentials conjugate to the        \n"
           "      parametric compositions. The number of parametric chemical \n"
           "      potentials provided must match the number of independent   \n"
           "      compositions.                                            \n\n"

           "    /\"comp\" (JSON object, \"canonical\" ensemble only, option "
           "1)\n"
           "      The parametric composition(s)                            \n\n"

           "      /\"a\", /\"b\", ...: (number)                              \n"
           "      The parametric composition. The number of entries provided \n"
           "      must match the number of independent compositions.       \n\n"

           "    /\"comp\" (JSON array, \"canonical\" ensemble only, option 2)\n"
           "      The parametric composition(s)                          \n\n"

           "      [number, ...]                                              \n"
           "      An array containing the parametric composition. The number \n"
           "      of entries provided must match the number of independent   \n"
           "      compositions.                                            \n\n"

           "    /\"comp_n\" (JSON object, \"canonical\" ensemble only, option "
           "3)\n"
           "      The mol composition per unitcell                         \n\n"

           "      /\"A\", /\"B\", ...: (number)                              \n"
           "      The mol composition per unitcell. The entries provided must\n"
           "      match occupant names in the 'prim.json' file. The values   \n"
           "      are summed, normalized, and then converted to parametric   \n"
           "      composition.\n\n"

           "    /\"comp_n\" (JSON array, \"canonical\" ensemble only, option "
           "4)\n"
           "      The mol composition per unitcell                         \n\n"

           "      [number, ...]                                              \n"
           "      An array containing the mol composition per unitcell. The  \n"
           "      number of entries provided must match the number           \n"
           "      components. The values are summed, normalized, and         \n"
           "      converted to parametric composition.  \n\n"

           "    /\"include_formation_energy\": (bool, optional)              \n"
           "      If true, include the formation energy in the potential     \n"
           "      energy. Default is true, except if corr_matching_pot or    \n"
           "      random_alloy_corr_matching_pot are included, in which case \n"
           "      the default value is false.                              \n\n"

           "    /\"param_comp_quad_pot_target\": (1d array, optional)        \n"
           "    /\"param_comp_quad_pot_vector\": (1d array, optional)        \n"
           "    /\"param_comp_quad_pot_matrix\": (2d array, optional)        \n"
           "      If \"_target\" is provided, include a variance-constrained \n"
           "      potential component with regards to the parameteric        \n"
           "      composition. If \"_vector\" is provided, the variance-     \n"
           "      constrained potential component is:                        \n"
           "                                                                 \n"
           "          Epot += sum_i v_i (x_i - x^0_i)^2                      \n"
           "                                                                 \n"
           "      where:                                                     \n"
           "      - x_i: parametric composition vector                       \n"
           "      - x^0_i: parametric composition component target           \n"
           "      - v_i: variance-constrained curvature (from \"_vector\")   \n"
           "                                                                 \n"
           "      If \"_matrix\" is provided, the variance-constrained       \n"
           "      potential component is:                                    \n"
           "                                                                 \n"
           "          Epot += (x - x^0)^T * M (x - x^0)                      \n"
           "                                                                 \n"
           "      where:                                                     \n"
           "      - x: parametric composition vector                         \n"
           "      - x^0: parametric composition vector target                \n"
           "      - M: variance-constrained curvature (from \"_matrix\")   \n\n"
           "                                                                 \n"
           "      Note that all inputs must be consistent in dimension with  \n"
           "      the parameter composition.                               \n\n"
           "                                                                 \n"
           "      All values may be set to vary in \"incremental_conditions\","
           "\n"
           "      or they should be set to 0.0 to remain fixed.            \n\n"

           "    /\"order_parameter_quad_pot_target\": (1d array, optional)   \n"
           "    /\"order_parameter_quad_pot_vector\": (1d array, optional)   \n"
           "    /\"order_parameter_quad_pot_matrix\": (2d array, optional)   \n"
           "      If \"_target\" is provided, include a variance-constrained \n"
           "      potential component with regards to the order parameters.  \n"
           "      If \"_vector\" is provided, the variance- constrained      \n"
           "      potential component is:                                    \n"
           "                                                                 \n"
           "          Epot += sum_i v_i (x_i - x^0_i)^2                      \n"
           "                                                                 \n"
           "      where:                                                     \n"
           "      - x_i: order parameter component                           \n"
           "      - x^0_i: order parameter component target                  \n"
           "      - v_i: variance-constrained curvature (from \"_vector\")   \n"
           "                                                                 \n"
           "      If \"_matrix\" is provided, the variance-constrained       \n"
           "      potential component is:                                    \n"
           "                                                                 \n"
           "          Epot += (x - x^0)^T * M (x - x^0)                      \n"
           "                                                                 \n"
           "      where:                                                     \n"
           "      - x: order parameter component                             \n"
           "      - x^0: order parameter component target                    \n"
           "      - M: variance-constrained curvature (from \"_matrix\")   \n\n"
           "                                                                 \n"
           "      Note that all inputs must be consistent in dimension with  \n"
           "      the order parameter basis.                                 \n"
           "                                                                 \n"
           "      All values may be set to vary in \"incremental_conditions\","
           "\n"
           "      or they should be set to 0.0 to remain fixed.            \n\n"

           "    /\"corr_matching_pot\": (JSON object, optional)              \n"
           "      If provided, include a correlation-matching potential      \n"
           "      component defined as:                                      \n"
           "                                                                 \n"
           "          Epot = -w_{exact}*N_{exact} +                          \n"
           "              \\sum_i v_i * | \\Gamma_{j_i} - "
           "\\Gamma^{target}_{j_i}) |, \n"
           "                                                                 \n"
           "      where:                                                     \n"
           "      - N_{exact} is that maximum value such that                \n"
           "        | \\Gamma_{j_i} - \\Gamma^{target}_{j_i}) | < tol for all "
           "i < N_{exact}\n"
           "      - exact_matching_weight = w_{exact},                       \n"
           "      - targets[i].index = w_i                                   \n"
           "      - targets[i].value = \\Gamma^{target}_{j_i}                \n"
           "      - targets[i].weight = v_i                                  \n"
           "                                                                 \n"
           "      This potential component follows the form described in:    \n"
           "      - van de Walle, et al. CALPHAD 42 (2013) 13–18             \n"
           "                                                                 \n"
           "      The expected input JSON format is:                         \n"
           "                                                                 \n"
           "         \"corr_matching_pot\": {                                \n"
           "           \"exact_matching_weight\": float  # default=0.0       \n"
           "           \"tol\": 1e-5                     # default=1e-5      \n"
           "           \"targets\": [ {                                      \n"
           "               \"index\": int,    # linear function index        \n"
           "               \"value\": number, # target correlation value     \n"
           "               \"weight\": float  # weight of diff in sum        \n"
           "             },                                                  \n"
           "             ...                                                 \n"
           "           ]                                                     \n"
           "         }                                                       \n"
           "                                                                 \n"
           "      The exact_matching_weight and targets value and weight may \n"
           "      be set to vary in \"incremental_conditions\", or they      \n"
           "      should be set to 0.0 to remain fixed.                      \n"
           "                                                                 \n"
           "      IMPORTANT: Must be used with an `eci.json` file where all  \n"
           "      basis functions have an a value for \"eci\", even if set to\n"
           "      0.0. Otherwise, changes in the potential will not be       \n"
           "      calculated correctly.                                    \n\n"

           "    /\"random_alloy_corr_matching_pot\": (JSON object, optional) \n"
           "      If provided, include a random alloy correlation-matching   \n"
           "      potential component of the same form as "
           "\"corr_matching_pot\",\n"
           "      but with target correlation values determined from         \n"
           "      sublattice composition probabilities in a random alloy. The\n"
           "      target weights are all set to 1.0. The expected input JSON \n"
           "      format is:                                                 \n"
           "                                                                 \n"
           "         \"random_alloy_corr_matching_pot\": {                   \n"
           "           \"exact_matching_weight\": float  # default=0.0       \n"
           "           \"tol\": 1e-5                     # default=1e-5      \n"
           "           \"sublattice_prob\": [ {                              \n"
           "             [number, number, ...],  # sublattice 0              \n"
           "             [number, number, ...],  # sublattice 1              \n"
           "             ...                                                 \n"
           "           ]                                                     \n"
           "         }                                                       \n"
           "                                                                 \n"
           "      The sublattice probabilities are given for the             \n"
           "      corresponding occupant on each sublattice in the prim      \n"
           "      structure definition. Requirements:                        \n"
           "      - There must be one sublattice probability array for every \n"
           "        sublattice in the prim, even if only a single occupant is\n"
           "        allowed.                                                 \n"
           "      - The number of components must match the number of        \n"
           "        occupants allowed on the corresponding sublattice in the \n"
           "        prim definition.                                         \n"
           "      - Symmetrically equivalent sublattices must be given the   \n"
           "        same probabilities.                                      \n"
           "      - On each sublattice, the probabilities must sum to 1.0    \n"
           "                                                                 \n"
           "      The exact_matching_weight and sublattice probabilities     \n"
           "      may all be set to vary in \"incremental_conditions\", or   \n"
           "      they should be set to 0.0 to remain fixed.                 \n"
           "                                                                 \n"
           "      IMPORTANT: Must be used with an `eci.json` file where all  \n"
           "      basis functions have an a value for \"eci\", even if set to\n"
           "      0.0. Otherwise, changes in the potential will not be       \n"
           "      calculated correctly.                                    \n\n"

           "    /\"tolerance\": (number)                                       "
           " \n"
           "      Specifies a numerical tolerance for comparing conditions.    "
           " \n\n"

           "  /\"custom_conditions\":\n"
           "    (JSON array of JSON objects) An array specifying a custom     "
           "\n"
           "    path of conditions.                                           "
           "\n\n"

           "  Restarts: Metropolis Monte Carlo calculations that are stopped   "
           "\n"
           "  before the entire path has been calculated can be restarted as   "
           "\n"
           "  long as the conditions of the existing calculations agree with   "
           "\n"
           "  the conditions specified in the input settings. This means that  "
           "\n"
           "  the \"final_conditions\" might be changed to increase the length "
           "\n"
           "  of a path, or additional \"custom_conditions\" might be added,   "
           "\n"
           "  but the \"incremental_conditions\" may not be changed. Upon      "
           "\n"
           "  restart, the results summary file is checked for the last        "
           "\n"
           "  finished conditions. Then the path is resumed from the next set  "
           "\n"
           "  of conditions. It is the responsibility of the user to ensure    "
           "\n"
           "  that other important settings, such as the \"model\" are not     "
           "\n"
           "  changed inappropriately.                                         "
           "\n\n\n"

           "\"debug\" (bool, default false):                                   "
           "\n\n"

           "  If true, will print as much information as possible to assist in "
           "\n"
           "  debugging input file settings.                                   "
           "\n\n\n";

    log()
        << "EXAMPLE: Settings for an incremental Metropolis calculation     \n"
           "with increasing temperature in automatic convergence mode.\n";
    log() << "-------\n";
    log() << "{\n"
             "  \"comment\" : \"This is a sample input file. Unrecognized "
             "attributes (like the ones prepended with '_' are ignored.\",\n"
             "  \"debug\" : false,\n"
             "  \"ensemble\" : \"grand_canonical\",\n"
             "  \"method\" : \"metropolis\",\n"
             "  \"model\" : {\n"
             "    \"formation_energy\" : \"formation_energy\"\n"
             "  },\n"
             "  \"supercell\" : [\n"
             "    [10, 0, 0],\n"
             "    [0, 10, 0],\n"
             "    [0, 0, 10]\n"
             "  ],\n"
             "  \"data\" : {\n"
             "    \"sample_by\" : \"pass\",\n"
             "    \"sample_period\" : 1,\n"
             "    \"_N_sample\" : 1000, \n"
             "    \"_N_pass\" : 1000,\n"
             "    \"_N_step\" : 1000,\n"
             "    \"_max_pass\" : 10000,\n"
             "    \"min_pass\" : 1000,\n"
             "    \"_max_step\" : 10000,\n"
             "    \"_max_sample\" : 500,\n"
             "    \"_min_sample\" : 100,\n"
             "    \"confidence\" : 0.95,\n"
             "    \"measurements\" : [ \n"
             "      { \n"
             "        \"quantity\" : \"formation_energy\"\n"
             "      },\n"
             "      { \n"
             "        \"quantity\" : \"potential_energy\"\n"
             "      },\n"
             "      { \n"
             "        \"quantity\" : \"atom_frac\"\n"
             "      },\n"
             "      { \n"
             "        \"quantity\" : \"site_frac\"\n"
             "      },\n"
             "      { \n"
             "        \"quantity\" : \"comp\",\n"
             "        \"precision\" : 1e-3\n"
             "      },\n"
             "      { \n"
             "        \"quantity\" : \"comp_n\"\n"
             "      }\n"
             "    ],\n"
             "    \"storage\" : {\n"
             "      \"write_observations\" : false,\n"
             "      \"write_trajectory\" : false,\n"
             "      \"output_format\" : [\"csv\", \"json\"]\n"
             "    }\n"
             "  },\n"
             "  \"driver\" : {\n"
             "    \"mode\" : \"incremental\", \n"
             "    \"motif\" : {\n"
             "      \"configname\" : \"auto\",\n"
             "      \"_configname\" : \"SCEL3_3_1_1_0_2_2/0\",\n"
             "      \"_configdof\" : \"path/to/final_state.json\"\n"
             "    },\n"
             "    \"initial_conditions\" : {\n"
             "      \"param_chem_pot\" : {\n"
             "        \"a\" : -1.75\n"
             "      },\n"
             "      \"temperature\" : 100.0,\n"
             "      \"tolerance\" : 0.001\n"
             "    },\n"
             "    \"final_conditions\" : {\n"
             "      \"param_chem_pot\" : {\n"
             "        \"a\" : -1.75\n"
             "      },\n"
             "      \"temperature\" : 1000.0,\n"
             "      \"tolerance\" : 0.001\n"
             "    },\n"
             "    \"incremental_conditions\" : {\n"
             "      \"param_chem_pot\" : {\n"
             "        \"a\" : 0.0\n"
             "      },\n"
             "      \"temperature\" : 10.0,\n"
             "      \"tolerance\" : 0.001\n"
             "    }\n"
             "  }\n"
             "}\n";
    log() << "-------\n\n";

    log()
        << "EXAMPLE: Settings for an custom drive mode LTE1 calculation with\n"
           "increasing temperature.\n";
    log() << "-------\n";
    log() << "{\n"
             "  \"comment\" : \"This is a sample input file. Unrecognized "
             "attributes (like the ones prepended with '_' are ignored.\",\n"
             "  \"debug\" : false,\n"
             "  \"ensemble\" : \"grand_canonical\",\n"
             "  \"method\" : \"lte1\",\n"
             "  \"model\" : {\n"
             "    \"formation_energy\" : \"formation_energy\"\n"
             "  },\n"
             "  \"supercell\" : [\n"
             "    [9, 0, 0],\n"
             "    [0, 9, 0],\n"
             "    [0, 0, 9]\n"
             "  ],\n"
             "  \"data\" : {\n"
             "    \"storage\" : {\n"
             "      \"write_observations\" : false,\n"
             "      \"write_trajectory\" : false,\n"
             "      \"output_format\" : [\"csv\", \"json\"]\n"
             "    }\n"
             "  },\n"
             "  \"driver\" : {\n"
             "    \"mode\" : \"incremental\", \n"
             "    \"motif\" : {\n"
             "      \"configname\" : \"auto\",\n"
             "      \"_configname\" : \"SCEL3_3_1_1_0_2_2/0\",\n"
             "      \"_configdof\" : \"path/to/final_state.json\"\n"
             "    },\n"
             "    \"custom_conditions\" : [\n"
             "      {\n"
             "        \"param_chem_pot\" : {\n"
             "          \"a\" : 0.0\n"
             "        },\n"
             "        \"temperature\" : 100.0,\n"
             "        \"tolerance\" : 0.001\n"
             "      },\n"
             "      {\n"
             "        \"param_chem_pot\" : {\n"
             "          \"a\" : 0.0\n"
             "        },\n"
             "        \"temperature\" : 200.0,\n"
             "        \"tolerance\" : 0.001\n"
             "      },\n"
             "      {\n"
             "        \"param_chem_pot\" : {\n"
             "          \"a\" : 0.0\n"
             "        },\n"
             "        \"temperature\" : 400.0,\n"
             "        \"tolerance\" : 0.001\n"
             "      },\n"
             "      {\n"
             "        \"param_chem_pot\" : {\n"
             "          \"a\" : 0.0\n"
             "        },\n"
             "        \"temperature\" : 800.0,\n"
             "        \"tolerance\" : 0.001\n"
             "      }\n"
             "    ]\n"
             "  }\n"
             "}\n";
    log() << "-------\n";
  }

  return 0;
}

}  // namespace CASM
