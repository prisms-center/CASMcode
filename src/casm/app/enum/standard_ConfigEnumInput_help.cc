#include "casm/app/enum/standard_ConfigEnumInput_help.hh"

namespace CASM {

  // At some point the help text may go in text files, but for now...

  std::string standard_ConfigEnumInput_help() {
    return
      "  scelnames: Array of strings (optional, override with --scelnames) \n"
      "    Names of supercells used as initial state of occupation enumeration. All\n"
      "    sites not being enumerated will be set to the first listed occupant, and all\n"
      "    other DoFs will be set to zero.\n"
      "    Ex: \"scelnames\" : [\"SCEL1_1_1_1_0_0_0\",\"SCEL2_2_1_1_0_0_0\"]\n\n"

      "  supercell_selection: string (optional) \n"
      "    Name of a selection of supercells to use as initial states for enumeration.\n\n"

      "  supercells: object, ScelEnum JSON settings (optional, override with --min, --max)\n"
      "    Indicate supercells to use as initial states of enumeration in terms of size\n"
      "    and unit cell via a JSON object conforming to the format of 'ScelEnum' JSON\n"
      "    settings \"min\", \"max\", \"dirs\", and \"unit_cell\". See 'ScelEnum' \n"
      "    description for more details.\n\n"

      "  confignames: Array of strings (optional, override with --confignames) \n"
      "    Names of configurations to be used as initial states for enumeration. All \n"
      "    specified sublattices or sites will be enumerated on and all other DoFs will\n"
      "    maintain the values of the initial state.\n"
      "    Ex: \"confignames\" : [\"SCEL1_1_1_1_0_0_0/1\",\"SCEL2_2_1_1_0_0_0/3\"]\n\n"

      "  config_selection: string (optional) \n"
      "    Name of a selection of configurations to use as initial states for enumeration.\n\n"

      "  sublats: array of integers (optional, default none) \n"
      "    Restricts enumeration to specified sublattices. Each sublattice index corresponds\n"
      "    to a basis site in prim.json, indexed from 0.\n"
      "    Ex: \"sublats\" : [0,2]\n\n"

      "  sites: array of 4-entry integer arrays (optional, default none) \n"
      "    Restricts enumeration to specified sites. Sites are specified in [b,i,j,k] convention,\n"
      "    where 'b' is sublattice index and [i,j,k] specifies linear combinations of primitive-\n"
      "    cell lattice vectors.\n"
      "    Ex: \"sites\" : [[0,0,0,0],\n"
      "                   [2,0,0,0]]\n\n"

      "  cluster_specs: object (optional) \n"
      "    JSON object specifying orbits of clusters to generate. Each orbit prototype is used \n"
      "    to select sites to enumerate on each selected supercell or configuration. If there are \n"
      "    4 supercells or configurations selected, and there are 10 orbits generated, then there \n"
      "    will be 4*10=40 initial states generated. The \"cluster_specs\" option cannot be used \n"
      "    with the \"sublats\" or \"sites\" options. Expect format is: \n\n"

      "    method: string (required) \n"
      "      Specify which cluster orbit generating method will be used. Supported: \n"
      "      - \"periodic_max_length\": Clusters differing by a lattice translation are considered \n"
      "        equivalent. Cluster generation is truncated by specifying the maximum distance \n"
      "        between sites in a cluster for 2-point, 3-point, etc. clusters. The point clusters \n"
      "        comprising the asymmetric unit of the prim structure are always included. For \n"
      "        enumeration input purposes, after the cluster orbits are generated using the prim \n"
      "        factor group symmetry, the orbits are broken into sub-orbits reflecting the \n"
      "        configuration factor group symmetry of the initial enumeration state. Any orbits that \n"
      "        are duplicated under periodic boundary conditions are removed. \n\n"

      "    params: object (required) \n"
      "      Specifies parameters for the method selected by `method`. Options depend on the \n"
      "      `method` chosen: \n\n"

      "      For method==\"periodic_max_length\": \n"
      "        orbit_branch_specs: object (optional)\n"
      "          Cluster generation is truncated by specifying the maximum distance\n"
      "          between sites in a cluster for each orbit branch (i.e. 2-point, 3-point, etc.\n"
      "          clusters). The 1-point clusters comprising the asymmetric unit of the prim\n"
      "          structure are always included. \n\n"

      "            Example: \n"
      "              \"orbit_branch_specs\": {\n"
      "                \"2\": { \"max_length\": 10.0 },\n"
      "                \"3\": { \"max_length\": 8.0 },\n"
      "                       ...\n"
      "              }\n\n"

      "        orbit_specs: array (optional) \n"
      "          An array of clusters which are used to generate and include orbits of clusters \n"
      "          whether or not they meet the `max_length` truncation criteria. See the \n"
      "          cluster input format below. Use the \"include_subclusters\" option to force \n"
      "          generation of orbits for all subclusters of the specified cluster. \n"

      "            Example cluster, with \"Direct\" coordinates: \n"
      "              { \n"
      "                \"coordinate_mode\" : \"Direct\", \n"
      "                \"sites\" : [ \n"
      "                  [ 0.000000000000, 0.000000000000, 0.000000000000 ], \n"
      "                  [ 1.000000000000, 0.000000000000, 0.000000000000 ], \n"
      "                  [ 2.000000000000, 0.000000000000, 0.000000000000 ], \n"
      "                  [ 3.000000000000, 0.000000000000, 0.000000000000 ]], \n"
      "                \"include_subclusters\" : true \n"
      "              } \n\n"

      "            Example cluster, with \"Integral\" coordinates: \n"
      "              { \n"
      "                \"coordinate_mode\" : \"Integral\", \n"
      "                \"sites\" : [ \n"
      "                  [ 0, 0, 0, 0 ], \n"
      "                  [ 0, 1, 0, 0 ], \n"
      "                  [ 1, 0, 0, 0 ]], \n"
      "                \"include_subclusters\" : true \n"
      "              } \n\n"

      "  filter: string (optional, default=None, override with --filter)\n"
      "    A query command to use to filter which Configurations are kept.          \n\n"

      "  dry_run: bool (optional, default=false, override with --dry-run)\n"
      "    Perform dry run.\n\n"

      "  primitive_only: bool (optional, default=true)\n"
      "    If true, only the primitive form of a configuration is saved in the      \n"
      "    configuration list. Otherwise, both primitive and non-primitive          \n"
      "    configurations are saved. \n\n"

      "  output_configurations: bool (optional, default=false)\n"
      "    If true, write formatted data for each enumerated configuration. Formatting options are \n"
      "    given by `output_options`. \n\n"

      "  output_configurations_options: object (optional) \n"
      "    Set output options for when `output_configurations==true`. \n\n"

      "      path: string (optional, default=\"enum.out\") Output file name.\n"
      "      json: bool (optional, default=false) If true, write JSON output files. Else CSV style.\n"
      "      json_arrays: bool (optional, default=false) If true, write data in JSON arrays. \n"
      "      compress: bool (optional, default=false) If true, compress data using gz. If `path` \n"
      "      does not end in '.gz' it will be appended. \n"
      "      output_filtered_configurations: bool (optional, default=false) If true, also include \n"
      "      output from configurations that excluded by the `filter` option.\n\n";
  };
}
