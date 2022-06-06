#include "casm/app/enum/standard_ConfigEnumInput_help.hh"

#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

// At some point the help text may go in text files, but for now...

std::string standard_ConfigEnumInput_help() {
  return parse_ConfigEnumInput_desc() +
         "  filter: string (optional, default=None, override with --filter)\n"
         "    A query command to use to filter which Configurations are kept.  "
         "        \n\n"

         "  dry_run: bool (optional, default=false, override with --dry-run)\n"
         "    Perform dry run.\n\n"

         "  primitive_only: bool (optional, default=true)\n"
         "    If true, only the primitive form of a configuration is saved in "
         "the      \n"
         "    configuration list. Otherwise, both primitive and non-primitive  "
         "        \n"
         "    configurations are saved. \n\n"

         "  output_configurations: bool (optional, default=false)\n"
         "    If true, write formatted data for each enumerated configuration. "
         "Formatting options are \n"
         "    given by `output_configurations_options`. \n\n"

         "  output_configurations_options: object (optional) \n"
         "    Set output options for when `output_configurations==true`. \n\n"

         "      path: string (optional, default=\"enum.out\") \n"
         "        Output file name.\n"
         "      json: bool (optional, default=false) \n"
         "        If true, write JSON output files. Else CSV style.\n"
         "      json_arrays: bool (optional, default=false) \n"
         "        If true, write data in JSON arrays. \n"
         "      compress: bool (optional, default=false) \n"
         "        If true, compress data using gz. If `path` does not end \n"
         "        in '.gz' it will be appended. \n"
         "      properties: array of strain (optional, default=[])        \n"
         "        Specify additional configuration properties to output,\n"
         "        using properties available from `casm query`.         \n"
         "      output_filtered_configurations: bool (optional, default=false) "
         "\n"
         "        If true, also include output from configurations that \n"
         "        excluded by the `filter` option.\n\n";
};
}  // namespace CASM
