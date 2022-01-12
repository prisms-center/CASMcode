#include "casm/app/enum/methods/ScelEnumInterface.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/dataformatter/ScelEnumIO_impl.hh"
#include "casm/app/enum/enumerate_supercells_impl.hh"
#include "casm/app/enum/io/enumerate_supercells_json_io.hh"
#include "casm/casm_io/dataformatter/DatumFormatterAdapter.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/SupercellIO.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

std::string ScelEnumInterface::desc() const {
  return

      "ScelEnum: \n\n"

      "  min: int, >0 (optional, default=1, override with --min)               "
      "\n"
      "    The minimum volume supercell to enumerate. The volume is measured   "
      "\n"
      "    relative the unit cell being used to generate supercells.           "
      "\n"
      "\n"
      "  max: int, >= min (optional, override with --max)                      "
      "\n"
      "    The maximum volume supercell to enumerate. The volume is measured   "
      "\n"
      "    relative the unit cell being used to generate supercells. One of    "
      "\n"
      "    \"max\" or \"all\" must be given.                                   "
      "\n"
      "\n"
      "  dirs: string (optional, default=\"abc\")\n"
      "    This option may be used to restrict the supercell enumeration to 1, "
      "\n"
      "    2 or 3 of the lattice vectors, to get 1-, 2-, or 3-dimensional      "
      "\n"
      "    supercells. By specifying combinations of 'a', 'b', and 'c', you    "
      "\n"
      "    determine which of the unit cell lattice vectors you want to        "
      "\n"
      "    enumerate over. For example, to enumerate 1-dimensional supercells  "
      "\n"
      "    along the 'c' use \"dirs\":\"c\". If you want 2-dimensional         "
      "\n"
      "    supercells along the 'a' and 'c' lattice vectors, specify           "
      "\n"
      "    \"dirs\":\"ac\". \n"
      "\n"
      "  unit_cell: 3x3 matrix of int, or string (default=identity matrix)     "
      "\n"
      "    This option may be used to specify the unit cell. It may be         "
      "\n"
      "    specified using a 3x3 matrix of int, representing the "
      "transformation\n"
      "    matrix, T, such that U = P*T, where P are the primitive lattice     "
      "\n"
      "    and U are the unit cell lattice vectors. For example, a unit cell   "
      "\n"
      "    that whose lattice vectors are (2*a+b, b, c) (with respect to the   "
      "\n"
      "    the primitive cell vectors) could be specified using:\n"
      "\n"
      "      \"unit_cell\" : [\n"
      "        [2, 0, 0],\n"
      "        [1, 1, 0],\n"
      "        [0, 0, 1]\n"
      "       ]\n"
      "\n"
      "    Or it may be specified with the name of an existing supercell to "
      "use as \n"
      "    the unit cell, for example: \n"
      "\n"
      "      \"unit_cell\" : \"SCEL2_1_1_2_0_0_0\"\n"
      "\n"
      "  filter: string (optional, default=None, override with --filter)\n"
      "    A query command to use to filter which Configurations are kept.     "
      "     \n\n"
      "\n"
      "  dry_run: bool (optional, default=false, override with --dry-run)\n"
      "    Perform dry run.\n"
      "\n"
      "  output_supercells: bool (optional, default=false)\n"
      "    If true, write formatted data for each enumerated supercell. "
      "Formatting options are \n"
      "    given by `output_supercells_options`. \n\n"

      "  output_supercells_options: object (optional) \n"
      "    Set output options for when `output_supercells==true`. \n\n"

      "      path: string (optional, default=\"enum.out\") Output file "
      "name.\n"
      "      json: bool (optional, default=false) If true, write JSON "
      "output files. Else CSV style.\n"
      "      json_arrays: bool (optional, default=false) If true, write "
      "data in JSON arrays. \n"
      "      compress: bool (optional, default=false) If true, compress "
      "data using gz. If `path` \n"
      "      does not end in '.gz' it will be appended. \n"
      "      output_filtered_supercells: bool (optional, default=false) "
      "If true, also include \n"
      "      output from supercells that were excluded by the `filter` "
      "option.\n"
      "\n"
      "Examples:\n"
      "\n"
      "    To enumerate supercells up to and including size 4:\n"
      "      casm enum --method ScelEnum -i '{\"max\": 4}' \n"
      "\n"
      "    To enumerate 2d supercells up to and including size 4:\n"
      "      casm enum --method ScelEnum -i '{\"max\": 4, \"dirs\": \"ab\"}' \n"
      "\n"
      "    If the prim is primitive FCC, two dimensional supercells of the \n"
      "    conventional FCC unit cell up to and including 4x the unit cell "
      "volume\n"
      "    could be enumerated using:\n"
      "\n"
      "     casm enum --method ScelEnum -i \n"
      "     '{\n"
      "        \"min\": 1,\n"
      "        \"max\": 4,\n"
      "        \"dirs\": \"ab\",\n"
      "        \"unit_cell\" : [\n"
      "          [-1,  1,  1],\n"
      "          [ 1, -1,  1],\n"
      "          [ 1,  1, -1]\n"
      "        ]\n"
      "      }'\n"
      "\n";
}

std::string ScelEnumInterface::name() const { return "ScelEnum"; }

void ScelEnumInterface::run(PrimClex &primclex, jsonParser const &json_options,
                            jsonParser const &cli_options_as_json) const {
  // combine JSON options and CLI options
  jsonParser json_combined =
      combine_supercell_enum_json_options(json_options, cli_options_as_json);

  // Read input data from JSON
  ParentInputParser parser{json_combined};
  std::runtime_error error_if_invalid{"Error reading ScelEnum JSON input"};

  auto options_parser_ptr = parser.parse_as<EnumerateSupercellsOptions>(
      ScelEnumByProps::enumerator_name, primclex,
      primclex.settings().query_handler<Supercell>().dict());
  report_and_throw_if_invalid(parser, log(), error_if_invalid);
  EnumerateSupercellsOptions const &options = *options_parser_ptr->value;
  log().set_verbosity(options.verbosity);

  auto scel_enum_props_parser_ptr =
      parser.parse_as<xtal::ScelEnumProps>(primclex.db<Supercell>());
  xtal::ScelEnumProps const &scel_enum_props =
      *scel_enum_props_parser_ptr->value;

  // `output_supercells` formatter
  typedef ScelEnumData<ScelEnumByProps> ScelEnumDataType;
  DataFormatter<ScelEnumDataType> formatter;
  formatter.push_back(
      ScelEnumIO::name<ScelEnumDataType>(),
      ScelEnumIO::selected<ScelEnumDataType>(),
      ScelEnumIO::is_new<ScelEnumDataType>(),
      ScelEnumIO::is_existing<ScelEnumDataType>(),
      ScelEnumIO::transformation_matrix_to_super<ScelEnumDataType>());
  if (options.filter) {
    formatter.push_back(ScelEnumIO::is_excluded_by_filter<ScelEnumDataType>());
  }

  ScelEnumByProps enumerator{primclex.shared_prim(), scel_enum_props};
  enumerate_supercells(options, enumerator, primclex.db<Supercell>(),
                       formatter);
}

}  // namespace CASM
