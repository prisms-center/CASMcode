#include "casm/app/sym/write_prim_symmetry.hh"

#include <boost/filesystem/fstream.hpp>

#include "casm/app/APICommand.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/io/json_io.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/json_io.hh"

namespace CASM {

class DirectoryStructure;
class Structure;

/// Write/print prim symmetry
///
/// Notes:
/// - The default `casm sym` action writes the following symmetry files:
///   - Lattice point group: <root>/symmetry/lattice_point_group.json
///   - Crystal factor group: <root>/symmetry/factor_group.json
///   - Crystal point group: <root>/symmetry/crystal_point_group.json
/// - Optionally, it also prints symmetry information to `log`:
///   - Print lattice point group with `print_lattice_point_group`
///   - Print crystal factor group with `print_factor_group`
///   - Print crystal point group with `print_crystal_point_group`
///   - Control coordinate printing mode (FRAC vs CART) with `coordtype`
/// - Use cases:
///   - Use this method to print symmetry info to screen
///   - Use this method to write symmetry files if they have been deleted
///   - Use this method to write symmetry files if the crystallography tolerance
///   in project
///     settings is modified
void write_prim_symmetry_impl(Structure const &prim,
                              DirectoryStructure const &dir,
                              COORD_TYPE coordtype, bool brief, Log &log,
                              bool print_lattice_point_group,
                              bool print_factor_group,
                              bool print_crystal_point_group) {
  // Selects coordinate printing mode (CART vs FRAC) == coordtype until "C" goes
  // out of scope
  xtal::COORD_MODE C(coordtype);
  SymGroup lattice_pg{SymGroup::lattice_point_group(prim.lattice())};

  if (!print_lattice_point_group && !print_factor_group &&
      !print_crystal_point_group) {
    log << "  Lattice point group size: " << lattice_pg.size() << std::endl;
    log << "  Lattice point group is: " << lattice_pg.get_name() << std::endl
        << std::endl;

    log << "  Factor group size: " << prim.factor_group().size() << std::endl;
    log << "  Crystal point group is: " << prim.point_group().get_name()
        << std::endl;
  }

  auto print_group = [&](std::string name, SymGroup const &group) {
    if (brief) {
      log << name << ": (" << xtal::COORD_MODE::NAME(coordtype)
          << " representation)" << std::endl;
      SymInfoOptions opt{coordtype};
      brief_description(log, group, prim.lattice(), opt);
    } else {
      log << name << ":" << std::endl;
      group.print(log, coordtype);
    }
  };

  if (print_lattice_point_group) {
    print_group("Lattice point group", lattice_pg);
  }

  if (print_factor_group) {
    print_group("Factor group", prim.factor_group());
  }

  if (print_crystal_point_group) {
    print_group("Crystal point group", prim.point_group());
  }

  // Write symmetry info files
  dir.new_symmetry_dir();

  auto write_group = [&](fs::path path, SymGroup const &group) {
    fs::ofstream outfile;
    jsonParser json;
    outfile.open(path);
    write_symgroup(group, json);
    json.print(outfile);
    outfile.close();
  };

  write_group(dir.lattice_point_group(), lattice_pg);
  write_group(dir.factor_group(), prim.factor_group());
  write_group(dir.crystal_point_group(), prim.point_group());
}

/// Describe the default `casm sym` option
std::string write_prim_symmetry_desc() {
  std::string description =

      "Write / print prim symmetry information (default behavior):             "
      " \n\n"

      "  The default `casm sym` action writes the following symmetry files:    "
      " \n"
      "  - Lattice point group: <root>/symmetry/lattice_point_group.json       "
      " \n"
      "  - Crystal factor group: <root>/symmetry/factor_group.json             "
      " \n"
      "  - Crystal point group: <root>/symmetry/crystal_point_group.json       "
      " \n\n"

      "  Optionally, it also prints symmetry information:                      "
      " \n"
      "  - Print lattice point group with --lattice-point-group                "
      " \n"
      "  - Print crystal factor group with --print-factor-group                "
      " \n"
      "  - Print crystal point group with --print-crystal-point-group          "
      " \n"
      "  - Print brief symmetry operation descriptions with --brief        "
      " \n"
      "  - Control coordinate printing mode (FRAC vs CART) with --coord        "
      " \n\n"

      "  Use cases:                                                            "
      " \n"
      "  - Use this method to print symmetry info to screen                    "
      " \n"
      "  - Use this method to write symmetry files if they have been deleted   "
      " \n"
      "  - Use this method to write symmetry files if the crystallography      "
      " \n"
      "    tolerance in project settings is modified                           "
      " \n\n\n";

  return description;
}

/// Write/print prim symmetry
void write_prim_symmetry(PrimClex &primclex, jsonParser const &json_options,
                         jsonParser const &cli_options_as_json) {
  jsonParser json_combined{json_options};
  std::map<std::string, std::string> cli_to_combined_keys{
      {"print_lattice_point_group",
       "print_lattice_point_group"},                 // --lattice-point-group
      {"print_factor_group", "print_factor_group"},  // --factor-group
      {"print_crystal_point_group",
       "print_crystal_point_group"},           // --crystal-point-group
      {"coordinate_mode", "coordinate_mode"},  // --coord
      {"brief", "brief"}                       // --brief
  };
  combine_json_options(cli_to_combined_keys, cli_options_as_json,
                       json_combined);

  ParentInputParser parser{json_combined};

  bool print_lattice_point_group;
  bool print_factor_group;
  bool print_crystal_point_group;
  COORD_TYPE coordtype;
  bool brief;

  parser.optional_else(print_lattice_point_group, "print_lattice_point_group",
                       false);
  parser.optional_else(print_factor_group, "print_factor_group", false);
  parser.optional_else(print_crystal_point_group, "print_crystal_point_group",
                       false);
  parser.optional_else(coordtype, "coordinate_mode", FRAC);
  std::runtime_error error_if_invalid{"Error reading `casm sym` input"};
  report_and_throw_if_invalid(parser, log(), error_if_invalid);

  parser.optional_else(brief, "brief", false);

  write_prim_symmetry_impl(primclex.prim(), primclex.dir(), coordtype, brief,
                           log(), print_lattice_point_group, print_factor_group,
                           print_crystal_point_group);
}

}  // namespace CASM
