#include "casm/app/sym/symmetrize.hh"

#include <boost/filesystem.hpp>

#include "casm/app/APICommand.hh"
#include "casm/app/io/json_io.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/symmetry/SymGroup.hh"

namespace symmetrize_impl {
using namespace CASM;

void _print_factor_group_convergence(const Structure &struc, double small_tol,
                                     double large_tol, double increment,
                                     std::ostream &print_stream) {
  std::vector<double> tols;
  std::vector<bool> is_group;
  std::vector<int> num_ops, num_enforced_ops;
  std::vector<std::string> name;

  xtal::Lattice lattice = struc.lattice();

  double orig_tol = lattice.tol();
  for (double i = small_tol; i < large_tol; i += increment) {
    tols.push_back(i);
    lattice.set_tol(i);

    xtal::SymOpVector factor_group_operations =
        xtal::make_factor_group(struc.structure());
    CASM::SymGroup factor_group =
        adapter::Adapter<SymGroup, xtal::SymOpVector>()(factor_group_operations,
                                                        lattice);

    factor_group.get_multi_table();
    num_ops.push_back(factor_group.size());
    is_group.push_back(factor_group.is_group(i));
    factor_group.enforce_group(i);
    num_enforced_ops.push_back(factor_group.size());
    name.push_back(factor_group.get_name());
  }
  lattice.set_tol(orig_tol);

  for (Index i = 0; i < tols.size(); i++) {
    print_stream << tols[i] << "\t" << num_ops[i] << "\t" << is_group[i] << "\t"
                 << num_enforced_ops[i] << "\t name: " << name[i] << "\n";
  }

  return;
}

/// Return SymGroup calculated for basic_structure with specified tolerance
std::vector<xtal::SymOp> make_enforced_factor_group(
    xtal::BasicStructure basic_structure, double enforced_tol) {
  // symmetrized lattice has tolerance set to "enforced_tol"
  Lattice symmetrized_lattice{
      xtal::symmetrize(basic_structure.lattice(), enforced_tol)
          .lat_column_mat(),
      enforced_tol};
  basic_structure.set_lattice(symmetrized_lattice, FRAC);

  return make_factor_group(basic_structure, enforced_tol);
}

/// Adjust a structure's lattice and basis to increase factor group symmetry
///
/// Method:
/// - Read a VASP POSCAR from "input_poscar_location" and, if not primitive,
/// make its primitive
///   structure.
/// - Construct a symmetrized lattice such that the point group of the lattice
/// is equal to that of
///   the input structure's lattice when using "enforced_tol" to check for
///   lattice equivalence. The symmetrized lattice vectors are obtaining through
///   a process of applying operations in the enforced point group, averaging,
///   and then rotating to match the original lattice orientation.
/// - Construct a symmetrized structure with factor group equal to that
/// generated for the input
///   structure when using "enforced_tol" to check for structure equivalence.
///   The symmetrized structure's basis is constructed by setting basis site
///   coordinates equal to the average coordinate positions found after applying
///   each operation in the enforced factor group to the original basis sites.
/// - Write the symmetrized structure to "output_poscar_location"
/// (default="POSCAR_sym")
///
void symmetrize_v1(fs::path poscar_path, double tol) {
  Log &log = CASM::log();
  // fs::path poscar_path = opt().poscar_path();
  // double tol = opt().tol();
  log << "\n***************************\n" << std::endl;
  log << "Symmetrizing: " << poscar_path << std::endl;
  log << "with tolerance: " << tol << std::endl;
  Structure struc(poscar_path);
  struc = Structure(xtal::make_primitive(struc));

  int biggest = struc.factor_group().size();
  xtal::BasicStructure basic_tmp = struc;
  // a) symmetrize the lattice vectors
  Lattice lat = basic_tmp.lattice();
  lat = xtal::symmetrize(lat, tol);
  lat.set_tol(tol);
  basic_tmp.set_lattice(lat, FRAC);

  Structure tmp(basic_tmp);

  tmp.factor_group();
  // b) find factor group with same tolerance
  symmetrize_impl::_print_factor_group_convergence(
      tmp, tmp.structure().lattice().tol(), tol,
      (tol - tmp.structure().lattice().tol()) / 10.0, log);
  // c) symmetrize the basis sites
  SymGroup g = tmp.factor_group();
  tmp = xtal::symmetrize(tmp, g);

  // TODO: Why are we doing this twice?
  g = tmp.factor_group();
  tmp = xtal::symmetrize(tmp, g);
  if (tmp.factor_group().is_group(tol) &&
      (tmp.factor_group().size() > biggest)) {
    struc = tmp;
  }
  struc = Structure(xtal::make_primitive(struc));
  fs::ofstream file_i;
  fs::path POSCARpath_i = "POSCAR_sym";
  file_i.open(POSCARpath_i);
  VaspIO::PrintPOSCAR p_i(xtal::make_simple_structure(struc),
                          struc.structure().title());
  p_i.print(file_i);
  file_i.close();
}

/// Adjust a structure's lattice and basis to increase factor group symmetry
///
/// Method:
/// - Read a VASP POSCAR from "input_poscar_location" and, if not primitive,
/// make its primitive structure.
/// - Construct a symmetrized lattice such that the point group of the lattice
/// is equal to that of the input structure's lattice when using "enforced_tol"
/// to check for lattice equivalence. The symmetrized lattice vectors are
/// obtaining through a process of applying operations in the enforced point
/// group, averaging, and then rotating to match the original lattice
/// orientation.
/// - Construct a symmetrized structure with factor group equal to that
/// generated for the input structure when using "enforced_tol" to check for
/// structure equivalence. The symmetrized structure's basis is constructed by
/// setting basis site coordinates equal to the average coordinate positions
/// found after applying each operation in the enforced factor group to the
/// original basis sites.
/// - Write the symmetrized structure to "POSCAR_sym" or "prim.sym.json"
///
/// \param input_poscar_location Location of VASP POSCAR style file to read
/// \param enforced_tol Tolerance used to generate the "enforced factor group"
/// \param input_tol Tolerance used to generate original input structure factor
/// group
///
void symmetrize_v2(fs::path input_poscar_location, double enforced_tol) {
  using namespace symmetrize_impl;

  Log &log = CASM::log();

  double primitive_tol = enforced_tol;

  log << "\n***************************\n" << std::endl;
  log << "Symmetrizing: " << input_poscar_location << std::endl;
  log << "with tolerance: " << enforced_tol << std::endl;

  // Read PRIM or prim.json:
  xtal::BasicStructure basic_structure;
  ParsingDictionary<AnisoValTraits> const *aniso_val_dict = nullptr;
  std::string prim_file_type;
  try {
    basic_structure = read_prim(input_poscar_location, enforced_tol,
                                aniso_val_dict, prim_file_type);
  } catch (std::exception &e) {
    log << "Error reading input structure" << std::endl;
    log << e.what() << std::endl;
    throw e;
  }

  fs::path output_poscar_location;
  if (prim_file_type == "vasp") {
    output_poscar_location = "POSCAR_sym";
  } else {
    output_poscar_location = "prim.sym.json";
  }

  // make primitive if necessary
  bool input_structure_is_primitive =
      xtal::is_primitive(basic_structure, primitive_tol);
  if (!input_structure_is_primitive) {
    log << "The input structure is not primitive." << std::endl;
    log << "The primitive structure will be generated and symmetrized."
        << std::endl;
    basic_structure = xtal::make_primitive(basic_structure, primitive_tol);
  }

  // symmetrize lattice using "enforced_tol"
  Lattice symmetrized_lattice{
      xtal::symmetrize(basic_structure.lattice(), enforced_tol)
          .lat_column_mat(),
      enforced_tol};
  basic_structure.set_lattice(symmetrized_lattice, FRAC);

  // make the factor group to be enforced
  std::vector<xtal::SymOp> enforced_factor_group =
      make_factor_group(basic_structure);

  // symmetrize the structure
  basic_structure = symmetrize(basic_structure, enforced_factor_group);

  Structure structure{basic_structure};

  // check factor group of the result and compare to original
  bool factor_group_valid = structure.factor_group().is_group(enforced_tol);
  Index factor_group_size = structure.factor_group().size();

  log << "Factor group size with enforced tolerance: " << factor_group_size
      << std::endl;

  // Print result
  if (prim_file_type == "vasp") {
    VaspIO::PrintPOSCAR printer{xtal::make_simple_structure(structure),
                                basic_structure.title()};

    fs::ofstream file{output_poscar_location};
    printer.print(file);
  } else {
    jsonParser json;
    bool include_va = true;
    write_prim(structure, json, FRAC, include_va);
    fs::ofstream file{output_poscar_location};
    file << json << std::endl;
  }
}
}  // namespace symmetrize_impl

namespace CASM {

/// Describe the symmetrize method
std::string symmetrize_desc() {
  std::string description =

      "Symmetrize a POSCAR (--symmetrize POSCAR_LOCATION --tol ENFORCED_TOL): "
      "\n\n"

      "  The --symmetrize option implements a method to adjust a          \n"
      "  structure's lattice and basis to increase factor group symmetry. \n\n"

      "  Method: \n"
      "  - Read a VASP POSCAR file from POSCAR_LOCATION and, if not       \n"
      "    primitive, make its primitive structure.                       \n"
      "  - Construct a symmetrized lattice such that the point group of   \n"
      "    the lattice is equal to that of the input structure's lattice  \n"
      "    when using ENFORCED_TOL to check for lattice equivalence. The  \n"
      "    symmetrized lattice vectors are obtaining through a process of \n"
      "    applying operations in the enforced point group, averaging, and\n"
      "    then rotating to match the original lattice orientation.       \n"
      "  - Construct a symmetrized structure with factor group equal to   \n"
      "    that generated for the input structure when using ENFORCED_TOL \n"
      "    to check for structure equivalence. The symmetrized structure's\n"
      "    basis is constructed by setting basis site coordinates equal to\n"
      "    the average coordinate positions found after applying each     \n"
      "    operation in the enforced factor group to the original basis   \n"
      "    sites.                                                         \n"
      "  - Write the symmetrized structure to \"POSCAR_sym\" or           \n"
      "    \"prim.sym.json\".                                             \n\n";

  return description;
}

/// Adjust a structure's lattice and basis to increase factor group symmetry
void symmetrize(jsonParser const &json_options,
                jsonParser const &cli_options_as_json) {
  Log &log = CASM::log();

  std::map<std::string, std::string> cli_to_combined_keys{
      {"tol", "tol"},               // --tol
      {"symmetrize", "symmetrize"}  // --symmetrize
  };

  jsonParser json_combined{json_options};
  combine_json_options(cli_to_combined_keys, cli_options_as_json,
                       json_combined);
  ParentInputParser parser{json_combined};

  std::string poscar_path;
  parser.require(poscar_path, "symmetrize");

  double enforced_tol;
  parser.require(enforced_tol, "tol");

  std::runtime_error error_if_invalid{
      "Error reading `casm sym --symmetrize` input"};
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // symmetrize_impl::symmetrize_v1(poscar_path, enforced_tol);

  symmetrize_impl::symmetrize_v2(poscar_path, enforced_tol);
}

}  // namespace CASM
