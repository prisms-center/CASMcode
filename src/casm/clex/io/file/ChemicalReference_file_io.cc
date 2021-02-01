#include "casm/clex/io/file/ChemicalReference_file_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/io/json/ChemicalReference_json_io.hh"

namespace CASM {

/// \brief Read chemical reference states from JSON file
///
/// See documentation with `jsonConstructor<ChemicalReference>` for expected
/// format.
ChemicalReference read_chemical_reference(fs::path filename,
                                          BasicStructure const &prim,
                                          double tol) {
  try {
    jsonParser json{filename};
    if (json.find("chemical_reference") == json.end()) {
      throw std::runtime_error(
          "Error reading chemical reference states: Expected "
          "\"chemical_reference\" entry");
    }
    return json["chemical_reference"].get<ChemicalReference>(prim, tol);
  } catch (...) {
    err_log() << "Error reading chemical reference states from " << filename
              << std::endl;
    /// re-throw exceptions
    throw;
  }
}

/// \brief Write chemical reference states to JSON file
///
/// See documentation with `to_json(ChemicalReference const&, jsonParser &)`
/// for JSON format.
void write_chemical_reference(const ChemicalReference &chem_ref,
                              fs::path filename) {
  SafeOfstream outfile;
  outfile.open(filename);

  jsonParser json;
  to_json(chem_ref, json["chemical_reference"]);
  json.print(outfile.ofstream());

  outfile.close();
}

}  // namespace CASM
