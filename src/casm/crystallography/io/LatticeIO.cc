#include "casm/crystallography/io/LatticeIO.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {
// write Lattice in json as array of vectors
jsonParser &to_json(const xtal::Lattice &lat, jsonParser &json) {
  json.put_array();
  json.push_back(lat[0]);
  json.push_back(lat[1]);
  json.push_back(lat[2]);
  return json;
}

// read Lattice from a json array of Eigen::Vector3d
void from_json(xtal::Lattice &lat, const jsonParser &json, double xtal_tol) {
  try {
    lat = xtal::Lattice(json[0].get<Eigen::Vector3d>(),
                        json[1].get<Eigen::Vector3d>(),
                        json[2].get<Eigen::Vector3d>(), xtal_tol);
  } catch (...) {
    /// re-throw exceptions
    throw;
  }
}
}  // namespace CASM
