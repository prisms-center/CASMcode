#ifndef LATTICEIO_HH
#define LATTICEIO_HH

namespace CASM {
namespace xtal {
class Lattice;
}

class jsonParser;

// write Lattice in json as array of vectors
jsonParser &to_json(const xtal::Lattice &lat, jsonParser &json);
void from_json(xtal::Lattice &lat, const jsonParser &json, double xtal_tol);
}  // namespace CASM

#endif
