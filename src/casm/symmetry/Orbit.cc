#include "casm/symmetry/Orbit_impl.hh"

namespace CASM {
namespace Orbit_impl {

std::ostream &operator<<(std::ostream &sout, const RelEqMap &map) {
  sout << "RelEqMap  a: " << map.a << std::endl;
  for (const auto &row : map.map) {
    sout << "  b: " << row.b << "  op: ";
    for (const auto &val : row.values) {
      sout << val << " ";
    }
    sout << std::endl;
  }
  return sout;
}
}  // namespace Orbit_impl
}  // namespace CASM
