#include "casm/container/io/PermutationIO.hh"

#include <vector>

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/container/Permutation.hh"

namespace CASM {

std::ostream &operator<<(std::ostream &out, const Permutation &perm) {
  out << perm.perm_array();
  return out;
}

//**************************************************************

jsonParser &to_json(const Permutation &perm, jsonParser &json) {
  return CASM::to_json(perm.perm_array(), json);
}

void from_json(Permutation &perm, const jsonParser &json) {
  try {
    std::vector<Index> permutation_vector;
    from_json(permutation_vector, json);
    perm = Permutation(permutation_vector);
  } catch (...) {
    /// re-throw exceptions
    throw;
  }
}
}  // namespace CASM
