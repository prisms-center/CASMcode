#ifndef SIMPLESTRUCTUREIO_HH
#define SIMPLESTRUCTUREIO_HH

#include <set>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/enum.hh"

namespace CASM {
namespace xtal {
class SimpleStructure;
}
/// Output SimpleStructure to JSON
jsonParser &to_json(xtal::SimpleStructure const &simple_structure,
                    jsonParser &json,
                    std::set<std::string> const &excluded_species = {"Va", "VA",
                                                                     "va"},
                    COORD_TYPE coordinate_mode = CART);

/// Read SimpleStructure from JSON
void from_json(xtal::SimpleStructure &simple_structure, const jsonParser &json);

}  // namespace CASM

#endif
