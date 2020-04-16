#ifndef SIMPLESTRUCTUREIO_HH
#define SIMPLESTRUCTUREIO_HH

#include <set>
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/enum.hh"

namespace CASM {
  namespace xtal {
    class SimpleStructure;
  }
  /// \brief Output to JSON, excluding any molecular or atomic species contained in 'excluded_species'
  jsonParser &to_json(xtal::SimpleStructure const &_struc,
                      jsonParser &json_supplement,
                      std::set<std::string> const &excluded_species = {"Va", "VA", "va"},
                      std::string prefix = "",
                      COORD_TYPE mode = CART);

  /// \brief Read from JSON
  void from_json(xtal::SimpleStructure &_struc, const jsonParser &json, std::string prefix = "");

}

#endif
