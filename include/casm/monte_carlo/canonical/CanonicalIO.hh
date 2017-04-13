#ifndef CASM_CanonicalIO_HH
#define CASM_CanonicalIO_HH

#include <string>
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/monte_carlo/canonical/Canonical.hh"
#include "casm/monte_carlo/MonteIO.hh"

namespace CASM {

  template<typename T>
  class DataFormatter;

  /// \brief Make a LTE results formatter
  DataFormatter<ConstMonteCarloPtr> make_results_formatter(const Canonical &mc);

  /// \brief Store CanonicalConditions in JSON format
  jsonParser &to_json(const CanonicalConditions &conditions, jsonParser &json);

  /// \brief Read CanonicalConditions from JSON format
  void from_json(CanonicalConditions &conditions, const CompositionConverter &comp_converter, const jsonParser &json);

}

#endif
