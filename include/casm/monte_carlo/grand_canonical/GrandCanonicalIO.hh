#ifndef CASM_GrandCanonicalIO_HH
#define CASM_GrandCanonicalIO_HH

#include <string>

namespace CASM {

  template<typename T, typename U> class GenericDatumFormatter;
  template<typename T> class DataFormatter;
  class PrimClex;
  class jsonParser;

  class MonteCarlo;
  typedef const MonteCarlo *ConstMonteCarloPtr;
  class GrandCanonical;
  class MonteSettings;
  class Log;
  class GrandCanonicalConditions;

  /// \brief Make a LTE results formatter
  DataFormatter<ConstMonteCarloPtr> make_results_formatter(const GrandCanonical &mc);

  /// \brief Make a results formatter
  DataFormatter<ConstMonteCarloPtr> make_lte_results_formatter(const GrandCanonical &mc, const double &phi_LTE1, const std::string &configname);


  /// \brief Store GrandCanonicalConditions in JSON format
  jsonParser &to_json(const GrandCanonicalConditions &conditions, jsonParser &json);

  /// \brief Read GrandCanonicalConditions from JSON format
  void from_json(GrandCanonicalConditions &conditions, const PrimClex &primclex, const jsonParser &json);


  /// \brief Print single spin flip LTE
  GenericDatumFormatter<double, ConstMonteCarloPtr> GrandCanonicalLTEFormatter(const double &phi_LTE1);

  /// \brief Will create new file or append to existing results file the results of the latest run
  void write_lte_results(const MonteSettings &settings, const GrandCanonical &mc, const double &phi_LTE1, const std::string &configname, Log &_log);

}

#endif
