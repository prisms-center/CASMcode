#ifndef CASM_CanonicalSettings
#define CASM_CanonicalSettings

#include "casm/enumerator/OrderParameter.hh"
#include "casm/monte_carlo/MonteSettings.hh"

namespace CASM {
namespace Monte {

class Canonical;
class CanonicalConditions;

class CanonicalSettings : public EquilibriumMonteSettings {
 public:
  /// \brief Default constructor
  CanonicalSettings() {}

  /// \brief Construct EquilibriumMonteSettings by reading a settings JSON file
  CanonicalSettings(const PrimClex &primclex, const fs::path &read_path);

  // --- CanonicalConditions settings ---------------------

  /// \brief Expects initial_conditions
  CanonicalConditions initial_conditions(Canonical const &mc) const;

  /// \brief Expects final_conditions
  CanonicalConditions final_conditions(Canonical const &mc) const;

  /// \brief Expects incremental_conditions
  CanonicalConditions incremental_conditions(Canonical const &mc) const;

  /// \brief Expects custom_conditions
  std::vector<CanonicalConditions> custom_conditions(Canonical const &mc) const;

  // --- Project settings ---------------------

  /// \brief Get formation energy cluster expansion
  ClexDescription formation_energy(const PrimClex &primclex) const;

  /// \brief Make order parameter calculator
  std::shared_ptr<OrderParameter> make_order_parameter(
      const PrimClex &primclex) const;

  // --- Sampler settings ---------------------

  /// \brief Construct MonteSamplers as specified in the MonteSettings
  template <typename SamplerInsertIterator>
  SamplerInsertIterator samplers(const PrimClex &primclex,
                                 SamplerInsertIterator result) const;

 private:
  CompositionConverter m_comp_converter;

  mutable bool m_order_parameter_checked = false;
  mutable std::shared_ptr<OrderParameter> m_order_parameter;

  CanonicalConditions _conditions(std::string name, Canonical const &mc,
                                  bool incremental = false) const;
  CanonicalConditions _conditions(const jsonParser &json, Canonical const &mc,
                                  bool incremental = false) const;

  template <typename jsonParserIteratorType>
  std::tuple<bool, double> _get_precision(jsonParserIteratorType it) const;

  template <typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator _make_non_zero_eci_correlations_samplers(
      const PrimClex &primclex, jsonParserIteratorType it,
      SamplerInsertIterator result) const;

  template <typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator _make_order_parameter_samplers(
      const PrimClex &primclex, jsonParserIteratorType it,
      SamplerInsertIterator result) const;

  template <typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator _make_query_samplers(
      const PrimClex &primclex, jsonParserIteratorType it,
      SamplerInsertIterator result) const;
};

}  // namespace Monte
}  // namespace CASM

#endif
