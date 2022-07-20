#ifndef CASM_GrandCanonicalSettings
#define CASM_GrandCanonicalSettings

#include "casm/enumerator/OrderParameter.hh"
#include "casm/monte_carlo/MonteSettings.hh"

namespace CASM {
namespace Monte {

class GrandCanonical;
class GrandCanonicalConditions;

class GrandCanonicalSettings : public EquilibriumMonteSettings {
 public:
  /// \brief Default constructor
  GrandCanonicalSettings() {}

  /// \brief Construct EquilibriumMonteSettings by reading a settings JSON file
  GrandCanonicalSettings(const PrimClex &primclex, const fs::path &read_path);

  // --- GrandCanonicalConditions settings ---------------------

  /// \brief Expects initial_conditions
  GrandCanonicalConditions initial_conditions(GrandCanonical const &mc) const;

  /// \brief Expects final_conditions
  GrandCanonicalConditions final_conditions(GrandCanonical const &mc) const;

  /// \brief Expects incremental_conditions
  GrandCanonicalConditions incremental_conditions(
      GrandCanonical const &mc) const;

  /// \brief Expects custom_conditions
  std::vector<GrandCanonicalConditions> custom_conditions(
      GrandCanonical const &mc) const;

  // --- Project settings ---------------------

  /// \brief Get formation energy cluster expansion
  ClexDescription formation_energy(const PrimClex &primclex) const;

  /// \brief Make order parameter calculator
  std::shared_ptr<OrderParameter> make_order_parameter(
      const PrimClex &primclex) const;

  // --- Sampler settings ---------------------

  /// \brief Construct Samplers as specified in the Settings
  template <typename SamplerInsertIterator>
  SamplerInsertIterator samplers(const PrimClex &primclex,
                                 SamplerInsertIterator result) const;

 private:
  mutable bool m_order_parameter_checked = false;
  mutable std::shared_ptr<OrderParameter> m_order_parameter;

  GrandCanonicalConditions _conditions(std::string name,
                                       GrandCanonical const &mc) const;
  GrandCanonicalConditions _conditions(const jsonParser &json,
                                       GrandCanonical const &mc) const;

  template <typename jsonParserIteratorType>
  std::tuple<bool, double> _get_precision(jsonParserIteratorType it) const;

  template <typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator _make_comp_samplers(const PrimClex &primclex,
                                            jsonParserIteratorType it,
                                            SamplerInsertIterator result) const;

  template <typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator _make_comp_n_samplers(
      const PrimClex &primclex, jsonParserIteratorType it,
      SamplerInsertIterator result) const;

  template <typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator _make_site_frac_samplers(
      const PrimClex &primclex, jsonParserIteratorType it,
      SamplerInsertIterator result) const;

  template <typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator _make_atom_frac_samplers(
      const PrimClex &primclex, jsonParserIteratorType it,
      SamplerInsertIterator result) const;

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
