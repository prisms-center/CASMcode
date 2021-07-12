#ifndef CASM_GrandCanonicalSettings
#define CASM_GrandCanonicalSettings

#include "casm/monte_carlo/MonteSettings.hh"

namespace CASM {
namespace Monte {

class GrandCanonicalConditions;

class GrandCanonicalSettings : public EquilibriumMonteSettings {
 public:
  /// \brief Default constructor
  GrandCanonicalSettings() {}

  /// \brief Construct EquilibriumMonteSettings by reading a settings JSON file
  GrandCanonicalSettings(const PrimClex &primclex, const fs::path &read_path);

  // --- GrandCanonicalConditions settings ---------------------

  /// \brief Expects initial_conditions
  GrandCanonicalConditions initial_conditions() const;

  /// \brief Expects final_conditions
  GrandCanonicalConditions final_conditions() const;

  /// \brief Expects incremental_conditions
  GrandCanonicalConditions incremental_conditions() const;

  /// \brief Expects custom_conditions
  std::vector<GrandCanonicalConditions> custom_conditions() const;

  // --- Project settings ---------------------

  /// \brief Get formation energy cluster expansion
  ClexDescription formation_energy(const PrimClex &primclex) const;

  // --- Sampler settings ---------------------

  /// \brief Construct Samplers as specified in the Settings
  template <typename SamplerInsertIterator>
  SamplerInsertIterator samplers(const PrimClex &primclex,
                                 SamplerInsertIterator result) const;

 private:
  GrandCanonicalConditions _conditions(std::string name) const;
  GrandCanonicalConditions _conditions(const jsonParser &json) const;

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
  SamplerInsertIterator _make_query_samplers(
      const PrimClex &primclex, jsonParserIteratorType it,
      SamplerInsertIterator result) const;
};

}  // namespace Monte
}  // namespace CASM

#endif
