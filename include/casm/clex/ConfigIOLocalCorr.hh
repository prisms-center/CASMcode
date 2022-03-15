#ifndef CASM_clex_ConfigIOLocalCorr
#define CASM_clex_ConfigIOLocalCorr

#include "casm/casm_io/dataformatter/DataFormatterTools.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clusterography/IntegralCluster.hh"

namespace CASM {

class Configuration;

namespace ConfigIO {

/// \brief Returns site-centric correlation values
///
/// Evaluated basis function values for each site, normalized by cluster orbit
/// multiplicity. Periodic cluster expansions only.
///
class LocalCorr : public BaseValueFormatter<jsonParser, Configuration> {
 public:
  static const std::string Name;

  static const std::string Desc;

  LocalCorr();

  // --- Required implementations -----------

  /// \brief Returns the atom fraction
  jsonParser evaluate(Configuration const &config) const override;

  /// \brief Clone using copy constructor
  std::unique_ptr<LocalCorr> clone() const {
    return std::unique_ptr<LocalCorr>(this->_clone());
  }

  // --- Specialized implementation -----------

  /// \brief Returns true if initialization is ok
  bool validate(const Configuration &config) const override;

  /// \brief If not yet initialized, use the global clexulator from the PrimClex
  bool init(const Configuration &_tmplt) const override;

  /// \brief Parse arguments
  bool parse_args(const std::string &args) override;

 private:
  /// \brief Clone using copy constructor
  LocalCorr *_clone() const override { return new LocalCorr(*this); }

  mutable std::vector<Clexulator> m_clexulator;
  mutable std::string m_clex_name;
  mutable std::vector<Index> m_sublattice_to_asymmetric_unit;
  mutable std::optional<IntegralCluster> m_prototype_phenom;
  mutable std::vector<SymOp> m_equivalents_generating_ops;
  mutable std::vector<IntegralCluster> m_phenom_orbit;
};

}  // namespace ConfigIO
}  // namespace CASM

#endif
