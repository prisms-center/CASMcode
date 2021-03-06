#include "TestConfiguration.hh"

#include "Common.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/FillSupercell.hh"

namespace test {

TestConfiguration::TestConfiguration(const PrimClex &primclex,
                                     const Configuration &_config)
    : TestSupercell(primclex, _config.ideal_lattice()), config(_config) {}

TestConfiguration::TestConfiguration(const PrimClex &primclex,
                                     const Eigen::Matrix3l &T,
                                     Eigen::VectorXi const &_occupation)
    : TestConfiguration(primclex,
                        xtal::make_superlattice(primclex.prim().lattice(), T),
                        _occupation) {}

TestConfiguration::TestConfiguration(const PrimClex &primclex,
                                     const Lattice &lat,
                                     Eigen::VectorXi const &_occupation)
    : TestSupercell(primclex, lat), config(Configuration::zeros(this->scel)) {
  config.set_occupation(_occupation);
}

namespace {
Configuration make_superconfig(const Configuration &unit,
                               const Eigen::Matrix3l &T) {
  std::shared_ptr<Supercell> scel = std::make_shared<Supercell>(
      &unit.primclex(), unit.supercell().transf_mat() * T);
  FillSupercell f{scel};
  return f(unit);
}
}  // namespace

TestConfiguration::TestConfiguration(const PrimClex &primclex,
                                     const Configuration &unit,
                                     const Eigen::Matrix3l &T, double tol)
    : TestConfiguration(primclex, make_superconfig(unit, T)) {}

const std::vector<PermuteIterator> &TestConfiguration::config_permute_fg()
    const {
  if (!m_config_permute_fg.size()) {
    m_config_permute_fg = config.factor_group();
  }
  return m_config_permute_fg;
}

const SymGroup &TestConfiguration::config_sym_fg() const {
  if (!m_config_sym_fg.size()) {
    m_config_sym_fg =
        make_sym_group(config_permute_fg().begin(), config_permute_fg().end(),
                       this->scel.sym_info().supercell_lattice());
  }
  return m_config_sym_fg;
}

}  // namespace test
