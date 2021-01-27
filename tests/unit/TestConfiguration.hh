#ifndef UNIT_TestConfiguration
#define UNIT_TestConfiguration

#include "TestSupercell.hh"

using namespace CASM;

namespace test {

struct TestConfiguration : public test::TestSupercell {
  TestConfiguration(const PrimClex &primclex, const Configuration &_config);

  TestConfiguration(const PrimClex &primclex, const Eigen::Matrix3l &T,
                    const std::vector<int> &_occupation);

  TestConfiguration(const PrimClex &primclex, const Lattice &lat,
                    const std::vector<int> &_occupation);

  TestConfiguration(const PrimClex &primclex, const Configuration &unit,
                    const Eigen::Matrix3l &T, double tol = TOL);

  Configuration config;
  const std::vector<PermuteIterator> &config_permute_fg() const;
  const SymGroup &config_sym_fg() const;

  mutable std::vector<PermuteIterator> m_config_permute_fg;
  mutable SymGroup m_config_sym_fg;
};
}  // namespace test

#endif
