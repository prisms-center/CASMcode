#ifndef UNIT_TestConfiguration
#define UNIT_TestConfiguration

#include "TestSupercell.hh"

using namespace CASM;

namespace test {

  struct TestConfiguration : public test::TestSupercell {

    TestConfiguration(
      const PrimClex &primclex,
      const Eigen::Matrix3i &T,
      const std::vector<int> &_occupation);

    TestConfiguration(
      const PrimClex &primclex,
      const Lattice &lat,
      const std::vector<int> &_occupation);

    Configuration config;
    std::vector<PermuteIterator> config_permute_fg;
    SymGroup config_sym_fg;
  };
}

#endif
