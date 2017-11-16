#ifndef UNIT_TestDiffTransConfigOccPerturbations
#define UNIT_TestDiffTransConfigOccPerturbations

#include "TestConfiguration.hh"

using namespace CASM;

namespace test {

  struct TestDiffTransConfigOccPerturbations : public test::TestConfiguration {

    TestDiffTransConfigOccPerturbations(
      const TestConfiguration &test_config,
      const DiffusionTransformation &diff_trans,
      const std::vector<int> &_occupation);

    TestDiffTransConfigOccPerturbations(
      const PrimClex &primclex,
      const Lattice &lat,
      const std::vector<int> &_occupation);

    Configuration config;
    std::vector<PermuteIterator> config_permute_fg;
    SymGroup config_sym_fg;
  };
}

#endif
