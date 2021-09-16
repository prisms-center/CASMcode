#include "casm/monte2/Sampling.hh"

#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte2/Conversions.hh"
#include "casm/monte2/events/OccCandidate.hh"
#include "casm/monte2/events/OccEventProposal.hh"
#include "casm/monte2/events/OccLocation.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

void run_case(std::shared_ptr<Structure const> shared_prim, MTRand &mtrand,
              Monte2::Sampler &sampler) {
  ScopedNullLogging logging;

  Eigen::Matrix3l T;
  T << 9, 0, 0, 0, 9, 0, 0, 0, 9;
  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, T);
  Monte2::Conversions convert(shared_supercell);

  // config with default occupation
  Configuration config(shared_supercell);
  Monte2::State state{config};

  // construct OccCandidateList
  Monte2::OccCandidateList cand_list(convert);

  // construct OccLocation
  Monte2::OccLocation occ_loc(convert, cand_list);
  occ_loc.initialize(state.configuration.occupation());

  Index count = 0;
  Monte2::OccEvent e;
  ConfigDoF &configdof = state.configuration.configdof();
  clexulator::ConfigDoFValues &dof_values = configdof.values();
  while (count < 1000000) {
    if (count % 1000 == 0) {
      sampler.sample(state);
    }
    propose_canonical_event(e, occ_loc, cand_list.canonical_swap(), mtrand);
    occ_loc.apply(e, dof_values.occupation);
    ++count;
  }
}

TEST(SamplingTest, CompNSamplingTest) {
  MTRand mtrand;
  auto shared_prim = std::make_shared<Structure const>(test::ZrO_prim());

  Monte2::StateSamplingFunction comp_n_sampling_f(
      "comp_n", "Composition per unit cell",
      [](Monte2::State const &state) { return comp_n(state.configuration); });

  Monte2::Sampler comp_n_sampler(comp_n_sampling_f);

  run_case(shared_prim, mtrand, comp_n_sampler);

  EXPECT_EQ(comp_n_sampler.n_samples(), 1000);
  EXPECT_EQ(comp_n_sampler.n_components(), 3);
}
