#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte2/Conversions.hh"
#include "casm/monte2/events/OccCandidate.hh"
#include "casm/monte2/events/OccEventProposal.hh"
#include "casm/monte2/events/OccLocation.hh"
#include "casm/monte2/state/State.hh"
#include "casm/monte2/state/StateSampler.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

void random_config(Configuration &config, Monte2::Conversions &convert,
                   MTRand &mtrand);

void run_case(std::shared_ptr<Structure const> shared_prim, Eigen::Matrix3l T,
              MTRand &mtrand,
              Monte2::StateSampler<CASM::Configuration> &sampler) {
  ScopedNullLogging logging;

  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, T);
  Monte2::Conversions convert(*shared_prim, T);

  // config with default occupation
  Configuration config(shared_supercell);
  Monte2::State<Configuration> state{config};
  random_config(state.configuration, convert, mtrand);
  Eigen::VectorXi &occupation =
      state.configuration.configdof().values().occupation;

  // construct OccCandidateList
  Monte2::OccCandidateList cand_list(convert);

  // construct OccLocation
  Monte2::OccLocation occ_loc(convert, cand_list);
  occ_loc.initialize(occupation);

  Index count = 0;
  Monte2::OccEvent e;
  while (count < 1000000) {
    if (count % 1000 == 0) {
      sampler.sample(state);
    }
    propose_canonical_event(e, occ_loc, cand_list.canonical_swap(), mtrand);
    occ_loc.apply(e, occupation);
    ++count;
  }
}

TEST(SamplingTest, CompNSamplingTest) {
  MTRand mtrand;
  auto shared_prim = std::make_shared<Structure const>(test::ZrO_prim());

  Eigen::Matrix3l T;
  T << 9, 0, 0, 0, 9, 0, 0, 0, 9;
  xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter(
      T, shared_prim->basis().size());
  std::vector<std::string> components = xtal::struc_molecule_name(*shared_prim);

  CompositionCalculator composition_calculator(
      unitcellcoord_index_converter, components,
      xtal::make_index_converter(*shared_prim, components));

  Monte2::StateSamplingFunction<Configuration> comp_n_sampling_f(
      "comp_n", "Composition per unit cell",
      [&](Monte2::State<Configuration> const &state) {
        return composition_calculator.mean_num_each_component(
            state.configuration.occupation());
      });

  Monte2::StateSampler<Configuration> comp_n_sampler(comp_n_sampling_f);

  std::map<std::string, std::shared_ptr<Monte2::Sampler>> samplers;
  samplers["comp_n"] = std::make_shared<Monte2::Sampler>();

  std::map<std::string, Monte2::StateSampler<Configuration>> state_samplers;
  state_samplers.emplace(
      std::piecewise_construct, std::forward_as_tuple("comp_n"),
      std::forward_as_tuple(comp_n_sampling_f, samplers["comp_n"]));

  run_case(shared_prim, T, mtrand, state_samplers.at("comp_n"));

  Monte2::Sampler const &sampler = *samplers.at("comp_n");
  EXPECT_EQ(sampler.n_samples(), 1000);
  EXPECT_EQ(sampler.n_components(), 3);
}
