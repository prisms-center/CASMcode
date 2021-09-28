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

/// Set config to random occupation
void random_config(Configuration &config, Monte2::Conversions &convert,
                   MTRand &mtrand) {
  config.init_occupation();
  for (Index l = 0; l < config.size(); ++l) {
    int Nocc = convert.occ_size(convert.l_to_asym(l));
    config.set_occ(l, mtrand.randInt(Nocc - 1));
  }
}

/// Set a single non-default occupant
void dilute_config(Configuration &config, Monte2::Conversions &convert,
                   MTRand &mtrand) {
  config.init_occupation();
  for (Index i = 0; i < convert.species_size(); ++i) {
    for (Index l = 0; l < config.size(); ++l) {
      Index asym = convert.l_to_asym(l);
      if (config.occ(l) == 0 && convert.species_allowed(asym, i)) {
        config.set_occ(l, convert.occ_index(asym, i));
        break;
      }
    }
  }
}

void check_occ_init(Configuration &config, Monte2::OccLocation &occ_loc,
                    Monte2::Conversions &convert,
                    Monte2::OccCandidateList &cand_list) {
  // check OccLocation initialization
  for (Index mol_id = 0; mol_id < occ_loc.mol_size(); ++mol_id) {
    ASSERT_EQ(mol_id, occ_loc.mol(mol_id).id);
  }

  for (Index l = 0; l < config.size(); ++l) {
    Index mol_id = occ_loc.l_to_mol_id(l);
    if (mol_id == occ_loc.mol_size()) {  // non-variable site
      continue;
    }
    auto &mol = occ_loc.mol(mol_id);

    ASSERT_EQ(mol.l, l);
    ASSERT_EQ(mol.id, mol_id);

    Index asym = convert.l_to_asym(l);
    ASSERT_EQ(asym, mol.asym);

    ASSERT_EQ(config.occ(l), convert.occ_index(asym, mol.species_index));
    ASSERT_EQ(convert.species_index(asym, config.occ(l)), mol.species_index);

    Index cand_index = cand_list.index(mol.asym, mol.species_index);
    ASSERT_EQ(mol.id, occ_loc.mol_id(cand_index, mol.loc));
  }
}

void check_occ(Configuration &config, Monte2::OccEvent &e,
               Monte2::OccLocation &occ_loc, Monte2::Conversions &convert,
               Monte2::OccCandidateList &cand_list) {
  // check that occ_loc / config / mol are consistent for initial state of
  // config
  for (const auto &occ : e.occ_transform) {
    Index l = occ.l;
    Index mol_id = occ_loc.l_to_mol_id(l);
    auto &mol = occ_loc.mol(mol_id);

    ASSERT_EQ(mol.l, l);
    ASSERT_EQ(mol.id, mol_id);
    ASSERT_EQ(occ.mol_id, mol_id);

    Index asym = convert.l_to_asym(l);
    ASSERT_EQ(asym, mol.asym);

    ASSERT_EQ(config.occ(l), convert.occ_index(asym, mol.species_index));
    ASSERT_EQ(convert.species_index(asym, config.occ(l)), mol.species_index);

    Index cand_index = cand_list.index(mol.asym, mol.species_index);
    ASSERT_EQ(mol.id, occ_loc.mol_id(cand_index, mol.loc));
  }
}

template <typename ConfigInit>
void run_case(std::shared_ptr<Structure const> shared_prim, ConfigInit f,
              MTRand &mtrand) {
  ScopedNullLogging logging;

  Eigen::Matrix3l T;
  T << 9, 0, 0, 0, 9, 0, 0, 0, 9;
  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, T);
  Monte2::Conversions convert(*shared_prim, T);

  // config with random occupation
  Configuration config(shared_supercell);
  f(config, convert, mtrand);

  // construct OccCandidateList
  Monte2::OccCandidateList cand_list(convert);

  // construct OccLocation
  Monte2::OccLocation occ_loc(convert, cand_list);
  occ_loc.initialize(config.occupation());

  check_occ_init(config, occ_loc, convert, cand_list);

  Index count = 0;
  Monte2::OccEvent e;
  ConfigDoF &configdof = config.configdof();
  clexulator::ConfigDoFValues &dof_values = configdof.values();
  while (count < 1000000) {
    // if(count % 100000 == 0) { std::cout << "count: " << count << std::endl; }
    propose_canonical_event(e, occ_loc, cand_list.canonical_swap(), mtrand);
    check_occ(config, e, occ_loc, convert, cand_list);
    occ_loc.apply(e, dof_values.occupation);
    check_occ(config, e, occ_loc, convert, cand_list);
    ++count;
  }
}

TEST(OccLocationTest, ZrO_RandomConfig) {
  MTRand mtrand;
  auto shared_prim = std::make_shared<Structure const>(test::ZrO_prim());
  run_case(shared_prim, random_config, mtrand);
}

TEST(OccLocationTest, ZrO_DiluteConfig) {
  MTRand mtrand;
  auto shared_prim = std::make_shared<Structure const>(test::ZrO_prim());
  run_case(shared_prim, dilute_config, mtrand);
}
TEST(OccLocationTest, FCCTernary_RandomConfig) {
  MTRand mtrand;
  auto shared_prim =
      std::make_shared<Structure const>(test::FCC_ternary_prim());
  run_case(shared_prim, random_config, mtrand);
}

TEST(OccLocationTest, FCCTernary_DiluteConfig) {
  MTRand mtrand;
  auto shared_prim =
      std::make_shared<Structure const>(test::FCC_ternary_prim());
  run_case(shared_prim, dilute_config, mtrand);
}
