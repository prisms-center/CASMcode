#include "gtest/gtest.h"

/// What is being tested:
///   the command line executable

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/app/casm_functions.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte_carlo/Conversions.hh"
#include "casm/monte_carlo/OccCandidate.hh"
#include "casm/monte_carlo/OccLocation.hh"

using namespace CASM;


void random_config(Configuration &config, Monte::Conversions &convert, MTRand &mtrand) {
  config.init_occupation();
  for(Index l = 0; l < config.size(); ++l) {
    int Nocc = convert.occ_size(convert.l_to_asym(l));
    config.set_occ(l, mtrand.randInt(Nocc - 1));
  }
}

void dilute_config(Configuration &config, Monte::Conversions &convert, MTRand &mtrand) {
  config.init_occupation();
  for(Index i = 0; i < convert.species_size(); ++i) {
    for(Index l = 0; l < config.size(); ++l) {
      Index asym = convert.l_to_asym(l);
      if(config.occ(l) == 0 &&
         convert.species_allowed(asym, i)) {
        config.set_occ(l, convert.occ_index(asym, i));
        break;
      }
    }
  }
}

void check_occ_init(Configuration &config, Monte::OccLocation &occ_loc, Monte::Conversions &convert, Monte::OccCandidateList &cand_list) {
  // check OccLocation initialization
  for(Index mol_id = 0; mol_id < occ_loc.size(); ++mol_id) {
    ASSERT_EQ(mol_id, occ_loc.mol(mol_id).id);
  }

  for(Index l = 0; l < config.size(); ++l) {
    Index mol_id = occ_loc.l_to_mol_id(l);
    if(mol_id == occ_loc.size()) { // non-variable site
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

void check_occ(Configuration &config, Monte::OccEvent &e, Monte::OccLocation &occ_loc, Monte::Conversions &convert, Monte::OccCandidateList &cand_list) {
  // check that occ_loc / config / mol are consistent for initial state of config
  for(const auto &occ : e.occ_transform) {
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

template<typename ProjType, typename ConfigInit>
void run_case(ProjType &proj, ConfigInit f, MTRand &mtrand) {
  proj.check_init();
  proj.check_composition();

  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);

  Eigen::Matrix3i T;
  T << 9, 0, 0,
  0, 9, 0,
  0, 0, 9;
  Supercell scel(&primclex, T);
  Monte::Conversions convert(scel);

  // config with random occupation
  Configuration config(scel);
  f(config, convert, mtrand);

  // construct OccCandidateList
  Monte::OccCandidateList cand_list(convert);
  //std::cout << std::make_pair(cand_list, convert) << std::endl;

  // construct OccLocation
  Monte::OccLocation occ_loc(convert, cand_list);
  occ_loc.initialize(config);

  check_occ_init(config, occ_loc, convert, cand_list);

  Index count = 0;
  Monte::OccEvent e;
  ConfigDoF &configdof = config.configdof();
  while(count < 1000000) {
    //if(count % 100000 == 0) { std::cout << "count: " << count << std::endl; }
    occ_loc.propose_canonical(e, cand_list.canonical_swap(), mtrand);
    check_occ(config, e, occ_loc, convert, cand_list);
    occ_loc.apply(e, configdof);
    check_occ(config, e, occ_loc, convert, cand_list);
    ++count;
  }
}


TEST(OccLocationTest, ZrO_RandomConfig) {

  MTRand mtrand;
  test::ZrOProj proj;
  run_case(proj, random_config, mtrand);
}

TEST(OccLocationTest, ZrO_DiluteConfig) {

  MTRand mtrand;
  test::ZrOProj proj;
  run_case(proj, dilute_config, mtrand);
}
TEST(OccLocationTest, FCCTernary_RandomConfig) {

  MTRand mtrand;
  test::FCCTernaryProj proj;
  run_case(proj, random_config, mtrand);
}

TEST(OccLocationTest, FCCTernary_DiluteConfig) {

  MTRand mtrand;
  test::FCCTernaryProj proj;
  run_case(proj, dilute_config, mtrand);
}
