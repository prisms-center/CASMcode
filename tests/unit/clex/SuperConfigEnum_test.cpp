#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/SuperConfigEnum.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(SuperConfigEnumTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

  Eigen::Matrix3d target_M;
  target_M << 8.0, 0.0, 0.0,
           0.0, 8.0, 0.0,
           0.0, 0.0, 8.0;

  Supercell target_scel(&primclex, Lattice {target_M});

  Eigen::Matrix3d motif_M;
  motif_M <<  4.0, 0.0, 0.0,
          0.0, 4.0, 0.0,
          0.0, 0.0, 4.0;

  Supercell motif_scel = Supercell(&primclex, Lattice {motif_M});

  auto sub_config = [&](std::initializer_list<int> occ) {
    Configuration tconfig(motif_scel);
    tconfig.set_occupation(occ);
    return tconfig;
  };

  std::vector<Configuration> sub_configs = {
    sub_config({0, 0, 0, 0}),
    sub_config({1, 1, 1, 1})
  };


  {
    SuperConfigEnum e(target_scel, sub_configs.begin(), sub_configs.end());
    BOOST_CHECK_EQUAL(pow(2, 8), std::distance(e.begin(), e.end()));
  }

  {
    SuperConfigEnum e(target_scel, sub_configs.begin(), sub_configs.end());
    std::map<Configuration, Array<int> > cmap;

    for(const auto &tconfig : e) {
      auto res = cmap.insert(std::make_pair(tconfig.primitive().in_canonical_supercell(), e.counter()));
      (void) res;
    }

    //for(const auto &val : cmap) {
    //  std::cout << "counter: " << val.second << "  occ: " << val.first.occupation() << std::endl;
    //}

    BOOST_CHECK_EQUAL(cmap.size(), 22);
  }

}

BOOST_AUTO_TEST_SUITE_END()
