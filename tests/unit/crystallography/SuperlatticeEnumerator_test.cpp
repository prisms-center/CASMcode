#include "gtest/gtest.h"
#include "autotools.hh"

#include<boost/filesystem.hpp>

/// What is being tested:
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"

/// What is being used to test it:
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/external/Eigen/Dense"
#include "casm/crystallography/Niggli.hh"

using namespace CASM;
boost::filesystem::path testdir(autotools::abs_srcdir() + "/tests/unit/crystallography");

void autofail() {
  EXPECT_EQ(1, 0);
  return;
}


jsonParser mat_test_case(const std::string &pos_filename, int minvol, int maxvol) {

  const Structure test_struc(testdir / pos_filename);
  const Lattice test_lat = test_struc.lattice();
  const SymGroup effective_pg = test_struc.factor_group();

  std::vector<Eigen::Matrix3i> enumerated_mats;

  ScelEnumProps enum_props(minvol, maxvol + 1);
  SuperlatticeEnumerator test_enumerator(effective_pg.begin(), effective_pg.end(), test_lat, enum_props);

  double tol = TOL;
  for(auto it = test_enumerator.begin(); it != test_enumerator.end(); ++it) {
    enumerated_mats.push_back(it.matrix());

    // -- check niggli generation

    Lattice niggli1 = niggli(*it, tol);
    Lattice niggli2 = niggli(niggli1, tol);
    bool check_niggli = almost_equal(
                          niggli1.lat_column_mat(),
                          niggli2.lat_column_mat(),
                          tol);

    EXPECT_EQ(check_niggli, true);

    // -- check canonical generation

    Lattice canon = xtal::canonical::equivalent(*it, effective_pg, tol);
    Lattice canon2 = xtal::canonical::equivalent(canon, effective_pg, tol);
    bool check = almost_equal(
                   canon.lat_column_mat(),
                   canon2.lat_column_mat(),
                   tol);

    EXPECT_EQ(check, true);

  }

  jsonParser mat_dump;
  mat_dump["input"]["min_vol"] = minvol;
  mat_dump["input"]["max_vol"] = maxvol;
  mat_dump["input"]["source"] = pos_filename;
  mat_dump["output"]["mats"] = enumerated_mats;

  return mat_dump;
}

jsonParser lat_test_case(const std::string &pos_filename, int minvol, int maxvol) {
  const Structure test_struc(testdir / pos_filename);
  const Lattice test_lat = test_struc.lattice();
  const SymGroup effective_pg = test_struc.factor_group();

  std::vector<Lattice> enumerated_lats;
  ScelEnumProps enum_props(minvol, maxvol + 1);
  test_lat.generate_supercells(enumerated_lats, effective_pg, enum_props);

  jsonParser lat_dump;
  lat_dump["input"]["min_vol"] = minvol;
  lat_dump["input"]["max_vol"] = maxvol;
  lat_dump["input"]["source"] = pos_filename;
  lat_dump["output"]["lats"] = enumerated_lats;

  return lat_dump;
}

jsonParser generate_all_test_cases() {
  jsonParser all_test_cases;

  //********************************************************************//

  std::vector<jsonParser> all_mat_tests;
  all_mat_tests.push_back(mat_test_case("POS1.txt", 1, 6));
  all_mat_tests.push_back(mat_test_case("PRIM1.txt", 2, 9));
  all_mat_tests.push_back(mat_test_case("PRIM2.txt", 4, 7));
  all_mat_tests.push_back(mat_test_case("PRIM4.txt", 1, 8));
  all_test_cases["mat_test_cases"] = all_mat_tests;

  //********************************************************************//

  std::vector<jsonParser> all_lat_tests;
  all_lat_tests.push_back(lat_test_case("POS1.txt", 2, 6));
  all_lat_tests.push_back(lat_test_case("PRIM1.txt", 2, 9));
  all_lat_tests.push_back(lat_test_case("PRIM2.txt", 3, 7));
  all_lat_tests.push_back(lat_test_case("PRIM4.txt", 1, 8));
  all_lat_tests.push_back(lat_test_case("PRIM5.txt", 1, 8));

  all_test_cases["lat_test_cases"] = all_lat_tests;

  //********************************************************************//

  jsonParser test1 = all_test_cases;
  jsonParser test2 = all_test_cases;

  return all_test_cases;
}


void trans_enum_test() {
  Lattice testlat(Lattice::fcc());
  SymGroup pg = SymGroup::lattice_point_group(testlat);

  // int dims = 3;
  Eigen::Matrix3i transmat;

  transmat << -1, 1, 1,
           1, -1, 1,
           1, 1, -1;

  Lattice bigunit = make_supercell(testlat, transmat);

  ScelEnumProps enum_props(1, 5 + 1, "abc", transmat);
  SuperlatticeEnumerator enumerator(pg.begin(), pg.end(), testlat, enum_props);

  std::vector<Lattice> enumerated_lat(enumerator.begin(), enumerator.end());

  for(Index i = 0; i > enumerated_lat.size(); i++) {
    EXPECT_TRUE(enumerated_lat[i].is_supercell_of(bigunit));
  }

  return;
}


void restricted_test() {
  std::vector<Lattice> all_test_lats;
  all_test_lats.push_back(Lattice::fcc());
  all_test_lats.push_back(Lattice::bcc());
  all_test_lats.push_back(Lattice::cubic());
  all_test_lats.push_back(Lattice::hexagonal());

  for(Index t = 0; t < all_test_lats.size(); t++) {
    Lattice testlat = all_test_lats[t];
    SymGroup pg = SymGroup::lattice_point_group(testlat);

    // int dims = 1;

    ScelEnumProps enum_props(1, 15 + 1, "a");
    SuperlatticeEnumerator enumerator(pg.begin(), pg.end(), testlat, enum_props);

    int l = 1;
    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      Eigen::Matrix3i comp_transmat;
      comp_transmat << (l), 0, 0,
                    0, 1, 0,
                    0, 0, 1;

      EXPECT_TRUE(it.matrix() == canonical_hnf(pg.begin(), pg.end(), comp_transmat, testlat));
      l++;
    }
  }

  return;
}


//Tests in here were created by first getting results from
//before HermiteCounter existed and then making sure the results
//didn't change after it was introduced

//Unfortunately, the old niggli routines weren't working as
//they should have been, so these hard coded examples to check
//had to be regenerated...

TEST(SuperlatticeEnumeratorTest, EnumeratorConsistency) {

  boost::filesystem::path old_test_path = testdir / "test_cases.json";
  boost::filesystem::path current_test_path = testdir / "current_test_results.json";

  jsonParser current_test_results = generate_all_test_cases();
  current_test_results.write(current_test_path);

  //Find out where things fail
  //Comparison will fail if you don't compare from pre-written files.
  jsonParser curr(current_test_path);
  jsonParser existing(old_test_path);

  boost::filesystem::path failure_point = find_diff(curr, existing, TOL);

  if(!failure_point.empty()) {
    std::cout << "Difference at: " << failure_point << "\n" << std::endl;

    auto &_existing = existing.at(failure_point);
    auto &_curr = curr.at(failure_point);
    if(_existing.type() != _curr.type()) {
      std::cout << "Different types\n" << std::endl;
    }
    else if(_existing.is_array() && (_existing.size() != _curr.type())) {
      std::cout << "Different array sizes" << std::endl;
      std::cout << "  Expected: " << _existing.size() << std::endl;
      std::cout << "  Found: " << _curr.size() << std::endl;

    }
    else if(_existing.is_obj() && (_existing.size() != _curr.type())) {
      std::cout << "Different object sizes\n" << std::endl;
      std::cout << "  Expected: " << _existing.size() << std::endl;
      std::cout << "  Found: " << _curr.size() << std::endl;

    }
    else {
      std::cout << "Different values\n" << std::endl;
    }

    std::cout << "Expected: \n" << existing.at(failure_point) << "\n"
              << "Found: \n" << curr.at(failure_point) << std::endl;
  }


  EXPECT_TRUE(failure_point.empty());

  current_test_results["WARNING"] = "This has been added as an inconvenience to anyone who is thinking of replacing the \
current test_results.json file. Do not replace anything unless you're certain the old \
results were incorrect, and these are an improvement. If you are sure you want to proceed, eliminate this key.";

  current_test_results.write(current_test_path);

}

TEST(SuperlatticeEnumeratorTest, RestrictedEnumeration) {
  trans_enum_test();
  restricted_test();
}
