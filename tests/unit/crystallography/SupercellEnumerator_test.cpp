#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include<boost/filesystem.hpp>

/// What is being tested:
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/container/LinearAlgebra.hh"

/// What is being used to test it:
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/external/Eigen/Dense"
#include "casm/crystallography/Niggli.hh"

using namespace CASM;
boost::filesystem::path testdir("tests/unit/crystallography");

void autofail() {
  BOOST_CHECK_EQUAL(1, 0);
  return;
}

void hermite_init() {
  int dims = 5;
  int det = 30;

  HermiteCounter hermit_test(det, dims);

  Eigen::VectorXi init_diagonal(Eigen::VectorXi::Ones(dims));
  init_diagonal(0) = det;

  BOOST_CHECK_EQUAL(init_diagonal, hermit_test.diagonal());
  BOOST_CHECK_EQUAL(0, hermit_test.position());

  auto tricounter = HermiteCounter_impl::_upper_tri_counter(hermit_test.diagonal());
  Eigen::VectorXi startcount(Eigen::VectorXi::Zero(HermiteCounter_impl::upper_size(dims)));
  BOOST_CHECK_EQUAL(tricounter.current(), startcount);


  Eigen::VectorXi endcount(Eigen::VectorXi::Zero(HermiteCounter_impl::upper_size(dims)));
  endcount(0) = (det - 1);
  endcount(1) = (det - 1);
  endcount(2) = (det - 1);
  endcount(3) = (det - 1);

  auto finalcounterstate = tricounter;

  for(; tricounter.valid(); ++tricounter) {
    finalcounterstate = tricounter;
  }

  BOOST_CHECK_EQUAL(finalcounterstate.current(), endcount);


  return;
}

void spill_test() {
  Eigen::VectorXi diagonal0(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi diagonal1(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi diagonal2(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi diagonal3(Eigen::VectorXi::Ones(5));


  int p = 0;
  diagonal0(p) = 2;
  int p0 = HermiteCounter_impl::_spill_factor(diagonal0, p, 2);
  BOOST_CHECK_EQUAL(p0, p + 1);
  BOOST_CHECK_EQUAL(diagonal0(p), 1);
  BOOST_CHECK_EQUAL(diagonal0(p + 1), 2);

  p = 3;
  diagonal1(p) = 6;
  int p1 = HermiteCounter_impl::_spill_factor(diagonal1, p, 2);
  BOOST_CHECK_EQUAL(p1, p + 1);
  BOOST_CHECK_EQUAL(diagonal1(p), 3);
  BOOST_CHECK_EQUAL(diagonal1(p + 1), 2);

  p = 3;
  diagonal2(p) = 6;
  int p2 = HermiteCounter_impl::_spill_factor(diagonal2, p, 4);
  BOOST_CHECK_EQUAL(p2, p + 1);
  BOOST_CHECK_EQUAL(diagonal2(p), 1);
  BOOST_CHECK_EQUAL(diagonal2(p + 1), 6);

  p = 2;
  diagonal3(p) = 8;
  int p3 = HermiteCounter_impl::_spill_factor(diagonal3, p, 4);
  BOOST_CHECK_EQUAL(p3, p + 1);
  BOOST_CHECK_EQUAL(diagonal3(p), 2);
  BOOST_CHECK_EQUAL(diagonal3(p + 1), 4);

  return;
}

void next_position_test() {
  //Example increment from one possible diagonal to the next
  Eigen::VectorXi diagonal(Eigen::VectorXi::Ones(5));
  Eigen::VectorXi next_diagonal(Eigen::VectorXi::Ones(5));
  diagonal(0) = 6;
  next_diagonal(0) = 3;
  next_diagonal(1) = 2;
  int p = 0;

  p = HermiteCounter_impl::next_spill_position(diagonal, p);

  BOOST_CHECK_EQUAL(diagonal, next_diagonal);
  BOOST_CHECK_EQUAL(p, 1);


  diagonal = Eigen::VectorXi::Ones(5);
  next_diagonal = Eigen::VectorXi::Ones(5);
  //[1 2 1 1 3]
  diagonal(1) = 2;
  diagonal(4) = 3;
  //[1 1 6 1 1]
  next_diagonal(2) = 6;

  p = 4;
  p = HermiteCounter_impl::next_spill_position(diagonal, p);

  BOOST_CHECK_EQUAL(diagonal, next_diagonal);
  BOOST_CHECK_EQUAL(p, 2);

  //*************/
  //Make sure every enumerated diagonal has the right determinant

  int det = 2 * 3 * 5 * 7;
  int dims = 5;

  Eigen::VectorXi diag = Eigen::VectorXi::Ones(dims);
  diag(0) = det;

  p = 0;
  while(p != diag.size()) {
    int testdet = 1;
    for(int i = 0; i < diag.size(); i++) {
      testdet = testdet * diag(i);
    }
    BOOST_CHECK_EQUAL(det, testdet);
    p = CASM::HermiteCounter_impl::next_spill_position(diag, p);
  }

  return;
}

void triangle_count_test() {
  HermiteCounter::Index totals = HermiteCounter_impl::upper_size(7);
  BOOST_CHECK_EQUAL(totals, -7 + 7 + 6 + 5 + 4 + 3 + 2 + 1);

  int dims = 5;
  // int det = 30;

  Eigen::VectorXi mid_diagonal(Eigen::VectorXi::Ones(dims));
  mid_diagonal(0) = 5;
  mid_diagonal(1) = 3;
  mid_diagonal(4) = 2;

  auto countertest = HermiteCounter_impl::_upper_tri_counter(mid_diagonal);
  auto finalcount = countertest;

  for(; countertest.valid(); countertest++) {
    finalcount = countertest;
  }

  //The initial matrix is 5x5 with diagonal [ 5 3 1 1 2 ], so it has determinant=30
  //The Hermite matrix with highest ranking for this determinant will therefore be:
  //    5 4 4 4 4
  //    0 3 2 2 2
  //    0 0 1 0 0
  //    0 0 0 1 0
  //    0 0 0 0 2
  //Which gives the upper triangular vector [ 4 4 4 4 2 2 2 0 0 0 ]

  Eigen::VectorXi end_count_value(Eigen::VectorXi::Zero(HermiteCounter_impl::upper_size(5)));
  end_count_value(0) = 4;
  end_count_value(1) = 4;
  end_count_value(2) = 4;
  end_count_value(3) = 4;
  end_count_value(4) = 2;
  end_count_value(5) = 2;
  end_count_value(6) = 2;
  //Rest of the values are zero

  BOOST_CHECK_EQUAL(finalcount.current(), end_count_value);

  return;
}

void matrix_construction_test() {
  Eigen::VectorXi diag;
  diag.resize(4);
  diag << 2, 4, 6, 8;
  Eigen::VectorXi upper;
  upper.resize(3 + 2 + 1);
  upper << 11, 12, 13,
        21, 22,
        33;

  Eigen::MatrixXi diagmat;
  diagmat.resize(4, 4);
  diagmat << 2, 11, 12, 13,
          0, 4, 21, 22,
          0, 0, 6, 33,
          0, 0, 0, 8;

  BOOST_CHECK_EQUAL(diagmat, HermiteCounter_impl::_zip_matrix(diag, upper));

  return;
}

void increment_test() {
  HermiteCounter hermit_test(6, 4);

  Eigen::MatrixXi hermmat;
  hermmat.resize(4, 4);

  //Test starting status
  hermmat << 6, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

  BOOST_CHECK_EQUAL(hermmat, hermit_test());

  //Test next status
  ++hermit_test;
  hermmat << 6, 1, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

  BOOST_CHECK_EQUAL(hermmat, hermit_test());

  //Jump to just before you need a new diagonal
  hermmat << 6, 5, 5, 5,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

  while(hermit_test() != hermmat) {
    ++hermit_test;
  }

  //Check diagonal jump
  ++hermit_test;
  hermmat << 3, 0, 0, 0,
          0, 2, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

  BOOST_CHECK_EQUAL(hermmat, hermit_test());

  //Check invalidation and last status
  auto lastherm = hermmat;
  while(hermit_test.determinant() != 7) {
    lastherm = hermit_test();
    ++hermit_test;
  }

  hermmat << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 6;

  BOOST_CHECK_EQUAL(hermmat, lastherm);

  //Check determinant jump
  hermit_test = HermiteCounter(3, 4);

  //Jump to just before you need a new determinant

  hermmat << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 3;

  while(hermit_test() != hermmat) {
    ++hermit_test;
  }

  //Check determinant jump
  ++hermit_test;

  hermmat << 4, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

  BOOST_CHECK_EQUAL(hermmat, hermit_test());
  return;
}

void reset_test() {
  HermiteCounter hermit_test(1, 3);

  Eigen::MatrixXi hermmat;
  hermmat.resize(3, 3);


  Eigen::MatrixXi startmat = hermit_test();

  //Skip to one of the bigger determinants
  hermmat << 2, 1, 1,
          0, 2, 1,
          0, 0, 1;

  while(hermit_test() != hermmat) {
    ++hermit_test;
  }

  hermmat << 4, 0, 0,
          0, 1, 0,
          0, 0, 1;

  hermit_test.reset_current();

  BOOST_CHECK_EQUAL(hermmat, hermit_test());

  hermit_test.jump_to_determinant(1);

  BOOST_CHECK_EQUAL(startmat, hermit_test());

  return;
}

void expand_dims_test() {
  Eigen::MatrixXi expandmat(Eigen::MatrixXi::Ones(5, 5));
  expandmat = expandmat * 3;

  Eigen::VectorXi expanddims(8);
  expanddims << 1, 1, 1, 0, 1, 0, 0, 1;

  Eigen::MatrixXi expandedmat(8, 8);
  expandmat << 3, 3, 3, 0, 3, 0, 0, 3,
            3, 3, 3, 0, 3, 0, 0, 3,
            3, 3, 3, 0, 3, 0, 0, 3,
            0, 0, 0, 1, 0, 0, 0, 0,
            3, 3, 3, 0, 3, 0, 0, 3,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            3, 3, 3, 0, 3, 0, 0, 1;

  BOOST_CHECK_EQUAL(expandmat, HermiteCounter_impl::_expand_dims_old(expandmat, expanddims));

  HermiteCounter minicount(1, 4);
  for(int i = 0; i < 12; i++) {
    ++minicount;
  }

  Eigen::MatrixXi endcount(4, 4);
  endcount << 1, 0, 0, 0,
           0, 2, 1, 1,
           0, 0, 1, 0,
           0, 0, 0, 1;

  BOOST_CHECK_EQUAL(endcount, minicount());

  Eigen::MatrixXi transmat(Eigen::MatrixXi::Identity(6, 6));

  Eigen::MatrixXi expanded = HermiteCounter_impl::_expand_dims(minicount(), transmat);
  Eigen::MatrixXi blockmat(6, 6);
  blockmat << 1, 0, 0, 0, 0, 0,
           0, 2, 1, 1, 0, 0,
           0, 0, 1, 0, 0, 0,
           0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 1, 0,
           0, 0, 0, 0, 0, 1;

  BOOST_CHECK_EQUAL(blockmat, expanded);

  Eigen::Matrix2Xi miniherm;
  miniherm << 2, 1,
           0, 3;

  Eigen::Matrix3i minitrans;
  minitrans << 1, 0, 0,
            0, 0, 1,
            0, 1, 0;

  Eigen::Matrix3i miniexpand;
  miniexpand << 2, 1, 0,
             0, 0, 1,
             0, 3, 0;

  BOOST_CHECK_EQUAL(HermiteCounter_impl::_expand_dims(miniherm, minitrans), miniexpand);


  return;
}

jsonParser mat_test_case(const std::string &pos_filename, int minvol, int maxvol) {

  const Structure test_struc(testdir / pos_filename);
  const Lattice test_lat = test_struc.lattice();
  const SymGroup effective_pg = test_struc.factor_group();

  Array<Eigen::Matrix3i> enumerated_mats;

  ScelEnumProps enum_props(minvol, maxvol + 1);
  SupercellEnumerator<Lattice> test_enumerator(test_lat, effective_pg, enum_props);

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

    BOOST_CHECK_EQUAL(check_niggli, true);

    // -- check canonical generation

    Lattice canon = canonical_equivalent_lattice(*it, effective_pg, tol);
    Lattice canon2 = canonical_equivalent_lattice(canon, effective_pg, tol);
    bool check = almost_equal(
                   canon.lat_column_mat(),
                   canon2.lat_column_mat(),
                   tol);

    BOOST_CHECK_EQUAL(check, true);

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

  Array<Lattice> enumerated_lats;
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
  all_mat_tests.push_back(mat_test_case("POS1", 1, 6));
  all_mat_tests.push_back(mat_test_case("PRIM1", 2, 9));
  all_mat_tests.push_back(mat_test_case("PRIM2", 4, 7));
  all_mat_tests.push_back(mat_test_case("PRIM4", 1, 8));

  all_test_cases["mat_test_cases"] = all_mat_tests;

  //********************************************************************//

  std::vector<jsonParser> all_lat_tests;
  all_lat_tests.push_back(lat_test_case("POS1", 2, 6));
  all_lat_tests.push_back(lat_test_case("PRIM1", 2, 9));
  all_lat_tests.push_back(lat_test_case("PRIM2", 3, 7));
  all_lat_tests.push_back(lat_test_case("PRIM4", 1, 8));
  all_lat_tests.push_back(lat_test_case("PRIM5", 1, 8));

  all_test_cases["lat_test_cases"] = all_lat_tests;

  //********************************************************************//

  jsonParser test1 = all_test_cases;
  jsonParser test2 = all_test_cases;

  return all_test_cases;
}

void unroll_test() {
  Eigen::MatrixXi mat5(5, 5);
  mat5 << 1, 12, 11, 10, 9,
       0, 2, 13, 15, 8,
       0, 0, 3, 14, 7,
       0, 0, 0, 4, 6,
       0, 0, 0, 0, 5;

  Eigen::VectorXi vec5(5 + 4 + 3 + 2 + 1);
  vec5 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  BOOST_CHECK_EQUAL(vec5, HermiteCounter_impl::_canonical_unroll(mat5));

  return;

  Eigen::Matrix3i mat3;
  mat3 << 1, 6, 5,
       0, 2, 4,
       0, 0, 3;

  Eigen::Vector3i vec3(3 + 2 + 1);
  vec3 << 1, 2, 3, 4, 5, 6;

  BOOST_CHECK_EQUAL(vec3, HermiteCounter_impl::_canonical_unroll(mat3));

  return;
}

void compare_test() {
  Eigen::Matrix3i low, high;

  low << 1, 9, 9,
      0, 9, 9,
      0, 9, 9;

  high << 2, 0, 0,
       0, 1, 0,
       0, 0, 1;

  BOOST_CHECK(HermiteCounter_impl::_canonical_compare(low, high));

  low << 1, 9, 9,
      0, 9, 9,
      0, 9, 9;

  high << 1, 10, 9,
       0, 9, 9,
       0, 9, 9;

  BOOST_CHECK(HermiteCounter_impl::_canonical_compare(low, high));

  return;
}

void trans_enum_test() {
  Lattice testlat(Lattice::fcc());
  SymGroup pg;
  testlat.generate_point_group(pg);
  // int dims = 3;
  Eigen::Matrix3i transmat;

  transmat << -1, 1, 1,
           1, -1, 1,
           1, 1, -1;

  Lattice bigunit = make_supercell(testlat, transmat);

  ScelEnumProps enum_props(1, 5 + 1, "abc", transmat);
  SupercellEnumerator<Lattice> enumerator(testlat, pg, enum_props);

  std::vector<Lattice> enumerated_lat(enumerator.begin(), enumerator.end());

  for(Index i = 0; i > enumerated_lat.size(); i++) {
    BOOST_CHECK(enumerated_lat[i].is_supercell_of(bigunit));
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
    SymGroup pg;
    testlat.generate_point_group(pg);
    // int dims = 1;

    ScelEnumProps enum_props(1, 15 + 1, "a");
    SupercellEnumerator<Lattice> enumerator(testlat, pg, enum_props);

    int l = 1;
    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      Eigen::Matrix3i comp_transmat;
      comp_transmat << (l), 0, 0,
                    0, 1, 0,
                    0, 0, 1;

      BOOST_CHECK(it.matrix() == canonical_hnf(comp_transmat, pg, testlat));
      l++;
    }
  }

  return;
}


BOOST_AUTO_TEST_SUITE(SupercellEnumeratorTest)

BOOST_AUTO_TEST_CASE(HermiteConstruction) {
  hermite_init();
}

BOOST_AUTO_TEST_CASE(HermiteImpl) {
  spill_test();
  next_position_test();
  triangle_count_test();
  matrix_construction_test();
  reset_test();
  unroll_test();
  compare_test();
}

BOOST_AUTO_TEST_CASE(HermiteCounting) {
  increment_test();
}

//Tests in here were created by first getting results from
//before HermiteCounter existed and then making sure the results
//didn't change after it was introduced

//Unfortunately, the old niggli routines weren't working as
//they should have been, so these hard coded examples to check
//had to be regenerated...

BOOST_AUTO_TEST_CASE(EnumeratorConsistency) {

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


  BOOST_CHECK(failure_point.empty());

  current_test_results["WARNING"] = "This has been added as an inconvenience to anyone who is thinking of replacing the \
current test_results.json file. Do not replace anything unless you're certain the old \
results were incorrect, and these are an improvement. If you are sure you want to proceed, eliminate this key.";

  current_test_results.write(current_test_path);

}

BOOST_AUTO_TEST_CASE(RestrictedEnumeration) {
  trans_enum_test();
  restricted_test();
}

BOOST_AUTO_TEST_SUITE_END()
