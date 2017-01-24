#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ConfigEnumInterpolation.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ConfigEnumInterpolationTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  double tol = primclex.crystallography_tol();

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.get_prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(c, a - b, a + b - c)};

  Configuration c_begin(scel);
  c_begin.init_displacement();
  c_begin.init_deformation();

  // include displacements & deformation
  Eigen::Vector3d zero(0., 0., 0.);
  Eigen::Vector3d dx(0.001, 0., 0.);
  Eigen::Vector3d dy(0., 0.001, 0.);
  Eigen::Vector3d dz(0., 0., 0.001);
  Eigen::Matrix3d dF;
  dF << 0.01, 0.01, 0.0,
  0.0, -0.01, 0.0,
  0.0, 0.0, 0.05;
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

  Configuration c_final {c_begin};
  c_final.set_disp(0, dx);
  c_final.set_disp(1, dy);
  c_final.set_deformation(I + dF);

  ConfigEnumInterpolation e(c_begin, c_final, 11);

  Index i = 0;
  for(const auto &tconfig : e) {
    //std::cout << "\n--\n" << tconfig << std::endl;
    double f = (1.*i) / (e.size() - 1);
    BOOST_CHECK_EQUAL(almost_equal(tconfig.disp(0), dx * f, tol), true);
    BOOST_CHECK_EQUAL(almost_equal(tconfig.disp(1), dy * f, tol), true);
    BOOST_CHECK_EQUAL(almost_equal(tconfig.deformation(), I + dF * f, tol), true);
    ++i;
  }

}

BOOST_AUTO_TEST_SUITE_END()
