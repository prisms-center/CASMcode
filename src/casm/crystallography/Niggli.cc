#include "casm/crystallography/Niggli.hh"

#include <set>

#include "casm/crystallography/Lattice.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace xtal {

std::vector<Eigen::Matrix3d> _cell_invariant_transforms() {
  std::vector<Eigen::Matrix3d> result(24, Eigen::Matrix3d::Zero());
  Index i = 0;
  result[i].setIdentity();

  ++i;
  result[i](0, 1) = result[i](1, 2) = result[i](2, 0) = 1;

  ++i;
  result[i](0, 2) = result[i](1, 0) = result[i](2, 1) = 1;

  ++i;
  result[i](0, 1) = result[i](1, 0) = result[i](2, 2) = -1;

  ++i;
  result[i](0, 2) = result[i](1, 1) = result[i](2, 0) = -1;

  ++i;
  result[i](0, 0) = result[i](1, 2) = result[i](2, 1) = -1;

  Eigen::Matrix3d tmat(-Eigen::Matrix3d::Identity());

  for (Index j = 0; j < 3; ++j) {
    tmat(j, j) = 1;
    for (Index k = 0; k < 6; ++k) {
      result[++i] = result[k] * tmat;
    }
    tmat(j, j) = -1;
  }

  return result;
}

std::vector<Eigen::Matrix3d> const &NiggliRep::cell_invariant_transforms() {
  static std::vector<Eigen::Matrix3d> result = _cell_invariant_transforms();
  return result;
}

NiggliRep::NiggliRep(const Eigen::Matrix3d &init_lat_col_mat)
    : m_metrical_matrix(init_lat_col_mat.transpose() * init_lat_col_mat),
      m_scale_factor(cuberoot(init_lat_col_mat.determinant())) {}

NiggliRep::NiggliRep(const Lattice &init_lat)
    : NiggliRep(init_lat.lat_column_mat()) {}

double NiggliRep::A() const { return m_metrical_matrix(0, 0); }

double NiggliRep::B() const { return m_metrical_matrix(1, 1); }

double NiggliRep::C() const { return m_metrical_matrix(2, 2); }

double NiggliRep::ksi() const { return m_metrical_matrix(2, 1) * 2; }

double NiggliRep::eta() const { return m_metrical_matrix(2, 0) * 2; }

double NiggliRep::zeta() const { return m_metrical_matrix(1, 0) * 2; }

const Eigen::Matrix3d &NiggliRep::metrical_matrix() const {
  return m_metrical_matrix;
}

bool NiggliRep::meets_criteria_1(double compare_tol) const {
  bool does_meet = true;
  compare_tol = m_scale_factor * compare_tol;
  // For niggli you need A<=B
  if (CASM::compare(B(), A(), compare_tol)) {
    does_meet = false;
  }

  // If A==B, then niggli requires |ksi|<=|eta|
  else if (almost_equal(A(), B(), compare_tol) &&
           compare(std::abs(eta()), std::abs(ksi()), compare_tol)) {
    does_meet = false;
  }

  return does_meet;
}

bool NiggliRep::meets_criteria_2(double compare_tol) const {
  bool does_meet = true;

  compare_tol = m_scale_factor * compare_tol;
  // For niggli you need B<=C
  if (compare(C(), B(), compare_tol)) {
    does_meet = false;
  }

  // If B==C, then niggli requires |eta|<=|zeta|
  else if (almost_equal(B(), C(), compare_tol) &&
           compare(std::abs(zeta()), std::abs(eta()), compare_tol)) {
    does_meet = false;
  }

  return does_meet;
}

bool NiggliRep::meets_criteria_3(double compare_tol) const {
  compare_tol = m_scale_factor * compare_tol;
  bool does_meet = false;
  // For type one niggli cells, the angles are all less than 90 degrees
  if (compare(0.0, ksi(), compare_tol) && compare(0.0, eta(), compare_tol) &&
      compare(0.0, zeta(), compare_tol)) {
    does_meet = true;
  }

  return does_meet;
}

bool NiggliRep::meets_criteria_4(double compare_tol) const {
  compare_tol = m_scale_factor * compare_tol;
  bool does_meet = false;
  // For type two niggli cells, the angles are all more than or equal to 90
  // degrees
  if (!compare(0.0, ksi(), compare_tol) && !compare(0.0, eta(), compare_tol) &&
      !compare(0.0, zeta(), compare_tol)) {
    does_meet = true;
  }

  return does_meet;
}

bool NiggliRep::meets_criteria_5(double compare_tol) const {
  compare_tol = m_scale_factor * compare_tol;
  bool does_meet = true;

  // For niggli you need |ksi|<=B
  if (compare(B(), std::abs(ksi()), compare_tol)) {
    does_meet = false;
  }

  // If ksi==B, then niggli requires zeta<=2*eta
  else if (almost_equal(ksi(), B(), compare_tol) &&
           compare(2 * eta(), zeta(), compare_tol)) {
  }

  // If ksi==-B, then niggli requires zeta==0
  else if (almost_equal(ksi(), (-1) * B(), compare_tol) &&
           !almost_zero(zeta(), compare_tol)) {
    does_meet = false;
  }

  return does_meet;
}

bool NiggliRep::meets_criteria_6(double compare_tol) const {
  compare_tol = m_scale_factor * compare_tol;
  bool does_meet = true;

  // For niggli you need |eta|<=A
  if (compare(A(), std::abs(eta()), compare_tol)) {
    does_meet = false;
  }

  // If ksi==A, then niggli requires zeta<=2*ksi
  else if (almost_equal(eta(), A(), compare_tol) &&
           compare(2 * ksi(), zeta(), compare_tol)) {
    does_meet = false;
  }

  // If ksi==-A, then niggli requires zeta==0
  else if (almost_equal(eta(), (-1) * A(), compare_tol) &&
           !almost_zero(zeta(), compare_tol)) {
    does_meet = false;
  }

  return does_meet;
}

bool NiggliRep::meets_criteria_7(double compare_tol) const {
  compare_tol = m_scale_factor * compare_tol;
  bool does_meet = true;

  // For niggli you need |zeta|<=A
  if (compare(A(), std::abs(zeta()), compare_tol)) {
    does_meet = false;
  }

  // If zeta==A, then you need eta<=2*ksi
  else if (almost_equal(zeta(), A(), compare_tol) &&
           compare(2 * ksi(), eta(), compare_tol)) {
    does_meet = false;
  }

  // If zeta==-A, then you need eta==0
  else if (almost_equal(zeta(), (-1) * A(), compare_tol) &&
           !almost_zero(eta(), compare_tol)) {
    does_meet = false;
  }

  return does_meet;
}

bool NiggliRep::meets_criteria_8(double compare_tol) const {
  compare_tol = m_scale_factor * compare_tol;
  bool does_meet = true;

  double clobber = A() + B() + C() + ksi() + eta() + zeta();

  // For niggli you need C<=A+B+C+ksi+eta+zeta
  if (compare(clobber, C(), compare_tol)) {
    does_meet = false;
  }

  // If C==A+B+C+ksi+eta+zeta, then 2*A+2*eta+zeta<=0
  else if (almost_equal(clobber, C(), compare_tol) &&
           compare(0.0, 2 * A() + 2 * eta() + zeta(), compare_tol)) {
    does_meet = false;
  }

  return does_meet;
}

bool NiggliRep::is_niggli_type1(double compare_tol) const {
  return (is_niggli(compare_tol) && meets_criteria_3(compare_tol));
}

bool NiggliRep::is_niggli_type2(double compare_tol) const {
  return (is_niggli(compare_tol) && meets_criteria_4(compare_tol));
}

bool NiggliRep::is_niggli(double compare_tol) const {
  bool totally_niggli = false;
  if (meets_criteria_1(compare_tol) && meets_criteria_2(compare_tol) &&
      (meets_criteria_3(compare_tol) != meets_criteria_4(compare_tol)) &&
      meets_criteria_5(compare_tol) && meets_criteria_6(compare_tol) &&
      meets_criteria_7(compare_tol) && meets_criteria_8(compare_tol)) {
    totally_niggli = true;
  }

  return totally_niggli;
}

void NiggliRep::debug_criteria(double compare_tol) const {
  std::cout << meets_criteria_1(compare_tol) << meets_criteria_2(compare_tol)
            << meets_criteria_3(compare_tol) << meets_criteria_4(compare_tol)
            << meets_criteria_5(compare_tol) << meets_criteria_6(compare_tol)
            << meets_criteria_7(compare_tol) << meets_criteria_8(compare_tol)
            << std::endl;

  std::cout << A() << " " << B() << " " << C() << " " << ksi() << " " << eta()
            << " " << zeta() << std::endl;

  return;
}

Index NiggliRep::niggli_index(double compare_tol) const {
  return Index(meets_criteria_1(compare_tol)) +
         Index(meets_criteria_2(compare_tol)) +
         Index(meets_criteria_3(compare_tol) || meets_criteria_4(compare_tol)) +
         Index(meets_criteria_5(compare_tol)) +
         Index(meets_criteria_6(compare_tol)) +
         Index(meets_criteria_7(compare_tol)) +
         Index(meets_criteria_8(compare_tol));
}

Index niggli_index(const Eigen::Matrix3d &test_lat_mat, double compare_tol) {
  return NiggliRep(test_lat_mat).niggli_index(compare_tol);
}

bool is_niggli(const Eigen::Matrix3d &test_lat_mat, double compare_tol) {
  NiggliRep ntest(test_lat_mat);

  // ntest.debug_criteria(compare_tol);
  return ntest.is_niggli(compare_tol);
}

bool is_niggli(const Lattice &test_lat, double compare_tol) {
  return is_niggli(test_lat.lat_column_mat(), compare_tol);
}

// Helper function, declared only in this file. Returns all niggli cells of
// 'in_lat'
std::set<Eigen::Matrix3d, StandardOrientationCompare> _niggli_set(
    const Lattice &in_lat, double compare_tol, bool keep_handedness) {
  Lattice reduced_in = in_lat.reduced_cell();

  if (!keep_handedness) {
    reduced_in.make_right_handed();
  }

  std::set<Eigen::Matrix3d, StandardOrientationCompare> result(
      (StandardOrientationCompare(compare_tol)));

  // The parallelepipiped of the reduced cell already has the shape of the
  // finally niggli cell, but we must select which edges will form the lattice
  // vectors (satisfying the Niggli property) and which orientation of the
  // lattice vectors to use (satisfying the CASM-canonical property) For now, we
  // optimize both properties simultaneously by looping over all transformation
  // matrices with determinant 1 and elements -1, 0 or 1. There are 3480 such
  // matrices We could vastly improve this to maximum 24 + 24 transformations if
  // we knew the (CRYSTAL) point-group There are 24 possible lattice vector
  // choices to screen for Niggli, starting from the reduced cell Once the
  // Niggli condition is achieved, there are maximum 24 proper point group
  // rotations to consider to screen for canonical orientation
  const std::vector<Eigen::Matrix3i> &candidate_trans_mats =
      positive_unimodular_matrices();
  Eigen::Matrix3d candidate_lat_mat;
  for (auto it = candidate_trans_mats.begin(); it != candidate_trans_mats.end();
       ++it) {
    candidate_lat_mat = reduced_in.lat_column_mat() * it->cast<double>();
    if (is_niggli(candidate_lat_mat, compare_tol)) {
      result.insert(candidate_lat_mat);
    }
  }
  return result;
}

/**
 * Converts a given Lattice into a new Lattice that satisfies
 * the Niggli criteria. Because there are several representations
 * of the Lattice that satisfy these conditions, the one with
 * the most standard orientation will be returned.
 *
 * The returned Niggli cell will be right handed regardless
 * of the inputted lattice, unless keep_handedness is set
 * to true.
 */

Lattice niggli(const Lattice &in_lat, double compare_tol,
               bool keep_handedness) {
  auto lat_set = _niggli_set(in_lat, compare_tol, keep_handedness);
  if (lat_set.empty()) {
    throw std::runtime_error(
        "In niggli(), did not find any lattices matching niggli criteria!");
  }
  return Lattice(*lat_set.rbegin(), in_lat.tol());
}

/**
 * Generates a vector of values from a Lattice column matrix, whose
 * lexicographical value determines how standard it's orientation is.
 *
 * First entry is the number of non-negative entries in the matrix.
 * Next entries are the diagonal values.
 * Next entries are the negative of the absolute value of the off diagonal
 * entries. Finally, the signs of the off diagonal entries.
 *
 * Maximizing the lexicographical value will favor matrices that
 * are diagonal if possible, followed by a preference for upper triangular
 * matrices. Off diagonal entries are minimized.
 */

Eigen::VectorXd spatial_unroll(const Eigen::Matrix3d &lat_mat,
                               double compare_tol) {
  // We want to give a preference to lattices with more positive values than
  // negative ones Count how many non-negative entries there are
  int non_negatives = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (float_sgn(lat_mat(i, j), compare_tol) != -1) {
        non_negatives++;
      }
    }
  }

  Eigen::VectorXd lat_spatial_descriptor(16);

  lat_spatial_descriptor <<
      // Favor positive values
      non_negatives,

      // Favor large diagonal values
      lat_mat(0, 0), lat_mat(1, 1), lat_mat(2, 2),

      // Favor small off diagonal values
      -std::abs(lat_mat(2, 1)), -std::abs(lat_mat(2, 0)),
      -std::abs(lat_mat(1, 0)), -std::abs(lat_mat(1, 2)),
      -std::abs(lat_mat(0, 2)), -std::abs(lat_mat(0, 1)),

      // Favor upper triangular
      float_sgn(lat_mat(2, 1), compare_tol),
      float_sgn(lat_mat(2, 0), compare_tol),
      float_sgn(lat_mat(1, 0), compare_tol),
      float_sgn(lat_mat(1, 2), compare_tol),
      float_sgn(lat_mat(0, 2), compare_tol),
      float_sgn(lat_mat(0, 1), compare_tol);

  return lat_spatial_descriptor;
}

/**
 * Given two lattice column matrices, determine which one of them is oriented
 * in a more standard manner in space relative to the Cartesian coordinates.
 *
 * The routine returns "low_score_lat_mat<high_score_lat_mat", which is true if
 * high_score_lat_mat has a more standard spatial orientation.
 */

bool standard_orientation_spatial_compare(
    const Eigen::Matrix3d &low_score_lat_mat,
    const Eigen::Matrix3d &high_score_lat_mat, double compare_tol) {
  Eigen::VectorXd low_score_lat_unroll =
      spatial_unroll(low_score_lat_mat, compare_tol);
  Eigen::VectorXd high_score_lat_unroll =
      spatial_unroll(high_score_lat_mat, compare_tol);

  return float_lexicographical_compare(low_score_lat_unroll,
                                       high_score_lat_unroll, compare_tol);
}

/**
 * Compares two lattice column matrices to each other and determines which
 * one has a more standard orientation based on whether said matrices are
 * (bi)symmetric and how aligned the vectors are along the Cartesian directions.
 *
 * A bisymmetric matrix will always have better standard orientation than one
 * that is simply symmetric. In turn, a symmetric matrix will always have a
 * better standard orientation than on that is not. The criteria with lowest
 * impact is the orientation of the lattice vectors in space, as defined in
 * ::standard_orientation_spatial_compare
 *
 * The routine returns "low_score_mat<high_score_mat", which is true if
 * high_score_mat has a more standard orientation.
 */

bool standard_orientation_compare(const Eigen::Matrix3d &low_score_lat_mat,
                                  const Eigen::Matrix3d &high_score_lat_mat,
                                  double compare_tol) {
  bool low_score_is_symmetric =
      Eigen::is_symmetric(low_score_lat_mat, compare_tol);
  bool low_score_is_bisymmetric =
      Eigen::is_bisymmetric(low_score_lat_mat, compare_tol);

  bool high_score_is_symmetric =
      Eigen::is_symmetric(high_score_lat_mat, compare_tol);
  bool high_score_is_bisymmetric =
      Eigen::is_bisymmetric(high_score_lat_mat, compare_tol);

  // The high_score is less symmetric than the low score, therefore high<low
  if ((low_score_is_bisymmetric && !high_score_is_bisymmetric) ||
      (low_score_is_symmetric && !high_score_is_symmetric)) {
    return false;
  }

  // The high_score is more symmetric than the low_score, therefore high>low
  // automatically (matrix symmetry trumps spatial orientation)
  else if ((!low_score_is_bisymmetric && high_score_is_bisymmetric) ||
           (!low_score_is_symmetric && high_score_is_symmetric)) {
    return true;
  }

  // If you made it here then the high_score and the low_score have the same
  // symmetry level, so we check for the spatial orientation
  else {
    return standard_orientation_spatial_compare(
        low_score_lat_mat, high_score_lat_mat, compare_tol);
  }
}

}  // namespace xtal
}  // namespace CASM
