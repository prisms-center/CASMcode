#include <iostream>
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/external/Eigen/Dense"
#include "casm/misc/CASM_math.hh"

using namespace CASM;

namespace Eigen {
  /**
   * Checks to see whether the given matrix is symmetric
   * by checking if its transpose is equal to itself.
   * Only works for square matrices.
   * (Reflected along 0,0 to n,n)
   */

  template <typename Derived>
  inline
  bool is_symmetric(const Eigen::MatrixBase<Derived> &test_mat, double test_tol = CASM::TOL) {
    return CASM::almost_zero(test_mat - test_mat.transpose(), test_tol);
  }

  /**
   * Checks to see if the given matrix is persymmetric, i.e.
   * whether it's symmetric along the cross diagonal.
   * Only works for square matrices.
   * (Reflected along 0,n to n,0)
   */

  template <typename Derived>
  inline
  bool is_persymmetric(const Eigen::MatrixBase<Derived> &test_mat, double test_tol = CASM::TOL) {
    //Reverse order of columns and rows
    auto rev_mat = test_mat.colwise().reverse().eval().rowwise().reverse().eval();
    return CASM::almost_zero(test_mat - rev_mat.transpose(), test_tol);
  }
}

namespace testing {

  std::vector<Eigen::Matrix3i> unimodular_matrices(bool positive, bool negative) {
    std::vector<Eigen::Matrix3i> uni_det_mats;
    int totalmats = 3480;

    if(positive && negative) {
      totalmats = totalmats * 2;
    }

    uni_det_mats.reserve(totalmats);

    EigenCounter<Eigen::Matrix3i> transmat_count(Eigen::Matrix3i::Constant(-1), Eigen::Matrix3i::Constant(1), Eigen::Matrix3i::Constant(1));

    for(; transmat_count.valid(); ++transmat_count) {
      if(positive && transmat_count.current().determinant() == 1) {
        uni_det_mats.push_back(transmat_count.current());
      }

      if(negative && transmat_count.current().determinant() == -1) {
        uni_det_mats.push_back(transmat_count.current());
      }
    }

    return uni_det_mats;
  }

  std::vector<Eigen::Matrix3i> positive_unimodular_matrices() {
    return unimodular_matrices(true, false);
  }

  std::vector<Eigen::Matrix3i> negative_unimodular_matrices() {
    return unimodular_matrices(false, true);
  }

  std::vector<Eigen::Matrix3i> unimodular_matrices() {
    return unimodular_matrices(true, true);
  }

  /**
   * Returns a values of A, B, C, ksi, eta and zeta.
   * If the initialization cell has vector lenghts a, b and c,
   * with angles alpha, beta, gamma, then the
   * niggli values are
   *
   * A=a*a
   * B=b*b
   * C=c*c
   * ksi=2bc*cos(alpha)
   * eta=2ac*cos(beta)
   * zeta=2ab*cos(gamma)
   *
   * Put another way, if the initialization lattice has
   * vectors a, b and c then
   *
   * A=aa
   * B=bb
   * C=cc
   * ksi=2*bc
   * eta=2*ac
   * zeta=2*ab
   *
   * <a href="http://dx.doi.org/10.1107/S0567739476000636">[doi:10.1107/S0567739476000636]</a>
   */

  class NiggliRep {
  public:

    NiggliRep(const Lattice &init_lat, double m_tolerance);

    ///Square of lattice length a
    double A() const;

    ///Square of lattice length b
    double B() const;

    ///Square of lattice length c
    double C() const;

    ///2bc*cos(alpha)
    double ksi() const;

    ///2ac*cos(beta)
    double eta() const;

    ///2ab*cos(gamma)
    double zeta() const;

    ///A>B OR (A==B, |ksi| > |eta|)
    bool meets_criteria_1() const;

    ///B>C OR (B==C, |eta| > |zeta|)
    bool meets_criteria_2() const;

    ///ksi*eta*ksi>0
    bool meets_criteria_3() const;

    ///ksi*eta*zeta<=0
    bool meets_criteria_4() const;

    ///|ksi|>B OR (ksi==B, 2*eta<zeta) OR (ksi==-B, zeta<0)
    bool meets_criteria_5() const;

    ///|eta|>A OR (eta==A, 2*ksi<zeta) OR (eta==-A, zeta<0)
    bool meets_criteria_6() const;

    ///|zeta|>A OR (zeta==A, 2*ksi<eta) OR (zeta==-A, eta<0)
    bool meets_criteria_7() const;

    ///ksi+eta+zeta+A+B<0 OR (ksi+eta+zeta+A+B==0, 2*A+2*eta+zeta>0)
    bool meets_criteria_8() const;

    ///True if all criteria return false, excluding criteria 3 and 4, which don't imply a non-Niggli form
    bool is_niggli() const;

  private:

    ///Transpose of initialization lattice dotted with itself
    const Eigen::Matrix3d m_self_dotted_lat;

    ///Tolerance value to be used for double comparisons
    const double m_tolerance;
  };

  NiggliRep::NiggliRep(const Lattice &init_lat, double init_tol):
    m_self_dotted_lat(init_lat.lat_column_mat().transpose() * init_lat.lat_column_mat()),
    m_tolerance(init_tol) {
  }

  double NiggliRep::A() const {
    return m_self_dotted_lat(0, 0);
  }

  double NiggliRep::B() const {
    return m_self_dotted_lat(1, 1);
  }

  double NiggliRep::C() const {
    return m_self_dotted_lat(2, 2);
  }

  double NiggliRep::ksi() const {
    return m_self_dotted_lat(2, 1);
  }

  double NiggliRep::eta() const {
    return m_self_dotted_lat(2, 0);
  }

  double NiggliRep::zeta() const {
    return m_self_dotted_lat(1, 0);
  }

  bool NiggliRep::meets_criteria_1() const {
    bool does_meet = false;

    if(CASM::compare(B(), A(), m_tolerance)) {
      does_meet = true;
    }

    else if(almost_equal(A(), B(), m_tolerance) && compare(std::abs(eta()), std::abs(ksi()), m_tolerance)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_2() const {
    bool does_meet = false;

    if(compare(C(), B(), m_tolerance)) {
      does_meet = true;
    }

    else if(almost_equal(B(), C(), m_tolerance) && compare(std::abs(zeta()), std::abs(eta()), m_tolerance)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_3() const {
    bool does_meet = false;

    if(compare(0.0, ksi()*eta()*ksi(), m_tolerance)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_4() const {
    bool does_meet = false;

    if(compare(ksi()*eta()*zeta(), 0.0, m_tolerance) || almost_zero(ksi()*eta()*zeta(), m_tolerance)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_5() const {
    bool does_meet = false;

    if(compare(B(), std::abs(ksi()), m_tolerance)) {
      does_meet = true;
    }

    else if(almost_equal(ksi(), B(), m_tolerance) && compare(2 * eta(), zeta(), m_tolerance)) {
      does_meet = true;
    }

    else if(almost_equal(ksi(), (-1)*B(), m_tolerance) && compare(zeta(), 0.0, m_tolerance)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_6() const {
    bool does_meet = false;

    if(compare(A(), std::abs(eta()), m_tolerance)) {
      does_meet = true;
    }

    else if(almost_equal(eta(), A(), m_tolerance) && compare(2 * ksi(), zeta(), m_tolerance)) {
      does_meet = true;
    }

    else if(almost_equal(eta(), (-1)*A(), m_tolerance) && compare(zeta(), 0.0, m_tolerance)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_7() const {
    bool does_meet = false;

    if(compare(A(), std::abs(zeta()), m_tolerance)) {
      does_meet = true;
    }

    else if(almost_equal(zeta(), A(), m_tolerance) && compare(2 * ksi(), eta(), m_tolerance)) {
      does_meet = true;
    }

    else if(almost_equal(zeta(), (-1)*A(), m_tolerance) && compare(eta(), 0.0, m_tolerance)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_8() const {
    bool does_meet = false;

    if(compare(ksi() + eta() + zeta() + A() + B(), 0.0, m_tolerance)) {
      does_meet = true;
    }

    else if(almost_zero(ksi() + eta() + zeta() + A() + B(), m_tolerance) && compare(0.0, 2 * (A() + eta()) + zeta(), m_tolerance)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::is_niggli() const {
    bool totally_niggli = true;

    if(meets_criteria_1() ||
       meets_criteria_2() ||
       //meets_criteria_3() ||
       //meets_criteria_4() ||
       meets_criteria_5() ||
       meets_criteria_6() ||
       meets_criteria_7() ||
       meets_criteria_8()) {
      totally_niggli = false;
    }

    return totally_niggli;
  }


  Lattice niggli(const Lattice &in_lat) {
    const Lattice reduced_in = in_lat.get_reduced_cell();

    Lattice best_lat = reduced_in;

    std::vector<Eigen::Matrix3i> canditate_trans_mats = positive_unimodular_matrices();

    for(auto it = canditate_trans_mats.begin(); it != canditate_trans_mats.end(); ++it) {
      Lattice candidate_lat(reduced_in.lat_column_mat() * (it->cast<double>()));
      candidate_lat.print(std::cout);
    }

    return best_lat;
  }

  void niggli_rep_test(const Lattice &in_lat) {
    int dims = 3;
    int minvol = 1;
    int maxvol = 10;

    SymGroup pg;
    in_lat.generate_point_group(pg);

    SupercellEnumerator<Lattice> latenumerator(in_lat, pg, minvol, maxvol + 1, dims);
    std::vector<Lattice> enumerated_lat(latenumerator.begin(), latenumerator.end());

    for(auto it = enumerated_lat.begin(); it != enumerated_lat.end(); ++it) {
      Lattice old_niggli = niggli(*it, pg, TOL);
      NiggliRep test_rep(old_niggli, TOL);
      NiggliRep raw_test_rep(*it, TOL);
      std::cout << "Your niggli test says " << test_rep.is_niggli() << " and " << raw_test_rep.is_niggli() << std::endl;
    }

    return;
  }

  void symmetric_testing() {
    Eigen::MatrixXd symmat(5, 5), persymmat(4, 4), nonsymmat;

    symmat << 1, 2, 3, 4, 5,
           2, 6, 7, 8, 9,
           3, 7, 10, 11, 12,
           4, 8, 11, 13, 14,
           5, 9, 12, 14, 15;

    std::cout << "symmetric matrix is symmetric? " << is_symmetric(symmat) << std::endl;
    std::cout << "symmetric matrix is symmetric? " << is_persymmetric(symmat) << std::endl;

    persymmat << 4, 3, 2, 1,
              7, 6, 5, 2,
              9, 8, 6, 3,
              10, 9, 7, 4;

    std::cout << "persymmetric matrix is symmetric? " << is_symmetric(persymmat) << std::endl;
    std::cout << "persymmetric matrix is symmetric? " << is_persymmetric(persymmat) << std::endl;

    return;
  }

}

int main() {
  Lattice testlat = Lattice::fcc();
  SymGroup pg;
  testlat.generate_point_group(pg);

  int dims = 1;
  int minvol = 1;
  int maxvol = 10;

  SupercellEnumerator<Lattice> latenumerator(testlat, pg, minvol, maxvol + 1, dims);
  std::vector<Lattice> enumerated_lat(latenumerator.begin(), latenumerator.end());

  std::cout << "Enumerated from " << minvol << " to " << maxvol << " and got " << enumerated_lat.size() << " lattices" << std::endl;


  int l = 1;
  for(auto it = enumerated_lat.begin(); it != enumerated_lat.end(); ++it) {
    Eigen::Matrix3i comp_transmat;
    comp_transmat << (l), 0, 0,
                  0, 1, 0,
                  0, 0, 1;

    Lattice comparelat = make_supercell(testlat, comp_transmat);

    Lattice nigglicompare = niggli(comparelat, pg, TOL);
    Lattice nigglitest = niggli(*it, pg, TOL);

    if(nigglicompare == nigglitest) {
      std::cout << "Lattice on iteration " << l << " matched." << std::endl;
    }
    else {
      std::cout << "Lattice on iteration " << l << " did NOT match." << std::endl;
    }

    testing::NiggliRep comparenig(comparelat, CASM::TOL);
    testing::NiggliRep nigglicomparenig(nigglicompare, CASM::TOL);
    testing::NiggliRep nigglitestnig(nigglitest, CASM::TOL);

    l++;
  }

  //testing::niggli(testlat);
  testing::niggli_rep_test(Lattice::fcc());
  testing::symmetric_testing();

  return 0;
}
