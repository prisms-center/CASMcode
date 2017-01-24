#ifndef NIGGLI_HH
#define NIGGLI_HH

#include "casm/external/Eigen/Core"

namespace CASM {
  class Lattice;
  class SymGroup;

  /** \ingroup Lattice
   *  @{
   */

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
   * \see
   * I. Krivy and B. Gruber, Acta Cryst. (1976). A32, 297.
   * <a href="http://dx.doi.org/10.1107/S0567739476000636">[doi:10.1107/S0567739476000636]</a>
   * R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, Acta Cryst. (2004). A60, 1.
   * <a href="http://dx.doi.org/10.1107/S010876730302186X"> [doi:10.1107/S010876730302186X]</a>
   *
   * <a href="https://www.phenix-online.org/papers/sh5006_reprint.pdf"> </a>
   */

  class NiggliRep {
  public:

    NiggliRep(const Lattice &init_lat);

    NiggliRep(const Eigen::Matrix3d &init_lat_col_mat);

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

    const Eigen::Matrix3d &metrical_matrix() const;

    ///A<=B OR (A==B, |ksi| <= |eta|)
    bool meets_criteria_1(double compare_tol) const;

    ///B<=C OR (B==C, |eta| <= |zeta|)
    bool meets_criteria_2(double compare_tol) const;

    ///For type I: ksi>0 && eta>0 && zeta>0 (all angles < 90)
    bool meets_criteria_3(double compare_tol) const;

    ///For type II: ksi<=0 && eta<=0 && zeta<=0 (all angles >= 90)
    bool meets_criteria_4(double compare_tol) const;

    ///|ksi|<=B OR (ksi==B, zeta<=2*eta) OR (ksi==-B, zeta==0)
    bool meets_criteria_5(double compare_tol) const;

    ///|eta|<=A OR (eta==A, zeta<=2*ksi) OR (eta==-A, zeta==0)
    bool meets_criteria_6(double compare_tol) const;

    ///|zeta|<=A OR (zeta==A, eta<=2*ksi) OR (zeta==-A, eta==0)
    bool meets_criteria_7(double compare_tol) const;

    ///ksi+eta+zeta+A+B<0 OR (ksi+eta+zeta+A+B==0, 2*A+2*eta+zeta>0)
    ///C<=A+B+C+ksi+eta+zeta OR (C==A+B+C+ksi+eta+zeta, 2*A+2*eta+zeta<=0)
    bool meets_criteria_8(double compare_tol) const;

    ///True if all conditions are true, and either 4 OR 3 is false
    bool is_niggli(double compare_tol) const;

    ///True if all conditions except 4 are true
    bool is_niggli_type1(double compare_tol) const;

    ///True if all conditions except 3 are true
    bool is_niggli_type2(double compare_tol) const;

    void debug_criteria(double compare_tol) const;

  private:

    ///Transpose of initialization lattice dotted with itself
    const Eigen::Matrix3d m_metrical_matrix;

  };


  ///Find the niggli, most standard oriented version of the given orbit (defined by the given SymGroup) of lattices
  Lattice canonical_equivalent_lattice(const Lattice &in_lat, const SymGroup &point_grp, double compare_tol);

  ///Convert the given lattice into it's niggli reduced form, with the most standard orientation possilbe
  Lattice niggli(const Lattice &in_lat, double compare_tol, bool keep_handedness = false);

  ///Check whether the given lattice (represented as a matrix) is in niggli TYPE ?? reduced form (does not check for orientation)
  bool is_niggli(const Eigen::Matrix3d &test_lat_mat, double compare_tol);

  ///Check whether the given lattice is primitive (does not check for orientation)
  bool is_niggli(const Lattice &test_lat, double compare_tol);

  /// \brief Generate a vector whose lexicographical value determines how well it's oriented in space
  Eigen::VectorXd spatial_unroll(const Eigen::Matrix3d &lat_mat, double compare_tol);

  /// \brief Compare the spatial orientation (ignoring matrix symmetry) and determine which one is oriented more standard. True if high is more standard.
  bool standard_orientation_spatial_compare(const Eigen::Matrix3d &low_score_lat_mat, Eigen::Matrix3d &high_score_lat_mat, double compare_tol);

  /// \brief Determine whether high_score has a more standard format than low_score
  bool standard_orientation_compare(const Eigen::Matrix3d &low_score_lat_mat, const Eigen::Matrix3d &high_score_lat_mat, double compare_tol);

  /** @} */
}

#endif
