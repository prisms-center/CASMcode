#ifndef NIGGLI_HH
#define NIGGLI_HH

#include "casm/external/Eigen/Core"

namespace CASM {
  class Lattice;
  class SymGroup;

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
   */

  class NiggliRep {
  public:

    NiggliRep(const Lattice &init_lat, double m_tolerance);

    NiggliRep(const Eigen::Matrix3d &init_lat_col_mat, double m_tolerance);

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


  ///Find the niggli, most standard oriented version of the given orbit (defined by the given SymGroup) of lattices
  Lattice canonical_equivalent_lattice(const Lattice &in_lat, const SymGroup &point_grp, double compare_tol);

  ///Convert the given lattice into it's niggli reduced form, with the most standard orientation possilbe
  Lattice niggli(const Lattice &in_lat, double compare_tol);

  ///Check whether the given lattice (represented as a matrix) is in niggli reduced form (does not check for orientation)
  bool is_niggli(const Eigen::Matrix3d &test_lat_mat, double compare_tol);

  ///Check whether the given lattice is primitive (does not check for orientation)
  bool is_niggli(const Lattice &test_lat, double compare_tol);

  /// \brief Compare the spatial orientation (ignoring matrix symmetry) and determine which one is oriented more standard. True if high is more standard.
  bool standard_orientation_spatial_compare(const Eigen::Matrix3d &low_score_lat_mat, Eigen::Matrix3d &high_score_lat_mat, double compare_tol);

  /// \brief Determine whether high_score has a more standard format than low_score
  bool standard_orientation_compare(const Eigen::Matrix3d &low_score_lat_mat, const Eigen::Matrix3d &high_score_lat_mat, double compare_tol);

  //************************************************************************************//

  /// \brief Returns an equivalent Lattice in Niggli form with a standard orientation
  Lattice niggli(const Lattice &lat, const SymGroup &point_grp, double tol);

  /// \brief Rotate the Lattice to a standard orientation using allowed point group operations
  Lattice standard_orientation(const Lattice &lat, const SymGroup &point_grp, double tol);

  namespace niggli_impl {

    /// \brief Returns an equivalent Lattice in Niggli form, but without setting standard orientation
    Lattice _niggli(const Lattice &lat, double compare_tol);

    /// \brief Same as ::_niggli but with Matrix3d type instead of Lattice
    Eigen::Matrix3d _niggli_mat(Eigen::Matrix3d lat_col_mat, double compare_tol);
  }

}

#endif
