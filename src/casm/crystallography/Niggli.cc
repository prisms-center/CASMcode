#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  NiggliRep::NiggliRep(const Eigen::Matrix3d &init_lat_col_mat):
    m_metrical_matrix(init_lat_col_mat.transpose() * init_lat_col_mat) {
  }

  NiggliRep::NiggliRep(const Lattice &init_lat):
    NiggliRep(init_lat.lat_column_mat()) {
  }

  double NiggliRep::A() const {
    return m_metrical_matrix(0, 0);
  }

  double NiggliRep::B() const {
    return m_metrical_matrix(1, 1);
  }

  double NiggliRep::C() const {
    return m_metrical_matrix(2, 2);
  }

  double NiggliRep::ksi() const {
    return m_metrical_matrix(2, 1) * 2;
  }

  double NiggliRep::eta() const {
    return m_metrical_matrix(2, 0) * 2;
  }

  double NiggliRep::zeta() const {
    return m_metrical_matrix(1, 0) * 2;
  }

  const Eigen::Matrix3d &NiggliRep::metrical_matrix() const {
    return m_metrical_matrix;
  }

  bool NiggliRep::meets_criteria_1(double compare_tol) const {
    bool does_meet = true;

    //For niggli you need A<=B
    if(CASM::compare(B(), A(), compare_tol)) {
      does_meet = false;
    }

    //If A==B, then niggli requires |ksi|<=|eta|
    else if(almost_equal(A(), B(), compare_tol) && compare(std::abs(eta()), std::abs(ksi()), compare_tol)) {
      does_meet = false;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_2(double compare_tol) const {
    bool does_meet = true;

    //For niggli you need B<=C
    if(compare(C(), B(), compare_tol)) {
      does_meet = false;
    }

    //If B==C, then niggli requires |eta|<=|zeta|
    else if(almost_equal(B(), C(), compare_tol) && compare(std::abs(zeta()), std::abs(eta()), compare_tol)) {
      does_meet = false;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_3(double compare_tol) const {
    bool does_meet = false;

    //For type one niggli cells, the angles are all less than 90 degrees
    if(compare(0.0, ksi(), compare_tol) && compare(0.0, eta(), compare_tol) && compare(0.0, zeta(), compare_tol)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_4(double compare_tol) const {
    bool does_meet = false;

    //For type two niggli cells, the angles are all more than or equal to 90 degrees
    if(!compare(0.0, ksi(), compare_tol) && !compare(0.0, eta(), compare_tol) && !compare(0.0, zeta(), compare_tol)) {
      does_meet = true;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_5(double compare_tol) const {
    bool does_meet = true;

    //For niggli you need |ksi|<=B
    if(compare(B(), std::abs(ksi()), compare_tol)) {
      does_meet = false;
    }

    //If ksi==B, then niggli requires zeta<=2*eta
    else if(almost_equal(ksi(), B(), compare_tol) && compare(2 * eta(), zeta(), compare_tol)) {
      does_meet = false;
    }

    //If ksi==-B, then niggli requires zeta==0
    else if(almost_equal(ksi(), (-1)*B(), compare_tol) && !almost_zero(zeta(), compare_tol)) {
      does_meet = false;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_6(double compare_tol) const {
    bool does_meet = true;

    //For niggli you need |eta|<=A
    if(compare(A(), std::abs(eta()), compare_tol)) {
      does_meet = false;
    }

    //If ksi==A, then niggli requires zeta<=2*ksi
    else if(almost_equal(eta(), A(), compare_tol) && compare(2 * ksi(), zeta(), compare_tol)) {
      does_meet = false;
    }

    //If ksi==-A, then niggli requires zeta==0
    else if(almost_equal(eta(), (-1)*A(), compare_tol) && !almost_zero(zeta(), compare_tol)) {
      does_meet = false;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_7(double compare_tol) const {
    bool does_meet = true;

    //For niggli you need |zeta|<=A
    if(compare(A(), std::abs(zeta()), compare_tol)) {
      does_meet = false;
    }

    //If zeta==A, then you need eta<=2*ksi
    else if(almost_equal(zeta(), A(), compare_tol) && compare(2 * ksi(), eta(), compare_tol)) {
      does_meet = false;
    }

    //If zeta==-A, then you need eta==0
    else if(almost_equal(zeta(), (-1)*A(), compare_tol) && !almost_zero(eta(), compare_tol)) {
      does_meet = false;
    }

    return does_meet;
  }

  bool NiggliRep::meets_criteria_8(double compare_tol) const {
    bool does_meet = true;

    double clobber = A() + B() + C() + ksi() + eta() + zeta();

    //For niggli you need C<=A+B+C+ksi+eta+zeta
    if(compare(clobber, C(), compare_tol)) {
      does_meet = false;
    }

    //If C==A+B+C+ksi+eta+zeta, then 2*A+2*eta+zeta<=0
    else if(almost_equal(clobber, C(), compare_tol) && compare(0.0, 2 * A() + 2 * eta() + zeta(), compare_tol)) {
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

    if(meets_criteria_1(compare_tol) &&
       meets_criteria_2(compare_tol) &&
       (meets_criteria_3(compare_tol) != meets_criteria_4(compare_tol)) &&
       meets_criteria_5(compare_tol) &&
       meets_criteria_6(compare_tol) &&
       meets_criteria_7(compare_tol) &&
       meets_criteria_8(compare_tol)) {
      totally_niggli = true;
    }

    return totally_niggli;
  }

  void NiggliRep::debug_criteria(double compare_tol) const {
    std::cout << meets_criteria_1(compare_tol) <<
              meets_criteria_2(compare_tol) <<
              meets_criteria_3(compare_tol) <<
              meets_criteria_4(compare_tol) <<
              meets_criteria_5(compare_tol) <<
              meets_criteria_6(compare_tol) <<
              meets_criteria_7(compare_tol) <<
              meets_criteria_8(compare_tol) << std::endl;

    return;
  }

  bool is_niggli(const Eigen::Matrix3d &test_lat_mat, double compare_tol) {
    NiggliRep niggtest(test_lat_mat);

    //niggtest.debug_criteria(compare_tol);
    return niggtest.is_niggli(compare_tol);
  }

  bool is_niggli(const Lattice &test_lat, double compare_tol) {
    return is_niggli(test_lat.lat_column_mat(), compare_tol);
  }

  /**
   * Converts a given Lattice into a new Lattice that satisfies
   * the Niggli criteria. Because there are several representations
   * of the Lattice that satisfy these conditions, the one with
   * the most standard orientation will be returned.
   */

  Lattice niggli(const Lattice &in_lat, double compare_tol) {
    const Lattice reduced_in = in_lat.get_reduced_cell();

    //const Lattice reduced_in = niggli_impl::_niggli(in_lat, compare_tol);
    Eigen::Matrix3d best_lat_mat = reduced_in.lat_column_mat();

    bool first_niggli = true;

    //Like the point group, but brute forcing for every possible transformation matrix ever with determinant 1 and elements -1, 0 or 1
    std::vector<Eigen::Matrix3i> canditate_trans_mats = positive_unimodular_matrices();

    Eigen::Matrix3i iknowthisoworks;
    iknowthisoworks << 1, 0, -1,
                    0, 1, 0,
                    0, 0, 1;

    for(auto it = canditate_trans_mats.begin(); it != canditate_trans_mats.end(); ++it) {
      Eigen::Matrix3d candidate_lat_mat = reduced_in.lat_column_mat() * it->cast<double>();

      if(is_niggli(candidate_lat_mat, compare_tol) && standard_orientation_compare(best_lat_mat, candidate_lat_mat, compare_tol)) {
        best_lat_mat = candidate_lat_mat;
      }

      //Always accept the first niggli cell regardless of the orientation
      else if(first_niggli && is_niggli(candidate_lat_mat, compare_tol)) {
        best_lat_mat = candidate_lat_mat;
        first_niggli = false;

      }
    }

    return Lattice(best_lat_mat);
  }

  /**
   * Return a canonical Lattice by first converting the
   * given Lattice into the standard Niggli form,
   * followed by applying the point group of the Lattice
   * so that the one oriented in the most standard manner
   * is selected.
   */

  Lattice canonical_equivalent_lattice(const Lattice &in_lat, const SymGroup &point_grp, double compare_tol) {
    //Ensure you at least get *something* back that's niggli
    Lattice most_canonical = niggli(in_lat, compare_tol);
    Eigen::Matrix3d most_canonical_lat_mat = most_canonical.lat_column_mat();

    for(auto it = point_grp.begin(); it != point_grp.end(); ++it) {
      Eigen::Matrix3d transformed_lat_mat = it->matrix() * in_lat.lat_column_mat();
      Lattice transformed_lat = Lattice(transformed_lat_mat);

      Eigen::Matrix3d candidate_lat_mat = niggli(transformed_lat, compare_tol).lat_column_mat();

      //Skip operations that change the handedness of the lattice
      if(candidate_lat_mat.determinant() < 0.0) {
        continue;
      }

      if(is_niggli(candidate_lat_mat, compare_tol) && standard_orientation_compare(most_canonical_lat_mat, candidate_lat_mat, compare_tol)) {
        most_canonical_lat_mat = candidate_lat_mat;
      }
    }

    return Lattice(most_canonical_lat_mat);
  }

  /**
   * Given two lattice column matrices, determine which one of them is oriented
   * in a more standard manner in space relative to the Cartesian coordinates.
   *
   * The routine returns "low_score_lat_mat<high_score_lat_mat", which is true if high_score_lat_mat
   * has a more standard spatial orientation.
   */

  bool standard_orientation_spatial_compare(const Eigen::Matrix3d &low_score_lat_mat, const Eigen::Matrix3d &high_score_lat_mat, double compare_tol) {
    //I have no idea who wrote this or why we decided this was the way to go.
    //This is simply a cut and paste from the old Lattice::standard_orientation.
    if(almost_equal(low_score_lat_mat(0, 0), high_score_lat_mat(0, 0), compare_tol)) {
      if(almost_equal(low_score_lat_mat(1, 0), high_score_lat_mat(1, 0), compare_tol)) {
        if(almost_equal(low_score_lat_mat(2, 0), high_score_lat_mat(2, 0), compare_tol)) {
          if(almost_equal(low_score_lat_mat(1, 1), high_score_lat_mat(1, 1), compare_tol)) {
            if(almost_equal(low_score_lat_mat(0, 1), high_score_lat_mat(0, 1), compare_tol)) {
              if(almost_equal(low_score_lat_mat(2, 1), high_score_lat_mat(2, 1), compare_tol)) {
                if(almost_equal(low_score_lat_mat(2, 2), high_score_lat_mat(2, 2), compare_tol)) {
                  if(almost_equal(low_score_lat_mat(0, 2), high_score_lat_mat(0, 2), compare_tol)) {
                    if(almost_equal(low_score_lat_mat(1, 2), high_score_lat_mat(1, 2), compare_tol)) {
                      return false;
                    }
                    else if(high_score_lat_mat(1, 2) > low_score_lat_mat(1, 2)) {
                      return  true;
                    }
                  }
                  else if(high_score_lat_mat(0, 2) > low_score_lat_mat(0, 2)) {
                    return  true;
                  }
                }
                else if(high_score_lat_mat(2, 2) > low_score_lat_mat(2, 2)) {
                  return  true;
                }
              }
              else if(high_score_lat_mat(2, 1) > low_score_lat_mat(2, 1)) {
                return  true;
              }
            }
            else if(high_score_lat_mat(0, 1) > low_score_lat_mat(0, 1)) {
              return  true;
            }
          }
          else if(high_score_lat_mat(1, 1) > low_score_lat_mat(1, 1)) {
            return  true;
          }
        }
        else if(high_score_lat_mat(2, 0) > low_score_lat_mat(2, 0)) {
          return  true;
        }
      }
      else if(high_score_lat_mat(1, 0) > low_score_lat_mat(1, 0)) {
        return  true;
      }
    }
    else if(high_score_lat_mat(0, 0) > low_score_lat_mat(0, 0)) {
      return  true;
    }

    return false;
  }

  /**
   * Compares two lattice column matrices to each other and determines which
   * one has a more standard orientation based on whether said matrices are
   * (bi)symmetric and how aligned the vectors are along the Cartesian directions.
   *
   * A bisymmetric matrix will always have better standard orientation than one that
   * is simply symmetric. In turn, a symmetric matrix will always have a better standard
   * orientation than on that is not. The criteria with lowest impact is the orientation
   * of the lattice vectors in space, as defined in ::standard_orientation_spatial_compare
   *
   * The routine returns "low_score_mat<high_score_mat", which is true if high_score_mat
   * has a more standard orientation.
   */

  bool standard_orientation_compare(const Eigen::Matrix3d &low_score_lat_mat, const Eigen::Matrix3d &high_score_lat_mat, double compare_tol) {
    bool low_score_is_symmetric = Eigen::is_symmetric(low_score_lat_mat, compare_tol);
    bool low_score_is_bisymmetric = Eigen::is_bisymmetric(low_score_lat_mat, compare_tol);

    bool high_score_is_symmetric = Eigen::is_symmetric(high_score_lat_mat, compare_tol);
    bool high_score_is_bisymmetric = Eigen::is_bisymmetric(high_score_lat_mat, compare_tol);

    //The high_score is less symmetric than the low score, therefore high<low
    if((low_score_is_bisymmetric && !high_score_is_bisymmetric) || (low_score_is_symmetric && !high_score_is_symmetric)) {
      return false;
    }

    //The high_score is more symmetric than the low_score, therefore high>low automatically (matrix symmetry trumps spatial orientation)
    else if((!low_score_is_bisymmetric && high_score_is_bisymmetric) || (!low_score_is_symmetric && high_score_is_symmetric)) {
      return true;
    }

    //If you made it here then the high_score and the low_score have the same symmetry level, so we check for the spatial orientation
    else {
      return standard_orientation_spatial_compare(low_score_lat_mat, high_score_lat_mat, compare_tol);
    }
  }

  //************************************************OLD CRAP GOES HERE************************************************************************//

  //    namespace niggli_impl
  //    {
  //
  //        /// Check that x < y, given some tolerance
  //        ///   Returns x < (y - tol)
  //        ///   Helper function for niggli
  //        bool _lt(double x, double y, double tol)
  //        {
  //            return x < (y - tol);
  //        }
  //
  //        /// Check that x > y, given some tolerance
  //        ///   returns y < x - tol
  //        ///   Helper function for niggli
  //        bool _gt(double x, double y, double tol)
  //        {
  //            return y < (x - tol);
  //        }
  //
  //        /// Check that x == y, given some tolerance
  //        ///   returns !(_lt(x,y,tol) || _gt(x,y,tol))
  //        ///   Helper function for niggli
  //        bool _eq(double x, double y, double tol)
  //        {
  //            return !(_lt(x, y, tol) || _gt(x, y, tol));
  //        }
  //
  //        /// Product of off-diagonal signs of S (= lat.transpose()*lat)
  //        ///   Helper function for niggli
  //        int _niggli_skew_product_step3(const Eigen::Matrix3d &S, double tol)
  //        {
  //            int S12 = 0;
  //            if(_gt(S(1, 2), 0.0, tol))
  //            {
  //                S12 = 1;
  //            }
  //            else if(_lt(S(1, 2), 0.0, tol))
  //            {
  //                S12 = -1;
  //            }
  //
  //            int S02 = 0;
  //            if(_gt(S(0, 2), 0.0, tol))
  //            {
  //                S02 = 1;
  //            }
  //            else if(_lt(S(0, 2), 0.0, tol))
  //            {
  //                S02 = -1;
  //            }
  //
  //            return S12 * S02 * S12;
  //        }
  //
  //        /// Product of off-diagonal signs of S (= lat.transpose()*lat)
  //        ///   Helper function for niggli
  //        int _niggli_skew_product(const Eigen::Matrix3d &S, double tol)
  //        {
  //            int S12 = 0;
  //            if(_gt(S(1, 2), 0.0, tol))
  //            {
  //                S12 = 1;
  //            }
  //            else if(_lt(S(1, 2), 0.0, tol))
  //            {
  //                S12 = -1;
  //            }
  //
  //            int S02 = 0;
  //            if(_gt(S(0, 2), 0.0, tol))
  //            {
  //                S02 = 1;
  //            }
  //            else if(_lt(S(0, 2), 0.0, tol))
  //            {
  //                S02 = -1;
  //            }
  //
  //            int S01 = 0;
  //            if(_gt(S(0, 1), 0.0, tol))
  //            {
  //                S01 = 1;
  //            }
  //            else if(_lt(S(0, 1), 0.0, tol))
  //            {
  //                S01 = -1;
  //            }
  //
  //            return S12 * S02 * S01;
  //        }
  //
  //        /// \brief Returns an equivalent \ref Lattice in Niggli form
  //        ///
  //        /// \returns an equivalent \ref Lattice in Niggli form
  //        ///
  //        /// \param lat a \ref Lattice
  //        /// \param tol tolerance for floating point comparisons
  //        ///
  //        /// The Niggli cell is a unique choice of lattice vectors for a particular lattice.
  //        /// It minimizes lattice vector lengths, and chooses a particular angular orientation.
  //        ///
  //        /// This implementation function does not set the standard spatial orientation.
  //        ///
  //        /// \see
  //        /// I. Krivy and B. Gruber, Acta Cryst. (1976). A32, 297.
  //        /// <a href="http://dx.doi.org/10.1107/S0567739476000636">[doi:10.1107/S0567739476000636]</a>
  //        /// R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, Acta Cryst. (2004). A60, 1.
  //        /// <a href="http://dx.doi.org/10.1107/S010876730302186X"> [doi:10.1107/S010876730302186X]</a>
  //        ///
  //        Eigen::Matrix3d _niggli_mat(Eigen::Matrix3d reduced, double tol)
  //        {
  //
  //            //std::cout << "begin _niggli(const Lattice &lat, double tol)" << std::endl;
  //
  //            // Get helper functions
  //            using namespace niggli_impl;
  //
  //            //  S = lat.transpose()*lat, is a matrix of lattice vector scalar products,
  //            //    this is sometimes called the metric tensor
  //            //
  //            //    S(0,0) = a*a, S(0,1) = a*b, etc.
  //            //
  //            //Eigen::Matrix3d reduced = lat.lat_column_mat();
  //            Eigen::Matrix3d S = reduced.transpose() * reduced;
  //
  //            while(true)
  //            {
  //
  //                // in notes: a = reduced.col(0), b = reduced.col(1), c = reduced.col(2)
  //                //   a*a is scalar product
  //
  //
  //                // 1)
  //                // if a*a > b*b or
  //                //    a*a == b*b and std::abs(b*c*2.0) > std::abs(a*c*2.0),
  //                // then permute a and b
  //                if(_gt(S(0, 0), S(1, 1), tol) ||
  //                   (_eq(S(0, 0), S(1, 1), tol) && _gt(std::abs(S(1, 2)), std::abs(S(0, 2)), tol)))
  //                {
  //                    reduced.col(0).swap(reduced.col(1));
  //                    reduced *= -1.0;
  //                    S = reduced.transpose() * reduced;
  //                }
  //
  //
  //                // 2)
  //                // if b*b > c*c or
  //                //    b*b == c*c and std::abs(a*c*2.0) > std::abs(a*b*2.0),
  //                // then permute b and c
  //                if(_gt(S(1, 1), S(2, 2), tol) ||
  //                   (_eq(S(1, 1), S(2, 2), tol) && _gt(std::abs(S(0, 2)), std::abs(S(0, 1)), tol)))
  //                {
  //                    reduced.col(1).swap(reduced.col(2));
  //                    reduced *= -1.0;
  //                    S = reduced.transpose() * reduced;
  //
  //                    continue;
  //                }
  //
  //                // 3)
  //                // if (2*b*c)*(2*a*c)*(2*a*b) > 0,  *** is this right? **
  //                if(_niggli_skew_product(S, tol) > 0)
  //                {
  //
  //                    // if (2*b*c)*(2*a*c)*(2*b*c) > 0,
  //                    // if(_niggli_skew_product_step3(S, tol) > 0) {
  //
  //                    // if (b*c) < 0.0, flip a
  //                    if(_lt(S(1, 2), 0.0, tol))
  //                    {
  //                        reduced.col(0) *= -1;
  //                    }
  //
  //                    // if (a*c) < 0.0, flip b
  //                    if(_lt(S(0, 2), 0.0, tol))
  //                    {
  //                        reduced.col(1) *= -1;
  //                    }
  //
  //                    // if (a*b) < 0.0, flip c
  //                    if(_lt(S(0, 1), 0.0, tol))
  //                    {
  //                        reduced.col(2) *= -1;
  //                    }
  //
  //                    S = reduced.transpose() * reduced;
  //
  //                }
  //
  //                // 4)
  //                // if (2*b*c)*(2*a*c)*(2*a*b) <= 0,
  //                if(_niggli_skew_product(S, tol) <= 0)
  //                {
  //
  //                    int i = 1, j = 1, k = 1;
  //                    int *p;
  //
  //
  //                    // !! paper says (a*b), but I think they mean (b*c)
  //                    // if (b*c) > 0.0, i = -1
  //                    if(_gt(S(1, 2), 0.0, tol))
  //                    {
  //                        i = -1;
  //                    }
  //                    // else if( !(b*c < 0)), p = &i;
  //                    else if(!(_lt(S(1, 2), 0.0, tol)))
  //                    {
  //                        p = &i;
  //                    }
  //
  //                    // if (a*c) > 0.0, j = -1
  //                    if(_gt(S(0, 2), 0.0, tol))
  //                    {
  //                        j = -1;
  //                    }
  //                    // else if( !(a*c < 0)), p = &j;
  //                    else if(!(_lt(S(0, 2), 0.0, tol)))
  //                    {
  //                        p = &j;
  //                    }
  //
  //                    // if (a*b) > 0.0, k = -1
  //                    if(_gt(S(0, 1), 0.0, tol))
  //                    {
  //                        k = -1;
  //                    }
  //                    // else if( !(a*b < 0)), p = &k;
  //                    else if(!(_lt(S(0, 1), 0.0, tol)))
  //                    {
  //                        p = &k;
  //                    }
  //
  //                    if(i * j * k < 0)
  //                    {
  //                        *p = -1;
  //                    }
  //
  //                    reduced.col(0) *= i;
  //                    reduced.col(1) *= j;
  //                    reduced.col(2) *= k;
  //
  //                    S = reduced.transpose() * reduced;
  //
  //                }
  //
  //                // 5)
  //                // if std::abs(2.0*S(1,2)) > S(1,1) or
  //                //    (2.0*S(1,2) == S(1,1) and 2.0*2.0*S(0,2) < 2.0*S(0,1)) or
  //                //    (2.0*S(1,2) == -S(1,1) and 2.0*S(0,1) < 0.0)
  //                if(_gt(std::abs(2.0 * S(1, 2)), S(1, 1), tol) ||
  //                   (_eq(2.0 * S(1, 2), S(1, 1), tol) && _lt(2.0 * S(0, 2), S(0, 1), tol)) ||
  //                   (_eq(2.0 * S(1, 2), -S(1, 1), tol) && _lt(S(0, 1), 0.0, tol)))
  //                {
  //
  //                    // if (b*c) > 0.0, subtract b from c
  //                    if(_gt(S(1, 2), 0.0, tol))
  //                    {
  //                        reduced.col(2) -= reduced.col(1);
  //                    }
  //                    // else if (b*c) < 0.0, add b to c
  //                    else if(_lt(S(1, 2), 0.0, tol))
  //                    {
  //                        reduced.col(2) += reduced.col(1);
  //                    }
  //
  //                    S = reduced.transpose() * reduced;
  //
  //                    continue;
  //                }
  //
  //                // 6)
  //                // if std::abs(2.0*S(0,2)) > S(0,0) or
  //                //    (2.0*S(0,2) == S(0,0) and 2.0*2.0*S(1,2) < 2.0*S(0,1)) or
  //                //    (2.0*S(0,2) == -S(0,0) and 2.0*S(0,1) < 0.0)
  //                if(_gt(std::abs(2.0 * S(0, 2)), S(0, 0), tol) ||
  //                   (_eq(2.0 * S(0, 2), S(0, 0), tol) && _lt(2.0 * S(1, 2), S(0, 1), tol)) ||
  //                   (_eq(2.0 * S(0, 2), -S(0, 0), tol) && _lt(S(0, 1), 0.0, tol)))
  //                {
  //
  //                    // if (a*c) > 0.0, subtract a from c
  //                    if(_gt(S(0, 2), 0.0, tol))
  //                    {
  //                        reduced.col(2) -= reduced.col(0);
  //                    }
  //                    // else if (a*c) < 0.0, add a to c
  //                    else if(_lt(S(0, 2), 0.0, tol))
  //                    {
  //                        reduced.col(2) += reduced.col(0);
  //                    }
  //
  //                    S = reduced.transpose() * reduced;
  //
  //                    continue;
  //                }
  //
  //                // 7)
  //                // if std::abs(2.0*S(0,1)) > S(0,0) or
  //                //    (2.0*S(0,1) == S(0,0) and 2.0*2.0*S(1,2) < 2.0*S(0,2)) or
  //                //    (2.0*S(0,1) == -S(0,0) and 2.0*S(0,2) < 0.0)
  //                if(_gt(std::abs(2.0 * S(0, 1)), S(0, 0), tol) ||
  //                   (_eq(2.0 * S(0, 1), S(0, 0), tol) && _lt(2.0 * S(1, 2), S(0, 2), tol)) ||
  //                   (_eq(2.0 * S(0, 1), -S(0, 0), tol) && _lt(S(0, 2), 0.0, tol)))
  //                {
  //
  //                    // if (a*b) > 0.0, subtract a from b
  //                    if(_gt(S(0, 1), 0.0, tol))
  //                    {
  //                        reduced.col(1) -= reduced.col(0);
  //                    }
  //                    // else if (a*b) < 0.0, add a to b
  //                    else if(_lt(S(0, 1), 0.0, tol))
  //                    {
  //                        reduced.col(1) += reduced.col(0);
  //                    }
  //
  //                    S = reduced.transpose() * reduced;
  //
  //                    continue;
  //                }
  //
  //                // 8)
  //                // let tmp = 2*b*c + 2*a*c + 2*a*b + a*a + b*b
  //                // if  tmp < 0.0 or
  //                //     tmp == 0 and 2*(a*a + 2*a*c) + 2*a*b > 0
  //                double tmp = 2.0 * S(1, 2) + 2.0 * S(0, 2) + 2.0 * S(0, 1) + S(0, 0) + S(1, 1);
  //                if(_lt(tmp, 0.0, tol) ||
  //                   (_eq(tmp, 0.0, tol) && _gt(2.0 * (S(0, 0) + 2.0 * S(0, 2)) + 2.0 * S(0, 1), 0.0, tol)))
  //                {
  //
  //                    // add a and b to c
  //                    reduced.col(2) += reduced.col(0);
  //                    reduced.col(2) += reduced.col(1);
  //
  //                    S = reduced.transpose() * reduced;
  //
  //                    continue;
  //                }
  //
  //                break;
  //
  //            } // end while
  //
  //            return reduced;
  //        }
  //
  //        Lattice _niggli(const Lattice &lat, double compare_tol)
  //        {
  //            return Lattice(_niggli_mat(lat.lat_column_mat(), compare_tol));
  //        }
  //
  //
  //    }
  //
  //    /// \brief Returns an equivalent \ref Lattice in Niggli form
  //    ///
  //    /// \returns an equivalent \ref Lattice in Niggli form
  //    ///
  //    /// \param lat a \ref Lattice
  //    /// \param tol tolerance for floating point comparisons
  //    ///
  //    /// The Niggli cell is a unique choice of lattice vectors for a particular lattice.
  //    /// It minimizes lattice vector lengths, and chooses a particular angular orientation.
  //    ///
  //    /// With the angular orientation fixed, the final spatial orientation of the lattice is
  //    /// set using standard_orientation.
  //    ///
  //    /// \see
  //    /// I. Krivy and B. Gruber, Acta Cryst. (1976). A32, 297.
  //    /// <a href="http://dx.doi.org/10.1107/S0567739476000636">[doi:10.1107/S0567739476000636]</a>
  //    /// R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, Acta Cryst. (2004). A60, 1.
  //    /// <a href="http://dx.doi.org/10.1107/S010876730302186X"> [doi:10.1107/S010876730302186X]</a>
  //    ///
  //    Lattice niggli(const Lattice &lat, const SymGroup &point_grp, double tol)
  //    {
  //        Lattice reduced = niggli_impl::_niggli(lat, tol);
  //        return standard_orientation(reduced, point_grp, tol);
  //    }
  //
  //    /**
  //     * Rotates the Lattice to a standard orientation using point group operations to
  //     * reorient the unit cell.
  //     *
  //     * The unit cell vectors of the lattice are chosen so that they are as aligned as possible
  //     * with the Cartesian coordinate system, but preference is always given to lattices
  //     * whose matrices are bisymmetric or symmetric.
  //     */
  //
  //    Lattice standard_orientation(const Lattice &in_lat, const SymGroup &point_grp, double compare_tol)
  //    {
  //
  //        Eigen::Matrix3d best_lat_mat = in_lat.lat_column_mat();
  //
  //        for(int i = 0; i < point_grp.size(); i++)
  //        {
  //
  //            const Eigen::Matrix3d candidate_lat_mat  = point_grp[i].matrix() * in_lat.lat_column_mat();
  //
  //            //Skip any operation that changes the handedness of the lattice
  //            if(candidate_lat_mat.determinant() < 0.0)
  //            {
  //                continue;
  //            }
  //
  //            if(standard_orientation_compare(best_lat_mat, candidate_lat_mat, compare_tol))
  //            {
  //                best_lat_mat = candidate_lat_mat;
  //            }
  //        }
  //
  //        return Lattice(best_lat_mat);
  //    }
}
