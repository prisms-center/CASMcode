#include "casm/crystallography/Lattice.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
  namespace xtal {

    Lattice::Lattice(Eigen::Ref<const Eigen::Vector3d> const &vec1,
                     Eigen::Ref<const Eigen::Vector3d> const &vec2,
                     Eigen::Ref<const Eigen::Vector3d> const &vec3,                     double xtal_tol,
                     bool force)
      : m_tol(xtal_tol) {
      m_lat_mat << vec1, vec2, vec3;
      if(!force && m_lat_mat.determinant() < 0) {
        // this->make_right_handed();
        // throw std::runtime_error("Attempted to construct a left-handed lattice. Try again or override if you know what you're doing");
      }
      m_inv_lat_mat = m_lat_mat.inverse();
    }

    std::vector<Eigen::Matrix3d> _skew_transforms() {
      // Creates 12 skew matrices + 12 "double skew" matrices
      std::vector<Eigen::Matrix3d> skew(24, Eigen::Matrix3d::Identity());
      Index i, j, k;
      Index l = 0;
      for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
          if(i != j) {
            // Checks to make sure we're not at a diagonal
            skew[l](i, j) = 1;
            skew[l + 12](i, j) = -1;
            ++l;
          }
        }
      }

      // the 12 "double skew" matrices corresponding to body diagonal substitution
      for(i = 0; i < 3; i++) {
        // column

        j = (i + 1) % 3;
        k = (i + 2) % 3;

        skew[l](j, i) = 1;
        skew[l](k, i) = 1;

        skew[l + 12](j, i) = -1;
        skew[l + 12](k, i) = -1;

        ++l;

        skew[l](j, i) = 1;
        skew[l](k, i) = -1;

        skew[l + 12](j, i) = -1;
        skew[l + 12](k, i) = 1;

        ++l;
      }
      return skew;
    }

    std::vector<Eigen::Matrix3d> const &Lattice::skew_transforms() {
      static std::vector<Eigen::Matrix3d> result = _skew_transforms();
      return result;
    }

    /// Construct Lattice from a matrix of lattice vectors, where lattice vectors are columns
    ///(e.g., lat_mat is equivalent to lat_column_mat())
    Lattice::Lattice(const Eigen::Ref<const Eigen::Matrix3d> &lat_mat, double xtal_tol, bool force)
      : m_lat_mat(lat_mat), m_inv_lat_mat(lat_mat.inverse()), m_tol(xtal_tol) {
      if(!force && m_lat_mat.determinant() < 0) {
        // this->make_right_handed();
        // throw std::runtime_error("Attempted to construct a left-handed lattice. Try again or override if you know what you're doing");
      }
    }

    Lattice Lattice::fcc(double tol) {
      Eigen::Matrix3d latmat;
      // clang-format off
      latmat << 0, 1, 1,
             1, 0, 1,
             1, 1, 0;
      // clang-format on
      latmat /= pow(latmat.determinant(), 1.0 / 3.0);
      return Lattice(latmat, tol);
    }

    Lattice Lattice::bcc(double tol) {
      Eigen::Matrix3d latmat;
      // clang-format off
      latmat << -1, 1, 1,
             1, -1, 1,
             1, 1, -1;
      // clang-format on
      latmat /= pow(latmat.determinant(), 1.0 / 3.0);
      return Lattice(latmat, tol);
    }

    Lattice Lattice::cubic(double tol) {
      return Lattice(Eigen::Matrix3d::Identity(), tol);
    }

    Lattice Lattice::hexagonal(double tol) {
      Eigen::Matrix3d latmat;
      // clang-format off
      latmat << 1, -1.0 / sqrt(3.0),
             0, 0, 2.0 / sqrt(3.0),
             0, 0, 0, sqrt(3.0);
      //clang-format off

      return Lattice(latmat.transpose(), tol);
    }

    Lattice Lattice::scaled_lattice(double scale) const {
      return Lattice(scale * lat_column_mat(), tol());
    }

    double Lattice::length(Index i) const {
      // Calculates Lengths
      return m_lat_mat.col(i).norm();
    }

    double Lattice::angle(Index i) const {
      return (180.0 / M_PI) * CASM::angle(m_lat_mat.col((i + 1) % 3), m_lat_mat.col((i + 2) % 3));
    }

    void Lattice::read(std::istream &stream) {
      double scale;
      stream >> scale;
      stream >> m_lat_mat;
      m_lat_mat *= scale;
      m_lat_mat.transposeInPlace();
      m_inv_lat_mat = m_lat_mat.inverse();

      return;
    }

    void Lattice::print(std::ostream &stream, int _prec) const {
      int tprec = stream.precision();
      std::ios::fmtflags tflags = stream.flags();
      stream.precision(_prec);
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      Eigen::IOFormat format(_prec);

      stream << 1.0 << '\n';
      stream << m_lat_mat.transpose().format(format);

      stream.precision(tprec);
      stream.flags(tflags);

      return;
    }

    Lattice Lattice::reciprocal() const {
      return Lattice(2 * M_PI * inv_lat_column_mat().transpose(), tol());
    }

    double Lattice::boxiness() const {
      return 1 / (this->inv_lat_column_mat().colwise().norm().sum());
    }

    std::vector<int> Lattice::calc_kpoints(std::vector<int> prim_kpoints, Lattice prim_lat) {
      std::vector<int> super_kpoints = prim_kpoints;
      //    Lattice prim_recip_lat = (lattice.primitive->reciprocal());
      Lattice prim_recip_lat = (prim_lat.reciprocal());
      Lattice recip_lat = (*this).reciprocal();
      double prim_density = (prim_kpoints[0] * prim_kpoints[1] * prim_kpoints[2]) / (prim_recip_lat.volume());
      double super_density = 0;

      std::vector<double> prim_vec_lengths;

      for(int i = 0; i < 3; i++) {
        prim_vec_lengths.push_back(prim_recip_lat.length(i));
      }

      double shortest = *std::min_element(prim_vec_lengths.begin(), prim_vec_lengths.end());
      int short_ind = find_index(prim_vec_lengths, shortest);

      double scale = (prim_kpoints[short_ind] / shortest);

      for(int i = 0; i < 3; i++) {
        super_kpoints[i] = int(ceil(scale * recip_lat.length(i)));
      }

      super_density = (super_kpoints[0] * super_kpoints[1] * super_kpoints[2]) / (recip_lat.volume());

      while(super_density < prim_density) {
        prim_kpoints[short_ind]++;

        scale = (prim_kpoints[short_ind] / shortest);

        for(int i = 0; i < 3; i++) {
          super_kpoints[i] = int(ceil(scale * recip_lat.length(i)));
        }

        super_density = (super_kpoints[0] * super_kpoints[1] * super_kpoints[2]) / (recip_lat.volume());
      }

      return super_kpoints;
    }

    /**
     * \brief This function finds the reduced cell from the given primitive cell.
     *
     * First, the translation vectors of the primitive cell are determined, then reduced.
     * The Bravais lattice is determined using Niggli's transformations. The method used to
     * determine the basis vectors are to find body and face diagonals across the cell.
     * The vectors satisfying these main conditions can be found by ensuring:
     * 		 i) the shortest face diagonal is not shorter than the longest edge of the same face
     * 		 ii) the shortest body diagonal is not shorter than the longest edge of the cell
     *
     * This is accomplished by enumerating the skew matrices (identity matrix with one of the
     * off-diagonal terms =1 or =-1).  This will replace one of the lattice vectors with the
     * face diagonal corresponding to that linear combination.  If the length of any of the
     * new skewed lattice vectors is shorter than the length of the original cell, it is replaced.
     * All skew matrices have determinants equal to 1, which preserves the volume of the original
     * cell.
     *
     */
    Lattice Lattice::reduced_cell2() const {
      std::vector<Eigen::Matrix3d> const &skew(skew_transforms());
      std::vector<Eigen::Matrix3d> const &ntrans(NiggliRep::cell_invariant_transforms());

      Eigen::Matrix3d trans;
      trans.setIdentity();

      Eigen::Matrix3d t_metric, reduced_metric, t_imetric, reduced_imetric;

      reduced_metric = lat_column_mat().transpose() * lat_column_mat();
      reduced_imetric = reduced_metric.inverse();

      bool minimized = false;
      while(!minimized) {
        minimized = true;
        for(Index i = 0; i < skew.size(); ++i) {
          // Lattice vectors times their transpose give length and angle info
          t_metric = skew[i].transpose() * reduced_metric * skew[i];
          t_imetric = skew[(i + 12) % 12].transpose() * reduced_imetric * skew[(i + 12) % 12];

          if(almost_equal(t_metric.trace(), reduced_metric.trace(), tol() * tol())) {
            if(almost_equal(t_imetric.trace(), reduced_imetric.trace(), tol() * tol())) {
              if(niggli_index(lat_column_mat() * trans * skew[i], tol()) > niggli_index(lat_column_mat() * trans, tol())) {
                /*if(fabs(t_metric(0,1)) < fabs(reduced_metric(0,1))
                 || fabs(t_metric(0,2)) < fabs(reduced_metric(0,2))
                 || fabs(t_metric(1,2)) < fabs(reduced_metric(1,2))) {*/

                reduced_metric = t_metric;
                reduced_imetric = t_imetric;
                trans *= skew[i];
                minimized = false;
              }
            }
            else if(t_imetric.trace() < reduced_imetric.trace()) {
              reduced_metric = t_metric;
              reduced_imetric = t_imetric;
              trans *= skew[i];
              minimized = false;
            }
          }
          else if(t_metric.trace() < reduced_metric.trace()) {
            reduced_metric = t_metric;
            reduced_imetric = t_imetric;
            trans *= skew[i];
            minimized = false;
          }
        }

        Index b_index = niggli_index(lat_column_mat() * trans, tol());
        for(Index i = 0; i < ntrans.size(); ++i) {
          Index tb_index = niggli_index(lat_column_mat() * trans * ntrans[i], tol());
          if(tb_index > b_index) {
            b_index = tb_index;
            reduced_metric = ntrans[i].transpose() * reduced_metric * ntrans[i];
            reduced_imetric = reduced_metric.inverse();
            trans *= ntrans[i];
            minimized = false;
          }
        }
      }

      return Lattice(lat_column_mat() * trans, tol());
    }

    /**
     * \brief This function finds the reduced cell from the given primitive cell.
     *
     * This implementation is the LLL algorithm as laid out by Hoffstein, Jeffrey; Pipher, Jill; Silverman, J.H. (2008).
     * An Introduction to Mathematical Cryptography. Springer. ISBN 978-0-387-77993-5.
     * Some corrections have been made from  Silverman, Joseph. "Introduction to Mathematical Cryptography Errata" (PDF).
     * Brown University Mathematics Dept. Retrieved 5 May 2015.
     *
     */
    Lattice Lattice::reduced_cell() const {
      Eigen::Matrix3d reduced_lat = lat_column_mat();
      bool right_handed = (reduced_lat.determinant() > 0);
      Eigen::HouseholderQR<Eigen::Matrix3d> qr(reduced_lat);
      Eigen::Matrix3d Q = qr.householderQ();
      Eigen::Matrix3d R = Q.inverse() * reduced_lat;
      Eigen::Matrix3d ortho = Q;
      ortho.col(0) = ortho.col(0) * R(0, 0);
      ortho.col(1) = ortho.col(1) * R(1, 1);
      ortho.col(2) = ortho.col(2) * R(2, 2);
      Index k = 1;
      while(k < 3) {
        for(int j = k - 1; j >= 0; j--) {
          double mu = reduced_lat.col(k).dot(ortho.col(j)) / ortho.col(j).squaredNorm();
          if(fabs(mu) > 0.5000001) {
            reduced_lat.col(k) = reduced_lat.col(k) - round(mu) * reduced_lat.col(j);
            Eigen::HouseholderQR<Eigen::Matrix3d> qr2(reduced_lat);
            Q = qr2.householderQ();
            R = Q.inverse() * reduced_lat;
            ortho = Q;
            ortho.col(0) = ortho.col(0) * R(0, 0);
            ortho.col(1) = ortho.col(1) * R(1, 1);
            ortho.col(2) = ortho.col(2) * R(2, 2);
          }
        }
        double mu2 = reduced_lat.col(k).dot(ortho.col(k - 1)) / ortho.col(k - 1).squaredNorm();
        if((ortho.col(k) + mu2 * ortho.col(k - 1)).squaredNorm() > (0.75) * ortho.col(k - 1).squaredNorm()) {
          k = k + 1;
        }
        else {
          Eigen::Vector3d tmp = reduced_lat.col(k);
          reduced_lat.col(k) = reduced_lat.col(k - 1);
          reduced_lat.col(k - 1) = tmp;
          Eigen::HouseholderQR<Eigen::Matrix3d> qr3(reduced_lat);
          Q = qr3.householderQ();
          R = Q.inverse() * reduced_lat;
          ortho = Q;
          ortho.col(0) = ortho.col(0) * R(0, 0);
          ortho.col(1) = ortho.col(1) * R(1, 1);
          ortho.col(2) = ortho.col(2) * R(2, 2);
          if(k > 1) {
            k = k - 1;
          }
        }
      }
      if(right_handed) {
        return Lattice(reduced_lat).make_right_handed();
      }
      return Lattice(reduced_lat);
    }

    double Lattice::max_voronoi_measure(const Eigen::Vector3d &pos, Eigen::Vector3d &lattice_trans) const {
      Eigen::MatrixXd::Index maxv;
      double maxproj = (voronoi_table() * pos).maxCoeff(&maxv);

      lattice_trans = (2. * floor(maxproj / 2. + (0.5 - TOL / 2.)) / m_voronoi_table.row(maxv).squaredNorm()) * m_voronoi_table.row(maxv);

      return maxproj;
    }

    int Lattice::voronoi_number(const Eigen::Vector3d &pos) const {

      int tnum = 0;
      double tproj = 0;

      Eigen::MatrixXd const &vtable = voronoi_table();

      for(Index nv = 0; nv < vtable.rows(); nv++) {
        tproj = vtable.row(nv) * pos;
        if(almost_equal(tproj, 1.0)) {
          tnum++;
        }
        else if(tproj > 1.0) {
          return -1;
        }
      }

      return tnum;
    }

    void Lattice::_generate_voronoi_table() const {
      // There are no fewer than 12 points in the voronoi table
      m_voronoi_table.resize(12, 3);

      m_inner_voronoi_radius = 1e20;

      Eigen::Vector3d tpoint;
      int i;
      int nrows = 1;
      Lattice tlat_reduced(reduced_cell());

      // Count over all lattice vectors, face diagonals, and body diagonals
      // originating from origin;
      EigenCounter<Eigen::Vector3i> combo_count(Eigen::Vector3i(-1, -1, -1), Eigen::Vector3i(1, 1, 1), Eigen::Vector3i(1, 1, 1));

      // For each linear combination, check to see if it is on a face, edge, or vertex of the voronoi cell
      for(; combo_count.valid(); ++combo_count) {
        if(combo_count().isZero())
          continue;
        // A linear combination does not fall on the voronoi boundary if the angle between
        // any two of the vectors forming that combination are acute
        for(i = 0; i < 3; i++) {
          if(combo_count[(i + 1) % 3] == 0 || combo_count[(i + 2) % 3] == 0)
            continue;
          if((180. / M_PI) * CASM::angle(combo_count[(i + 1) % 3] * tlat_reduced[(i + 1) % 3],
                                         combo_count[(i + 2) % 3] * tlat_reduced[(i + 2) % 3]) +
             TOL <
             90.) {
            break;
          }
        }

        if(i == 3) {
          if(nrows > m_voronoi_table.rows())
            m_voronoi_table.conservativeResize(nrows, Eigen::NoChange);

          tpoint = tlat_reduced.lat_column_mat() * combo_count().cast<double>();

          double t_rad = tpoint.norm();
          if((t_rad / 2.) < m_inner_voronoi_radius)
            m_inner_voronoi_radius = t_rad / 2.;

          m_voronoi_table.row(nrows - 1) = (2. / (t_rad * t_rad)) * tpoint;
          nrows++;
        }
      }
      return;
    }

    /**
     * This function, given a linearly independent set of lattice vectors, finds the
     * dimensions along the unit cell vectors such that a sphere of given radius
     * fits within a uniform grid of 2dim[1]x2dim[2]x2dim[3] lattice points
     * centered at the origin.
     *
     * The algorithm works by getting the normal (e.g. n1) to each pair of lattice
     * vectors (e.g. a2, a3), scaling this normal to have length radius and
     * then projecting this normal parallel to the a2,a3 plane onto the
     * remaining lattice vector a1. This will tell us the number of a1 vectors
     * needed to make a grid to encompass the sphere.
     */

    Eigen::Vector3i Lattice::enclose_sphere(double radius) const {

      // reciprocal vectors
      Eigen::Matrix3d recip(inv_lat_column_mat().transpose());

      for(int i = 0; i < 3; i++) {
        recip.col(i) *= radius / recip.col(i).norm();
      }
      // recip contains three column vectors of length 'radius' pointed along plane normals
      auto lambda = [](double val) {
        return ceil(val);
      };
      return (inv_lat_column_mat() * recip).cwiseAbs().unaryExpr(lambda).colwise().maxCoeff().cast<int>();
    }

    /**
     * A lattice is considered right handed when the
     * determinant of the lattice vector matrix is positive.
     * This routine will flip the sign of all the lattice
     * vectors if it finds that the determinant is negative
     */

    Lattice &Lattice::make_right_handed() {

      if(lat_column_mat().determinant() < 0) {
        m_lat_mat = -m_lat_mat;
        m_inv_lat_mat = m_lat_mat.inverse();
      }

      return *this;
    }

    Eigen::Vector3i Lattice::millers(Eigen::Vector3d plane_normal) const {
      // Get fractional coordinates of plane_normal in recip_lattice
      // These are h, k, l
      // For miller indeces h, k and l    plane_normal[CART]=h*a.recip+k*b.recip+l*c.recip
      return scale_to_int(lat_column_mat().transpose() * plane_normal, tol());
    }

    /**
     *  Using miller indeces, the intercept of the plane is calculated. All intercepts
     *  are multiplied by a factor such that their values are integers. This corresponds
     *  to the plane being shifted so that it lands on three lattice sites. Using these
     *  lattice sites new vectors A and B that lie on the desired plane can be constructed.
     *  The choice of vector C is somewhat arbitrary. lattice_in_plane will first
     *  construct C to be the shortest most orthogonal vector possible. The user is then
     *  given the option to make C more orthogonal to A and B in exchange for a greater
     *  length.
     *  The returned lattice will be the one containing the most orthogonal C
     *  that is still smaller than max_vol
     */

    Lattice Lattice::lattice_in_plane(Eigen::Vector3i millers, int max_vol) const {
      // John G 121030
      // Hold new lattice vectors in these. Then at the end we make an actual Lattice out of it
      Eigen::Matrix3d surface_cell, last_surface_cell; // Holds new lattice vectors, two of which are in the surface plane

      // Miller indeces of 100, 010 or 001 mean you don't need a new cell to expose the plane, however
      // you may want to reorient the vectors so that the ab plane is exposed (useful for Structure::stitch)

      if(millers == Eigen::Vector3i(0, 1, 0)) {
        return Lattice(lat_column_mat().col(2), lat_column_mat().col(0), lat_column_mat().col(1), tol());
      }

      else if(millers == Eigen::Vector3i(1, 0, 0)) {
        return Lattice(lat_column_mat().col(1), lat_column_mat().col(2), lat_column_mat().col(0), tol());
      }

      else if(millers == Eigen::Vector3i(0, 0, 1)) {
        // || millers==Eigen::Vector3i(0,1,0) || millers==Eigen::Vector3i(1,0,0))
        return *this;
      }

      // Miller indeces of xx0, 0xx or x0x mean one of the lattice vectors is the same as
      // one of the primitive lattice vectors. This means we only need to get TWO
      // points that are on the plane, which we connect to get the second lattice vector
      else if(millers[0] == 0 || millers[1] == 0 || millers[2] == 0) {
        int zero = 0;

        // Find out which miller index is 0
        while(true) {
          if(millers[zero] == 0) {
            break;
          }

          else {
            zero++;
          }
        }

        surface_cell.col(0) = lat_column_mat().col(zero);

        Eigen::Vector3d H_miller_point, K_miller_point;
        Eigen::Vector3d HK;
        Eigen::Vector3d millers_dubs;

        // Turn integer millers into doubles for mathematical purposes (inverse)
        millers_dubs = millers.cast<double>();

        Eigen::Vector3d inv_miller_dubs;
        Eigen::Vector3i inv_miller;

        // In actualility, the inverse miller of 0 is infinity, we set it to 0 here so that we can use
        // scale_to_int without going crazy. Since we won't be using this inverse miller index it's OK
        inv_miller_dubs[zero] = 0;
        inv_miller_dubs[(zero + 1) % 3] = 1.0 / millers_dubs[(zero + 1) % 3];
        inv_miller_dubs[(zero + 2) % 3] = 1.0 / millers_dubs[(zero + 2) % 3];

        inv_miller = scale_to_int(inv_miller_dubs, TOL);
        H_miller_point = inv_miller[(zero + 1) % 3] * lat_column_mat().col((zero + 1) % 3);
        K_miller_point = inv_miller[(zero + 2) % 3] * lat_column_mat().col((zero + 2) % 3);

        HK = K_miller_point - H_miller_point;
        surface_cell.col(1) = HK;
      }

      else {
        // Get three points that lie on the plane
        // We'll want to find points that lie on the plane AND land of lattice points. In order to do
        // this we need the miller inverses multiplied by a factor that makes them integers
        Eigen::Vector3d H_miller_point, K_miller_point, L_miller_point;
        Eigen::Vector3d inv_miller_dubs;
        Eigen::Vector3i inv_miller;
        // Turn integer millers into doubles for mathematical purposes (inverse)
        Eigen::Vector3d millers_dubs;
        millers_dubs = millers.cast<double>();

        inv_miller_dubs[0] = 1.0 / millers_dubs[0];
        inv_miller_dubs[1] = 1.0 / millers_dubs[1];
        inv_miller_dubs[2] = 1.0 / millers_dubs[2];

        inv_miller = scale_to_int(inv_miller_dubs, TOL);

        H_miller_point = inv_miller[0] * lat_column_mat().col(0);
        K_miller_point = inv_miller[1] * lat_column_mat().col(1);
        L_miller_point = inv_miller[2] * lat_column_mat().col(2);

        // Get three vectors that connect the three points on the plane to each other. Any two of the following
        // vectors could be used for constructing the new lattice, but it's convenient to pick the two
        // most orthogonal vectors
        Eigen::Vector3d HK, KL, LH;
        Eigen::Vector3d tangles;

        HK = K_miller_point - H_miller_point;
        KL = L_miller_point - K_miller_point;
        LH = H_miller_point - L_miller_point;

        // John G 121212
        // The vectors that we got at this point are valid, but sometimes larger than they need to be.
        Eigen::Matrix3d templat;
        templat.col(0) = HK;
        templat.col(1) = KL;
        templat.col(2) = LH;

        // Find shortest vector
        int s = 0;
        for(int i = 1; i < 3; i++) {
          if(lat_column_mat().col(i).norm() < lat_column_mat().col(s).norm()) {
            s = i;
          }
        }

        // Try dividing by integers and see if they're still lattice combinations. If they are
        // shorten and continue
        for(int i = 0; i < 3; i++) {
          int maxval = round(templat.col(i).norm() / lat_column_mat().col(i).norm() + 1);
          for(int j = maxval; j > 1; j--) {
            Eigen::Vector3d shortened = inv_lat_column_mat() * templat.col(i) / j; // Shorten and convert to fractional
            bool combo = true;

            for(int k = 0; k < 3; k++) {
              if(!almost_zero(std::abs(round(shortened[k]) - shortened[k]))) {
                combo = false;
                break;
              }
            }

            if(combo) {
              templat.col(i) = lat_column_mat() * shortened;
              break;
            }
          }
        }

        HK = templat.col(0);
        KL = templat.col(1);
        LH = templat.col(2);

        // We select the two vectors that spawn the smallest area

        double HKKL, KLLH, LHHK;
        HKKL = HK.cross(KL).norm();
        KLLH = KL.cross(LH).norm();
        LHHK = LH.cross(HK).norm();

        if(HKKL <= KLLH && HKKL <= LHHK) {
          surface_cell.col(0) = HK;
          surface_cell.col(1) = KL;
        }

        else if(KLLH <= HKKL && KLLH <= LHHK) {
          surface_cell.col(0) = KL;
          surface_cell.col(1) = LH;
        }

        else {
          surface_cell.col(0) = LH;
          surface_cell.col(1) = HK;
        }
      }
      //\John G 121212
      // We now have lattice vectors a and b. The c vector can be any vector that is a linear combination
      // of the original primitive cell lattice vectors. Ideally the vector will be short and orthogonal
      // to the plane.
      Eigen::Vector3d normal;
      Eigen::Vector3i L_combination;
      int factor;

      normal = inv_lat_column_mat() * surface_cell.col(0).cross(surface_cell.col(1)); // 101112
      factor = 1;

      // Divide by largest value in normal vector. We'll do something similar to when finding the miller
      // indeces. We'll approximate the real normal vector with an integer normal vector, increasing
      // the "resolution" of the integer vector in steps, until the desired accuary is reached.
      if(fabs(normal[0]) >= fabs(normal[1]) && fabs(normal[0]) >= fabs(normal[2]) && normal[0] != 0) {
        normal = normal / normal[0];
      }

      else if(fabs(normal[1]) >= fabs(normal[2]) && fabs(normal[1]) >= fabs(normal[0]) && fabs(normal[1]) != 0) {
        normal = normal / normal[1];
      }

      else {
        normal = normal / normal[2];
      }

      Eigen::Vector3d tnormal;
      // orthoscore represents how close the linear combination is to the plane normal. 1 is perfect, 0 is stupid.
      double orthoscore = 1;
      double torthoscore = 0;

      tnormal = normal;

      double new_vol = 0;
      bool one_found = false;

      do {
        normal = tnormal * factor;

        // Get linear combinations by rounding the normal to integers
        L_combination[0] = int(floor(normal[0] + 0.5));
        L_combination[1] = int(floor(normal[1] + 0.5));
        L_combination[2] = int(floor(normal[2] + 0.5));

        // After getting the normal vector in terms of integers, make a linear combination of the initial
        // lattice vectors to get third new lattice vector
        surface_cell.col(2) = lat_column_mat() * L_combination.cast<double>();
        orthoscore = fabs(cos(CASM::angle(lat_column_mat() * normal, surface_cell.col(2))));
        // Only use new linear combination if it's more orthogonal than the previous one
        if(orthoscore > torthoscore + TOL) {
          torthoscore = orthoscore;

          last_surface_cell = surface_cell; // Remember currect generated cell, in case we can't find anything better later
          one_found = true;

          new_vol = fabs(surface_cell.col(2).dot(surface_cell.col(0).cross(surface_cell.col(1))));
        }

        factor++;

        if(factor == 100) {
          if(!one_found) {
            throw std::runtime_error("failed get_lattice_in_plane");
          }
          std::cerr << "Reached an outrageous size. Returning last generated cell" << std::endl;
          surface_cell = last_surface_cell;
          break;
        }

      }
      while(new_vol / this->volume() < max_vol && orthoscore < 1);    // John G 121030

      Lattice surface_lat(surface_cell.col(0), surface_cell.col(1), surface_cell.col(2));
      surface_lat.make_right_handed();

      return surface_lat;
    }

    /// \brief Compare two Lattice
    ///
    /// - First compares is_niggli(*this, TOL) with is_niggli(RHS, TOL)
    /// - Then compares via standard_orientation_compare(this->lat_column_mat(), RHS.lat_column_mat(), TOL)
    bool Lattice::operator<(const Lattice &RHS) const {
      bool A_is_niggli = is_niggli(*this, TOL);
      bool B_is_niggli = is_niggli(RHS, TOL);
      if(A_is_niggli != B_is_niggli) {
        return B_is_niggli;
      }

      return standard_orientation_compare(lat_column_mat(), RHS.lat_column_mat(), TOL);
    }

    /// Are lattice vectors identical for two lattices
    bool Lattice::_eq(const Lattice &RHS) const {
      return almost_equal(RHS.lat_column_mat(), lat_column_mat());
    }

    bool Lattice::is_right_handed() const {
      if(this->volume() < 0) {
        return false;
      }
      else {
        return true;
      }
    }

    ///\brief returns Lattice that is smallest possible superlattice of both input Lattice
    ///
    //*******************************************************************************************
    //
    //  Finds "superduper" Lattice L_{sd} (represented as a matrix with lattice vectors as its columns
    //  such that L_{sd} satisfies
    //        L_{sd} = L_1 * N_1 = L_2 * N_2,     (*1*)
    //  where N_1 and N_2 are integer matrices such that Eq.(*1*) is satisfied and det(N_1) and det(N_2) are minimized.
    //
    //  It is assumed that L_1 = L * M_1 and L_2 = L * M_2  (i.e., L_1 and L_2 are superlattices of PRIM lattice L having
    //  integer transformation matrices M_1 and M_2, respectively).
    //
    //  Algorithm proceeds by noting inv(L_2)*L_1 = N_2*inv(N_1) = inv(M_2)*inv(M_1) = A/n, where A is an integer matrix and n = det(M_2).
    //  Assuming that 'n' is small (n<10000), we can attempt to find 'n' and 'A'.
    //
    //  Solution: N_2 = A*N_1/n, s.t. det(N_1) is minimized and N_2 is integer
    //
    //  To minimize det(N_1), find smith normal form A = U*S*V, where U and V have det(1), S is diagonal,
    //  and all entries are integer.  Then choose N_1 = inv(V)*R, where R is a diagonal integer matrix with entries
    //             R(i,i)=n/gcf(S(i,i),n)
    //  The resulting solution will have det(M_1*N_1)>=lcm(det(M_1),det(M_2))
    //
    //*******************************************************************************************
    Lattice make_superduperlattice(const Lattice &lat1, const Lattice &lat2) {

      Eigen::Matrix3d dA(lat2.inv_lat_column_mat() * lat1.lat_column_mat());
      long N = 1, num, denom;
      for(Index i = 0; i < 3; i++) {
        for(Index j = 0; j < 3; j++) {
          nearest_rational_number(dA(i, j), num, denom);
          dA *= double(denom);
          N *= denom;
        }
      }

      Eigen::Matrix3l U, S, V;
      smith_normal_form(lround(dA), U, S, V);

      denom = N;

      // reuse matrix S for matrix 'R', as above
      for(Index i = 0; i < 3; i++) {
        S(i, i) = N / gcf(S(i, i), N);
      }
      // Matrix 'N_1', as above is now equal to inverse(V) * S

      Lattice tlat(lat1.lat_column_mat() * (inverse(V) * S).cast<double>());
      return tlat.reduced_cell();
    }

    //*******************************************************************************************

    Lattice replace_vector(const Lattice &lat, const Eigen::Vector3d &new_vector, double tol) {

      // replace a lattice vector with translation
      Lattice new_lat{lat};
      double min_vol = std::abs(volume(new_lat));

      for(int i = 0; i < 3; i++) {

        Lattice tmp_lat = lat;
        tmp_lat[i] = new_vector;
        double vol = std::abs(volume(tmp_lat));

        if(vol < min_vol && vol > tol) {
          min_vol = vol;
          new_lat = tmp_lat;
        }
      }

      return new_lat;
    }

    /// Check if scel is a superlattice of unitcell unit and some integer transformation matrix T
    // First, find unit*Matrix3i approximation of 'scel', then check if reconstructing 'unit' from this approximation
    // results in residual vectors less than length 'tol'
    std::pair<bool, Eigen::Matrix3d> is_superlattice(const Lattice &scel, const Lattice &unit, double tol) {
      // check scel = unit*T, with integer T
      std::pair<bool, Eigen::Matrix3d> result(std::make_pair(true, unit.inv_lat_column_mat() * scel.lat_column_mat()));

      Eigen::Matrix3d diff = unit.lat_column_mat() - scel.lat_column_mat() * iround(result.second).cast<double>().inverse();

      result.first = almost_zero((diff.transpose() * diff).diagonal(), tol * tol) && !almost_zero(result.second, tol);

      return result;
    }


    Eigen::Matrix3l make_transformation_matrix_to_super(const Lattice &tiling_unit, const Lattice &superlattice, double tol) {

      Eigen::Matrix3d direct_transformation_matrix;
      bool is_integer_transformation;

      //TODO: convention is usually "prim" always goes first, but this is contradicte by is_superlattice. Which should change?
      std::tie(is_integer_transformation, direct_transformation_matrix) = is_superlattice(superlattice, tiling_unit, tol);

      if(!is_integer_transformation) {
        throw std::runtime_error("The provided tiling unit and superlattice are not related by a non-singular integer transformation.");
      }

      Eigen::Matrix3l rounded_transformation_matrix = round(direct_transformation_matrix).cast<long>();
      return rounded_transformation_matrix;
    }
  } // namespace xtal

} // namespace CASM
