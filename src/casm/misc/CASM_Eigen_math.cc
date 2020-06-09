#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/container/Counter.hh"
#include <iostream>

namespace CASM {
  /*
  Eigen::MatrixXd reduced_column_echelon(Eigen::Ref<const Eigen::MatrixXd> const &M, double _tol) {
    Eigen::MatrixXd R(M);
    Index col = 0;
    for(Index row = 0; row < R.rows(); ++row) {
      Index i = 0;
      for(i = col; i < R.cols(); ++i) {
        if(!almost_zero(R(row, i), _tol)) {
          if(i != col)
            R.col(i).swap(R.col(col));
          break;
        }
      }
      if(i == R.cols())
        continue;
      R.col(col) /= R(row, col);
      for(i = 0; i < R.cols(); ++i) {
        if(i != col)
          R.col(i) -= R(row, i) * R.col(col);
      }
      ++col;
    }
    return R.leftCols(col);
  }
  */

  //Given a string and the size expected in the array. Converts it to an Eigen::VectorXd
  Eigen::VectorXd eigen_vector_from_string(const std::string &tstr, const int &size) {
    std::stringstream tstream;
    tstream << tstr;
    Eigen::VectorXd tvec(size);
    int i = 0;
    while(!tstream.eof()) {
      if(i == (size)) {
        std::cerr << "WARNING in eigen_vector_from_string. You have more elements in your array, than the value of size you passed. Ignoring any extraneous elements" << std::endl;
        break;
      }
      tstream >> tvec(i);
      i++;
    }
    if(i < (size - 1)) {
      std::cerr << "ERROR in eigen_vector_from_string. Fewer elements in your array than the value of size you passed. Can't handle this. Quitting." << std::endl;
      exit(666);
    }
    return tvec;
  }



  namespace HungarianMethod_impl {

    //*******************************************************************************************
    /**
     * Implement Hungarian algorithm to find optimal mapping given a
     * cost matrix. Returns vector whose value at index l specifies the
     * assignment of relaxed basis $value to POS basis $l.
     *
     *
     */
    //*******************************************************************************************

    void reduce_cost(Eigen::MatrixXd &cost_matrix, double _infinity) {

      // cost matrix dimension
      int dim = cost_matrix.rows();

      // Step 1. For every row subtract the minimum element from every
      // entry in that row.
      for(int i = 0; i < dim; i++) {
        double row_min = cost_matrix.row(i).minCoeff();
        if(row_min > _infinity)
          continue;
        for(int j = 0; j < dim; j++) {
          cost_matrix(i, j) -= row_min;
        }
      }
    }

    //*******************************************************************************************
    void find_zeros(const Eigen::MatrixXd &cost_matrix, Eigen::MatrixXi &zero_marks, const double _tol) {
      // cost matrix dimension
      int dim = cost_matrix.rows();
      int star = -1;
      bool is_star;
      // Step 2. Find a zero in the matrix
      // find zeros and star
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          if(almost_zero(cost_matrix(i, j), _tol)) {
            is_star = false;
            for(int k = 0; k < dim; k++) {
              if(zero_marks(i, k) == star) {
                is_star = true;
              }
            }
            for(int l = 0; l < dim; l++) {
              if(zero_marks(l, j) == star) {
                is_star = true;
              }
            }
            if(!is_star) {
              zero_marks(i, j) = star;
            }
          }
        }
      }
    }

    //*******************************************************************************************

    bool check_assignment(const Eigen::MatrixXi &zero_marks, Eigen::VectorXi &col_covered) {


      // cover columns with starred zeros and check assignment. If N colmn
      // covered we are done else return false.
      for(int i = 0; i < zero_marks.rows(); i++) {
        for(int j = 0; j < zero_marks.cols(); j++) {
          if(zero_marks(i, j) == -1) {
            col_covered(j) = 1;
          }
        }
      }

      if(col_covered.sum() == zero_marks.rows()) {
        return true;
      }
      else return false;
    }

    //*******************************************************************************************

    int prime_zeros(const Eigen::MatrixXd &cost_matrix, Eigen::VectorXi &row_covered, Eigen::VectorXi &col_covered, Eigen::MatrixXi &zero_marks, double &min, Eigen::VectorXi &first_prime_zero, const double _tol, const double _infinity) {

      int prime = 1;
      int covered = 1;
      int uncovered = 0;

      bool DONE = false;


      // In step 4 we are looking for a noncovered zero.
      // If one is found it is primed. If no star is found
      // in its row, then step 5 is called. If a star is
      // found, the row containing the prime is covered,
      // and the column containing the star is uncovered.
      // This is continued until no uncovered zeros exist.
      while(!DONE) {
        DONE = true;
        for(int i = 0; i < cost_matrix.rows(); i++) {
          for(int j = 0; j < cost_matrix.cols(); j++) {
            if(almost_zero(cost_matrix(i, j), _tol) && col_covered(j) == 0 && row_covered(i) == 0) {
              //find an uncovered zero and prime it
              zero_marks(i, j) = prime;
              // Determine if there is a starred zero in the row, if so, cover row i
              // and uncover the column of the starred zero, k.
              bool is_star = false;
              for(int k = 0; k < zero_marks.cols(); k++) {
                if(zero_marks(i, k) == -1 && !is_star) {
                  row_covered(i) = covered;
                  col_covered(k) = uncovered;
                  is_star = true;
                }
              }

              // If no starred zero, save the last primed zero and proceed to STEP 5.
              if(!is_star) {
                first_prime_zero(0) = i;
                first_prime_zero(1) = j;
                //alternating_path(cost_matrix, zero_marks, row_covered, col_covered, first_prime_zero);
                return 5;
              }
              //else {
              //    row_covered(i) = 1;
              //    col_covered(j) = 0;
              //}


            }
          }
        }
        for(int i = 0; i < cost_matrix.rows(); i++) {
          for(int j = 0; j < cost_matrix.cols(); j++) {
            if(row_covered(i) == 0 && col_covered(j) == 0 && almost_zero(cost_matrix(i, j), _tol)) {
              DONE = false;
            }
          }
        }
      }

      // loop through all uncovered elements
      // to find the minimum value remaining.
      // This will be used in step 6.
      //min = numeric_limits<double>::infinity();
      min = _infinity;
      int return_code(-1);
      for(int i = 0; i < cost_matrix.rows(); i++) {
        for(int j = 0; j < cost_matrix.cols(); j++) {
          if(row_covered(i) == 0 && col_covered(j) == 0 && cost_matrix(i, j) < min) {
            min = cost_matrix(i, j);
            return_code = 6;
          }
        }
      }
      //update_costs(cost_matrix, zero_marks, row_covered, col_covered, min);
      return return_code;
    }

    //*******************************************************************************************

    int alternating_path(const Eigen::MatrixXd &cost_matrix, const Eigen::VectorXi &first_prime_zero, Eigen::MatrixXi &zero_marks, Eigen::VectorXi &row_covered, Eigen::VectorXi &col_covered) {


      bool done = false;

      // Initialize vectors to hold the locations
      // of the alternating path. The largest path
      // is at most 2*dimension of cost matrix -1.
      Eigen::VectorXi row_path(2 * cost_matrix.rows() - 1);
      Eigen::VectorXi col_path(2 * cost_matrix.cols() - 1);
      //
      //fill the paths initially with -1's so they
      //are not mistaken for points in path.
      row_path.fill(-1);
      col_path.fill(-1);
      //
      // Keep track of how long the path is.
      int path_counter = 1;

      // Set the first prime zero in the path
      // as the location specified from step 4.
      row_path(0) = first_prime_zero(0);
      col_path(0) = first_prime_zero(1);

      // These bools help test if there
      // are stars found in the prime's column.
      // Or if there are any primes found in
      // star's column.
      bool star_found = false;
      bool prime_found = false;

      // Build alternating path beginning with:
      //      uncovered primed zero (Z0)
      // ---> starred zero in same col as Z0 (Z1)
      // ---> primed zero in same row as Z1 (Z2)
      // ---> continue until a starred zero cannot be found
      //
      // The path is stored in row_path and col_path
      while(!done) {

        for(int i = 0; i < zero_marks.rows(); i++) {
          if(zero_marks(i, col_path(path_counter - 1)) == -1 && !star_found)  {
            path_counter += 1;
            row_path(path_counter - 1) = i;
            col_path(path_counter - 1) = col_path(path_counter - 2);
            star_found = true;
            prime_found = false;
          }
          if(i == zero_marks.rows() - 1 && !star_found) {
            done = true;
          }

        }
        if(star_found && !done) {
          for(int j = 0; j < zero_marks.cols(); j++) {
            if(zero_marks(row_path(path_counter - 1), j) == 1 && !prime_found) {
              path_counter += 1;
              row_path(path_counter - 1) = row_path(path_counter - 2);
              col_path(path_counter - 1) = j;
              prime_found = true;
              star_found = false;
            }
          }
        }
      }




      // Unstar each starred zero and star each primed zero
      for(int i = 0; i < path_counter; i ++)
        // Star primed zeros. I.e. even elements in path.
        if(i % 2 == 0 && zero_marks(row_path(i), col_path(i)) == 1) {
          zero_marks(row_path(i), col_path(i)) = -1;
        }
      // Unstar all starred zeros. I.e. odd elements.
        else if(zero_marks(row_path(i), col_path(i)) == -1) {
          zero_marks(row_path(i), col_path(i)) = 0;
        }
      for(int i = 0; i < zero_marks.rows(); i++) {
        for(int j = 0; j < zero_marks.cols(); j ++) {
          if(zero_marks(i, j) == 1) {
            zero_marks(i, j) = 0;
          }
        }
      }

      // Uncover all lines
      row_covered.fill(0);
      col_covered.fill(0);
      return 3;
      // Step 3 should be run next by the while loop in main.
    }

    //*******************************************************************************************
    int update_costs(const Eigen::VectorXi &row_covered, const Eigen::VectorXi &col_covered, const double min, Eigen::MatrixXd &cost_matrix) {

      // Step 6. Add value from 4 to all covered rows, and subtract it from uncovered columns.
      // DO NOT alter stars, primes or covers.
      // Return to 3.

      for(int i = 0; i < cost_matrix.rows(); i++) {
        for(int j = 0; j < cost_matrix.cols(); j++) {
          if(row_covered(i) == 1) {
            cost_matrix(i, j) += min;
          }
          if(col_covered(j) == 0) {
            cost_matrix(i, j) -= min;
          }
        }
      }
      //prime_zeros(cost_matrix, row_covered, col_covered, zero_marks);
      return 4;
    }

    // *******************************************************************************************
    /* Hungarian Algorithm Routines
     * Step 1: reduce the rows by smallest element
     * Step 2: star zeros
     * Step 3: cover columns with starred zeros and check assignement
     *         if K columns covered DONE. Else goto 4.
     * Step 4: Find uncovered zero and prime. if no starred zero in row
     *         goto 5. Otherwise, cover row, uncover column with star zero.
     *         Continue until all zeros are covered. Store smalles
     *         uncovered value goto 6.
     * Step 5: Build alternating prime and star zeros. Goto 3.
     * Step 6: Add value from 4 to all covered rows, and subtract it
     *         from uncovered columns. DO NOT alter stars, primes or covers.
     *         Return to 3.
     *
     */
    // *******************************************************************************************
    void hungarian_method(const Eigen::MatrixXd &cost_matrix_arg, std::vector<Index> &optimal_assignment, const double _tol) {
      double _infinity = 1E10;
      Eigen::MatrixXd cost_matrix = cost_matrix_arg;

      // cost matrix dimension
      int dim = cost_matrix.rows();


      // initialize matrix to carry stars and primes.
      Eigen::MatrixXi zero_marks = Eigen::MatrixXi::Zero(dim, dim);
      // initialize vectors to track (un)covered rows or columns
      Eigen::VectorXi row_covered(dim);
      Eigen::VectorXi col_covered(dim);


      int uncovered = 0;

      // begin uncovered
      row_covered.fill(uncovered);
      col_covered.fill(uncovered);

      //initialized by step 4
      double min = 0;

      // Vector to track the first_prime_zero found in step 4 and used in step 5.
      Eigen::VectorXi first_prime_zero(2);


      // Calling munkres step 1 to
      // reduce the rows of the cost matrix
      // by the smallest element in each row
      reduce_cost(cost_matrix, _infinity);

      // Calling munkres step 2 to find
      // a zero in the reduced matrix and
      // if there is no starred zero in
      // its row or column, star it.
      // Repeat for all elements.
      //
      find_zeros(cost_matrix, zero_marks, _tol);


      // Set DONE to false in order to control the while loop
      // which runs the main portion of the algorithm.
      bool DONE = false;

      // loop_counter could be used to set a maxumum number of iterations.
      int loop_counter = 0;


      // Main control loop for hungarian algorithm.
      // Can only exit the while loop if munkres step 3
      // returns true which means that there exists
      // an optimal assignment.
      //
      // Begin with step three.
      int next_step = 3;

      while(!DONE) {
        if(next_step == 3) {
          // Step 3 returns a boolean. TRUE if the assignment has been found
          // FALSE, if step 4 is required.
          if(check_assignment(zero_marks, col_covered) == true) {
            DONE = true;
          }
          // Set max number of iterations if desired.
          else if(loop_counter > 15) {
            DONE = true;
          }
          else next_step = 4;
        }
        // Step 4 returns an int, either 5 or 6 (or -1 if failure detected)
        if(next_step == 4) {
          next_step = prime_zeros(cost_matrix, row_covered, col_covered, zero_marks, min, first_prime_zero, _tol, _infinity);
        }
        // Step 5 returns an int, either 3 or 4.
        if(next_step == 5) {
          next_step = alternating_path(cost_matrix, first_prime_zero, zero_marks, row_covered, col_covered);
        }
        // Step 6 returns an int, always 4.
        if(next_step == 6) {
          next_step = update_costs(row_covered, col_covered, min, cost_matrix);
        }
        if(next_step == -1) {
          optimal_assignment.clear();
          return;
        }
        // Use loop counter if desired.
        // loop_counter++;
      }

      // Once the main control loop finishes, an
      // optimal assignment has been found. It
      // is denoted by the locations of the
      // starred zeros in the zero_marks matrix.
      // optimal_assignment vector contains the
      // results of the assignment as follows:
      // the index corresponds to the index of
      // the atom in the ideal POSCAR, and the
      // value at an index corresponds to the
      // atom of the relaxed structure.
      // i.e. optimal_assignemnt(1) = 0
      // means that the atom with index '1' in
      // the ideal POSCAR is mapped onto the atom
      // with index '0' in the relaxed CONTCAR.
      optimal_assignment.assign(cost_matrix.rows(), -1);

      for(int i = 0; i < zero_marks.rows(); i++) {
        for(int j = 0; j < zero_marks.cols(); j++) {
          if(zero_marks(i, j) == -1) {
            if(!valid_index(optimal_assignment[i])) {
              optimal_assignment[i] = j;
            }
            else {
              optimal_assignment.clear();
              break;
            }
          }
        }
        if(optimal_assignment.size() == 0 || !valid_index(optimal_assignment[i])) {
          optimal_assignment.clear();
          break;
        }

      }
      //std::cout << cost_matrix << std::endl;
      return;
    }

  }

  //*******************************************************************************************

  double hungarian_method(const Eigen::MatrixXd &cost_matrix, std::vector<Index> &optimal_assignment, const double _tol) {
    // Run hungarian algorithm on original cost_matrix
    HungarianMethod_impl::hungarian_method(cost_matrix, optimal_assignment, _tol);

    double tot_cost = 0.0;
    // Find the costs associated with the correct assignments
    for(int i = 0; i < optimal_assignment.size(); i++) {
      tot_cost += cost_matrix(i, optimal_assignment[i]);
    }
    if(optimal_assignment.size() == 0)
      tot_cost = 1e20;
    return tot_cost;
  }

  //*******************************************************************************************
  /**
   * Given a matrix, original_mat, this function calculates
   * both the Hermitian and anti-Hermitian parts.
   *
   * @param[in] original_mat The matrix we want to get the
   *                         Hermitian and anti-Hermitian
   *                         parts of
   * @param[in,out] hermitian_mat The Hermitian matrix
   * @param[in,out] antihermitian_mat The anti-Hermitian matrix
   */
  //*******************************************************************************************
  //Added by Ivy 11/28/12
  void get_Hermitian(Eigen::MatrixXcd &original_mat, Eigen::MatrixXcd &hermitian_mat, Eigen::MatrixXcd &antihermitian_mat) {

    Eigen::MatrixXcd conj_trans = original_mat.conjugate().transpose();

    hermitian_mat = 0.5 * (original_mat + conj_trans);
    antihermitian_mat = 0.5 * (original_mat - conj_trans);
  }

  //*******************************************************************************************
  /**
   * Returns true if the specified matrix is Hermitian.
   *
   * Checks if the matrix is the same as the conjugate
   */
  //*******************************************************************************************
  //Added by Ivy 11/28/12
  bool is_Hermitian(Eigen::MatrixXcd &mat) {

    if((mat - mat.conjugate().transpose()).isZero(TOL)) {
      return true;
    }
    return false;
  }



  /// \brief Return the hermite normal form, M == H*V
  ///
  /// \param M a square integer matrix
  ///
  /// \throws std::runtime_error if M is not full rank
  /// \returns std::pair<H,V>
  ///
  /// H set to upper triangular matrix H with:
  /// - H(i,i) > 0
  /// - 0 <= H(r,i) < H(i,i) for r > i
  /// - H(r,i) == 0 for r < i
  /// V set to unimodular matrix V
  /// - V.determinant() = +/- 1
  ///
  std::pair<Eigen::MatrixXi, Eigen::MatrixXi> hermite_normal_form(const Eigen::MatrixXi &M) {

    if(M.rows() != M.cols()) {
      throw std::runtime_error(
        std::string("Error in hermite_normal_form(const Eigen::MatrixXi &M)\n") +
        "  M must be square.");
    }

    if(boost::math::iround(M.cast<double>().determinant()) == 0) {
      throw std::runtime_error(
        std::string("Error in hermite_normal_form(const Eigen::MatrixXi &M)\n") +
        "  M must be full rank.");
    }

    Eigen::MatrixXi V = Eigen::MatrixXi::Identity(M.rows(), M.cols());
    Eigen::MatrixXi H = M;



    int N = M.rows();
    int i, r, c;


    // get non-zero entries along H diagonal
    for(i = N - 1; i >= 0; i--) {

      // swap columns if necessary to get non-zero entry at H(i,i)
      if(H(i, i) == 0) {
        for(c = i - 1; c >= 0; c--) {
          if(H(i, c) != 0) {
            H.col(c).swap(H.col(i));
            V.row(c).swap(V.row(i));
            break;
          }
        }
      }

      // get positive entry at H(i,i)
      if(H(i, i) < 0) {
        H.col(i) *= -1;
        V.row(i) *= -1;
      }

      // do column operations to get zeros to left of H(i,i)
      for(c = i - 1; c >= 0; c--) {
        if(H(i, c) == 0)
          continue;

        while(H(i, c) != 0) {

          while(H(i, c) < 0) {
            H.col(c) += H.col(i);
            V.row(i) -= V.row(c);
          }

          while(H(i, i) <= H(i, c) && H(i, c) != 0) {
            H.col(c) -= H.col(i);
            V.row(i) += V.row(c);
          }

          while(H(i, c) < H(i, i) && H(i, c) != 0) {
            H.col(i) -= H.col(c);
            V.row(c) += V.row(i);
          }

        }
      }
    }

    // now we have M = H*V, where H is upper triangular
    // so use column operations to enforce 0 <= H(i,c) < H(i,i) for c > i

    // columns left to right
    for(c = 1; c < H.cols(); c++) {
      // rows above the diagonal to the top
      for(r = c - 1; r >= 0; r--) {
        while(H(r, c) < 0) {
          H.col(c) += H.col(r);
          V.row(r) -= V.row(c);
        }
        while(H(r, c) >= H(r, r)) {
          H.col(c) -= H.col(r);
          V.row(r) += V.row(c);
        }
      }
    }

    // now we have M = H*V
    return std::make_pair(H, V);
  }


  /// \brief Get angle, in radians, between two vectors on range [0,pi]
  double angle(const Eigen::Ref<const Eigen::Vector3d> &a, const Eigen::Ref<const Eigen::Vector3d> &b) {
    return acos(a.dot(b) / (a.norm() * b.norm()));
  }

  ///return signed angle, in radians, between -pi and pi that describe separation in direction of two vectors
  double signed_angle(const Eigen::Ref<const Eigen::Vector3d> &a,
                      const Eigen::Ref<const Eigen::Vector3d> &b,
                      const Eigen::Ref<const Eigen::Vector3d> &pos_ref) {
    if(pos_ref.dot(a.cross(b)) < 0) {
      return -angle(a, b);
    }
    else
      return angle(a, b);
  }

  /// \brief Round entries that are within tol of being integer to that integer value
  Eigen::MatrixXd pretty(const Eigen::MatrixXd &M, double tol) {
    Eigen::MatrixXd Mp(M);
    for(int i = 0; i < M.rows(); i++) {
      for(int j = 0; j < M.cols(); j++) {
        if(std::abs(std::round(M(i, j)) - M(i, j)) < tol) {
          Mp(i, j) = std::round(M(i, j));
        }
      }
    }
    return Mp;
  }

  Eigen::Matrix3d polar_decomposition(Eigen::Matrix3d const &mat) {
    return Eigen::SelfAdjointEigenSolver <Eigen::Matrix3d>(mat.transpose() * mat).operatorSqrt();
  }

  std::vector<Eigen::Matrix3i> _unimodular_matrices(bool positive, bool negative, int range) {
    std::vector<Eigen::Matrix3i> uni_det_mats;
    int totalmats = 3480;

    if(positive && negative) {
      totalmats = totalmats * 2;
    }

    uni_det_mats.reserve(totalmats);

    EigenCounter<Eigen::Matrix3i> transmat_count(Eigen::Matrix3i::Constant(-range), Eigen::Matrix3i::Constant(range), Eigen::Matrix3i::Constant(1));

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

  const std::vector<Eigen::Matrix3i> &positive_unimodular_matrices() {
    static std::vector<Eigen::Matrix3i> static_positive(_unimodular_matrices(true, false));
    return static_positive;
  }

  const std::vector<Eigen::Matrix3i> &negative_unimodular_matrices() {
    static std::vector<Eigen::Matrix3i> static_negative(_unimodular_matrices(false, true));
    return static_negative;
  }


}

