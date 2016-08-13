#include "casm/misc/CASM_math.hh"

namespace CASM {
  //*******************************************************************************************

  int round(double val) {
    return int(val < 0 ? floor(val + 0.5) : ceil(val - 0.5));
  };

  //*******************************************************************************************

  double ran0(int &idum) {
    int IA = 16807;
    int IM = 2147483647;
    int IQ = 127773;
    int IR = 2836;
    int MASK = 123459876;
    double AM = 1.0 / IM;

    //minimal random number generator of Park and Miller
    //returns uniform random deviate between 0.0 and 1.0
    //set or rest idum to any  integer value (except the
    //unlikely value MASK) to initialize the sequence: idum must
    //not be altered between calls for successive deviates
    //in a sequence

    int k;
    idum = idum ^ MASK; //XOR the two integers
    k = idum / IQ;
    idum = IA * (idum - k * IQ) - IR * k;
    if(idum < 0) {
      idum = idum + IM;
    }
    double ran = AM * idum;
    idum = idum ^ MASK;
    return ran;

  }

  //*******************************************************************************************

  int dl_string_dist(const std::string &a, const std::string &b) {
    // "infinite" distance is just the max possible distance
    int max_val = a.size() + b.size();

    // make and initialize the character array indices
    std::map<char, int> DA;

    // make the distance matrix H
    Eigen::MatrixXi H(a.size() + 2, b.size() + 2);

    // initialize the left and top edges of H
    H(0, 0) = max_val;
    for(int i = 0; i <= a.size(); ++i) {
      DA[a[i]] = 0;
      H(i + 1, 0) = max_val;
      H(i + 1, 1) = i;
    }
    for(int j = 0; j <= b.size(); ++j) {
      DA[b[j]] = 0;
      H(0, j + 1) = max_val;
      H(1, j + 1) = j;
    }

    // fill in the distance matrix H
    // look at each character in a
    for(int i = 1; i <= a.size(); ++i) {
      int DB = 0;
      // look at each character in b
      for(int j = 1; j <= b.size(); ++j) {
        int i1 = DA[b[j - 1]];
        int j1 = DB;
        int cost;
        if(a[i - 1] == b[j - 1]) {
          cost = 0;
          DB   = j;
        }
        else
          cost = 1;
        H(i + 1, j + 1) = min(min(H(i, j) + cost, // substitution
                                  H(i + 1, j) + 1),  // insertion
                              min(H(i, j + 1) + 1,  // deletion
                                  H(i1, j1) + (i - i1 - 1) + 1 + (j - j1 - 1)));
      }
      DA[a[i - 1]] = i;
    }
    return H(a.size() + 1, b.size() + 1);
  }

  //*******************************************************************************************
  /// Find greatest common factor
  int gcf(int i1, int i2) {
    i1 = std::abs(i1);
    i2 = std::abs(i2);
    while(i1 != i2 && i1 != 1 && i2 != 1) {
      if(i1 < i2) {
        i2 -= i1;
      }
      else {
        i1 -= i2;
      }
    }
    if(i1 == 1) {
      return i1;
    }
    else {
      return i2;
    }

  }

  //*******************************************************************************************
  /// Find least common multiple
  int lcm(int i1, int i2) {
    return std::abs(i1 * (i2 / gcf(i1, i2)));
  }

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

  //*******************************************************************************************
  // This evaluates Gaussians using the formula:
  // f(x) = a*e^(-(x-b)^2/(c^2))
  //*******************************************************************************************

  double gaussian(double a, double x, double b, double c) {
    return a * exp(-((x - b) * (x - b)) / (c * c));
  }

  //*******************************************************************************************
  // This calculates Gaussian moments given by the integral:
  // m = \int_{\infty}^{\infty} dx x^pow*exp[-x^2/(2*sigma^2)]/(\sqrt(2*\pi)*sigma)
  //*******************************************************************************************

  double gaussian_moment(int expon, double sigma) {

    if(expon % 2) return 0.0;

    double m = pow(sigma, expon);

    expon -= 1;
    while(expon - 2 > 0) {
      m *= double(expon);
      expon -= 2;
    }
    return m;
  }

  //*******************************************************************************************
  // This calculates Gaussian moments given by the integral:
  // m = \int_{\infty}^{\infty} dx x^pow*exp[-(x-x0)^2/(2*sigma^2)]/(\sqrt(2*\pi)*sigma)
  //*******************************************************************************************

  double gaussian_moment(int expon, double sigma, double x0) {
    double m = 0;
    for(int i = 0; i <= expon; i++) {
      m += nchoosek(expon, i) * gaussian_moment(i, sigma) * pow(x0, expon - i);
    }

    return m;
  }


  //*******************************************************************************************
  // finds rational number that approximates 'val' to within 'tol' --
  //      --  almost_int(double(denominator)*val, tol) == true OR almost_int(val/double(numerator), tol) == true OR
  //      --  denominator is always positive
  //      --  sign(numerator)=sign(val)
  //      --  if(almost_zero(val,tol)) --> numerator = 0, denominator = 1
  void nearest_rational_number(double val, long &numerator, long &denominator, double tol) {
    if(almost_zero(val, tol)) {
      numerator = 0;
      denominator = 1;
      return;
    }
    long sgn(val < 0 ? -1 : 1);
    val = std::abs(val);
    double tdenom, tnum;
    long lim(max(long(100), long(1 / (10 * tol))));
    for(long i = 1; i < lim + 1; i++) {
      tdenom = double(i) / val;
      tnum = val / double(i);
      if(tdenom > 1 && almost_zero(tdenom - round(tdenom), tol)) {
        numerator = sgn * i;
        denominator = round(tdenom);
        return;
      }
      else if(tnum > 1 && almost_zero(tnum - round(tnum), tol)) {
        denominator = i;
        numerator = sgn * round(tnum);
        return;
      }
    }


  }

  //*******************************************************************************************
  /* Finds best irrational number approximation of double 'val' and
   * returns tex-formated string that contains irrational approximation
   * searches numbers of the form (x/y)^(1/z), where x and y range from 1 to 'lim'
   * z ranges from 1 to 'max_pow'
   */
  //*******************************************************************************************

  std::string irrational_to_tex_string(double val, int lim, int max_pow) {
    std::stringstream tstr;
    if(almost_zero(round(val) - val)) {
      tstr << round(val);
      return tstr.str();
    }
    if(val < 0) {
      tstr << '-';
      val = std::abs(val);
    }
    double tval(val), tdenom, tnum;
    int idenom, inum;
    for(int ipow = 1; ipow < max_pow + 1; ipow++) {
      for(int i = 1; i < lim + 1; i++) {
        tdenom = double(i) / tval;
        tnum = tval / double(i);
        if(tdenom > 1 && almost_zero(std::abs(tdenom - round(tdenom)))) {
          inum = i;
          idenom = round(tdenom);
        }
        else if(tnum > 1 && almost_zero(std::abs(tnum - round(tnum)))) {
          idenom = i;
          inum = round(tnum);
        }
        else {
          continue;
        }

        if(ipow == 1) {
          tstr << inum << '/' << idenom;
          return tstr.str();
        }
        if(ipow == 2) {
          tstr << "\\sqrt{" << inum;
          if(idenom != 1)
            tstr << '/' << idenom;
          tstr << '}';
          return tstr.str();
        }
        else {
          tstr << '(' << inum;
          if(idenom != 1)
            tstr << '/' << idenom;
          tstr << ")^{1/" << ipow << "}";
          return tstr.str();
        }
      }
      tval *= val;
    }
    tstr << val;
    return tstr.str();
  }


  //*******************************************************************************************
  //John G 010413
  int mod(int a, int b) {
    if(b < 0)
      return mod(-a, -b);

    int ret = a % b;
    if(ret < 0)
      ret += b;
    return ret;
  }

  //*******************************************************************************************

  double cuberoot(double number) {
    if(number < 0) {
      return -pow(-number, 1.0 / 3);
    }

    else {
      return pow(number, 1.0 / 3);
    }
  }

  //*******************************************************************************************
  /**
   *
   */
  //*******************************************************************************************

  void poly_fit(Eigen::VectorXcd &xvec, Eigen::VectorXcd &yvec, Eigen::VectorXcd &coeffs, int degree) {

    // Check that the dimensions are correct
    if(xvec.rows() != yvec.rows()) {
      std::cout << "******************************************\n"
                << "ERROR in poly_fit: Dimensions of xvec and \n"
                << "yvec do not match up! Cannot perform a \n"
                << "polynomial fit.\n"
                << "******************************************\n";
      exit(1);
    }

    Eigen::MatrixXcd polymat = Eigen::MatrixXcd::Ones(xvec.rows(), degree + 1);

    // Construct polymatrix
    for(int d = degree; d > 0; d--) {

      polymat.col(degree - d) = xvec;

      for(int m = 0; m < d - 1; m++) {
        polymat.col(degree - d) = polymat.col(degree - d).cwiseProduct(xvec);
      }
    }

    // SVD least squares fit
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd_solver(polymat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    coeffs = svd_solver.solve(yvec);

    return;

  }


  //************************************************************

  // returns 'i' if 'input' is equivalent to 'unique[i]', w.r.t. permutation of the equivalent elements of 'input'.
  // equivalent elements are specified by 'ind_equiv'
  // if 'input' specifies a new combination of integers, unique.size() is returned
  Index which_unique_combination(const Array<Index> &input, const Array<Array<Index> > &unique, const Array<Array<Index> > &ind_equiv) {
    Index tval, tcount;
    Index i, j, k, l;
    for(i = 0; i < unique.size(); i++) {
      //Is input equivalent to unique[i] w.r.t. ind_equiv?
      //Loop over groups of equivalent indices
      if(unique[i].size() != input.size()) continue;
      for(j = 0; j < ind_equiv.size(); j++) {

        //Loop over indices that are equivalent
        for(k = 0; k < ind_equiv[j].size(); k++) {
          tval = input[ind_equiv[j][k]];
          tcount = 0;
          for(l = 0; l < ind_equiv[j].size(); l++)
            tcount += int(unique[i][ind_equiv[j][l]] == tval) - int(input[ind_equiv[j][l]] == tval);

          if(tcount != 0)
            break;
        }
        if(k != ind_equiv[j].size())
          break;
      }
      if(j == ind_equiv.size())
        return i;
    }
    return i;
  }

  //****************************************

  Index which_unique_combination(const Array<Index> &input, const Array<Array<Index> > &unique) {
    Index tval, tcount;
    Index i, j, l;
    for(i = 0; i < unique.size(); i++) {
      //Is input equivalent to unique[i] w.r.t. ind_equiv?
      //Loop over groups of equivalent indices
      if(unique[i].size() != input.size()) continue;
      for(j = 0; j < input.size(); j++) {
        tval = input[j];
        tcount = 0;
        for(l = 0; l < input.size(); l++)
          tcount += int(unique[i][l] == tval) - int(input[l] == tval);

        if(tcount != 0)
          break;
      }

      if(j == input.size())
        return i;
    }
    return i;

  }

  //************************************************************
  /// Find least common multiple
  int lcm(const Array<int> &series) {
    if(!series.size()) return 0;
    int lcm_val(series[0]);
    for(Index i = 1; i < series.size(); i++)
      lcm_val = lcm(lcm_val, series[i]);

    return lcm_val;
  }

  //************************************************************

  ReturnArray<Array<int> > get_prime_factors(int target) {
    Array<Array<int> > factors_array;
    Array<int> factor_list;

    if(target <= 1) {
      std::cerr << "WARNING in CASM_Global_Definitions::get_prime_factors" << std::endl;
      std::cerr << "You're asking for prime factors of " << target << ". Returning empty array." << std::endl << std::endl;
      return factors_array;
    }
    int factor = 2;

    while(target != 1) {
      while(target % factor == 0) {
        factor_list.push_back(factor);
        target = target / factor;
      }

      factor++;

      if(factor_list.size() > 0) {
        factors_array.push_back(factor_list);
      }
      factor_list.clear();
    }

    return factors_array;
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

    //*******************************************************************************************

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
      double min;

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
}

