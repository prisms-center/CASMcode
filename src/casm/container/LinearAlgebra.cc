#include "casm/container/LinearAlgebra.hh"

#include <boost/math/special_functions/round.hpp>

namespace CASM {

  //John G 020613

  //*******************************************************************************************
  /**
   * Adds number of nearest neighbors together if they fall
   * within a certain given distance (range) and then
   * averages the distance.
   */
  //*******************************************************************************************
  Array<Array<double> > oneNN_blur(const Array<Array<double> > &oneNN, double range) {
    if(oneNN.size() != 2 || oneNN[0].size() != oneNN[1].size() || oneNN[0].size() < 2) {
      std::cerr << "ERROR in oneNN_blur" << std::endl;
      std::cerr << "Expected a double array with equally sized inner arrays of size greater than 1." << std::endl;
      exit(70);
    }

    Array<Array<double> > blurred_oneNN;
    blurred_oneNN.resize(2);

    blurred_oneNN[0].push_back(oneNN[0][0]);
    blurred_oneNN[1].push_back(oneNN[1][0]);

    for(Index i = 1; i < oneNN[0].size() - 1; i++) {
      double diff = std::abs(blurred_oneNN[1].back() - oneNN[1][i]);
      if(diff < range) {
        blurred_oneNN[1].back() = (blurred_oneNN[0].back() * blurred_oneNN[1].back() + oneNN[0][i] * oneNN[1][i]) / (blurred_oneNN[0].back() + oneNN[0][i]);
        blurred_oneNN[0].back() += oneNN[0][i];
      }

      else {
        blurred_oneNN[0].push_back(oneNN[0][i]);
        blurred_oneNN[1].push_back(oneNN[1][i]);
      }
    }
    return blurred_oneNN;
    //I have not tested whether the average lengths make sense!
  }

  //*******************************************************************************************
  /**
   * Goes through the full nearest neighbor table and blurs
   * each one
   */
  //*******************************************************************************************
  Array<Array<Array<double> > > NN_blur(const Array<Array<Array<double> > > &NN, double range) {
    Array<Array<Array<double> > > blurred_NN;
    for(Index i = 0; i < NN.size(); i++) {
      blurred_NN.push_back(oneNN_blur(NN[i], range));
    }
    return blurred_NN;
  }

  //\John G 020613

  //*******************************************************************************************
  /**
   * Adds number of nearest neighbors together if they fall
   * within a certain given distance (range) and then
   * averages the distance.
   */
  //*******************************************************************************************
  Array<Array<double> > NN_blur(const Array<Array<double> > &NN, double range) {
    if(NN.size() != 2 || NN[0].size() != NN[1].size() || NN[0].size() < 2) {
      std::cerr << "ERROR in NN_blur" << std::endl;
      std::cerr << "Expected a double array with equally sized inner arrays of size greater than 1." << std::endl;
      exit(70);
    }

    Array<Array<double> > blurred_NN;
    blurred_NN.resize(2);

    blurred_NN[0].push_back(NN[0][0]);
    blurred_NN[1].push_back(NN[1][0]);

    for(Index i = 0; i < NN[0].size() - 1; i++) {
      double diff = std::abs(NN[1][i] - NN[1][i + 1]);
      if(diff < range) {
        blurred_NN[1].back() = (blurred_NN[0].back() * blurred_NN[1].back() + NN[0][i + 1] * NN[1][i + 1]) / (blurred_NN[0].back() + NN[0][i + 1]);
        blurred_NN[0].back() += NN[0][i];
      }

      else {
        blurred_NN[0].push_back(NN[0][i + 1]);
        blurred_NN[1].push_back(NN[1][i + 1]);
      }
    }
    return blurred_NN;
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

  /// \brief Check if Eigen::MatrixXd is integer
  ///
  /// \returns true if every coefficient of M is integer, with precision specified by tol
  ///
  /// \param M Eigen::Matrix begin checked
  /// \param tol tolerance for checking
  ///
  /// For each coefficient, M(i,j), checks that \code std::abs(std::round(M(i, j)) - M(i, j)) > tol \endcode
  ///
  bool is_integer(const Eigen::MatrixXd& M, double tol) {
    for(int i = 0; i < M.rows(); i++) {
      for(int j = 0; j < M.cols(); j++) {
        if(std::abs(std::round(M(i, j)) - M(i, j)) > tol) {
          return false;
        }
      }
    }
    return true;
  }
  
  /// \brief Check if Eigen::MatrixXd is unimodular
  ///
  /// \returns true if M is unimodular, false otherwise
  ///
  /// \param M Eigen::MatrixXd to check
  /// \param tol tolerance
  ///
  /// Equivalent to \code almost_equal(std::abs(M.determinant()), 1.0, tol) \endcode
  ///
  bool is_unimodular(const Eigen::MatrixXd& M, double tol) {
    return almost_equal(std::abs(M.determinant()), 1.0, tol);
  }
  
  /// \brief Round Eigen::Matrix3d to Eigen::Matrix3i
  ///
  /// \returns an Eigen:Matrix3i
  ///
  /// \param M Eigen::Matrix3d to be rounded to integer
  ///
  /// For each coefficient, sets \code Mint(i,j) = boost::math::iround(Mdouble(i, j)) \endcode
  ///
  Eigen::Matrix3i iround(const Eigen::Matrix3d &M) {
    Eigen::Matrix3i Mint(M.rows(), M.cols());
    for(int i = 0; i < M.rows(); i++) {
      for(int j = 0; j < M.cols(); j++) {
        Mint(i, j) = boost::math::iround(M(i, j));
      }
    }
    return Mint;
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
  
  /// \brief Check if Eigen::MatrixXd is diagonal
  ///
  /// \returns true if only diagonal elements are non-zero, with precision specified by tol
  ///
  /// \param M Eigen::Matrix begin checked
  /// \param tol tolerance for checking
  ///
  bool is_diagonal(const Eigen::MatrixXd& M, double tol) {
    
    for( int i=0; i<M.rows(); i++) {
      for( int j=0; j<M.cols(); j++) {
        if( i != j) {
          if( std::abs(M(i,j)) > tol)
            return false;
        }
      }
    }
    return true;
  }
  
  /// \brief Check if Eigen::MatrixXi is diagonal
  ///
  /// \returns true if only diagonal elements are non-zero, with precision specified by tol
  ///
  /// \param M Eigen::Matrix begin checked
  ///
  bool is_diagonal(const Eigen::MatrixXi& M) {
    
    for( int i=0; i<M.rows(); i++) {
      for( int j=0; j<M.cols(); j++) {
        if( i != j) {
          if( M(i,j) != 0)
            return false;
        }
      }
    }
    return true;
  }
  
  /// \brief Return the minor of integer Matrix M element row, col
  ///
  /// \returns the minor of element row, col
  ///
  /// \param M matrix
  /// \param row,col element row and column
  ///
  /// The minor of element row, col of a matrix is the determinant of the submatrix of
  /// M which does not include any elements from row 'row' or column 'col'
  ///
  int minor(const Eigen::MatrixXi& M, int row, int col) {
    
    // create the submatrix of M which doesn't include any elements from M in 'row' or 'col'
    Eigen::MatrixXi subM(M.rows()-1, M.cols()-1);
    int _i, _j;
    for( int i=0; i<M.rows(); i++) {
      for( int j=0; j<M.cols(); j++) {
        if(i==row || j ==col)
          continue;
        (i<row) ? _i = i : _i = i-1;
        (j<col) ? _j = j : _j = j-1;
        subM(_i,_j) = M(i,j);
      }
    }
    
    // return determinant of subM
    return boost::math::iround(subM.cast<double>().determinant());
  }
  
  /// \brief Return cofactor matrix
  ///
  /// \returns cofactor matrix
  ///
  /// \param M matrix
  ///
  /// The cofactor matrix is \code C(i,j) = pow(-1,i+j)*minor(M,i,j) \endcode
  /// 
  Eigen::MatrixXi cofactor(const Eigen::MatrixXi& M) {
    
    Eigen::MatrixXi Mcof(M.rows(), M.cols());
    for( int i=0; i<Mcof.rows(); i++) {
      for( int j=0; j<Mcof.cols(); j++) {
        (((i+j) % 2) == 0) ? Mcof(i,j) = minor(M, i, j) : Mcof(i,j) = -minor(M, i, j);
      }
    }
    return Mcof;
  }
  
  /// \brief Return adjugate matrix
  ///
  /// \returns adjugate matrix
  ///
  /// \param M matrix
  ///
  /// The adjugate matrix is the transpose of the cofactor matrix.
  /// The adjugate matrix is related to the inverse of a matrix through
  /// \code M.inverse() == adjugate(M)/M.determinant() \endcode
  ///
  Eigen::MatrixXi adjugate(const Eigen::MatrixXi& M) {
    return cofactor(M).transpose();
  }
  
  namespace normal_form_impl {
  
    /// \brief Calculate greatest common factor of two integers, and bezout coefficients
    ///
    /// \returns greatest common factor
    ///
    /// \param i1,i2 two integers for which to find greatest common factor
    /// \param[out] p1, p2 bezout coefficients such that p1*i1 + p2*i2 = gcf(abs(i1),abs(i2));
    ///
    /// \see smith_normal_form
    ///
    int _extended_gcf(int i1, int i2, int &p1, int &p2) {
      
      int s1 = sgn(i1);
      int s2 = sgn(i2);

      int a = i1, b = i2;
      p1 = 1;
      p2 = 0;
      int tp1(0), tp2(1), quotient;

      i1 = std::abs(i1);
      i2 = std::abs(i2);

      while(i2 != 0) {

        quotient = i1 / i2;
        i1 = i1 % i2;
        std::swap(i1, i2);
        p1 = p1 - quotient * tp1;
        std::swap(p1, tp1);
        p2 = p2 - quotient * tp2;
        std::swap(p2, tp2);
      }

      p1 *= s1;
      p2 *= s2;
      return p1 * a + p2 * b;
    }
    
    /// \brief Return an elementary hermite transformation
    ///
    /// \returns an integer matrix, E, with determinant 1 such that E * [a;b] = [g;0]
    ///
    /// \param a,b integers 
    ///
    /// \see smith_normal_form
    ///
    Eigen::Matrix3i _elementary_hermite_op(int a, int b, int i, int j) {
      Eigen::Matrix3i tmat = Eigen::Matrix3i::Identity();
      int tgcf = _extended_gcf(a, b, tmat(i, i), tmat(i, j));
      if(!tgcf) {
        tmat = Eigen::Matrix3i::Identity();
        return tmat;
      }
      tmat(j, i) = -b / tgcf;
      tmat(j, j) = a / tgcf;
      return tmat;
    }
    
  }
  
  /// \brief Return the integer inverse matrix of an invertible integer matrix
  ///
  /// \returns a 3x3 integer matrix inverse
  ///
  /// \param M an 3x3 invertible integer matrix
  ///
  Eigen::Matrix3i inverse(const Eigen::Matrix3i& M) {
    Eigen::Matrix3i Minv = adjugate(M)/boost::math::iround(M.cast<double>().determinant());
    return Minv;
  }
  
  /// \brief Return the smith normal form, M == U*S*V
  ///
  /// \param M 3x3 integer matrix
  /// \param[out] U set to U: is integer, U.determinant() == 1
  /// \param[out] S set to S: is diagonal, integer, and nonnegative, S(i,i) divides S(i+1,i+1) for all i
  /// \param[out] V set to V: is integer, V.determinant() == 1
  ///
  /// Adapted from Matlab implementation written by John Gilbert (gilbert@parc.xerox.com):
  /// - http://www.mathworks.com/matlabcentral/newsreader/view_thread/13728
  ///
  void smith_normal_form(const Eigen::Matrix3i& M, 
                         Eigen::Matrix3i& U, 
                         Eigen::Matrix3i& S, 
                         Eigen::Matrix3i& V) {
    
    using namespace normal_form_impl;
    
    U = V = Eigen::Matrix3i::Identity();
    S = M;
    Eigen::Matrix3i tmat = U;

    int i, j;

    // Bidiagonalize S with elementary Hermite transforms.
    for(j = 0; j < 3; j++) {
      //take j as column index and zero elements of j below the diagonal
      for(i = j + 1; i < 3; i++) {
        if(S(i, j) == 0) continue;
        tmat = _elementary_hermite_op(S(j, j), S(i, j), j, i);
        S = tmat * S;
        U = U * inverse(tmat);
      }

      //take j as row index and zero elements of j after superdiagonal
      for(i = j + 2; i < 3; i++) {
        if(S(j, i) == 0) continue;
        tmat = _elementary_hermite_op(S(j, j + 1), S(j, i), j + 1, i);
        S = S * tmat.transpose();
        V = inverse(tmat.transpose()) * V;
      }
    }

    // std::cout << "About to hunt zeros, SNF decomposition is currently:\n" << U << "\n\n" << S << "\n\n" << V << "\n\n";
    //while(!S.is_diagonal()) {
    while(!is_diagonal(S)) {
      //find off-diagonal element
      int b = 0;
      for(b = 0; b < 2; b++) {
        if(S(b, b + 1))
          break;
      }

      // To guarantee reduction in S(b,b), first make S(b,b) positive
      // and make S(b,b+1) nonnegative and less than S(b,b).
      if(S(b, b) < 0) {
        for(i = 0; i < 3; i++) {
          S(b, i) = -S(b, i);
          U(i, b) = -U(i, b);
        }
      }

      //std::cout << "S before q:\n" << S << '\n';

      if(S(b, b)) {
        int q = S(b, b + 1) / S(b, b);
        if(S(b, b + 1) % S(b, b) < 0) q -= 1;
        tmat = Eigen::Matrix3i::Identity();
        tmat(b + 1, b) = -q;
        S = S * tmat.transpose();
        V = inverse(tmat.transpose()) * V;
        //std::cout << "tmat for q:\n" << tmat << '\n';
        //std::cout << "S after q:\n" << S << '\n';
      }
      else {
        tmat = Eigen::Matrix3i::Identity();
        tmat(b, b) = 0;
        tmat(b, b + 1) = 1;
        tmat(b + 1, b + 1) = 0;
        tmat(b + 1, b) = 1;
        S = S * tmat.transpose();
        V = inverse(tmat.transpose()) * V;

      }

      if(!S(b, b + 1)) continue;
      tmat = _elementary_hermite_op(S(b, b), S(b, b + 1), b, b + 1);
      S = S * tmat.transpose();
      V = inverse(tmat.transpose()) * V;
      for(j = 0; j < 2; j++) {
        tmat = _elementary_hermite_op(S(j, j), S(j + 1, j), j, j + 1);
        S = tmat * S;
        U = U * inverse(tmat);
        if(j + 2 >= 3) continue;
        tmat = _elementary_hermite_op(S(j, j + 1), S(j, j + 2), j + 1, j + 2);
        S = S * tmat.transpose();
        V = inverse(tmat.transpose()) * V;
      }
    }

    //Now it's diagonal -- make it non-negative
    for(j = 0; j < 3; j++) {
      if(S(j, j) < 0) {
        for(i = 0; i < 3; i++) {
          S(j, i) = -S(j, i);
          U(i, j) = -U(i, j);
        }
      }
    }

    //sort diagonal elements
    for(i = 0; i < 3; i++) {
      for(j = i + 1; j < 3; j++) {
        if(S(i, i) > S(j, j)) {
          S.row(i).swap(S.row(j));
          S.col(i).swap(S.col(j));
          U.col(i).swap(U.col(j));
          V.row(i).swap(V.row(j));
        }
      }
    }
    //enforce divisibility condition:
    for(i = 0; i < 3; i++) {
      if(S(i, i) == 0) continue;
      for(j = i + 1; j < 3; j++) {
        if(S(j, j) % S(i, i) == 0) continue;
        //Replace S(i,i), S(j,j) by their gcd and lcm respectively.
        tmat = Eigen::Matrix3i::Identity();
        Eigen::Matrix3i tmat2(tmat);
        int a(S(i, i)), b(S(j, j)), c, d, tgcf;
        tgcf = _extended_gcf(a, b, c, d);
        tmat(i, i) = 1;
        tmat(i, j) = d;
        tmat(j, i) = -b / tgcf;
        tmat(j, j) = (a * c) / tgcf;

        tmat2(i, i) = c;
        tmat2(i, j) = 1;
        tmat2(j, i) = -(b * d) / tgcf;
        tmat2(j, j) = a / tgcf;
        S = tmat * S * tmat2.transpose();
        U = U * inverse(tmat);
        V = inverse(tmat2.transpose()) * V;
      }
    }
    return;
  }
  
  /// \brief Round entries that are within tol of being integer to that integer value
  Eigen::MatrixXd pretty(const Eigen::MatrixXd& M, double tol) {
    Eigen::MatrixXd Mp(M);
    for(int i = 0; i < M.rows(); i++) {
      for(int j = 0; j < M.cols(); j++) {
        if(std::abs(std::round(M(i, j)) - M(i, j)) < tol) {
          Mp(i,j) = std::round(M(i,j));
        }
      }
    }
    return Mp;
  }

}

