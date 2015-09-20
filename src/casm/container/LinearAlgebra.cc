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

}

