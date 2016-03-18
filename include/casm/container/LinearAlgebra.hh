//
//  LinearAlgebra.hh
//  CASM
//
//

#ifndef LINEARALGEBRA_HH
#define LINEARALGEBRA_HH

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <cassert>

#include "casm/external/Eigen/Dense"
#include "casm/external/Eigen/Eigenvalues"

#include "casm/container/Array.hh"
#include "casm/misc/CASM_math.hh"
namespace CASM {

  ///Introduces a tolerance for a given nearest neighbor table
  //Array<Array<double> > one_NN_blur(const Array<Array<double> > &one_NN, double range);
  //Array<Array<Array<double> > > NN_blur(const Array<Array<Array<double> > > &NN, double range);

  void get_Hermitian(Eigen::MatrixXcd &original_mat, Eigen::MatrixXcd &hermitian_mat, Eigen::MatrixXcd &antihermitian_mat); //Ivy
  bool is_Hermitian(Eigen::MatrixXcd &mat); //Ivy
  void poly_fit(Eigen::VectorXcd &xvec, Eigen::VectorXcd &yvec, Eigen::VectorXcd &coeffs, int degree);


  /*
  template <typename T>
  double triple_prod(const Vector3<T> &vec0, const Vector3<T> &vec1, const Vector3<T> &vec2) {
    return vec0.dot(vec1.cross(vec2));
  }


  template <typename T>
  Vector3<T> cross_prod(const Vector3<T> &vec0, const Vector3<T> &vec1) {
    Vector3<T> cross_vec;
    cross_vec[0] = vec0[1] * vec1[2] - vec0[2] * vec1[1];
    cross_vec[1] = vec0[2] * vec1[0] - vec0[0] * vec1[2];
    cross_vec[2] = vec0[0] * vec1[1] - vec0[1] * vec1[0];
    return cross_vec;
  }


  template <typename T>
  Matrix3<int> round(const Matrix3<T> &A) {
    Matrix3<int> iA;
    iA[0] = round(A[0]);
    iA[1] = round(A[1]);
    iA[2] = round(A[2]);
    iA[3] = round(A[3]);
    iA[4] = round(A[4]);
    iA[5] = round(A[5]);
    iA[6] = round(A[6]);
    iA[7] = round(A[7]);
    iA[8] = round(A[8]);
    return iA;
  }

  template <typename T>
  Vector3<int> round(const Vector3<T> &A) {
    Vector3<int> iA;
    iA[0] = round(A[0]);
    iA[1] = round(A[1]);
    iA[2] = round(A[2]);
    return iA;
  }
  */
  /// \brief Check if Eigen::MatrixXd is integer
  bool is_integer(const Eigen::MatrixXd &M, double tol);

  /// \brief Check if Eigen::MatrixXd is unimodular
  bool is_unimodular(const Eigen::MatrixXd &M, double tol);

  /// \brief Round Eigen::MatrixXd to Eigen::MatrixXi
  Eigen::MatrixXi iround(const Eigen::MatrixXd &M);

  /// \brief Round Eigen::MatrixXd to Eigen::MatrixXl
  Eigen::MatrixXl lround(const Eigen::MatrixXd &M);

  /// \brief Return the hermite normal form, M == H*V
  std::pair<Eigen::MatrixXi, Eigen::MatrixXi> hermite_normal_form(const Eigen::MatrixXi &M);

  /// \brief Check if Eigen::MatrixXd is diagonal
  bool is_diagonal(const Eigen::MatrixXd &M, double tol);

  /// \brief Check if Eigen::MatrixXi is diagonal
  bool is_diagonal(const Eigen::MatrixXi &M);

  /// \brief Return the minor of integer Matrix M element row, col
  int matrix_minor(const Eigen::MatrixXi &M, int row, int col);

  /// \brief Return cofactor matrix
  Eigen::MatrixXi cofactor(const Eigen::MatrixXi &M);

  /// \brief Return adjugate matrix
  Eigen::MatrixXi adjugate(const Eigen::MatrixXi &M);

  /// \brief Return the integer inverse matrix of an invertible integer matrix
  Eigen::Matrix3i inverse(const Eigen::Matrix3i &M);

  /// \brief Return the smith normal form, M == U*S*V
  void smith_normal_form(const Eigen::Matrix3i &M,
                         Eigen::Matrix3i &U,
                         Eigen::Matrix3i &S,
                         Eigen::Matrix3i &V);

  /// \brief Round entries that are within tol of being integer to that integer value
  Eigen::MatrixXd pretty(const Eigen::MatrixXd &M, double tol);

}
#endif
