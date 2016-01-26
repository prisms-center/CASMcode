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
  Array<Array<double> > one_NN_blur(const Array<Array<double> > &one_NN, double range);
  Array<Array<Array<double> > > NN_blur(const Array<Array<Array<double> > > &NN, double range);

  void get_Hermitian(Eigen::MatrixXcd &original_mat, Eigen::MatrixXcd &hermitian_mat, Eigen::MatrixXcd &antihermitian_mat); //Ivy
  bool is_Hermitian(Eigen::MatrixXcd &mat); //Ivy
  void poly_fit(Eigen::VectorXcd &xvec, Eigen::VectorXcd &yvec, Eigen::VectorXcd &coeffs, int degree);

  template <typename T, Index N> class Vector;
  template <typename T> class Vector3;
  template <typename T, Index N> class Matrix;
  template <typename T> class Matrix3;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template <typename T, Index N>
  class Matrix {

  private:
    T Vals[N * N];

    //Static members hold important matrices for a particular class.
    //These matrices can be returned by reference for speedy initialization
    //ones_mat is a matrix of all ones, zeros_mat is a matrix of all zeros and
    //identity_mat is the NxN identity matrix.
    static Matrix ones_mat;
    static Matrix zeros_mat;
    static Matrix identity_mat;
    static bool static_is_initialized;
    static void initialize_static_members();

  public:
    typedef T value_type;
    typedef Index size_type;

    Matrix() { };

    Matrix(T init_val) {
      for(Index i = 0; i < size(); i++) at(i) = init_val;
    }

    // Accessors
    T &operator()(Index i, Index j);
    const T &operator()(Index i, Index j) const;

    const T &operator[](Index i) const;
    T &operator[](Index i);
    T &at(Index i, Index j);
    T &at(Index i);
    const T &at(Index i, Index j) const;
    const T &at(Index i) const;

    void fill(T fill_val);

    // sizes
    Index num_rows() const;
    Index num_cols() const;
    Index size() const;


    // Essential Matrix Operations
    T determinant() const;
    T trace() const;
    Matrix inverse() const;
    Matrix transpose() const;

    // Static Matrix generating routines;
    static const Matrix &identity() {
      initialize_static_members();
      return identity_mat;
    }
    static const Matrix &ones() {
      initialize_static_members();
      return ones_mat;
    }
    static const Matrix &zeros() {
      initialize_static_members();
      return zeros_mat;
    }


    // Algebraic Operators
    Matrix operator+(const Matrix &RHS) const;
    Matrix operator-(const Matrix &RHS) const;

    Matrix  operator*(const T &RHS) const;
    Matrix operator/(const T &RHS) const;


    // Self-modifying Algebraic Operatos
    Matrix &operator+=(const Matrix &RHS);
    Matrix &operator-=(const Matrix &RHS);
    Matrix &operator*=(const Matrix &RHS);
    Matrix &operator*=(const T &RHS);
    Matrix &operator/=(const T &RHS);


    // Comparison Operators
    bool operator==(const Matrix &RHS) const;
    bool is_equal(const Matrix &RHS, double tol_val = TOL) const;
    bool is_unitary(double tol = TOL) const;
    bool is_identity(double tol = TOL) const;
    bool is_integer(double tol = TOL) const;
    bool is_zero(double zero_tol = TOL) const;

    Vector<T, N> eigen(Matrix &eigenvecs);  //eigenvectors are columns of Matrix eigenvectors
    Vector<T, N> eigen_values();
    void eigen_vectors(Matrix &eigenvecs);  //eigenvectors are columns of Matrix eigenvectors

    bool SVD(Matrix &U, Matrix &S, Matrix &V);
    void LU(Matrix &L, Matrix &U);
    Matrix minors();

    Matrix gauss_elimination();

    Index column_rank();
    Index row_rank();

    Matrix unitary_transform(const Matrix &TransformMatrix);

  };

  template <typename T, Index N>
  jsonParser &to_json(const Matrix<T, N> &value, jsonParser &json);

  template <typename T, Index N>
  void from_json(const Matrix<T, N> &value, const jsonParser &json);

  template <typename T, Index N>
  bool Matrix<T, N>::static_is_initialized = false;

  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::ones_mat(1);

  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::zeros_mat(0);

  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::identity_mat(0);

  template <typename T, Index N>
  std::ostream &operator<<(std::ostream &stream, const Matrix<T, N> &out_matrix);

  template <typename T, Index N>
  std::istream &operator>>(std::istream &stream, Matrix<T, N> &in_matrix);

  template <typename T, class U, Index N>
  Matrix<U, N> operator*(const Matrix<T, N> &LHS, const Matrix<U, N> &RHS);

  template <typename T, Index N>
  Matrix<T, N> operator*(const Matrix<T, N> &LHS, const Matrix<int, N> &RHS);

  template <typename T, Index N>
  Matrix<T, N> operator*(T LHS, const Matrix<T, N> &RHS);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  template < typename T >
  class Matrix3 : public Matrix<T, 3> {

  public:
    using Matrix<T, 3>::at;
    using Matrix<T, 3>::operator();
    using Matrix<T, 3>::operator==;
    using Matrix<T, 3>::num_cols;
    using Matrix<T, 3>::num_rows;
    using Matrix<T, 3>::is_integer;
    using Matrix<T, 3>::is_zero;

    static Matrix3 ones() {
      return Matrix3<T>(1);
    }
    static Matrix3 zeros() {
      return Matrix3<T>(0);
    }
    static Matrix3 identity() {
      Matrix3<T> tMat(0);
      tMat[0] = tMat[4] = tMat[8] = 1;
      return tMat;
    }

    Matrix3() { }
    Matrix3(T init_val) {
      for(int i = 0; i < 9; i++) at(i) = init_val;
    }
    Matrix3(const Matrix<T, 3> &init_mat) : Matrix<T, 3>(init_mat) { }

    template<int dim1, int dim2, int flag1>
    Matrix3(const Eigen::Matrix<T, dim1, dim2, flag1> &_mat);

    T determinant() const;
    T trace() const;

    //Square of Frobenius norm
    double square_norm()const;

    //Frobenius norm
    double norm()const;

    Matrix3 inverse() const;
    Matrix3 adjugate() const;
    Matrix3 transpose() const;

    void smith_normal_form(Matrix3<T> &U, Matrix3<T> &S, Matrix3<T> &V) const;

    //void hermite_normal_form(Matrix3<T> &U, Matrix3<T> &S) const;
    //void hermite_normal_form(Matrix3<T> &S) const;

    bool is_diagonal(double compare_tol = TOL) const;
    bool is_unitary(double compare_tol = TOL) const;
    //bool is_almost_hermite_normal_form(double compare_tol = TOL) const;

    Matrix3 &operator=(const Eigen::Matrix<T, 3, 3> &RHS);

    Matrix3 &operator-=(const Matrix3<T> &RHS);
    Matrix3 &operator+=(const Matrix3<T> &RHS);
    Matrix3 &operator*=(T RHS);

    Matrix3 operator-(const Matrix3<T> &RHS) const;
    Matrix3 operator+(const Matrix3<T> &RHS) const;

    operator Eigen::Matrix<T, 3, 3>() const;
    operator Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>() const;

    Vector3< std::complex<double> > eigen() const;
    Vector3< std::complex<double> > eigen(Vector3<Vector3< std::complex<double> > > &eigen_vecs) const;

    void swap_rows(Index i, Index j);
    void swap_cols(Index i, Index j);

  private:
    static Matrix3 elementary_hermite_op(T a, T b, T i, T j);
  };

  template <typename T, class U>
  Matrix3<U> operator*(const Matrix3<T> &LHS, const Matrix3<U> &RHS);

  template <typename T>
  Matrix3<T> operator*(const Matrix3<T> &LHS, const Matrix3<int> &RHS);

  template <typename T>
  Matrix3<T> operator*(T LHS, const Matrix3<T> &RHS);



  template <typename T, Index N>
  class Vector {
  protected:
    T Vals[N];
  public:
    typedef T value_type;
    typedef Index size_type;

    Vector() {};

    Vector(const T &init_val);
    T &operator[](Index ind) {
      assert(ind >= 0 && ind < size());
      return Vals[ind];
    }
    const T &operator[](Index ind) const {
      assert(ind >= 0 && ind < size());
      return Vals[ind];
    }
    const T &at(Index ind) const {
      assert(ind >= 0 && ind < size());
      return Vals[ind];
    }
    T &at(Index ind) {
      assert(ind >= 0 && ind < size());
      return Vals[ind];
    }
    Index size() const {
      return N;
    }

    // Vector(); //Do we need default constructor to initialize Vals?
    //    T operator*(const Vector &RHS) const;  // Dot product of two vectors; e.g., T val; val=vec1*vec2;

    //Vector& operator=(const T &init_val);

    Vector operator+(const Vector &RHS) const;
    Vector operator-(const Vector &RHS) const;
    Vector operator/(const T &divisor) const;

    Vector operator-() const;

    Vector &operator+=(const Vector &RHS);
    Vector &operator-=(const Vector &RHS);

    bool operator==(const Vector &RHS) const;

    bool is_equal(const Vector &RHS, double eq_tol = TOL) const;
    bool is_zero(double zero_tol = TOL) const;

    void scale(double scale_factor);
    double length() const;

  };

  template <typename T, Index N>
  jsonParser &to_json(const Vector<T, N> &value, jsonParser &json);

  template <typename T, Index N>
  void from_json(const Vector<T, N> &value, const jsonParser &json);

  template <typename T, class U, Index N>
  Vector<U, N> operator*(const Matrix<T, N> &LHS, const Vector<U, N> &RHS);

  template <typename T, Index N>
  Vector<T, N> operator*(const Matrix<T, N> &LHS, const Vector<int, N> &RHS);

  template <typename T, class U, Index N>
  Vector<T, N> operator*(const U &LHS, const Vector<T, N> &RHS);

  template <typename T, Index N>
  std::ostream &operator<<(std::ostream &stream, const Vector<T, N> &vec_out);


  template <typename T, Index N>
  std::istream &operator>>(std::istream &stream, Vector<T, N> &vec_in);



  //*******************************************************************************************

  template < typename T >
  class Vector3 : public Vector<T, 3> {

  public:
    using Vector<T, 3>::at;
    using Vector<T, 3>::scale;


    Vector3() : Vector<T, 3>() {}
    Vector3(const T &init_val) : Vector<T, 3>(init_val) {}
    Vector3(const T &elem0, const T &elem1, const T &elem2) {
      at(0) = elem0;
      at(1) = elem1;
      at(2) = elem2;
    }

    template<typename Derived>
    Vector3(const Eigen::MatrixBase<Derived> &vec);

    bool is_zero(double zero_tol = TOL) const;
    Vector3 cross(const Vector3 &RHS) const;

    template<class U>
    U dot(const Vector3<U> &RHS) const;

    template<class U>
    Vector3 operator*(const U &RHS) const;

    template<class U>
    operator Vector3<U>() const;

    void normalize();

    double norm() const;
    //John G 120921
    ///Calculate angle between two vectors
    double get_angle(const Vector3 &RHS) const;
    T max() const;
    T min() const;

    void scale_by_max_mag();

    ///Take a vector of doubles, and multiply by some factor that turns it into a vector of integers (within a tolerance)
    Vector3< int > scale_to_int(double tolerance = TOL);

    //John G 121106
    double get_signed_angle(const Vector3 &RHS, const Vector3 &pos_ref) const;
    Matrix3<T> get_cross_product_mat() const;
    Matrix3<T> get_tensor_product(const Vector3<T> &tensvec) const;
    Matrix3<T> get_tensor_product() const;
    Matrix3<T> get_rotation_mat(double angle) const;
    //\John G

    void rand_ball(int &rand_seed);
    void rand(int &rand_seed);

    template<int dim1, int dim2, int flag1>
    void insert_as_col(Eigen::Matrix<T, dim1, dim2, flag1> &mat, int ind);

    template<int dim1, int dim2, int flag1>
    void insert_as_row(Eigen::Matrix<T, dim1, dim2, flag1> &mat, int ind);

    Vector3 &operator=(const Vector3<T> &RHS);

    Vector3 &operator/=(const T &RHS);
    Vector3 &operator+=(const Vector3 &RHS);
    Vector3 &operator-=(const Vector3 &RHS);

    template<class U>
    Vector3 &operator*=(const U &RHS);

    Vector3 operator/(const T &RHS) const;
    Vector3 operator+(const Vector3 &RHS) const;
    Vector3 operator-(const Vector3 &RHS) const;
    Vector3 operator-() const;

    template<int dim1, int dim2, int flag1>
    operator Eigen::Matrix<T, dim1, dim2, flag1>() const;
  };

  template <typename T, class U>
  Vector3<U> operator*(const Matrix3<T> &LHS, const Vector3<U> &RHS);

  template <typename T>
  Vector3<T> operator*(const Matrix3<T> &LHS, const Vector3<int> &RHS);

  template <typename T, class U>
  Vector3<T> operator*(const U &LHS, const Vector3<T> &RHS);

  template<typename T, int dim1, int dim2, int flag1>
  Vector3<T> operator*(const Eigen::Matrix<T, dim1, dim2, flag1> &LHS, const Vector3<T> &RHS);

  template <typename T>
  double triple_prod(const Vector3<T> &vec0, const Vector3<T> &vec1, const Vector3<T> &vec2);

  template <typename T>
  Vector3<T> cross_prod(const Vector3<T> &vec0, const Vector3<T> &vec1);

  template <typename T>
  std::ostream &operator<<(std::ostream &out, const Vector3<T> &vec_out);


  template <typename T>
  std::istream &operator>>(std::istream &in, Vector3<T> &RHS);



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  template <typename T, Index N>
  void Matrix<T, N>::initialize_static_members() {
    if(static_is_initialized)
      return;

    static_is_initialized = true;

    for(Index i = 0; i < (identity_mat.num_rows() < identity_mat.num_cols() ? identity_mat.num_rows() : identity_mat.num_cols()); i++)
      identity_mat(i, i) = 1.0;


  }

  // Accessors
  template <typename T, Index N>
  T &Matrix<T, N>::operator()(Index i, Index j) {
    assert((i * N + j) > -1 && (i * N + j) < size());
    return Vals[i * N + j];
  }

  // Accessors
  template <typename T, Index N>
  const T &Matrix<T, N>::operator()(Index i, Index j) const {
    assert((i * N + j) > -1 && (i * N + j) < size());
    return Vals[i * N + j];
  }

  template <typename T, Index N>
  const T &Matrix<T, N>::operator[](Index i) const {
    assert(i > -1 && i < size());
    return Vals[i];
  }

  template <typename T, Index N>
  T &Matrix<T, N>::operator[](Index i) {
    assert(i > -1 && i < size());
    return Vals[i];
  }

  template <typename T, Index N>
  T &Matrix<T, N>::at(Index i, Index j) {
    assert((i * N + j) > -1 && (i * N + j) < size());
    return Vals[i * N + j];
  }

  template <typename T, Index N>
  T &Matrix<T, N>::at(Index i) {
    assert(i > -1 && i < size());
    return Vals[i];
  }

  template <typename T, Index N>
  const T &Matrix<T, N>::at(Index i, Index j) const {
    assert((i * N + j) > -1 && (i * N + j) < size());
    return Vals[i * N + j];
  }

  template <typename T, Index N>
  const T &Matrix<T, N>::at(Index i) const {
    assert(i > -1 && i < size());
    return Vals[i];
  }

  template <typename T, Index N>
  void Matrix<T, N>::fill(T fill_val) {
    for(Index i = 0; i < size(); i++)
      at(i) = fill_val;
  }


  // sizes
  template <typename T, Index N>
  Index Matrix<T, N>::num_rows() const {
    return N;
  }

  template <typename T, Index N>
  Index Matrix<T, N>::num_cols() const {
    return N;
  }

  template <typename T, Index N>
  Index Matrix<T, N>::size() const {
    return N * N;
  }


  // Self-modifying Algebraic Operators
  template <typename T, Index N>
  Matrix<T, N> &Matrix<T, N>::operator+=(const Matrix<T, N> &RHS) {
    for(Index i = 0; i < num_rows(); i++) {
      for(Index j = 0; j < num_cols(); j++) {
        at(i, j) += RHS.at(i, j);
      }
    }
    return *this;
  }

  template <typename T, Index N>
  Matrix<T, N> &Matrix<T, N>::operator-=(const Matrix<T, N> &RHS) {
    for(Index i = 0; i < num_rows(); i++) {
      for(Index j = 0; j < num_cols(); j++) {
        at(i, j) -= RHS.at(i, j);
      }
    }
    return *this;
  }
  template <typename T, Index N>
  Matrix<T, N> &Matrix<T, N>::operator*=(const Matrix<T, N> &RHS) {

    //// IS THERE ANY WAY TO DO THIS IN PLACE !?!?!?!


  }

  template <typename T, Index N>
  Matrix<T, N> &Matrix<T, N>::operator*=(const T &RHS) {
    for(Index i = 0; i < num_rows()*num_cols(); i++)
      at(i) *= RHS;
    return *this;

  }

  template <typename T, Index N>
  Matrix<T, N> &Matrix<T, N>::operator/=(const T &RHS) {
    for(Index i = 0; i < num_rows()*num_cols(); i++)
      at(i) /= RHS;
    return *this;
  }


  // Algebraic Operators
  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::operator+(const Matrix &RHS) const {
    Matrix<T, N> tMat = *this;
    return tMat += RHS;

  }

  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::operator-(const Matrix &RHS) const {
    Matrix<T, N> tMat = *this;
    return tMat -= RHS;

  }

  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::operator*(const T &RHS) const {
    Matrix<T, N> tMat = *this;
    return tMat *= RHS;

  }

  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::operator/(const T &RHS) const {
    Matrix<T, N> tMat = *this;
    return tMat /= RHS;

  }

  //Comparison Operators
  //*******************************************************************************************
  /**
   * Overloads == operator for Matrix comparisons
   *
   * The == operator for Matrices checks that each
   * of the elements in the two Matrices being
   * compared are equal.
   * @return true if matrices are equal
   */
  //*******************************************************************************************
  template <typename T, Index N>
  bool Matrix<T, N>::operator==(const Matrix &RHS) const {
    for(Index i = 0; i < size(); i++) {
      if(!(at(i) == RHS.at(i)))
        return false;
    }

    return true;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  bool Matrix<T, N>::is_equal(const Matrix &RHS, double tol_val) const {
    for(Index i = 0; i < size(); i++) {
      if(std::abs(at(i) - RHS.at(i)) > tol_val)
        return false;
    }

    return true;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  bool Matrix<T, N>::is_unitary(double tol) const {
    if(almost_equal(determinant(), 1.0, tol)) {
      //      std::cout << "checking for unitarity.  Matrix times transpose is: \n" << transpose()*(*this);
      return (transpose() * (*this)).is_identity(tol);
    }
    else
      return false;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  bool Matrix<T, N>::is_identity(double tol) const {

    for(Index i = 0; i < num_rows(); i++) {
      for(Index j = 0; j < num_cols(); j++) {
        if((i == j && almost_equal(at(i, j), 1.0, tol)) || (i != j && almost_zero(at(i, j), tol))) return false;
      }
    }

    return true;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  bool Matrix<T, N>::is_integer(double tol) const {                 //better to take a tolerance isntead of always using TOL?
    for(Index i = 0; i < size(); i++) {
      if(std::abs(round(at(i)) - at(i)) > tol)
        return false;
    }
    return true;

  }
  //*******************************************************************************************

  template <typename T, Index N>
  bool Matrix<T, N>::is_zero(double tol) const {     //zero_tol never used
    for(Index i = 0; i < size(); i++) {
      if(std::abs(at(i)) > tol)
        return false;
    }
    return true;

  }

  //*******************************************************************************************

  // Essential Matrix Operations
  template <typename T, Index N>
  T Matrix<T, N>::determinant() const {

  }

  template <typename T, Index N>
  T Matrix<T, N>::trace() const {
    T temp = 0;

    for(Index i = 0; i < size(); i++) {
      temp += at(i, i);
    }

    return temp;
  }

  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::inverse() const {

  }

  template <typename T, Index N>
  Matrix<T, N> Matrix<T, N>::transpose() const {
    Matrix<T, N> tMat;

    for(Index i = 0; i < num_rows(); i++) {
      for(Index j = 0; j < num_cols(); j++) {
        tMat(i, j) = at(j, i);
      }
    }

    return tMat;
  }

  // ****************************************************************************************************

  template <typename T, Index N>
  jsonParser &to_json(const Matrix<T, N> &value, jsonParser &json) {
    json.put_array();
    for(Index i = 0; i < value.num_rows(); i++) {
      jsonParser json_row;
      json_row.put_array();
      for(Index j = 0; j < value.num_cols(); j++) {
        json_row.push_back(value(i, j));
      }
      json.push_back(json_row);
    }
    return json;
  }

  // ****************************************************************************************************

  template <typename T, Index N>
  void from_json(Matrix<T, N> &value, const jsonParser &json) {
    if(json.size() != N) {
      std::cerr << "Matrix size: " << N << std::endl;
      std::cerr << "JSON: " << json << std::endl;
      throw std::runtime_error("Error reading Matrix from JSON");
    }
    for(Index i = 0; i < value.num_rows(); i++) {
      if(json[i].size() != N) {
        std::cerr << "Matrix size: " << N << std::endl;
        std::cerr << "JSON: " << json << std::endl;
        throw std::runtime_error("Error reading Matrix from JSON");
      }
      for(Index j = 0; j < value.num_cols(); j++) {
        from_json(value.at(i, j), json[i][j]);
      }
    }
  }

  template <typename T, Index N>
  std::ostream &operator<<(std::ostream &stream, const Matrix<T, N> &out_matrix) {

    std::ios::fmtflags old_flags = stream.flags();
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    int twide = stream.width();

    for(Index i = 0; i < out_matrix.num_rows(); i++) {
      for(Index j = 0; j < out_matrix.num_cols(); j++) {
        stream << " " << std::setw(twide) << out_matrix.at(i, j);
      }
      stream << std::endl;
    }
    stream.flags(old_flags);
    return stream;
  }


  template <typename T, Index N>
  std::istream &operator>>(std::istream &stream, Matrix<T, N> &in_matrix) {
    for(Index i = 0; i < in_matrix.num_rows(); i++) {
      for(Index j = 0; j < in_matrix.num_cols(); j++) {
        stream >> in_matrix.at(i, j);
      }
    }
    return stream;
  }

  template <typename T, class U, Index N>
  Matrix<U, N> operator*(const Matrix<T, N> &LHS, const Matrix<U, N> &RHS) {
    Matrix<U, N> tMat(0);
    for(Index i = 0; i < LHS.num_rows(); i++) {
      for(Index j = 0; j < RHS.num_cols(); j++) {
        for(Index k = 0; k < LHS.num_cols(); k++) {
          tMat(i, j) += LHS.at(i, k) * RHS.at(k, j);
        }
      }
    }
    return tMat;
  }

  template <typename T, Index N>
  Matrix<T, N> operator*(const Matrix<T, N> &LHS, const Matrix<int, N> &RHS) {
    Matrix<T, N> tMat(0);
    for(Index i = 0; i < LHS.num_rows(); i++) {
      for(Index j = 0; j < RHS.num_cols(); j++) {
        for(Index k = 0; k < LHS.num_cols(); k++) {
          tMat(i, j) += LHS.at(i, k) * RHS.at(k, j);
        }
      }
    }
    return tMat;
  }

  template <typename T, Index N>
  Matrix<T, N> operator*(T LHS, const Matrix<T, N> &RHS) {
    Matrix<T, N> tMat(RHS);
    return tMat *= LHS;

  }


  template<typename T> template<int dim1, int dim2, int flag1>
  Matrix3<T>::Matrix3(const Eigen::Matrix<T, dim1, dim2, flag1> &_mat) {
    if(_mat.cols() != 3 && _mat.rows() != 3) {
      std::cerr << "CRITICAL ERROR: Attempting to construct a Matrix3<T> from an Eigen::Matrix of size (" << _mat.rows() << ", " << _mat.cols() << ").\n"
                << "                Exiting...\n";
      exit(1);
    }
    for(Index i = 0; i < 3; i++)
      for(Index j = 0; j < 3; j++)
        at(i, j) = _mat(i, j);
  }

  //*******************************************************************************************
  // Calculates and returns the determinant of a matrix.
  //*******************************************************************************************

  template <typename T>
  T Matrix3<T>::determinant() const {

    return (at(0, 0) * (at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)) -
            at(0, 1) * (at(1, 0) * at(2, 2) - at(1, 2) * at(2, 0)) +
            at(0, 2) * (at(1, 0) * at(2, 1) - at(1, 1) * at(2, 0)));
  }

  //*******************************************************************************************
  // Calculates and returns the trace of a matrix.
  //*******************************************************************************************

  template <typename T>
  T Matrix3<T>::trace() const {
    return at(0, 0) + at(1, 1) + at(2, 2);
  }

  //*******************************************************************************************
  // Returns the square of the Frobenius norm of a matrix.
  //*******************************************************************************************

  template <typename T>
  double Matrix3<T>::square_norm() const {
    return at(0, 0) * at(0, 0) + at(0, 1) * at(0, 1) + at(0, 2) * at(0, 2)
           + at(1, 0) * at(1, 0) + at(1, 1) * at(1, 1) + at(1, 2) * at(1, 2)
           + at(2, 0) * at(2, 0) + at(2, 1) * at(2, 1) + at(2, 2) * at(2, 2);
  }

  //*******************************************************************************************
  // Returns the Frobenius norm of a matrix.
  //*******************************************************************************************

  template <typename T>
  double Matrix3<T>::norm() const {
    return sqrt(square_norm());
  }

  //*******************************************************************************************
  // Calculates and returns the inverse of a matrix.
  //*******************************************************************************************

  template <typename T>
  Matrix3<T> Matrix3<T>::inverse() const {
    Matrix3<T> invMatrix;

    double det = determinant();
    invMatrix(0, 0) = (at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)) / det;
    invMatrix(0, 1) = (at(0, 2) * at(2, 1) - at(0, 1) * at(2, 2)) / det;
    invMatrix(0, 2) = (at(0, 1) * at(1, 2) - at(0, 2) * at(1, 1)) / det;
    invMatrix(1, 0) = (at(1, 2) * at(2, 0) - at(1, 0) * at(2, 2)) / det;
    invMatrix(1, 1) = (at(0, 0) * at(2, 2) - at(0, 2) * at(2, 0)) / det;
    invMatrix(1, 2) = (at(0, 2) * at(1, 0) - at(0, 0) * at(1, 2)) / det;
    invMatrix(2, 0) = (at(1, 0) * at(2, 1) - at(1, 1) * at(2, 0)) / det;
    invMatrix(2, 1) = (at(0, 1) * at(2, 0) - at(0, 0) * at(2, 1)) / det;
    invMatrix(2, 2) = (at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0)) / det;

    return invMatrix;

  }

  //*******************************************************************************************
  // Returns the adjugate matrix, defined adj(M)=det(M)*inv(M) -- useful for integer matrices, since it is always integer
  template <typename T>
  Matrix3<T> Matrix3<T>::adjugate() const {
    Matrix3<T> adjMatrix;

    adjMatrix(0, 0) = (+at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1));
    adjMatrix(0, 1) = (-at(0, 1) * at(2, 2) + at(0, 2) * at(2, 1));
    adjMatrix(0, 2) = (+at(0, 1) * at(1, 2) - at(0, 2) * at(1, 1));
    adjMatrix(1, 0) = (-at(1, 0) * at(2, 2) + at(1, 2) * at(2, 0));
    adjMatrix(1, 1) = (+at(0, 0) * at(2, 2) - at(0, 2) * at(2, 0));
    adjMatrix(1, 2) = (-at(0, 0) * at(1, 2) + at(0, 2) * at(1, 0));
    adjMatrix(2, 0) = (+at(1, 0) * at(2, 1) - at(1, 1) * at(2, 0));
    adjMatrix(2, 1) = (-at(0, 0) * at(2, 1) + at(0, 1) * at(2, 0));
    adjMatrix(2, 2) = (+at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0));

    return adjMatrix;

  }

  //*******************************************************************************************

  template <typename T>
  Matrix3<T> Matrix3<T>::transpose() const {
    Matrix3<T> tMatrix;
    tMatrix(0, 0) = at(0, 0);
    tMatrix(0, 1) = at(1, 0);
    tMatrix(0, 2) = at(2, 0);
    tMatrix(1, 0) = at(0, 1);
    tMatrix(1, 1) = at(1, 1);
    tMatrix(1, 2) = at(2, 1);
    tMatrix(2, 0) = at(0, 2);
    tMatrix(2, 1) = at(1, 2);
    tMatrix(2, 2) = at(2, 2);
    return tMatrix;
  }

  //*******************************************************************************************
  template <typename T>
  Matrix3<T> Matrix3<T>::elementary_hermite_op(T a, T b, T i, T j) {
    Matrix3<T> tmat(Matrix3<T>::identity());
    T tgcf = extended_gcf(a, b, tmat(i, i), tmat(i, j));
    if(!tgcf) {
      tmat = Matrix3<T>::identity();
      return tmat;
    }
    tmat(j, i) = -b / tgcf;
    tmat(j, j) = a / tgcf;
    return tmat;
  }

  //*******************************************************************************************
  // Get smith normal form of a 3x3 integer matrix.  Currently only works for full rank matrices
  // Decomposes *this into U*S*V, U,S,V have integer elements and det(U)=1, det(V)=1, det(S)=det(*this)
  // matrix S is diagonal and (S(i+1,i+1)%S(i,i) == 0); U and V have integer matrix inverses
  // Adapted from Matlab implementation written by John Gilbert (gilbert@parc.xerox.com)
  //   -> http://www.mathworks.com/matlabcentral/newsreader/view_thread/13728

  template <typename T>
  void Matrix3<T>::smith_normal_form(Matrix3<T> &U, Matrix3<T> &S, Matrix3<T> &V) const {
    assert(std::abs(at(0, 1)) < 100);
    U = V = Matrix3<T>::identity();
    S = *this;
    Matrix3<T> tmat = U;

    int i, j;

    // Bidiagonalize S with elementary Hermite transforms.
    for(j = 0; j < 3; j++) {
      //take j as column index and zero elements of j below the diagonal
      for(i = j + 1; i < 3; i++) {
        if(S(i, j) == 0) continue;
        tmat = elementary_hermite_op(S(j, j), S(i, j), j, i);
        S = tmat * S;
        U = U * tmat.inverse();
      }

      //take j as row index and zero elements of j after superdiagonal
      for(i = j + 2; i < 3; i++) {
        if(S(j, i) == 0) continue;
        tmat = elementary_hermite_op(S(j, j + 1), S(j, i), j + 1, i);
        S = S * tmat.transpose();
        V = tmat.transpose().inverse() * V;
      }
    }

    // std::cout << "About to hunt zeros, SNF decomposition is currently:\n" << U << "\n\n" << S << "\n\n" << V << "\n\n";
    while(!S.is_diagonal()) {
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
        T q = S(b, b + 1) / S(b, b);
        if(S(b, b + 1) % S(b, b) < 0) q -= 1;
        tmat = Matrix3<T>::identity();
        tmat(b + 1, b) = -q;
        S = S * tmat.transpose();
        V = tmat.transpose().inverse() * V;
        //std::cout << "tmat for q:\n" << tmat << '\n';
        //std::cout << "S after q:\n" << S << '\n';
      }
      else {
        tmat = Matrix3<T>::identity();
        tmat(b, b) = 0;
        tmat(b, b + 1) = 1;
        tmat(b + 1, b + 1) = 0;
        tmat(b + 1, b) = 1;
        S = S * tmat.transpose();
        V = tmat.transpose().inverse() * V;

      }

      if(!S(b, b + 1)) continue;
      tmat = elementary_hermite_op(S(b, b), S(b, b + 1), b, b + 1);
      S = S * tmat.transpose();
      V = tmat.transpose().inverse() * V;
      for(j = 0; j < 2; j++) {
        tmat = elementary_hermite_op(S(j, j), S(j + 1, j), j, j + 1);
        S = tmat * S;
        U = U * tmat.inverse();
        if(j + 2 >= 3) continue;
        tmat = elementary_hermite_op(S(j, j + 1), S(j, j + 2), j + 1, j + 2);
        S = S * tmat.transpose();
        V = tmat.transpose().inverse() * V;
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
          S.swap_rows(i, j);
          S.swap_cols(i, j);
          U.swap_cols(i, j);
          V.swap_rows(i, j);
        }
      }
    }
    //enforce divisibility condition:
    for(i = 0; i < 3; i++) {
      if(S(i, i) == 0) continue;
      for(j = i + 1; j < 3; j++) {
        if(S(j, j) % S(i, i) == 0) continue;
        //Replace S(i,i), S(j,j) by their gcd and lcm respectively.
        tmat = Matrix3<T>::identity();
        Matrix3<T> tmat2(tmat);
        T a(S(i, i)), b(S(j, j)), c, d, tgcf;
        tgcf = extended_gcf(a, b, c, d);
        tmat(i, i) = 1;
        tmat(i, j) = d;
        tmat(j, i) = -b / tgcf;
        tmat(j, j) = (a * c) / tgcf;

        tmat2(i, i) = c;
        tmat2(i, j) = 1;
        tmat2(j, i) = -(b * d) / tgcf;
        tmat2(j, j) = a / tgcf;
        S = tmat * S * tmat2.transpose();
        U = U * tmat.inverse();
        V = tmat2.transpose().inverse() * V;
      }
    }
    return;
  }

  //*******************************************************************************************

  template <typename T>
  Matrix3<T> &Matrix3<T>::operator*=(T RHS) {
    for(int i = 0; i < 9; i++)
      at(i) *= RHS;
    return *this;
  }

  //*******************************************************************************************

  template <typename T>
  Matrix3<T> &Matrix3<T>::operator=(const Eigen::Matrix<T, 3, 3> &RHS) {
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
        at(i, j) = RHS(i, j);

    return *this;
  }

  //*******************************************************************************************
  template <typename T>
  Matrix3<T> &Matrix3<T>::operator-=(const Matrix3<T> &RHS) {
    for(int i = 0; i < 9; i++)
      at(i) -= RHS[i];
    return *this;
  }

  //*******************************************************************************************

  template <typename T>
  Matrix3<T> &Matrix3<T>::operator+=(const Matrix3<T> &RHS) {
    for(int i = 0; i < 9; i++)
      at(i) += RHS[i];
    return *this;
  }

  //*******************************************************************************************

  template <typename T>
  Matrix3<T> Matrix3<T>::operator-(const Matrix3<T> &RHS) const {
    return Matrix3<T>(*this) -= RHS;
  }

  //*******************************************************************************************

  template <typename T>
  Matrix3<T> Matrix3<T>::operator+(const Matrix3<T> &RHS) const {
    return Matrix3<T>(*this) += RHS;
  }

  //*******************************************************************************************
  template <typename T>
  Matrix3<T>::operator Eigen::Matrix<T, 3, 3>() const {
    Eigen::Matrix<T, 3, 3> tEMat;;
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
        tEMat(i, j) = at(i, j);
    return tEMat;
  }

  //*******************************************************************************************
  template <typename T>
  Matrix3<T>::operator Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>() const {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tEMat(3, 3);
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
        tEMat(i, j) = at(i, j);
    return tEMat;
  }

  //*******************************************************************************************

  template <typename T>
  bool Matrix3<T>::is_unitary(double compare_tol) const {
    if(compare_tol < std::abs(std::abs(determinant()) - 1))
      return false;

    return transpose().is_equal(inverse(), compare_tol);
  }

  //*******************************************************************************************

  template <typename T>
  bool Matrix3<T>::is_diagonal(double compare_tol) const {
    return (almost_zero(at(0, 1), compare_tol)
            && almost_zero(at(1, 0), compare_tol)
            && almost_zero(at(2, 1), compare_tol)
            && almost_zero(at(1, 2), compare_tol)
            && almost_zero(at(2, 0), compare_tol)
            && almost_zero(at(0, 2), compare_tol));
  }

  //*******************************************************************************************

  template <typename T>
  Vector3< std::complex<double> > Matrix3<T>::eigen() const {
    Vector3< std::complex<double> > t_eigenvals;
    Eigen::EigenSolver<Eigen::Matrix<T, 3, 3> > tES((Eigen::Matrix<T, 3, 3>) *this);
    t_eigenvals[0] = tES.eigenvalues()[0];
    t_eigenvals[1] = tES.eigenvalues()[1];
    t_eigenvals[2] = tES.eigenvalues()[2];
    return t_eigenvals;
  }

  //*******************************************************************************************

  template <typename T>
  Vector3< std::complex<double> > Matrix3<T>::eigen(Vector3< Vector3< std::complex<double > > > &eigen_vecs) const {

    Vector3< std::complex<double> > t_eigenvals;
    Eigen::EigenSolver<Eigen::Matrix<T, 3, 3> > tES((Eigen::Matrix<T, 3, 3>) *this);

    for(int i = 0; i < 3; i++) {
      t_eigenvals[i] = tES.eigenvalues()[i];
      for(int j = 0; j < 3; j++)
        eigen_vecs[i][j] = tES.eigenvectors()(j, i);
    }
    return t_eigenvals;

  }

  //*******************************************************************************************

  template <typename T>
  void Matrix3<T>::swap_rows(Index i, Index j) {
    T tval(at(i, 0));
    at(i, 0) = at(j, 0);
    at(j, 0) = tval;

    tval = at(i, 1);
    at(i, 1) = at(j, 1);
    at(j, 1) = tval;

    tval = at(i, 2);
    at(i, 2) = at(j, 2);
    at(j, 2) = tval;
    return;
  }
  //*******************************************************************************************

  template <typename T>
  void Matrix3<T>::swap_cols(Index i, Index j) {
    T tval(at(0, i));
    at(0, i) = at(0, j);
    at(0, j) = tval;

    tval = at(1, i);
    at(1, i) = at(1, j);
    at(1, j) = tval;

    tval = at(2, i);
    at(2, i) = at(2, j);
    at(2, j) = tval;
    return;
  }

  //*******************************************************************************************

  template <typename T, class U>
  Matrix3<U> operator*(const Matrix3<T> &LHS, const Matrix3<U> &RHS) {
    Matrix3<U> tMat(0);
    for(Index i = 0; i < LHS.num_rows(); i++) {
      for(Index j = 0; j < RHS.num_cols(); j++) {
        for(Index k = 0; k < LHS.num_cols(); k++) {
          tMat(i, j) += LHS.at(i, k) * RHS.at(k, j);
        }
      }
    }
    return tMat;
  }

  template <typename T>
  Matrix3<T> operator*(const Matrix3<T> &LHS, const Matrix3<int> &RHS) {
    Matrix3<T> tMat(0);
    for(Index i = 0; i < LHS.num_rows(); i++) {
      for(Index j = 0; j < RHS.num_cols(); j++) {
        for(Index k = 0; k < LHS.num_cols(); k++) {
          tMat(i, j) += LHS.at(i, k) * T(RHS.at(k, j));
        }
      }
    }
    return tMat;
  }

  template <typename T>
  Matrix3<T> operator*(T LHS, const Matrix3<T> &RHS) {
    return Matrix3<T>(RHS) *= LHS;
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  //*******************************************************************************************
  /*
  //THIS IS A DOT PRODUCT
  template< typename T, Index N>
  T Vector<T,N>::operator*(const Vector &RHS) const{
  T Prod=0;
  for(int i=0; i<N; i++)
  Prod+=at(i)*RHS.at(i);
  return Prod;
  }
  */
  //***********************************define operator /, *, -,+=,-=,==,is_equal*********************
  template< typename T, Index N>
  Vector<T, N>::Vector(const T &init_val) {
    for(Index i = 0; i < size(); i++)
      at(i) = init_val;
    return;
  }

  template< typename T, Index N>
  Vector<T, N> Vector<T, N>::operator/(const T &divisor) const {
    Vector<T, N> tVec;
    for(Index i = 0; i < size(); i++)
      tVec[i] = at(i) / divisor;
    return tVec;
  }

  //*******************************************************************************************
  /*
     template <typename T, Index N>
     Vector<T,N> Vector<T,N>::operator*(const T &scale) const{
     Vector<T,N> tVec;
     for(Index i=0; i<N; i++)
     tVec[i] = at(i)*scale;
     return tVec;
     }
     */
  //*******************************************************************************************
  /*
     template< typename T, Index N>
     Vector<T,N>& Vector<T,N>::operator=(const T &init_val){
     for(Index i=0; i<size(); i++)
     at(i)=init_val;
     return *this;
     }
     */
  //*******************************************************************************************

  template< typename T, Index N>
  Vector<T, N> Vector<T, N>::operator+(const Vector<T, N> &RHS) const {
    Vector<T, N> tVec;
    for(Index i = 0; i < size(); i++)
      tVec[i] = at(i) + RHS.at(i);
    return tVec;
  }

  //*******************************************************************************************

  template< typename T, Index N>
  Vector<T, N> Vector<T, N> :: operator-(const Vector<T, N> &RHS) const {
    Vector<T, N> tVec;
    for(Index i = 0; i < size(); i++)
      tVec[i] = at(i) - RHS.at(i);
    return tVec;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  Vector<T, N> &Vector<T, N>::operator+=(const Vector<T, N> &RHS) {
    for(Index i = 0; i < size(); i++)
      at(i) += RHS.at(i);
    return *this;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  Vector<T, N> &Vector<T, N>::operator-=(const Vector<T, N> &RHS) {
    for(Index i = 0; i < size(); i++)
      at(i) -= RHS.at(i);
    return *this;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  bool Vector<T, N>::operator==(const Vector<T, N> &RHS) const {
    for(Index i = 0; i < size(); i++)
      if(at(i) != RHS.at(i)) return false;
    return true;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  bool Vector<T, N>::is_equal(const Vector<T, N> &RHS, double eq_tol) const {
    for(Index i = 0; i < size(); i++)
      if(std::abs(at(i) - RHS.at(i)) > eq_tol) return false;
    return true;
  }

  //*******************************************************************************************
  template <typename T, Index N>
  bool Vector<T, N>::is_zero(double zero_tol) const {
    for(Index i = 0; i < size(); i++) {
      if(std::abs(at(i)) > zero_tol)
        return false;
    }
    return true;
  }

  //*******************************************************************************************
  /**
   * Scales Vector by length len
   *
   */
  //*******************************************************************************************

  template <typename T, Index N>
  void Vector<T, N>::scale(double scale_factor) {
    for(Index i = 0; i < size(); i++)
      at(i) *= scale_factor;
    return;

  }


  //*******************************************************************************************

  template <typename T, Index N>
  double Vector<T, N>::length() const {
    double tlength = 0.0;
    for(Index i = 0; i < size(); i++) {
      tlength += at(i) * at(i);
    }
    tlength = sqrt(tlength);
    return tlength;
  }

  // ****************************************************************************************************
  /*
     template <typename T, class U, Index N>
     Vector<U,N> operator*(const Matrix<T,N> & LHS, const Vector<U,N> & RHS){

     Vector<U,N> tVec;
     Index i, j;
     for(i=0; i<LHS.num_rows(); i++){
     tVec[i]=0;
     for(j=0; j<LHS.num_cols(); j++){
     tVec[i]+=LHS.at(i,j)*RHS.at(j);
     }
     }
     return tVec;
     }

  // ****************************************************************************************************

  template <typename T, class U, Index N>
  Vector<T,N> operator*(const Matrix<T,N> & LHS, const Vector<int,N> & RHS){

  Vector<T,N> tVec;
  Index i, j;
  for(i=0; i<LHS.num_rows(); i++){
  tVec[i]=0;
  for(j=0; j<LHS.num_cols(); j++){
  tVec[i]+=LHS.at(i,j)*RHS.at(j);
  }
  }
  return tVec;
  }

  // ****************************************************************************************************

  template <typename T, class U, Index N>
  Vector<T,N> operator*(const U &LHS, const Vector<T,N> & RHS){
  Vector<T,N> tVec;
  for(Index i=0; i<RHS.size(); i++){
  tVec[i]=LHS*RHS.at(i);
  }
  return tVec;
  }
  */
  // ****************************************************************************************************

  template <typename T, Index N>
  jsonParser &to_json(const Vector<T, N> &value, jsonParser &json) {
    json.put_array();
    for(Index i = 0; i < value.size(); i++)
      json.push_back(value[i]);
    return json;
  }

  // ****************************************************************************************************

  template <typename T, Index N>
  void from_json(Vector<T, N> &value, const jsonParser &json) {
    if(json.size() != N) {
      std::cerr << "Vector size: " << N << std::endl;
      std::cerr << "JSON: " << json << std::endl;
      throw std::runtime_error("Error reading Vector from JSON");
    }
    for(Index i = 0; i < value.size(); i++) {
      from_json(value.at(i), json[i]);
    }
  }

  // ****************************************************************************************************

  template <typename T, Index N>
  std::ostream &operator<<(std::ostream &stream, const Vector<T, N> &vec_out) {
    std::ios::fmtflags old_flags = stream.flags();
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

    int twide = stream.width();
    if(twide < stream.precision())
      twide = stream.precision() + 3;
    for(Index i = 0; i < vec_out.size(); i++)
      stream << " " << std::setw(twide) << vec_out.at(i);

    stream.flags(old_flags);
    return stream;
  }

  //*******************************************************************************************

  template <typename T, Index N>
  std::istream &operator>>(std::istream &stream, Vector<T, N> &vec_in) {
    for(Index i = 0; i < vec_in.size(); i++)
      stream >> vec_in.at(i);
    return stream;
  }

  //*******************************************************************************************
  template<typename T> template<typename Derived>
  Vector3<T>::Vector3(const Eigen::MatrixBase<Derived> &_vec) {
    assert((_vec.cols() == 3 && _vec.rows() == 1) || (_vec.cols() == 1 && _vec.rows() == 3));
    for(Index i = 0; i < 3; i++)
      at(i) = _vec[i];
    return;
  }
  //*******************************************************************************************

  template <typename T>
  bool Vector3<T>::is_zero(double zero_tol) const {
    for(int i = 0; i < 3; i++) {
      if(!almost_zero(at(i), zero_tol))
        return false;
    }
    return true;
  }

  //*******************************************************************************************

  template <typename T>
  Vector3<T> Vector3<T>::cross(const Vector3 &RHS) const {
    Vector3<T> cross_vec;
    cross_vec[0] = at(1) * RHS[2] - at(2) * RHS[1];
    cross_vec[1] = at(2) * RHS.at(0) - at(0) * RHS[2];
    cross_vec[2] = at(0) * RHS[1] - at(1) * RHS.at(0);
    return cross_vec;
  }

  //*******************************************************************************************

  template <typename T> template <class U>
  U Vector3<T>::dot(const Vector3<U> &RHS) const {
    return at(0) * RHS[0] + at(1) * RHS[1] + at(2) * RHS[2];
  }

  //*******************************************************************************************

  /*
     template <typename T> template <class U>
     Vector3<T>::operator Vector3<U>() {
     return Vector3<U>(dynamic_cast<U>(at(0)),
     dynamic_cast<U>(at(1)),
     dynamic_cast<U>(at(2)));
     }
     */

  template <typename T> template <class U>
  Vector3<T>::operator Vector3<U>() const {
    return Vector3<U>(static_cast<U>(at(0)),
                      static_cast<U>(at(1)),
                      static_cast<U>(at(2)));
  }


  template <typename T> template <class U>
  Vector3<T> Vector3<T>::operator*(const U &RHS) const {
    return Vector3<T>(*this) *= RHS;
  }

  //*******************************************************************************************

  template <typename T>
  double Vector3<T>::norm() const {
    return sqrt(dot(*this));
  }

  //*******************************************************************************************

  ///return angle, in radians, between 0 and pi that describe separation in direction of two vectors
  template <typename T>
  double Vector3<T>::get_angle(const Vector3 &RHS) const {
    T tempval = (dot(RHS) / (norm() * RHS.norm()));

    if(almost_zero(tempval - 1))        //John G 121212     Rounding disaster averted
      return 0;

    else if(almost_zero(tempval + 1))
      return M_PI;

    else
      return acos(tempval);
  }

  //*******************************************************************************************

  ///return signed angle, in radians, between -pi and pi that describe separation in direction of two vectors
  template <typename T>
  double Vector3<T>::get_signed_angle(const Vector3 &RHS, const Vector3 &pos_ref) const {
    double angle = get_angle(RHS);

    if(pos_ref.dot(cross(RHS)) < 0) {
      angle = -angle;
    }

    return angle;
  }

  //*******************************************************************************************
  template <typename T>
  T Vector3<T>::max() const {
    T tmax = CASM::max(at(0), at(1));
    return CASM::max(tmax, at(2));
  }

  //*******************************************************************************************

  template <typename T>
  T Vector3<T>::min() const {
    T tmin = CASM::min(at(0), at(1));
    return CASM::min(tmin, at(2));
  }

  //*******************************************************************************************
  template <typename T>
  void Vector3<T>::scale_by_max_mag() {
    //Determine what the largest value in the vector is and divide everything by it.
    int largest = -1;
    T tlarge = 0.0;

    for(int i = 0; i < 3; i++) {
      if(almost_zero(at(i)) || std::abs(tlarge) > std::abs(at(i))) continue;
      tlarge = at(i);
      largest = i;
    }
    if(largest == -1) return;

    (*this) /= tlarge;
    return;

  }

  //*******************************************************************************************
  ///Take a vector of doubles, and multiply by some factor that turns it into a vector of integers (within a tolerance)
  template <typename T>
  Vector3< int > Vector3<T>::scale_to_int(double tolerance) {
    Vector3< int > ints;
    Vector3< double > dubs(*this);

    dubs.scale_by_max_mag();

    for(int i = 0; i < 3; i++) {
      if(almost_zero(dubs[i])) {
        dubs[i] = 0.0;
        ints[i] = 0;
      }
    }

    //Determine what the smallest non zero value is and divide all indeces by that value
    //so that one miller index=1 and the others are abs()>1
    int smallest = -1;
    double tsmall = 2; // all values are >=1
    for(int i = 0; i < 3; i++) {
      if(almost_zero(dubs[i]) || fabs(tsmall) < fabs(dubs[i])) continue;

      tsmall = dubs[i];
      smallest = i;
    }

    if(smallest != -1)
      dubs = dubs / tsmall;


    //We want to multiply the miller indeces by some factor such that all indeces become integers.
    //In order to do this we pick a tolerance to work with and round the miller indeces if they are close
    //enough to the integer value (e.g. 2.95 becomes 3). Choosing a tolerance that is too small will
    //result in the "primitive-slab" blowing up.

    //Begin choosing a factor and multiply all indeces by it (starting with 1). Then round the non-smallest
    //miller indeces (smallest index requires no rounding, since it will always be a perfect
    //integer thanks to the previous division).
    //Next take absolute value of difference between rounded indeces and actual values (int_diff 1 & 2).
    //If the difference for both indeces is smaller than the tolerance then you've reached the desired
    //accuracy and the rounded indeces can be used to construct the "primitive-slab" cell. If not, increase the
    //factor by 1 and try again, until the tolerance is met.
    bool within_tol = false;
    int factor = 1;

    Vector3< double > tdubs(dubs);
    int i;
    while(!within_tol) {
      dubs = tdubs * factor;

      for(i = 0; i < 3; i++) {
        if(tolerance < std::abs(round(dubs[i]) - dubs[i])) break;
      }

      if(i < 3) {
        factor++;
        continue;
      }

      for(i = 0; i < 3; i++)
        ints[i] = round(dubs[i]);

      within_tol = true;
    }

    return ints;
  }

  //John G 121106
  //*******************************************************************************************
  //Return cross product matrix of vector
  template <typename T>
  Matrix3<T> Vector3<T>::get_cross_product_mat() const {
    Matrix3<T> crossprodmat(0);

    crossprodmat.at(0, 1) = -at(2);
    crossprodmat.at(0, 2) = at(1);
    crossprodmat.at(1, 0) = at(2);
    crossprodmat.at(1, 2) = -at(0);
    crossprodmat.at(2, 0) = -at(1);
    crossprodmat.at(2, 1) = at(0);

    return crossprodmat;
  }

  //*******************************************************************************************
  //Return tensor product of two vectors
  template <typename T>
  Matrix3<T> Vector3<T>::get_tensor_product(const Vector3<T> &tensvec) const {
    Matrix3<T> tensorprodmat;

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        tensorprodmat.at(i, j) = at(i) * tensvec[j];
      }
    }

    return tensorprodmat;
  }

  //*******************************************************************************************
  //Return tesnsor product of itself
  template <typename T>
  Matrix3<T> Vector3<T>::get_tensor_product() const {
    Vector3<T> tvec = *this;

    return get_tensor_product(tvec);
  }

  //*******************************************************************************************
  //http://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
  //Return rotation matrix given an angle in radians
  //This might be better as type double instead of type T
  template <typename T>
  Matrix3<T> Vector3<T>::get_rotation_mat(double angle) const {
    Vector3<T> tvec = (*this);
    tvec.normalize();

    Matrix3<T> crossmat, tensormat;
    tensormat = tvec.get_tensor_product();
    crossmat = tvec.get_cross_product_mat();
    return cos(angle) * Matrix3<T>::identity() + sin(angle) * crossmat + (1 - cos(angle)) * tensormat;
  }

  //*******************************************************************************************
  //\John G

  //get random vector, evenly distributed within unit cube.
  template <typename T>
  void Vector3<T>::rand(int &rand_seed) {
    at(0) = 2.0 * ran0(rand_seed) - 1;
    at(1) = 2.0 * ran0(rand_seed) - 1;
    at(2) = 2.0 * ran0(rand_seed) - 1;
  }
  //*******************************************************************************************

  //get random vector, evenly distributed within unit sphere.
  //uses a try and discard approach, which is potentially faster than a constrained picking approach
  template <typename T>
  void Vector3<T>::rand_ball(int &rand_seed) {
    rand(rand_seed);
    while(dot(*this) > 1) {
      rand(rand_seed);
    }
    return;
  }

  //*******************************************************************************************

  template <typename T>
  void Vector3<T>::normalize() {
    T mag = norm();

    if(mag == 0) {
      return;
    }

    at(0) /= mag;
    at(1) /= mag;
    at(2) /= mag;
    return;
  }

  //*******************************************************************************************

  template <typename T>
  template <int dim1, int dim2, int flag1>
  void Vector3<T>::insert_as_col(Eigen::Matrix<T, dim1, dim2, flag1> &mat, int ind) {
    for(int i = 0; i < 3; i++)
      mat(i, ind) = at(i);
    return;
  }

  //*******************************************************************************************

  template <typename T>
  template <int dim1, int dim2, int flag1>
  void Vector3<T>::insert_as_row(Eigen::Matrix<T, dim1, dim2, flag1> &mat, int ind) {
    for(int i = 0; i < 3; i++)
      mat(ind, i) = at(i);
    return;
  }

  //*******************************************************************************************

  template <typename T>
  Vector3<T> &Vector3<T>::operator=(const Vector3<T> &RHS) {
    at(0) = RHS[0];
    at(1) = RHS[1];
    at(2) = RHS[2];
    return *this;
  }

  //*******************************************************************************************

  template <typename T>
  Vector3<T> &Vector3<T>::operator/=(const T &RHS) {
    at(0) /= RHS;
    at(1) /= RHS;
    at(2) /= RHS;
    return *this;
  }
  //*******************************************************************************************

  template <typename T>
  inline
  Vector3<T> &Vector3<T>::operator+=(const Vector3 &RHS) {
    at(0) += RHS[0];
    at(1) += RHS[1];
    at(2) += RHS[2];
    return *this;
  }
  //*******************************************************************************************

  template <typename T>
  inline
  Vector3<T> &Vector3<T>::operator-=(const Vector3 &RHS) {
    at(0) -= RHS[0];
    at(1) -= RHS[1];
    at(2) -= RHS[2];
    return *this;
  }
  //*******************************************************************************************

  template <typename T> template <class U>
  Vector3<T> &Vector3<T>::operator*=(const U &RHS) {
    at(0) *= RHS;
    at(1) *= RHS;
    at(2) *= RHS;
    return *this;
  }

  //*******************************************************************************************

  template <typename T>
  Vector3<T> Vector3<T>::operator/(const T &RHS) const {
    Vector3<T> tVec(*this);
    return tVec /= RHS;
  }

  //*******************************************************************************************

  template <typename T>
  Vector3<T> Vector3<T>::operator+(const Vector3 &RHS) const {
    return Vector3<T>(*this) += RHS;
  }

  //*******************************************************************************************

  template <typename T>
  Vector3<T> Vector3<T>::operator-(const Vector3 &RHS) const {
    return Vector3<T>(*this) -= RHS;
  }

  //*******************************************************************************************

  template <typename T>
  Vector3<T> Vector3<T>::operator-() const {
    Vector3<T> tvec(*this);

    tvec[0] = -tvec[0];
    tvec[1] = -tvec[1];
    tvec[2] = -tvec[2];
    return tvec;
  }

  //*******************************************************************************************

  template < typename T> template<int dim1, int dim2, int flag1>
  Vector3<T>::operator Eigen::Matrix<T, dim1, dim2, flag1>() const {
    Eigen::Matrix<T, dim1, dim2, flag1> tvec(3, 1);
    tvec(0, 0) = at(0);
    tvec(1, 0) = at(1);
    tvec(2, 0) = at(2);
    return tvec;
  }

  //****************************************************************************************************

  template <typename T>
  Vector3<T> operator*(const Matrix3<T> &LHS, const Vector3<int> &RHS) {

    Vector3<T> tVec;
    Index i, j;
    for(i = 0; i < LHS.num_rows(); i++) {
      tVec[i] = 0;
      for(j = 0; j < LHS.num_cols(); j++) {
        tVec[i] += LHS.at(i, j) * RHS.at(j);
      }
    }
    return tVec;
  }

  //*******************************************************************************************

  template <typename T, class U>
  Vector3<U> operator*(const Matrix3<T> &LHS, const Vector3<U> &RHS) {

    Vector3<U> tVec;

    Index i, j;
    for(i = 0; i < LHS.num_rows(); i++) {
      tVec[i] = 0;
      for(j = 0; j < LHS.num_cols(); j++) {
        tVec[i] += LHS.at(i, j) * RHS.at(j);
      }
    }
    return tVec;
  }


  //*******************************************************************************************

  template <typename T, class U>
  Vector3<T> operator*(const U &LHS, const Vector3<T> &RHS) {
    return Vector3<T>(RHS) *= LHS;
  }

  //*******************************************************************************************

  template<typename T, int dim1, int dim2, int flag1>
  Vector3<T> operator*(const Eigen::Matrix<T, dim1, dim2, flag1> &LHS, const Vector3<T> &RHS) {

    return Vector3<T>(LHS(0, 0) * RHS[0] + LHS(0, 1) * RHS[1] + LHS(0, 2) * RHS[2],
                      LHS(1, 0) * RHS[0] + LHS(1, 1) * RHS[1] + LHS(1, 2) * RHS[2],
                      LHS(2, 0) * RHS[0] + LHS(2, 1) * RHS[1] + LHS(2, 2) * RHS[2]);
  }
  //*******************************************************************************************

  template <typename T>
  std::ostream &operator<<(std::ostream &stream, const Vector3<T> &vec_out) {
    int twide = stream.width();
    if(twide < stream.precision())
      twide = stream.precision() + 3;
    for(Index i = 0; i < vec_out.size(); i++)
      stream << std::setw(twide) << vec_out.at(i);
    return stream;
  }

  //*******************************************************************************************

  template <typename T>
  std::istream &operator>>(std::istream &stream, Vector3<T> &vec_in) {
    for(Index i = 0; i < vec_in.size(); i++) {
      stream >> vec_in.at(i);
    }
    return stream;
  }

  //*******************************************************************************************

  template <typename T>
  double triple_prod(const Vector3<T> &vec0, const Vector3<T> &vec1, const Vector3<T> &vec2) {
    return vec0.dot(vec1.cross(vec2));
  }

  //*******************************************************************************************
  template <typename T>
  Vector3<T> cross_prod(const Vector3<T> &vec0, const Vector3<T> &vec1) {
    Vector3<T> cross_vec;
    cross_vec[0] = vec0[1] * vec1[2] - vec0[2] * vec1[1];
    cross_vec[1] = vec0[2] * vec1[0] - vec0[0] * vec1[2];
    cross_vec[2] = vec0[0] * vec1[1] - vec0[1] * vec1[0];
    return cross_vec;
  }

  //*******************************************************************************************
  template <typename T>
  bool almost_equal(const Vector3<T> &A, const Vector3<T> &B, double tol = TOL) {
    return almost_equal(A[0], B[0], tol)
           && almost_equal(A[1], B[1], tol)
           && almost_equal(A[2], B[2], tol);
  }
  //*******************************************************************************************
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
  //*******************************************************************************************
  template <typename T>
  Vector3<int> round(const Vector3<T> &A) {
    Vector3<int> iA;
    iA[0] = round(A[0]);
    iA[1] = round(A[1]);
    iA[2] = round(A[2]);
    return iA;
  }

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
