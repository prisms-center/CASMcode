#ifndef CASM_GLOBAL_DEFINTIONS_HH
#define CASM_GLOBAL_DEFINTIONS_HH

#include <iostream>
#include <cmath>
#include <cstddef>
#include <complex>
#include <string>
#include <sstream>
#include <vector>


#include "casm/external/Eigen/Dense"
#include "casm/external/boost.hh"

namespace CASM {

  namespace fs = boost::filesystem;
  namespace po = boost::program_options;

  using std::swap;

  typedef unsigned int uint;
  typedef unsigned long int ulint;
  typedef long int lint;

  ///For integer indexing:
  typedef  Eigen::MatrixXd::Index EigenIndex;

  //tolerance
  const double TOL = 0.00001;

  //Boltzmann Constant
  const double KB = 8.6173423E-05; //eV/K

  //Planck's Constant
  const double PLANCK = 4.135667516E-15; //eV-s

  template<class T>
  std::istream &operator>>(std::istream &_in, std::vector<T> &vec) {
    std::string line;
    std::getline(_in, line, '\n');
    std::stringstream tss(line);
    T tval;
    while(_in) {
      _in >> tval;
      vec.push_back(tval);
    }

    return _in;
  }

  ///For long integer indexing:
  typedef EigenIndex Index;
  bool valid_index(Index i);

  enum COORD_TYPE {FRAC = 0, CART = 1, COORD_DEFAULT = 2};
  std::istream &operator>>(std::istream &sin, COORD_TYPE &coord);



  enum PERIODICITY_TYPE {PERIODIC = 0, LOCAL = 1, PERIODICITY_DEFAULT = 2};

  enum CELL_TYPE {PRIM = 0, SCEL = 1};

  enum COMPLEX_OUTPUT_TYPE {REAL = 0, IMAG = 1, COMPLEX = 2}; // Added by Ivy

  enum CASMfileTypes {TYPEFILE, TYPEDIR, TYPEOTHER, IOERR}; // Moved from FileSystemInterface.cc by Ivy 01/16/14

  //************************************************************

  void print_splash(std::ostream &out);

  /// Apply a transformation, in place
  /// - Default is equivalent to \code obj.apply_sym(f, args...) \endcode
  template<typename Object, typename Transform, typename...Args>
  Object &apply(const Transform &f, Object &obj, Args &&...args) {
    return obj.apply_sym(f, std::forward<Args>(args)...);
  }

  /// Copy and apply a transformation
  /// - We can also include a specialization for cloneable objects, via SFINAE
  template<typename Object, typename Transform, typename...Args>
  Object copy_apply(const Transform &f, Object obj, Args &&...args) {
    return apply(f, obj, std::forward<Args>(args)...);
  }

};

namespace Eigen {

  typedef Matrix<long int, 3, 3> Matrix3l;
  typedef Matrix<long int, 3, 1> Vector3l;
  typedef Matrix<long int, Dynamic, Dynamic> MatrixXl;
  typedef Matrix<long int, Dynamic, 1> VectorXl;

  template<typename Derived>
  std::istream &operator >> (std::istream &s, MatrixBase<Derived> &m) {
    for(int i = 0; i < m.rows(); ++i)
      for(int j = 0; j < m.cols(); j++)
        s >> m(i, j);
    return s;
  }
}


namespace std {
  template<class T>
  ostream &operator<<(ostream &out, const vector<T> &vec) {
    if(vec.size() == 0)
      out << "[empty]  ";
    for(auto it = vec.cbegin(); it != vec.cend(); ++it) {
      out << *it << "  ";
    }
    return out;
  }
}


#endif
