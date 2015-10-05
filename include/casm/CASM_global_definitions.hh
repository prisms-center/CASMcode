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

#define BOOST_NO_SCOPED_ENUMS
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/container/stable_vector.hpp>

class MTRand;

namespace CASM {

  class jsonParser;
  class MonteCarloConditions;

  namespace fs = boost::filesystem;
  namespace po = boost::program_options;

  using std::swap;

  typedef unsigned int uint;
  typedef unsigned long int ulint;
  typedef long int lint;

  ///For unsigned integer indexing:
  //typedef std::size_t Index;
  //const Index max_index = -4;
  //bool valid_index(Index i) {
  //  return i < max_index;
  //};


  typedef  Eigen::MatrixXd::Index EigenIndex;

  //tolerance
  const double TOL = 0.00001;

  //Boltzmann Constant
  const double KB = 8.6173423E-05; //eV/K

  //Planck's Constant
  const double PLANCK = 4.135667516E-15; //eV-s

  namespace multivector_impl {
    template<typename T, size_t N>
    struct multivector_tmp;

    template<typename T>
    struct multivector_tmp<T, 0> {
      using type = T;
    };
    template<typename T, size_t N>
    struct multivector_tmp {
      using type = std::vector < typename multivector_tmp < T, N - 1 >::type >;
    };
  }

  template<typename T>
  struct multivector {
    template<size_t N>
    using X = typename multivector_impl::multivector_tmp<T, N>::type;
  };


  template<class T>
  std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec) {
    if(vec.size() == 0)
      out << "[empty]  ";
    for(auto it = vec.cbegin(); it != vec.cend(); ++it) {
      out << *it << "  ";
    }
    return out;
  }

  ///For long integer indexing:
  typedef EigenIndex Index;
  bool valid_index(Index i);

  enum COORD_TYPE {FRAC = 0, CART = 1, COORD_DEFAULT = 2};
  std::istream &operator>>(std::istream &sin, COORD_TYPE &coord);
  jsonParser &to_json(const COORD_TYPE &value, jsonParser &json);
  void from_json(COORD_TYPE &value, const jsonParser &json);

  enum PERIODICITY_TYPE {PERIODIC = 0, LOCAL = 1, PERIODICITY_DEFAULT = 2};
  jsonParser &to_json(const PERIODICITY_TYPE &value, jsonParser &json);
  void from_json(PERIODICITY_TYPE &value, const jsonParser &json);

  enum CELL_TYPE {PRIM = 0, SCEL = 1};
  jsonParser &to_json(const CELL_TYPE &value, jsonParser &json);
  void from_json(CELL_TYPE &value, const jsonParser &json);

  jsonParser &to_json(const MTRand &twister, jsonParser &json);

  //jsonParser &to_json(const MonteCarloConditions &ref_conditions, jsonParser &fill_json);

  //void from_json(MTRand &twister, const jsonParser &json);

  enum COMPLEX_OUTPUT_TYPE {REAL = 0, IMAG = 1, COMPLEX = 2}; // Added by Ivy

  enum CASMfileTypes {TYPEFILE, TYPEDIR, TYPEOTHER, IOERR}; // Moved from FileSystemInterface.cc by Ivy 01/16/14

  //************************************************************

  void print_splash(std::ostream &out);

};


//*****************************************************************************************************//

//Templates need the implementation with the declaration!!
//http://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor
#include "casm/casm_io/jsonParser.hh"
namespace Eigen {
  template <typename Derived>
  CASM::jsonParser &to_json(const Eigen::MatrixBase<Derived> &value, CASM::jsonParser &json) {
    json.put_array();
    for(int i = 0; i < value.rows(); i++) {
      CASM::jsonParser json_row;
      json_row.put_array();
      for(int j = 0; j < value.cols(); j++) {
        json_row.push_back(value(i, j));
      }
      json.push_back(json_row);
    }
    return json;
  }

  template <typename Derived>
  void from_json(Eigen::MatrixBase<Derived>  &value, const CASM::jsonParser &json) {
    try {
      //Eigen::MatrixBase<Derived> &value = const_cast< MatrixBase<Derived>& >(value_);
      value.derived().resize(json.size(), json[0].size());
      for(int i = 0; i < value.rows(); i++) {
        for(int j = 0; j < value.cols(); j++) {
          from_json(value(i, j), json[i][j]);
        }
      }
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }
}

#endif
