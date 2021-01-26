#ifndef UNIT_COMMON_HH
#define UNIT_COMMON_HH

#include <iostream>
#include <string>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/eigen.hh"
// #include "FCCTernaryProj.hh"
// #include "ZrOProj.hh"

using namespace CASM;

namespace test {
/// \brief Check expected JSON vs calculated JSON using BOOST_CHECK_EQUAL
void check(std::string test, const jsonParser &expected,
           const jsonParser &calculated, fs::path test_cases_path, bool quiet,
           double tol = 0.0);

/// \brief Create a new project directory, appending ".(#)" to ensure
/// it is a new project
fs::path proj_dir(fs::path init);

template <typename Container1DType>
void print_computed_result(std::ostream &sout, std::string name,
                           const Container1DType &vec) {
  sout << name << " = {";
  auto it = vec.begin();
  auto end = vec.end();
  while (it != end) {
    sout << *it;
    ++it;
    if (it != end) {
      sout << ", ";
    }
  }
  sout << "};" << std::endl;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> eigen_vector(std::initializer_list<T> v) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> eigen_v =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(v.size());
  int i = 0;
  for (int value : v) {
    eigen_v[i++] = value;
  }
  return eigen_v;
}

}  // namespace test

#endif
