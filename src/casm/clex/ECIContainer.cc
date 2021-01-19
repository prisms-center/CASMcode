#include "casm/clex/ECIContainer.hh"

#include <fstream>

#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {

/// \brief Evaluate property given an ECIContainer and correlations
double operator*(const ECIContainer &_eci, const Eigen::VectorXd &_corr) {
  double result(0);
  auto ind_it(_eci.index().cbegin()), ind_end(_eci.index().cend());
  auto eci_it(_eci.value().cbegin());
  for (; ind_it != ind_end; ++ind_it, ++eci_it)
    result += (*eci_it) * _corr[*ind_it];
  return result;
}

/// \brief Evaluate property given an ECIContainer and pointer to beginning of
/// range of correlation
double operator*(const ECIContainer &_eci, double const *_corr_begin) {
  double result(0);
  auto ind_it(_eci.index().cbegin()), ind_end(_eci.index().cend());
  auto eci_it(_eci.value().cbegin());
  while (ind_it != ind_end) {
    result += (*eci_it) * (*(_corr_begin + *ind_it));
    ++ind_it;
    ++eci_it;
  }
  return result;
}

}  // namespace CASM
